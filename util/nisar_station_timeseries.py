"""
Shared pipeline for building per-station NISAR GCOV backscatter time series
against ISMN soil moisture stations, used by the moistx/02_get_NISAR_*_data.ipynb
notebooks. Each notebook stays responsible for its own AOI-specific bits (loading
its own stations.geojson, picking its own output directory, and the
earthaccess/fsspec session setup) and calls compute_station_terrain() followed by
process_station() per station.
"""

import os
import queue
import re
import sys
import threading
import time
from typing import Any, Callable

import aiohttp
import earthaccess
import fsspec
import h5py
import numpy as np
import pandas as pd
import pyproj
import rasterio
from rasterio.windows import Window

# pyproj's own bundled PROJ data directory can be missing/broken in this env
# (CRSError: "no database context specified" on any pyproj.CRS(...)/
# Transformer.from_crs(...) call) -- the PROJ_DATA env var alone doesn't fix
# this, since pyproj prefers its own bundled copy over the env var unless
# explicitly overridden. Point it at the conda-forge env's own canonical
# proj.db instead (same one rasterio/GDAL use via their own activation
# scripts), which is intact. Done at import time, here, so every caller of
# this module (both the 02_get_NISAR_*_data.ipynb notebooks and, via
# nisar_gslc_windows.py, notebook 03) is covered without needing its own
# fix -- rather than relying on rasterio/GDAL happening to initialize the
# shared native PROJ context first, which is a matter of import order.
pyproj.datadir.set_data_dir(os.path.join(sys.prefix, "Library", "share", "proj"))


GCOV_SHORT_NAME = "NISAR_L2_GCOV_PROVISIONAL_V1"

CRID_RE = re.compile(r"^(NISAR_L2_PR_GCOV(?:_[^_]+){9})_([A-Z]\d+)_")
GRANULE_ID_RE = re.compile(
    r"^NISAR_L2_PR_GCOV_\d+_(?P<track>\d+)_[AD]_\d+_\d+_\w+_[A-Z]_(?P<start>\d{8}T\d{6})_"
)

ADJACENT_PASS_GAP_SECONDS = 60
MIN_ACQUISITIONS_PER_TRACK = 3
REFERENCE_INCIDENCE_ANGLE_DEG = 40.0

MAX_FETCH_RETRIES = 5
RETRY_BACKOFF_SECONDS = 5
RETRYABLE_STATUSES = {500, 502, 503, 504}

COPERNICUS_DEM_BASE_URL = "https://copernicus-dem-30m.s3.amazonaws.com"

BACKSCATTER_WINDOW_HALF = 1  # 3x3 window

FETCH_TOTAL_TIMEOUT_SECONDS = 60
FETCH_CONNECT_TIMEOUT_SECONDS = 15


def get_fsspec_https_session_with_timeout(
    total_seconds: int = FETCH_TOTAL_TIMEOUT_SECONDS,
    sock_connect_seconds: int = FETCH_CONNECT_TIMEOUT_SECONDS,
) -> fsspec.AbstractFileSystem:
    """
    Same bearer-token HTTPS session as earthaccess.get_fsspec_https_session(), but
    with an explicit aiohttp timeout. That function (earthaccess/store.py's
    get_fsspec_session) builds its aiohttp.ClientSession with only headers/trust_env
    in client_kwargs -- no timeout -- so a stalled fetch (a connection accepted but
    that stops sending data mid-stream; the DAAC's presigned S3/CloudFront URLs have
    been observed doing this under load) falls back to aiohttp's own default
    (ClientTimeout(total=300, sock_connect=30)) at best, and in practice can hang far
    longer once fsspec's own background block-prefetch/retry queuing is factored in --
    observed hanging 17+ minutes on a single HVHV window read before being manually
    interrupted. A short explicit timeout here makes a stalled fetch fail predictably
    instead of hanging indefinitely; paired with retrying on aiohttp.ClientError (not
    just ClientResponseError) below, so a timeout gets the same retry-with-backoff
    treatment as an HTTP error status instead of crashing outright.

    skip_instance_cache=True: HTTPFileSystem is `cachable` (fsspec's `_Cached`
    metaclass), so calling this function twice with the same arguments would
    otherwise hand back the *same* underlying instance -- including its one lazily-
    created aiohttp.ClientSession -- shared across every caller. Combined with
    run_with_timeout's documented behavior of leaving a timed-out extraction's
    thread running in the background rather than cancelling it, that shared instance
    is what let two genuinely concurrent granule reads run on the same session:
    diagnosed after saved GSLC window files turned up with a `native_id` attr naming
    a different granule than that same file's own `acquisition_datetime_utc` /
    `track_number` / `orbit_pass_direction` attrs (see process_station /
    process_station_gslc below, which now open a fresh instance per granule instead
    of sharing one across a station's whole granule loop).
    """
    token = earthaccess.__auth__.token["access_token"]
    client_kwargs = {
        "headers": {"Authorization": f"Bearer {token}"},
        "trust_env": False,
        "timeout": aiohttp.ClientTimeout(total=total_seconds, sock_connect=sock_connect_seconds),
    }
    return fsspec.filesystem("https", client_kwargs=client_kwargs, skip_instance_cache=True)


# Confirmed (against a real GCOV granule) that the per-request timeout above works
# fine on its own -- a raw small read completes in ~3s. The actual cost is that
# h5py's file-open needs many small, *sequential* metadata reads (superblock, group
# headers, chunk B-trees, ...) to build its view of the file, each individually
# fast but adding up to several minutes for some granules, with no single request
# ever slow enough to trip the per-request timeout. A per-request timeout can't
# bound that; only an external wall-clock deadline around the whole call can.
PER_GRANULE_TIMEOUT_SECONDS = 180


def run_with_timeout(func: Callable, args: tuple, timeout_seconds: float) -> tuple[bool, Any]:
    """
    Run func(*args) in a fresh daemon thread, waiting up to timeout_seconds for a
    result. There's no SIGALRM on Windows, so a background thread is the standard
    way to bound a call (like extract_granule_record / extract_and_save_gslc_window)
    that has no way to be cancelled from outside once started.

    A *fresh* thread per call, not a shared bounded pool (the first version of this
    used ThreadPoolExecutor(max_workers=4) + future.result(timeout=...)): a timed-out
    call's thread keeps running in the background, since Python can't forcibly kill a
    thread. With a bounded pool, enough timeouts permanently exhaust all workers --
    observed in practice: exactly 4 granules timed out (matching max_workers=4), and
    every single granule submitted after that also "timed out" after 180s without
    ever actually starting, since none could get a free worker. A fresh thread per
    call has no shared resource to exhaust. daemon=True means even if many calls
    time out over a long run, the leaked threads don't block the interpreter/kernel
    from exiting normally (a non-daemon thread would keep it alive indefinitely,
    since concurrent.futures.ThreadPoolExecutor's own threads are non-daemon and
    joined at exit).

    func, args: the call to run
    timeout_seconds: how long to wait before giving up on it
    returns: (True, func(*args)) on success within the deadline, or (False, None) on
        timeout. An exception raised by func itself propagates from this call (not
        swallowed), same as calling func(*args) directly, as long as it's raised
        within the deadline.
    """
    result_queue: queue.Queue = queue.Queue(maxsize=1)

    def _target():
        try:
            result_queue.put(("ok", func(*args)))
        except Exception as e:  # noqa: BLE001 -- re-raised on the calling thread below
            result_queue.put(("error", e))

    threading.Thread(target=_target, daemon=True).start()
    try:
        status, value = result_queue.get(timeout=timeout_seconds)
    except queue.Empty:
        return False, None
    if status == "error":
        raise value
    return True, value


OUTPUT_COLUMNS = [
    "datetime_utc",
    "NISAR_L_HH_db_mean3x3", "NISAR_L_HH_db_std3x3",
    "NISAR_L_HV_db_mean3x3", "NISAR_L_HV_db_std3x3",
    "NISAR_local_incidence_angle_deg", "NISAR_orb_dir", "NISAR_orb_number",
]


class _SkipGranule(Exception):
    """Raised to abandon a granule that isn't usable, distinct from a transient fetch error."""


# --- Product version dedup / track filtering / same-pass merging -----------------

CRID_PARTS_RE = re.compile(r"^([A-Za-z]*)(\d*)$")


def _crid_sort_key(crid: str) -> tuple[str, int]:
    """
    Numeric-aware sort key for a CRID like "D00340" (letter prefix + digits): a
    plain string comparison breaks once the digit run's width differs (e.g.
    "D100" would sort *before* "D99", since '1' < '9' lexicographically) even
    though 100 > 99 numerically.
    """
    m = CRID_PARTS_RE.match(crid)
    letters, digits = m.groups() if m else (crid, "")
    return (letters, int(digits) if digits else -1)


def keep_latest_versions(results: list, crid_re: re.Pattern = CRID_RE) -> list:
    """
    Keep only the highest-CRID (Composite Release Identifier) granule per unique
    acquisition, parsed from each result's native-id, so reprocessed granules that
    CMR returns under multiple CRIDs don't create duplicate rows.

    results: earthaccess search_data() results
    crid_re: native-id regex capturing (product, crid); defaults to GCOV's, pass
        a product-specific pattern (e.g. GSLC) for other product lines that share
        the same native-id token structure
    returns: results with older-CRID duplicates removed
    """
    best = {}
    for result in results:
        native_id = result["meta"]["native-id"]
        m = crid_re.match(native_id)
        product, crid = m.groups() if m else (native_id, "")
        if product not in best or _crid_sort_key(crid) > _crid_sort_key(best[product][0]):
            best[product] = (crid, result)
    return [crid_result[1] for crid_result in best.values()]


def filter_tracks_with_min_acquisitions(
    results: list,
    min_acquisitions: int = MIN_ACQUISITIONS_PER_TRACK,
    gap_seconds: int = ADJACENT_PASS_GAP_SECONDS,
    granule_id_re: re.Pattern = GRANULE_ID_RE,
) -> list:
    """
    Keep only results whose relative orbit/track covers the station at least
    min_acquisitions times in the search window, parsed from native-id's track
    number and start time (no download needed). Acquisitions from the same
    overpass split across adjacent frames (within gap_seconds) count as one.

    results: earthaccess search_data() results (post version-dedup)
    min_acquisitions: minimum acquisitions a track must have to be kept
    gap_seconds: max gap between acquisitions still counted as the same overpass
    granule_id_re: native-id regex capturing (track, start); defaults to GCOV's,
        pass a product-specific pattern for other product lines
    returns: results belonging to tracks with enough acquisitions
    """
    rows = []
    for i, result in enumerate(results):
        m = granule_id_re.match(result["meta"]["native-id"])
        if not m:
            continue
        rows.append({
            "i": i,
            "track": int(m.group("track")),
            "start": pd.to_datetime(m.group("start"), format="%Y%m%dT%H%M%S"),
        })
    if not rows:
        return []  # no native-id matched (e.g. no granules at all for this station)
    meta_df = pd.DataFrame(rows)

    kept_indices = set()
    for _, group in meta_df.groupby("track"):
        group = group.sort_values("start")
        n_acquisitions = 1 + (group["start"].diff().dt.total_seconds() > gap_seconds).sum()
        if n_acquisitions >= min_acquisitions:
            kept_indices.update(group["i"])

    return [results[i] for i in sorted(kept_indices)]


def merge_adjacent_acquisitions(
    station_df: pd.DataFrame, gap_seconds: int = ADJACENT_PASS_GAP_SECONDS
) -> pd.DataFrame:
    """
    Collapse same-pass adjacent-frame rows (gap to the previous acquisition <=
    gap_seconds) into one row per physical overpass. NISAR_L_*_db_mean3x3 is
    averaged in linear power (dB values can't be arithmetically averaged
    directly); NISAR_L_*_db_std3x3 is averaged directly (a simple approximation —
    same-pass boundary cases are rare enough that a statistically rigorous pooled
    combination isn't worth the extra complexity here).

    station_df: rows for one station, with columns matching OUTPUT_COLUMNS
    gap_seconds: max gap between rows still counted as the same overpass
    returns: one row per physical overpass
    """
    df = station_df.copy()
    df["_dt"] = pd.to_datetime(df["datetime_utc"])
    df = df.sort_values("_dt").reset_index(drop=True)

    group_id = (df["_dt"].diff().dt.total_seconds() > gap_seconds).cumsum()

    def _mean_db(db_values):
        power = 10 ** (db_values.dropna() / 10)
        return float(10 * np.log10(power.mean())) if len(power) else np.nan

    def _mean_plain(values):
        valid = values.dropna()
        return float(valid.mean()) if len(valid) else np.nan

    merged = df.groupby(group_id).agg(
        datetime_utc=("_dt", "mean"),
        NISAR_L_HH_db_mean3x3=("NISAR_L_HH_db_mean3x3", _mean_db),
        NISAR_L_HH_db_std3x3=("NISAR_L_HH_db_std3x3", _mean_plain),
        NISAR_L_HV_db_mean3x3=("NISAR_L_HV_db_mean3x3", _mean_db),
        NISAR_L_HV_db_std3x3=("NISAR_L_HV_db_std3x3", _mean_plain),
        NISAR_local_incidence_angle_deg=("NISAR_local_incidence_angle_deg", "mean"),
        NISAR_orb_dir=("NISAR_orb_dir", "first"),
        NISAR_orb_number=("NISAR_orb_number", "first"),
    ).reset_index(drop=True)

    merged["datetime_utc"] = merged["datetime_utc"].dt.strftime("%Y-%m-%dT%H:%M:%S.%f")
    return merged


def already_in_csv(
    result, existing_df: pd.DataFrame | None, gap_seconds: int = ADJACENT_PASS_GAP_SECONDS
) -> bool:
    """
    Identify a granule from its native-id alone (track + start time) to check
    whether it's already saved in a station's CSV, without opening/streaming
    anything. A CSV row matches if it's the same track and within gap_seconds of
    the granule's start time (same tolerance used by merge_adjacent_acquisitions,
    since a merged row's datetime_utc is the mean of its constituent times).

    result: a single earthaccess search_data() result
    existing_df: the station's existing CSV, with a parsed "_dt" datetime column, or None
    gap_seconds: same-overpass tolerance, see merge_adjacent_acquisitions
    returns: True if this granule is already represented in existing_df
    """
    if existing_df is None or existing_df.empty:
        return False
    m = GRANULE_ID_RE.match(result["meta"]["native-id"])
    if not m:
        return False
    track = int(m.group("track"))
    start = pd.to_datetime(m.group("start"), format="%Y%m%dT%H%M%S")
    same_track = existing_df.loc[existing_df["NISAR_orb_number"] == track, "_dt"]
    if same_track.empty:
        return False
    return (same_track - start).abs().dt.total_seconds().le(gap_seconds).any()


# --- DEM-derived terrain (local incidence angle inputs) ---------------------------

def _copernicus_dem_tile_url(lat: float, lon: float) -> str:
    """Copernicus GLO-30 tiles are named by their SW-corner integer degree."""
    lat_tile = int(np.floor(lat))
    lon_tile = int(np.floor(lon))
    ns = "N" if lat_tile >= 0 else "S"
    ew = "E" if lon_tile >= 0 else "W"
    tile = f"Copernicus_DSM_COG_10_{ns}{abs(lat_tile):02d}_00_{ew}{abs(lon_tile):03d}_00_DEM"
    return f"/vsicurl/{COPERNICUS_DEM_BASE_URL}/{tile}/{tile}.tif"


def _open_dem_tile(lat: float, lon: float, tile_cache: dict[str, Any]) -> Any:
    url = _copernicus_dem_tile_url(lat, lon)
    if url not in tile_cache:
        tile_cache[url] = rasterio.open(url)
    return tile_cache[url]


def slope_aspect_from_dem(lon: float, lat: float, tile_cache: dict[str, Any]) -> tuple[float, float, float]:
    """
    Local terrain slope/aspect (Horn's method) and elevation from a 3x3 Copernicus
    GLO-30 window centered on (lon, lat).

    lon, lat: station coordinates
    tile_cache: dict of open DEM tile datasets, reused across calls sharing a tile
    returns: (slope_deg, aspect_deg, elevation_m); aspect is the downhill compass
        bearing, clockwise from north; elevation_m is the window's center pixel
    """
    src = _open_dem_tile(lat, lon, tile_cache)
    row, col = src.index(lon, lat)
    row = min(max(row, 1), src.height - 2)  # keep the 3x3 window inside this tile
    col = min(max(col, 1), src.width - 2)
    window = Window(col_off=col - 1, row_off=row - 1, width=3, height=3)
    z = src.read(1, window=window).astype(float)

    # GLO-30 is posted in geographic coordinates (EPSG:4326); convert pixel spacing
    # from degrees to meters (longitude spacing shrinks with cos(latitude)) so the
    # slope gradient comes out in real-world meters/meter, not degrees/pixel.
    dx_deg, dy_deg = src.res
    dx_m = dx_deg * 111_320.0 * np.cos(np.radians(lat))
    dy_m = dy_deg * 111_320.0

    a, b, c = z[0]
    d, e, f = z[1]
    g, h, i = z[2]
    dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * dx_m)
    dzdy = ((a + 2 * b + c) - (g + 2 * h + i)) / (8 * dy_m)

    slope_deg = float(np.degrees(np.arctan(np.hypot(dzdx, dzdy))))
    # Aspect is the downhill bearing, i.e. the compass direction of -gradient =
    # (-dzdx, -dzdy) in (East, North); both components must be negated here, not
    # just the East one.
    aspect_deg = float(np.degrees(np.arctan2(-dzdx, -dzdy)) % 360)
    return slope_deg, aspect_deg, float(e)  # e is the window's center pixel


def surface_normal_enu(slope_deg: float, aspect_deg: float) -> np.ndarray:
    """
    Unit surface normal in local East-North-Up coordinates. Horizontal component
    points in the *downhill* (aspect) direction (e.g. a south-facing cliff's
    normal points due south).

    slope_deg, aspect_deg: from slope_aspect_from_dem
    returns: 3-element unit vector [East, North, Up]
    """
    b = np.radians(slope_deg)
    a = np.radians(aspect_deg)
    return np.array([np.sin(b) * np.sin(a), np.sin(b) * np.cos(a), np.cos(b)])


def compute_station_terrain(stations: list[dict]) -> None:
    """
    Compute and attach DEM-derived terrain to each station in place:
    dem_slope_deg, dem_aspect_deg, dem_elevation_m, surface_normal_enu. On a DEM
    lookup failure for a given station these are set to NaN (surface_normal_enu to
    an all-NaN vector), which propagates cleanly through downstream local
    incidence angle computation rather than raising.

    stations: list of station dicts, each with "name", "lon", "lat"; mutated in place
    """
    tile_cache: dict[str, Any] = {}
    print("Computing local terrain slope/aspect from Copernicus GLO-30 for each station...")
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
        for station in stations:
            try:
                slope_deg, aspect_deg, dem_elevation_m = slope_aspect_from_dem(
                    station["lon"], station["lat"], tile_cache
                )
                normal = surface_normal_enu(slope_deg, aspect_deg)
            except Exception as e:
                print(f"  WARNING: DEM lookup failed for {station['name']}: {e}")
                slope_deg = aspect_deg = dem_elevation_m = np.nan
                normal = np.array([np.nan, np.nan, np.nan])
            station["dem_slope_deg"] = slope_deg
            station["dem_aspect_deg"] = aspect_deg
            station["dem_elevation_m"] = dem_elevation_m
            station["surface_normal_enu"] = normal

        for tile_ds in tile_cache.values():
            tile_ds.close()

    slopes = [s["dem_slope_deg"] for s in stations if np.isfinite(s["dem_slope_deg"])]
    n_failed = len(stations) - len(slopes)
    print(f"Done. {len(tile_cache)} unique DEM tiles used for {len(stations)} stations"
          f" ({n_failed} failed lookups)." if n_failed else
          f"Done. {len(tile_cache)} unique DEM tiles used for {len(stations)} stations.")
    print(f"Slope range across stations: {min(slopes):.2f} - {max(slopes):.2f} deg (median {np.median(slopes):.2f} deg)")


# --- Incidence-angle calibration ---------------------------------------------------

def calibrate_to_reference_angle(
    sigma0_db: float, theta_local_deg: float, theta_ref_deg: float = REFERENCE_INCIDENCE_ANGLE_DEG
) -> float:
    """
    Cosine power law (n=1, the standard surface/Bragg-scattering normalization)
    against the *local* incidence angle, referenced to theta_ref_deg. Only used to
    build the *_calibrated columns; the raw NISAR_L_HH_db/NISAR_L_HV_db columns are
    left as measured so this can be recomputed with a different model or reverted.

    sigma0_db: raw backscatter in dB
    theta_local_deg: local (terrain-corrected) incidence angle in degrees
    theta_ref_deg: reference angle to normalize to
    returns: calibrated backscatter in dB, or NaN if either input is NaN
    """
    if not (np.isfinite(sigma0_db) and np.isfinite(theta_local_deg)):
        return np.nan
    return sigma0_db + 10 * np.log10(
        np.cos(np.radians(theta_ref_deg)) / np.cos(np.radians(theta_local_deg))
    )


# --- Per-granule extraction ---------------------------------------------------------

def _interp_radar_grid_field(
    rg, field_name: str, rg_h: np.ndarray, yi_rg: int, xi_rg: int, elevation_m: float
) -> float:
    """
    Nearest (y, x) cell, linear interpolation over height only, for a radarGrid
    (height, y, x) field. Reads only the (at most two) height levels needed
    instead of the full ~43 MB/field cube. Returns NaN without touching rg_h if
    elevation_m is NaN (e.g. from a failed DEM lookup for this station).
    """
    if not np.isfinite(elevation_m):
        return np.nan
    if elevation_m <= rg_h[0]:
        return float(rg[field_name][0, yi_rg, xi_rg])
    if elevation_m >= rg_h[-1]:
        return float(rg[field_name][-1, yi_rg, xi_rg])
    hi = int(np.searchsorted(rg_h, elevation_m))
    lo = hi - 1
    v_lo = float(rg[field_name][lo, yi_rg, xi_rg])
    v_hi = float(rg[field_name][hi, yi_rg, xi_rg])
    frac = (elevation_m - rg_h[lo]) / (rg_h[hi] - rg_h[lo])
    return v_lo + frac * (v_hi - v_lo)


def _window_mean_db(power_window: np.ndarray, valid_mask: np.ndarray) -> float:
    """
    Mean backscatter in dB over a window of linear power values, averaging in
    linear power first (dB values can't be arithmetically averaged directly),
    over only the pixels flagged valid.

    power_window: 2D array of raw (linear) backscatter power
    valid_mask: same-shape boolean array, True where the pixel should be included
    returns: mean backscatter in dB, or NaN if no pixel in the window is valid
    """
    valid_power = power_window[valid_mask]
    return float(10 * np.log10(valid_power.mean())) if valid_power.size else np.nan


def _window_std_db(power_window: np.ndarray, valid_mask: np.ndarray) -> float:
    """
    Sample standard deviation (ddof=1) of the per-pixel dB-converted values over
    a window, computed directly in the dB domain (unlike the mean, there's no
    clean power-domain equivalent of a spread statistic). NaN if fewer than 2
    valid pixels remain in the window.

    power_window: 2D array of raw (linear) backscatter power
    valid_mask: same-shape boolean array, True where the pixel should be included
    returns: sample std of dB values over the valid pixels, or NaN
    """
    valid_power = power_window[valid_mask]
    if valid_power.size < 2:
        return np.nan
    return float(np.std(10 * np.log10(valid_power), ddof=1))


def _extract_window_stats(
    group, yi: int, xi: int, ny: int, nx: int, half: int, has_hh_hv: bool
) -> tuple[float, float, float, float]:
    """
    Mean and standard deviation of HH/HV backscatter in dB over a
    (2*half+1)x(2*half+1) window centered on grid cell (yi, xi), clipped to
    whatever part of the window falls inside the frequency A grid extent (a
    station near an edge may get fewer valid pixels). Within that clipped
    window, individual pixels with an invalid mask value (0 or 255) or
    non-positive power are excluded. The mean is computed by averaging in linear
    power then converting to dB; the std is computed directly on the per-pixel
    dB values (see _window_std_db). Returns all-NaN if the granule isn't
    HH/HV-polarized.

    group: the granule's /science/LSAR/GCOV/grids/frequencyA h5py group
    yi, xi: nearest-pixel grid indices of the station (row, col)
    ny, nx: frequency A grid dimensions
    half: window half-size (1 for a 3x3 window)
    has_hh_hv: whether this granule has HHHH/HVHV datasets
    returns: (hh_db_mean, hh_db_std, hv_db_mean, hv_db_std)
    """
    if not has_hh_hv:
        return np.nan, np.nan, np.nan, np.nan
    row0, row1 = max(yi - half, 0), min(yi + half + 1, ny)
    col0, col1 = max(xi - half, 0), min(xi + half + 1, nx)
    mask_window = group["mask"][row0:row1, col0:col1]
    valid = ~np.isin(mask_window, (0, 255))
    hh_window = group["HHHH"][row0:row1, col0:col1]
    hv_window = group["HVHV"][row0:row1, col0:col1]
    hh_valid = valid & (hh_window > 0)
    hv_valid = valid & (hv_window > 0)
    hh_mean_db = _window_mean_db(hh_window, hh_valid)
    hh_std_db = _window_std_db(hh_window, hh_valid)
    hv_mean_db = _window_mean_db(hv_window, hv_valid)
    hv_std_db = _window_std_db(hv_window, hv_valid)
    return hh_mean_db, hh_std_db, hv_mean_db, hv_std_db


def extract_granule_record(
    result,
    station: dict,
    station_normal: np.ndarray,
    fs,
    fsspec_config: dict,
    max_retries: int = MAX_FETCH_RETRIES,
    retry_backoff_seconds: int = RETRY_BACKOFF_SECONDS,
    retryable_statuses: set[int] = RETRYABLE_STATUSES,
) -> dict | None:
    """
    Open one GCOV granule remotely (fsspec + h5py, no download) and extract one
    row for `station`: HH/HV backscatter mean and standard deviation in dB over a
    3x3 pixel window centered on the nearest frequency A pixel (excluding
    individually invalid pixels — see _extract_window_stats), and the true local
    (terrain-corrected) incidence angle at the center pixel (combining GCOV's
    per-acquisition LOS vector with the station's DEM-derived surface normal).

    result: one earthaccess search_data() result
    station: station dict with "lon", "lat", "dem_elevation_m" (see compute_station_terrain)
    station_normal: station["surface_normal_enu"], passed explicitly to avoid
        recomputing per granule
    fs: an earthaccess fsspec HTTPS session
    fsspec_config: kwargs for fs.open (cache_type, block_size, ...)
    returns: a record dict (matching OUTPUT_COLUMNS minus datetime_utc handling
        already done), or None if the granule isn't usable for this station
        (missing frequencyA group, or station outside the grid extent) — distinct
        from a record with NaN backscatter (e.g. non-HH/HV granule, invalid pixel
        mask), which is still returned so orbit/angle metadata isn't lost
    """
    native_id = result["meta"]["native-id"]
    url = result.data_links(access="external")[0]

    for attempt in range(max_retries):
        try:
            with fs.open(url, **fsspec_config) as fh, h5py.File(fh, "r") as f:
                if "/science/LSAR/GCOV/grids/frequencyA" not in f:
                    raise _SkipGranule
                group = f["/science/LSAR/GCOV/grids/frequencyA"]

                ident = f["/science/LSAR/identification"]
                acq_time = ident["zeroDopplerStartTime"][()].decode()
                orb_dir = ident["orbitPassDirection"][()].decode()
                orb_number = int(ident["trackNumber"][()])  # relative orbit/track number, not absoluteOrbitNumber

                epsg = int(group["projection"][()])
                has_hh_hv = "HHHH" in group and "HVHV" in group

                # xCoordinates/yCoordinates are stored contiguously (not chunked) and are
                # perfectly regularly spaced, so the nearest-pixel index can be computed
                # arithmetically from the origin, spacing, and array length instead of a
                # full-array read.
                x_spacing = float(group["xCoordinateSpacing"][()])
                y_spacing = float(group["yCoordinateSpacing"][()])
                x_origin = float(group["xCoordinates"][0])
                y_origin = float(group["yCoordinates"][0])
                nx = group["xCoordinates"].shape[0]
                ny = group["yCoordinates"].shape[0]

                transformer = pyproj.Transformer.from_crs(4326, epsg, always_xy=True)
                x0, y0 = transformer.transform(station["lon"], station["lat"])

                xi = int(round((x0 - x_origin) / x_spacing))
                yi = int(round((y0 - y_origin) / y_spacing))

                if not (0 <= xi < nx and 0 <= yi < ny):
                    raise _SkipGranule  # station falls outside frequency A's grid extent

                # incidenceAngle/losUnitVectorX/losUnitVectorY are each a (height, y, x)
                # cube, chunked (1, 512, 512) — ~43 MB/field if read in full. Instead:
                # nearest grid cell in (y, x), and linear interpolation between the two
                # nearest height levels, via _interp_radar_grid_field.
                rg = f["/science/LSAR/GCOV/metadata/radarGrid"]
                rg_h = rg["heightAboveEllipsoid"][:]
                rg_x0, rg_x1 = rg["xCoordinates"][0:2]
                rg_y0, rg_y1 = rg["yCoordinates"][0:2]
                nx_rg = rg["xCoordinates"].shape[0]
                ny_rg = rg["yCoordinates"].shape[0]

                xi_rg = int(round((x0 - rg_x0) / (rg_x1 - rg_x0)))
                yi_rg = int(round((y0 - rg_y0) / (rg_y1 - rg_y0)))
                xi_rg = min(max(xi_rg, 0), nx_rg - 1)
                yi_rg = min(max(yi_rg, 0), ny_rg - 1)

                # DEM-derived elevation (compute_station_terrain), not the ISMN-reported
                # one: same pixel driving the surface normal, and more robust than ISMN
                # metadata. Can be NaN if that station's DEM lookup failed.
                elevation_m = station["dem_elevation_m"]

                # Reconstruct the LOS unit vector (target -> sensor) in local ENU:
                # East/North come straight from GCOV; Up = cos(ellipsoidal incidence
                # angle) by definition (incidenceAngle is measured from the ellipsoid
                # normal, which *is* "Up" in a local ENU frame), then renormalize since
                # independently height-interpolating 3 fields can leave tiny norm error.
                ellip_incidence = _interp_radar_grid_field(rg, "incidenceAngle", rg_h, yi_rg, xi_rg, elevation_m)
                los_x = _interp_radar_grid_field(rg, "losUnitVectorX", rg_h, yi_rg, xi_rg, elevation_m)
                los_y = _interp_radar_grid_field(rg, "losUnitVectorY", rg_h, yi_rg, xi_rg, elevation_m)
                los_z = np.cos(np.radians(ellip_incidence))
                los_vec = np.array([los_x, los_y, los_z])
                los_vec = los_vec / np.linalg.norm(los_vec)

                # True local (terrain-corrected) incidence angle: angle between the LOS
                # vector and the station's DEM-derived surface normal, not the ellipsoid
                # normal that ellip_incidence above is measured against.
                local_incidence_deg = float(np.degrees(np.arccos(
                    np.clip(np.dot(los_vec, station_normal), -1.0, 1.0)
                )))

                hh_mean_db, hh_std_db, hv_mean_db, hv_std_db = _extract_window_stats(
                    group, yi, xi, ny, nx, BACKSCATTER_WINDOW_HALF, has_hh_hv
                )

                return {
                    "datetime_utc": acq_time,
                    "NISAR_L_HH_db_mean3x3": hh_mean_db,
                    "NISAR_L_HH_db_std3x3": hh_std_db,
                    "NISAR_L_HV_db_mean3x3": hv_mean_db,
                    "NISAR_L_HV_db_std3x3": hv_std_db,
                    "NISAR_local_incidence_angle_deg": local_incidence_deg,
                    "NISAR_orb_dir": orb_dir,
                    "NISAR_orb_number": orb_number,
                }
        except _SkipGranule:
            return None
        except aiohttp.ClientError as e:
            # The CloudFront/S3 distribution behind these presigned URLs occasionally
            # answers a byte-range fetch with a transient 503 (seen mid-stream, not
            # just on open) under sustained request volume; retrying the whole granule
            # with backoff clears it. Widened from just ClientResponseError (an HTTP
            # error status) to the broader ClientError so a stalled/timed-out fetch --
            # which has no .status at all, and previously just hung indefinitely since
            # the session had no timeout configured -- gets the same retry treatment,
            # not just an HTTP error response.
            status = getattr(e, "status", None)
            is_retryable = status in retryable_statuses if status is not None else True
            if not is_retryable or attempt == max_retries - 1:
                raise
            wait_s = retry_backoff_seconds * (2 ** attempt)
            reason = status if status is not None else type(e).__name__
            print(f"  {native_id}: transient {reason}, retrying in {wait_s}s (attempt {attempt + 1}/{max_retries})")
            time.sleep(wait_s)

    return None


# --- Per-station pipeline -----------------------------------------------------------

def process_station(
    station: dict,
    output_dir,
    temporal: tuple[str, str],
    fsspec_config: dict,
    all_gpolygons: dict | None = None,
) -> None:
    """
    Full per-station pipeline: search, dedup, track-filter, skip granules already
    saved, extract each new granule, merge same-pass frames, and append/save to
    `{station['name']}_NISAR_raw.csv` under output_dir. Prints a status line per
    station, matching the notebooks' original inline behavior.

    station: station dict with "name", "lon", "lat", "surface_normal_enu",
        "dem_elevation_m" (see compute_station_terrain)
    output_dir: directory to write/append the station's CSV in (pathlib.Path)
    temporal: (start, end) date strings for the CMR search
    fsspec_config: kwargs for fs.open (cache_type, block_size, ...)
    all_gpolygons: optional dict to collect native-id -> GPolygon into, for an
        overview map across all stations; left untouched if None
    """
    output_path = output_dir / f"{station['name']}_NISAR_raw.csv"

    results = earthaccess.search_data(
        short_name=GCOV_SHORT_NAME,
        point=(station["lon"], station["lat"]),
        temporal=temporal,
    )
    results = keep_latest_versions(results)
    results = filter_tracks_with_min_acquisitions(results)

    if all_gpolygons is not None:
        for result in results:
            native_id = result["meta"]["native-id"]
            all_gpolygons[native_id] = result["umm"]["SpatialExtent"]["HorizontalSpatialDomain"]["Geometry"]["GPolygons"][0]

    if output_path.exists():
        existing_df = pd.read_csv(output_path)
        existing_df["_dt"] = pd.to_datetime(existing_df["datetime_utc"])
    else:
        existing_df = None

    new_results = [r for r in results if not already_in_csv(r, existing_df)]

    if output_path.exists() and not new_results:
        print(f"{station['name']}: {len(results)} granule(s) after filtering, all already saved — skipping")
        return

    print(f"{station['name']}: {len(results)} granule(s) after version dedup and track filtering, {len(new_results)} new")

    station_normal = station["surface_normal_enu"]
    n_new = len(new_results)
    records = []
    for i, result in enumerate(new_results, start=1):
        native_id = result["meta"]["native-id"]
        print(f"  [{i}/{n_new}] opening {native_id} ...", end="", flush=True)
        # A fresh session per granule, not one shared across the whole loop --
        # see get_fsspec_https_session_with_timeout's skip_instance_cache note.
        per_granule_fs = get_fsspec_https_session_with_timeout()
        ok, record = run_with_timeout(
            extract_granule_record,
            (result, station, station_normal, per_granule_fs, fsspec_config),
            PER_GRANULE_TIMEOUT_SECONDS,
        )
        if not ok:
            # The call keeps running in its own daemon thread in the background
            # (can't be forcibly cancelled) -- we just stop waiting on it. Not
            # marked as saved, so it's retried on the next run.
            print(f" timed out after {PER_GRANULE_TIMEOUT_SECONDS}s, skipping (will retry on next run)")
            continue
        print(" ok" if record is not None else " no usable data")
        if record is not None:
            records.append(record)

    # Merged only within this run's new granules — a granule merging with one saved
    # in an earlier run never happens in practice, since same-pass adjacent frames
    # are only seconds apart and therefore always land in the same CMR search/run.
    new_df = pd.DataFrame(records, columns=OUTPUT_COLUMNS).sort_values("datetime_utc")
    n_before = len(new_df)
    new_df = merge_adjacent_acquisitions(new_df)[OUTPUT_COLUMNS]

    write_header = not output_path.exists()
    new_df.to_csv(output_path, mode="a", header=write_header, index=False)
    print(f"  {'Saved' if write_header else 'Appended'} {len(new_df)} row(s) to {output_path} ({n_before} before merging same-pass frames)")
