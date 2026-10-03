"""
Per-station NDVI from NASA's HLS (Harmonized Landsat Sentinel-2) collections, used by
notebooks/05_get_HLS_ndvi.ipynb to sanity-check SM-inversion outliers against
vegetation cover -- not a modeling input, just a quick visual QC signal alongside each
SAR acquisition date.

Unlike the NISAR GCOV/GSLC pipelines, no HDF5 remote-metadata-walk is needed: each HLS
band is its own Cloud-Optimized GeoTIFF, so a windowed `rasterio` read over HTTPS pulls
only the needed pixels via a couple of range requests, fast enough that (unlike the
GSLC pipeline) no PER_GRANULE_TIMEOUT_SECONDS wall-clock deadline / background-thread
abandonment is needed here -- reads are synchronous and a stalled one just raises
through the retry loop below like a normal HTTP error, so one shared `fs` session can
safely be reused across every call in a run (no analog of the GSLC pipeline's
concurrent-abandoned-thread bug is possible here, since nothing here is ever abandoned
mid-read).

:::{note} Verified against a live CMR search (2026-08-30, station Coldfoot)
HLS_L30_SHORT_NAME / HLS_S30_SHORT_NAME both returned real results, and the
HLS_BANDS suffixes (B04/B05/B8A/Fmask) match real asset URLs exactly. One thing
that did NOT match the original assumption: earthaccess's result `umm` has no
`CollectionReference.ShortName` for HLS granules (only `EntryTitle`) -- the
collection is identified via the `HLS.L30.`/`HLS.S30.` prefix on `native-id`
instead (see _short_name below).

Still unverified: the Fmask bit layout (FMASK_BIT_*) below, which is the HLS v2.0
User Guide's documented QA bit-packing but hasn't been checked against real decoded
values. Section 2 of the notebook prints the raw Fmask value distribution for one
window so the decoded cloud/shadow/snow flags can be visually sanity-checked before
relying on them to mask every station's data.
:::
"""

import time
from pathlib import Path
from typing import Any

import aiohttp
import numpy as np
import pandas as pd
import pyproj
import rasterio
from rasterio.windows import Window

from util.nisar_station_timeseries import (
    MAX_FETCH_RETRIES,
    RETRY_BACKOFF_SECONDS,
    RETRYABLE_STATUSES,
)


HLS_L30_SHORT_NAME = "HLSL30"  # Landsat 8/9 -- see module note above, unverified
HLS_S30_SHORT_NAME = "HLSS30"  # Sentinel-2 -- see module note above, unverified

# Red/NIR band letters differ between the two collections; NIR uses the "narrow"
# band (S30's B8A, not the native 10m B08) since that's the one HLS harmonizes
# against L30's B05 -- using B08 instead would reintroduce the cross-sensor
# spectral-response difference HLS is meant to remove.
HLS_BANDS = {
    HLS_L30_SHORT_NAME: {"red": "B04", "nir": "B05", "fmask": "Fmask"},
    HLS_S30_SHORT_NAME: {"red": "B04", "nir": "B8A", "fmask": "Fmask"},
}

REFLECTANCE_FILL_VALUE = -9999
REFLECTANCE_SCALE_FACTOR = 0.0001

# HLS v2.0 Fmask bit layout (LP DAAC User Guide) -- see module note above, unverified
# against a real granule.
FMASK_BIT_CIRRUS = 0
FMASK_BIT_CLOUD = 1
FMASK_BIT_ADJACENT_CLOUD_SHADOW = 2
FMASK_BIT_CLOUD_SHADOW = 3
FMASK_BIT_SNOW_ICE = 4
FMASK_BIT_WATER = 5

DEFAULT_WINDOW_FOOTPRINT_M = 80.0  # matches the NISAR GSLC / Sentinel-1 window
DEFAULT_MAX_GAP_DAYS = 5  # widest acceptable gap between a target date and an HLS scene


def _short_name(result) -> str | None:
    """
    HLSL30 or HLSS30, parsed from native-id's "HLS.L30."/"HLS.S30." prefix --
    result["umm"]["CollectionReference"] has no ShortName field for HLS granules
    (only EntryTitle), confirmed against a live search (see module note above).
    """
    native_id = result["meta"]["native-id"]
    if native_id.startswith("HLS.L30."):
        return HLS_L30_SHORT_NAME
    if native_id.startswith("HLS.S30."):
        return HLS_S30_SHORT_NAME
    return None


# --- Structure sanity-check (run once, manually, before batch processing) ---------

def inspect_hls_structure(result, fs) -> None:
    """
    List every asset URL for one HLS granule and open its Red/NIR/Fmask bands to
    print their shape/dtype/CRS, so the assumptions documented at the top of this
    module (band-name suffixes, Fmask bit layout) can be checked against a real
    file before running the full per-station batch below.

    result: one earthaccess search_data() result (HLSL30 or HLSS30)
    fs: an earthaccess fsspec HTTPS session
    """
    short_name = _short_name(result)
    bands = HLS_BANDS.get(short_name)
    urls = result.data_links(access="external")
    print(f"short_name: {short_name}")
    print(f"native-id: {result['meta']['native-id']}")
    print(f"\n{len(urls)} asset URL(s):")
    for u in urls:
        print(f"  {u}")

    if bands is None:
        print(f"\nMISSING: no HLS_BANDS entry for short_name '{short_name}'")
        return

    for role, band in bands.items():
        url = next((u for u in urls if u.endswith(f"{band}.tif")), None)
        if url is None:
            print(f"\nMISSING: no asset URL ending in '{band}.tif' ({role})")
            continue
        with fs.open(url) as fh, rasterio.open(fh) as src:
            print(f"\n{role} ({band}): shape={src.shape} dtype={src.dtypes[0]} crs={src.crs} "
                  f"res={src.res} nodata={src.nodata}")
            if band == bands["fmask"]:
                sample = src.read(1, window=Window(0, 0, min(200, src.width), min(200, src.height)))
                values, counts = np.unique(sample, return_counts=True)
                print(f"  Fmask value distribution in a 200x200 corner sample: "
                      f"{dict(zip(values.tolist(), counts.tolist()))}")


# --- Per-granule extraction ---------------------------------------------------------

def _fetch_band_window(
    url: str,
    lon: float,
    lat: float,
    window_footprint_m: float,
    fs,
    max_retries: int = MAX_FETCH_RETRIES,
    retry_backoff_seconds: int = RETRY_BACKOFF_SECONDS,
    retryable_statuses: set[int] = RETRYABLE_STATUSES,
) -> np.ndarray:
    """
    Open one HLS band COG remotely and read a window of pixels covering
    approximately window_footprint_m x window_footprint_m on the ground, centered
    on the pixel nearest (lon, lat). No download of the full (~120 MB) file --
    rasterio/GDAL issues HTTP range requests for just the needed blocks.

    Retries on transient errors the same way the NISAR pipelines do (the same
    CloudFront/S3 distribution occasionally answers a range request with a
    transient 503 under load) -- but with no PER_GRANULE_TIMEOUT_SECONDS wall-clock
    wrapper, since a COG windowed read is a couple of range requests, not the many
    small sequential HDF5 metadata reads that made the GSLC pipeline need one.

    returns: 2D array (native band dtype, typically int16) of the window's raw
        (unscaled) values
    """
    for attempt in range(max_retries):
        try:
            with fs.open(url) as fh, rasterio.open(fh) as src:
                transformer = pyproj.Transformer.from_crs(4326, src.crs, always_xy=True)
                x0, y0 = transformer.transform(lon, lat)
                row, col = src.index(x0, y0)

                window_size = max(1, round(window_footprint_m / abs(src.res[0])))
                half_before = window_size // 2
                half_after = window_size - half_before - 1
                row0, row1 = max(row - half_before, 0), min(row + half_after + 1, src.height)
                col0, col1 = max(col - half_before, 0), min(col + half_after + 1, src.width)

                return src.read(1, window=Window(col0, row0, col1 - col0, row1 - row0))
        except aiohttp.ClientError as e:
            status = getattr(e, "status", None)
            is_retryable = status in retryable_statuses if status is not None else True
            if not is_retryable or attempt == max_retries - 1:
                raise
            wait_s = retry_backoff_seconds * (2 ** attempt)
            reason = status if status is not None else type(e).__name__
            print(f"    transient {reason} reading {url}, retrying in {wait_s}s "
                  f"(attempt {attempt + 1}/{max_retries})")
            time.sleep(wait_s)
    raise RuntimeError(f"exhausted retries reading {url}")  # pragma: no cover


def _valid_fmask(fmask_window: np.ndarray) -> np.ndarray:
    """
    True where a pixel is NOT flagged cloud, cloud-adjacent, cloud-shadow, or
    snow/ice (see FMASK_BIT_* above) -- water is deliberately not excluded here,
    since a station's window shouldn't normally contain open water and leaving it
    in makes an unexpectedly-wet window visible instead of silently dropped.
    """
    cloud = (fmask_window >> FMASK_BIT_CLOUD) & 1
    adjacent = (fmask_window >> FMASK_BIT_ADJACENT_CLOUD_SHADOW) & 1
    shadow = (fmask_window >> FMASK_BIT_CLOUD_SHADOW) & 1
    snow = (fmask_window >> FMASK_BIT_SNOW_ICE) & 1
    return (cloud == 0) & (adjacent == 0) & (shadow == 0) & (snow == 0)


def compute_window_ndvi(
    result,
    lon: float,
    lat: float,
    fs,
    window_footprint_m: float = DEFAULT_WINDOW_FOOTPRINT_M,
) -> dict[str, Any]:
    """
    NDVI over a station window for one HLS granule: mean/std of (NIR-Red)/(NIR+Red)
    across whichever window pixels aren't flagged cloud/shadow/snow by Fmask.

    result: one earthaccess search_data() result (HLSL30 or HLSS30)
    lon, lat: station coordinates
    fs: an earthaccess fsspec HTTPS session
    window_footprint_m: target window size in meters on the ground
    returns: dict with ndvi_mean, ndvi_std, n_valid_px, n_total_px (NaN/0 if the
        granule has no usable band for this station, e.g. missing asset)
    """
    short_name = _short_name(result)
    bands = HLS_BANDS.get(short_name)
    empty = {"ndvi_mean": np.nan, "ndvi_std": np.nan, "n_valid_px": 0, "n_total_px": 0}
    if bands is None:
        return empty

    urls = result.data_links(access="external")
    band_urls = {
        role: next((u for u in urls if u.endswith(f"{band}.tif")), None)
        for role, band in bands.items()
    }
    if any(u is None for u in band_urls.values()):
        return empty

    red = _fetch_band_window(band_urls["red"], lon, lat, window_footprint_m, fs).astype(np.float32)
    nir = _fetch_band_window(band_urls["nir"], lon, lat, window_footprint_m, fs).astype(np.float32)
    fmask = _fetch_band_window(band_urls["fmask"], lon, lat, window_footprint_m, fs)

    valid = (
        _valid_fmask(fmask)
        & (red != REFLECTANCE_FILL_VALUE) & (nir != REFLECTANCE_FILL_VALUE)
    )
    n_total_px = valid.size
    n_valid_px = int(valid.sum())
    if n_valid_px == 0:
        return {"ndvi_mean": np.nan, "ndvi_std": np.nan, "n_valid_px": 0, "n_total_px": n_total_px}

    red_r = red[valid] * REFLECTANCE_SCALE_FACTOR
    nir_r = nir[valid] * REFLECTANCE_SCALE_FACTOR
    denom = nir_r + red_r
    ndvi = np.where(denom != 0, (nir_r - red_r) / denom, np.nan)

    return {
        "ndvi_mean": float(np.nanmean(ndvi)),
        "ndvi_std": float(np.nanstd(ndvi, ddof=1)) if n_valid_px > 1 else np.nan,
        "n_valid_px": n_valid_px,
        "n_total_px": n_total_px,
    }


# --- Per-station pipeline -----------------------------------------------------------

def _granule_date(result) -> pd.Timestamp:
    # HLS's BeginningDateTime is an ISO8601 string with an explicit "Z" (UTC) suffix,
    # which pandas parses as tz-aware -- but the NISAR side's acquisition_datetime_utc
    # (zeroDopplerStartTime, no tz marker) parses tz-naive, and process_station_ndvi
    # subtracts one from the other to find the nearest match. Both already represent
    # UTC wall-clock time, so stripping tz here (rather than adding one to the NISAR
    # side) keeps every date in this module consistently tz-naive.
    dt = pd.to_datetime(result["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"])
    return dt.tz_localize(None) if dt.tzinfo is not None else dt


def search_hls_for_station(station: dict, temporal: tuple[str, str]) -> list:
    """
    Search both HLSL30 and HLSS30 for one station's point, combined and sorted by
    acquisition date. No version dedup / track filtering needed here (unlike the
    NISAR pipelines) -- HLS has no relative-orbit concept, and CMR doesn't return
    duplicate CRID versions for it the way it does for NISAR GCOV/GSLC.

    station: station dict with "name", "lon", "lat"
    temporal: (start, end) date strings for the CMR search
    returns: combined results, sorted by acquisition date
    """
    import earthaccess

    results = []
    for short_name in (HLS_L30_SHORT_NAME, HLS_S30_SHORT_NAME):
        results.extend(earthaccess.search_data(
            short_name=short_name,
            point=(station["lon"], station["lat"]),
            temporal=temporal,
        ))
    return sorted(results, key=_granule_date)


def process_station_ndvi(
    station: dict,
    target_dates: list,
    search_results: list,
    output_dir: Path,
    fs,
    window_footprint_m: float = DEFAULT_WINDOW_FOOTPRINT_M,
    max_gap_days: int = DEFAULT_MAX_GAP_DAYS,
) -> Path:
    """
    For one station, match each of target_dates to the nearest available HLS scene
    (from search_results, already searched for this station -- see
    search_hls_for_station) within max_gap_days, compute window NDVI for each
    match, and save one row per target date to
    `{station['name']}_HLS_ndvi.csv` under output_dir.

    A target date with no HLS scene within max_gap_days still gets a row (NaN
    ndvi_mean, gap_days beyond the threshold) rather than being silently dropped,
    so a station with no coverage is visible instead of just missing from the
    output entirely.

    station: station dict with "name", "lon", "lat"
    target_dates: SAR acquisition dates to match NDVI against (datetime-like)
    search_results: combined HLSL30+HLSS30 results for this station, sorted by
        date (see search_hls_for_station)
    output_dir: directory to write the CSV into
    fs: an earthaccess fsspec HTTPS session
    window_footprint_m: see compute_window_ndvi
    max_gap_days: a match beyond this many days from its target date is still
        recorded (for visibility) but flagged via gap_days rather than trusted
    returns: path to the saved CSV
    """
    granule_dates = [_granule_date(r) for r in search_results]

    rows = []
    for target in sorted(pd.to_datetime(d) for d in target_dates):
        if not search_results:
            rows.append({"target_date": target, "hls_date": pd.NaT, "sensor": None,
                         "gap_days": np.nan, "ndvi_mean": np.nan, "ndvi_std": np.nan,
                         "n_valid_px": 0, "n_total_px": 0})
            continue

        gaps = [abs((gd - target).total_seconds()) for gd in granule_dates]
        i = int(np.argmin(gaps))
        result = search_results[i]
        hls_date = granule_dates[i]
        gap_days = abs((hls_date - target).total_seconds()) / 86400.0
        short_name = _short_name(result)

        if gap_days > max_gap_days:
            print(f"  {station['name']} {target.date()}: nearest HLS scene is "
                  f"{gap_days:.1f}d away (> {max_gap_days}d), recording anyway")

        ndvi = compute_window_ndvi(result, station["lon"], station["lat"], fs, window_footprint_m)
        rows.append({
            "target_date": target, "hls_date": hls_date, "sensor": short_name,
            "gap_days": gap_days, **ndvi,
        })

    df = pd.DataFrame(rows)
    for col in ("target_date", "hls_date"):
        df[col] = pd.to_datetime(df[col]).dt.strftime("%Y-%m-%dT%H:%M:%S")

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{station['name']}_HLS_ndvi.csv"
    df.to_csv(out_path, index=False)
    return out_path
