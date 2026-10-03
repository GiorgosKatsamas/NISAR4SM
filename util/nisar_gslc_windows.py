"""
Shared pipeline for saving per-acquisition NISAR GSLC image chips (complex SLC
values + geometry metadata) around each ISMN station, used by
01b_get_nisar_gslc_windows.ipynb.

Unlike util/nisar_station_timeseries.py (one scalar row per acquisition, appended
to a per-station CSV), this module writes one NetCDF file per acquisition per
station, containing a window of pixels around the station rather than a single
value. Product-version dedup, track-count filtering, and DEM-derived terrain are
reused directly from nisar_station_timeseries.py, since GCOV and GSLC share the
same native-id token structure and the same terrain/local-incidence-angle math.

:::{important} Unverified assumptions about GSLC's internal HDF5 layout
This module was written by analogy to the verified GCOV extraction pipeline and
to notebooks/GSLC_Interferogram.ipynb (which opens GSLC via xarray/h5netcdf, not
raw h5py), since no direct h5py-based GSLC read exists elsewhere in this repo to
copy from. Specifically unverified against a real granule:
- That `/science/LSAR/GSLC/grids/frequencyA/<pol>` datasets are complex-valued
  and readable by h5py as numpy complex64/complex128 arrays directly (NISAR/ISCE3
  products conventionally store complex data as an HDF5 compound type with `r`/`i`
  fields, which this module handles defensively, but this hasn't been checked
  against an actual file).
- That a `mask` dataset exists in frequencyA at all (GCOV has one; GSLC may not —
  handled gracefully by omitting the mask variable if absent).
- That `/science/LSAR/GSLC/metadata/radarGrid` exists with the same
  `incidenceAngle`/`losUnitVectorX`/`losUnitVectorY`/`heightAboveEllipsoid` fields
  as GCOV's radarGrid cube.
- That `projection`'s EPSG code is the dataset's raw scalar value (GCOV's
  convention) rather than only available via an `epsg_code` attribute (as
  GSLC_Interferogram.ipynb reads it through xarray) — both are tried.
- `GSLC_SHORT_NAME` below, inferred from ASF's `NISAR_L2_GCOV_PROVISIONAL_V1`
  naming convention — not confirmed via a live CMR collection search.

Run `inspect_gslc_structure()` (Section 2 of the notebook) against one real
granule before batch-processing all stations, so a wrong assumption here is
caught immediately rather than silently producing empty/wrong output.
:::
"""

import re
import time
from pathlib import Path

import aiohttp
import earthaccess
import h5py
import numpy as np
import pandas as pd
import pyproj
import xarray as xr

from util.nisar_station_timeseries import (
    MAX_FETCH_RETRIES,
    PER_GRANULE_TIMEOUT_SECONDS,
    RETRY_BACKOFF_SECONDS,
    RETRYABLE_STATUSES,
    _SkipGranule,
    _interp_radar_grid_field,
    calibrate_to_reference_angle,  # noqa: F401  (re-exported for notebook convenience)
    filter_tracks_with_min_acquisitions,
    get_fsspec_https_session_with_timeout,
    keep_latest_versions,
    run_with_timeout,
)


GSLC_SHORT_NAME = "NISAR_L2_GSLC_PROVISIONAL_V1"  # see module note above — unverified

GSLC_CRID_RE = re.compile(r"^(NISAR_L2_PR_GSLC(?:_[^_]+){9})_([A-Z]\d+)_")
GSLC_GRANULE_ID_RE = re.compile(
    r"^NISAR_L2_PR_GSLC_\d+_(?P<track>\d+)_[AD]_\d+_\d+_\w+_[A-Z]_(?P<start>\d{8}T\d{6})_"
)

CANDIDATE_POLARIZATIONS = ["HH", "HV", "VH", "VV"]

# Meters, not a fixed pixel count -- matches the HyP3 Sentinel-1 burst product's
# ground pixel size (the "20x4" looks default, 80m; see
# hyp3_isce2_burst_stack_for_ts_analysis.ipynb) so the two sensors' windows cover
# comparable ground footprints. GSLC's native grid is anisotropic (~10m range x
# ~5m azimuth, confirmed against a real granule) and every granule reports its own
# xCoordinateSpacing/yCoordinateSpacing, so extract_and_save_gslc_window derives
# separate x/y pixel counts per granule from this target rather than assuming a
# fixed spacing -- update this single constant if the HyP3 stack ever gets
# resubmitted at a different `looks` value (40m/"10x2" or 20m/"5x1").
DEFAULT_WINDOW_FOOTPRINT_M = 80.0

FREQUENCY_A_GROUP = "/science/LSAR/GSLC/grids/frequencyA"
RADAR_GRID_GROUP = "/science/LSAR/GSLC/metadata/radarGrid"
IDENTIFICATION_GROUP = "/science/LSAR/identification"


# --- Structure sanity-check (run once, manually, before batch processing) ---------

def inspect_gslc_structure(result, fs, fsspec_config: dict) -> None:
    """
    Open one GSLC granule remotely and print its frequencyA/radarGrid/
    identification layout, to check this module's layout assumptions
    (see module docstring) against a real file.

    Args:
        result: One earthaccess search_data() result.
        fs: An earthaccess fsspec HTTPS session.
        fsspec_config (dict): kwargs for fs.open (cache_type, block_size, ...).
    """
    url = result.data_links(access="external")[0]
    with fs.open(url, **fsspec_config) as fh, h5py.File(fh, "r") as f:
        print(f"native-id: {result['meta']['native-id']}")
        print(f"\nTop-level groups: {list(f['/science/LSAR'].keys())}")

        if FREQUENCY_A_GROUP not in f:
            print(f"\nMISSING: {FREQUENCY_A_GROUP}")
        else:
            group = f[FREQUENCY_A_GROUP]
            print(f"\n{FREQUENCY_A_GROUP} contents:")
            for k in group.keys():
                ds = group[k]
                if isinstance(ds, h5py.Dataset):
                    print(f"  {k}: shape={ds.shape} dtype={ds.dtype}")
                else:
                    print(f"  {k}/ (group)")
            if "projection" in group:
                proj_ds = group["projection"]
                print(f"  projection raw value: {proj_ds[()]}, attrs: {dict(proj_ds.attrs)}")

        if RADAR_GRID_GROUP not in f:
            print(f"\nMISSING: {RADAR_GRID_GROUP}")
        else:
            rg = f[RADAR_GRID_GROUP]
            print(f"\n{RADAR_GRID_GROUP} contents:")
            for k in rg.keys():
                ds = rg[k]
                if isinstance(ds, h5py.Dataset):
                    print(f"  {k}: shape={ds.shape} dtype={ds.dtype}")

        if IDENTIFICATION_GROUP in f:
            ident = f[IDENTIFICATION_GROUP]
            print(f"\n{IDENTIFICATION_GROUP} contents: {list(ident.keys())}")


# --- Per-granule extraction ---------------------------------------------------------

def _read_epsg(group) -> int:
    """
    EPSG code from a grids/frequencyA `projection` dataset.

    Args:
        group: An open h5py frequencyA group.

    Returns:
        int: EPSG code.

    Notes:
        Tried both as the raw scalar value (GCOV's convention) and via an
        `epsg_code` attribute (how GSLC_Interferogram.ipynb reads it through
        xarray), since it's unconfirmed which applies when read via h5py.
    """
    proj_ds = group["projection"]
    if "epsg_code" in proj_ds.attrs:
        return int(proj_ds.attrs["epsg_code"])
    return int(proj_ds[()])


def _read_complex_window(dataset, row0: int, row1: int, col0: int, col1: int) -> np.ndarray:
    """
    Read a window from an SLC polarization dataset as a numpy complex array.

    Args:
        dataset: An open h5py dataset (a polarization under frequencyA).
        row0 (int): Window row start.
        row1 (int): Window row end.
        col0 (int): Window column start.
        col1 (int): Window column end.

    Returns:
        np.ndarray: Complex-valued window.

    Notes:
        Handles both a native complex dtype (h5py's usual mapping for HDF5
        files that store complex data as an `r`/`i` compound type) and,
        defensively, an explicit `r`/`i`-named compound dtype read out as a
        structured array.
    """
    raw = dataset[row0:row1, col0:col1]
    if np.iscomplexobj(raw):
        return raw
    if raw.dtype.names and {"r", "i"}.issubset(raw.dtype.names):
        return raw["r"].astype(np.float32) + 1j * raw["i"].astype(np.float32)
    return raw.astype(np.complex64)


def extract_and_save_gslc_window(
    result,
    station: dict,
    station_normal: np.ndarray,
    fs,
    fsspec_config: dict,
    output_dir: Path,
    window_footprint_m: float = DEFAULT_WINDOW_FOOTPRINT_M,
    max_retries: int = MAX_FETCH_RETRIES,
    retry_backoff_seconds: int = RETRY_BACKOFF_SECONDS,
    retryable_statuses: set[int] = RETRYABLE_STATUSES,
) -> Path | None:
    """
    Open one GSLC granule remotely (fsspec + h5py, no download), extract a
    window of pixels centered on the frequency A pixel nearest `station`,
    and save it as a NetCDF file under `output_dir`.

    Args:
        result: One earthaccess search_data() result.
        station (dict): Station dict with "name", "lon", "lat",
            "dem_elevation_m" (see compute_station_terrain in
            nisar_station_timeseries.py).
        station_normal (np.ndarray): station["surface_normal_enu"].
        fs: An earthaccess fsspec HTTPS session.
        fsspec_config (dict): kwargs for fs.open (cache_type, block_size, ...).
        output_dir (Path): Directory to write the NetCDF file into.
        window_footprint_m (float): Target window size in meters on the
            ground (not pixels); see Notes.
        max_retries (int): Retry attempts for a transient fetch error.
        retry_backoff_seconds (int): Base backoff between retries (doubles
            each attempt).
        retryable_statuses (set[int]): HTTP statuses treated as transient.

    Returns:
        Path | None: Path to the saved file, or None if the granule wasn't
            usable for this station (missing frequencyA group, no recognized
            polarizations, or station outside the grid extent).

    Notes:
        The saved file contains, for every polarization present (HH/HV/VH/VV,
        whichever apply): `<pol>_real`/`<pol>_imag` (NetCDF has no native
        complex type). Also per-pixel geometry cubes covering the whole
        window -- `incidence_angle_deg` (ellipsoidal),
        `los_unit_vector_x`/`los_unit_vector_y`, and
        `local_incidence_angle_deg` (DEM-terrain-corrected, using the
        station's single precomputed surface normal) -- plus `mask` if
        present, and station/orbit/acquisition metadata as file attrs.
        `window_footprint_m` is converted to separate x/y pixel counts per
        granule from that granule's own xCoordinateSpacing/
        yCoordinateSpacing (GSLC's grid is anisotropic); the actual achieved
        footprint (not exactly `window_footprint_m`, due to pixel rounding)
        is written to the saved file's attrs.
    """
    native_id = result["meta"]["native-id"]
    url = result.data_links(access="external")[0]

    for attempt in range(max_retries):
        try:
            with fs.open(url, **fsspec_config) as fh, h5py.File(fh, "r") as f:
                if FREQUENCY_A_GROUP not in f:
                    raise _SkipGranule
                group = f[FREQUENCY_A_GROUP]

                available_pols = [p for p in CANDIDATE_POLARIZATIONS if p in group]
                if not available_pols:
                    raise _SkipGranule

                ident = f[IDENTIFICATION_GROUP]
                acq_time = ident["zeroDopplerStartTime"][()].decode()
                orb_dir = ident["orbitPassDirection"][()].decode()
                track_number = int(ident["trackNumber"][()])

                epsg = _read_epsg(group)
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

                # Pixel counts derived from this granule's own reported spacing,
                # not a fixed pixel count -- GSLC's grid is anisotropic (~10m
                # range x ~5m azimuth), so window_footprint_m/x_spacing and
                # window_footprint_m/y_spacing differ. x and y are each rounded
                # and centered independently; a window can't be evenly centered
                # on a single pixel with an even side length either way.
                window_size_x = max(1, round(window_footprint_m / abs(x_spacing)))
                window_size_y = max(1, round(window_footprint_m / abs(y_spacing)))
                half_before_x = window_size_x // 2
                half_after_x = window_size_x - half_before_x - 1
                half_before_y = window_size_y // 2
                half_after_y = window_size_y - half_before_y - 1
                row0, row1 = max(yi - half_before_y, 0), min(yi + half_after_y + 1, ny)
                col0, col1 = max(xi - half_before_x, 0), min(xi + half_after_x + 1, nx)

                win_x = group["xCoordinates"][col0:col1]
                win_y = group["yCoordinates"][row0:row1]

                pol_arrays = {
                    pol: _read_complex_window(group[pol], row0, row1, col0, col1)
                    for pol in available_pols
                }

                mask_window = group["mask"][row0:row1, col0:col1] if "mask" in group else None

                # Per-pixel local incidence angle across the window: reuses
                # _interp_radar_grid_field (nearest (y, x) radarGrid cell, linear
                # height interpolation) per pixel, memoized per unique radarGrid
                # cell since the window (a few hundred m at most) typically maps
                # to only one or two cells on radarGrid's much coarser spacing.
                rg = f[RADAR_GRID_GROUP]
                rg_h = rg["heightAboveEllipsoid"][:]
                rg_x0, rg_x1 = rg["xCoordinates"][0:2]
                rg_y0, rg_y1 = rg["yCoordinates"][0:2]
                nx_rg = rg["xCoordinates"].shape[0]
                ny_rg = rg["yCoordinates"].shape[0]
                elevation_m = station["dem_elevation_m"]

                geom_cache: dict[tuple[int, int], tuple[float, float, float]] = {}

                def _geom_at(y: float, x: float) -> tuple[float, float, float]:
                    xi_rg = min(max(int(round((x - rg_x0) / (rg_x1 - rg_x0))), 0), nx_rg - 1)
                    yi_rg = min(max(int(round((y - rg_y0) / (rg_y1 - rg_y0))), 0), ny_rg - 1)
                    key = (yi_rg, xi_rg)
                    if key not in geom_cache:
                        ellip = _interp_radar_grid_field(rg, "incidenceAngle", rg_h, yi_rg, xi_rg, elevation_m)
                        lx = _interp_radar_grid_field(rg, "losUnitVectorX", rg_h, yi_rg, xi_rg, elevation_m)
                        ly = _interp_radar_grid_field(rg, "losUnitVectorY", rg_h, yi_rg, xi_rg, elevation_m)
                        geom_cache[key] = (ellip, lx, ly)
                    return geom_cache[key]

                ny_win, nx_win = len(win_y), len(win_x)
                incidence_cube = np.full((ny_win, nx_win), np.nan, dtype=np.float32)
                los_x_cube = np.full((ny_win, nx_win), np.nan, dtype=np.float32)
                los_y_cube = np.full((ny_win, nx_win), np.nan, dtype=np.float32)
                local_incidence_cube = np.full((ny_win, nx_win), np.nan, dtype=np.float32)

                for r, y in enumerate(win_y):
                    for c, x in enumerate(win_x):
                        ellip, lx, ly = _geom_at(y, x)
                        incidence_cube[r, c] = ellip
                        los_x_cube[r, c] = lx
                        los_y_cube[r, c] = ly
                        if np.isfinite(ellip) and np.isfinite(lx) and np.isfinite(ly):
                            lz = np.cos(np.radians(ellip))
                            vec = np.array([lx, ly, lz])
                            vec = vec / np.linalg.norm(vec)
                            local_incidence_cube[r, c] = np.degrees(
                                np.arccos(np.clip(np.dot(vec, station_normal), -1.0, 1.0))
                            )

                data_vars = {}
                for pol, arr in pol_arrays.items():
                    data_vars[f"{pol}_real"] = (("y", "x"), arr.real.astype(np.float32))
                    data_vars[f"{pol}_imag"] = (("y", "x"), arr.imag.astype(np.float32))
                if mask_window is not None:
                    data_vars["mask"] = (("y", "x"), mask_window)
                data_vars["incidence_angle_deg"] = (("y", "x"), incidence_cube)
                data_vars["los_unit_vector_x"] = (("y", "x"), los_x_cube)
                data_vars["los_unit_vector_y"] = (("y", "x"), los_y_cube)
                data_vars["local_incidence_angle_deg"] = (("y", "x"), local_incidence_cube)

                ds = xr.Dataset(
                    data_vars=data_vars,
                    coords={"y": win_y, "x": win_x},
                    attrs={
                        "station_name": station["name"],
                        "station_lon": station["lon"],
                        "station_lat": station["lat"],
                        "epsg": epsg,
                        "acquisition_datetime_utc": acq_time,
                        "orbit_pass_direction": orb_dir,
                        "track_number": track_number,
                        "native_id": native_id,
                        "polarizations": ",".join(available_pols),
                        "window_footprint_m_requested": window_footprint_m,
                        "window_size_x_px": window_size_x,
                        "window_size_y_px": window_size_y,
                        "window_footprint_x_m_actual": window_size_x * abs(x_spacing),
                        "window_footprint_y_m_actual": window_size_y * abs(y_spacing),
                        "center_row_in_grid": yi,
                        "center_col_in_grid": xi,
                    },
                )

                dt_str = pd.to_datetime(acq_time).strftime("%Y%m%dT%H%M%S")
                dir_tag = "ASC" if orb_dir[:1].upper() == "A" else "DESC"
                fname = f"{station['name']}_{dt_str}_track{track_number:03d}_{dir_tag}_{native_id}.nc"
                output_dir.mkdir(parents=True, exist_ok=True)
                out_path = output_dir / fname
                ds.to_netcdf(out_path, engine="h5netcdf")
                return out_path
        except _SkipGranule:
            return None
        except aiohttp.ClientError as e:
            # Widened from just ClientResponseError (an HTTP error status) to the
            # broader ClientError so a stalled/timed-out fetch -- which has no
            # .status at all -- gets the same retry treatment instead of hanging
            # indefinitely (see get_fsspec_https_session_with_timeout in
            # nisar_station_timeseries.py, which is what actually bounds how long a
            # stalled fetch can hang before raising this in the first place).
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

def process_station_gslc(
    station: dict,
    output_dir: Path,
    temporal: tuple[str, str],
    fsspec_config: dict,
    window_footprint_m: float = DEFAULT_WINDOW_FOOTPRINT_M,
    all_gpolygons: dict | None = None,
) -> None:
    """
    Full per-station GSLC pipeline: search, dedup (CRID), track-filter, skip
    granules already saved, and extract+save each new granule as its own
    NetCDF window file under `output_dir`. Prints one status line per station.

    Args:
        station (dict): Station dict with "name", "lon", "lat",
            "surface_normal_enu", "dem_elevation_m" (see
            compute_station_terrain).
        output_dir (Path): Directory to write NetCDF files into.
        temporal (tuple[str, str]): (start, end) date strings for the CMR search.
        fsspec_config (dict): kwargs for fs.open (cache_type, block_size, ...).
        window_footprint_m (float): See extract_and_save_gslc_window.
        all_gpolygons (dict | None): If given, collects native-id -> GPolygon
            for an overview map across all stations.

    Notes:
        Already-saved granules are matched by native-id substring in existing
        filenames -- no network calls beyond the CMR search for already-covered
        stations.
    """
    results = earthaccess.search_data(
        short_name=GSLC_SHORT_NAME,
        point=(station["lon"], station["lat"]),
        temporal=temporal,
    )
    results = keep_latest_versions(results, crid_re=GSLC_CRID_RE)
    results = filter_tracks_with_min_acquisitions(results, granule_id_re=GSLC_GRANULE_ID_RE)

    if all_gpolygons is not None:
        for result in results:
            native_id = result["meta"]["native-id"]
            all_gpolygons[native_id] = result["umm"]["SpatialExtent"]["HorizontalSpatialDomain"]["Geometry"]["GPolygons"][0]

    existing_filenames = (
        {p.name for p in output_dir.glob(f"{station['name']}_*.nc")} if output_dir.exists() else set()
    )

    station_normal = station["surface_normal_enu"]
    new_results = [
        r for r in results
        if not any(r["meta"]["native-id"] in fname for fname in existing_filenames)
    ]
    n_skipped = len(results) - len(new_results)
    n_new = len(new_results)

    # Per-granule wall-clock deadline: opening a remote GSLC file needs many small,
    # sequential HDF5 metadata reads (superblock, group headers, chunk B-trees), each
    # individually fast but sometimes adding up to minutes for one granule -- no
    # single request is ever slow enough to trip fs's own per-request timeout, so
    # only an external deadline around the whole call can bound it. See
    # PER_GRANULE_TIMEOUT_SECONDS / run_with_timeout in nisar_station_timeseries.py --
    # a *fresh* daemon thread per call, not a shared bounded pool (that first version
    # broke in practice: once enough granules timed out to exhaust a small pool, every
    # later granule "timed out" immediately without ever starting, since none could
    # get a free worker).
    n_saved = n_unusable = n_timed_out = 0
    for i, result in enumerate(new_results, start=1):
        native_id = result["meta"]["native-id"]
        print(f"  [{i}/{n_new}] opening {native_id} ...", end="", flush=True)
        # A fresh session per granule, not one shared across the whole loop --
        # see get_fsspec_https_session_with_timeout's skip_instance_cache note.
        per_granule_fs = get_fsspec_https_session_with_timeout()
        ok, saved_path = run_with_timeout(
            extract_and_save_gslc_window,
            (result, station, station_normal, per_granule_fs, fsspec_config, output_dir, window_footprint_m),
            PER_GRANULE_TIMEOUT_SECONDS,
        )
        if not ok:
            # The call keeps running in its own daemon thread in the background
            # (can't be forcibly cancelled) -- we just stop waiting on it. Not
            # marked as saved, so it's retried on the next run.
            print(f" timed out after {PER_GRANULE_TIMEOUT_SECONDS}s, skipping (will retry on next run)")
            n_timed_out += 1
            continue
        if saved_path is not None:
            n_saved += 1
            print(" ok")
        else:
            n_unusable += 1
            print(" no usable data")

    print(
        f"{station['name']}: {len(results)} granule(s) after version dedup and track filtering — "
        f"{n_saved} saved, {n_skipped} already saved, {n_unusable} skipped (no usable data for this station), "
        f"{n_timed_out} timed out"
    )
