"""
Pairwise interferometric coherence and closure phase from NISAR GSLC window chips.

Reads the per-acquisition window chips saved by 01b_get_nisar_gslc_windows.ipynb
(util/nisar_gslc_windows.py) and writes two per-station CSVs, used by
02b_calc_nisar_coherence_closure.ipynb.

For each station/track: coherence + phase for every acquisition pair (the full
N*(N-1)/2 matrix, matching 03a_estimate_sm_nisar.ipynb's full-matrix inversion
input) and closure phase for every triplet (C(N,3)), over a window centered on
the station pixel (COHERENCE_WINDOW_HALF, default 3x3). Both HH and HV are
computed; pixels flagged invalid by the saved mask (0/255) are excluded.

Columns:
    coherence_<pol>, phase_<pol>: magnitude and angle of sum(s_i * conj(s_j)).
    closure_phase_<pol>, min_coherence_<pol>: triplet closure angle and its
        minimum pairwise leg coherence.
    temporal_baseline_days / _12_days / _23_days / _13_days: exact elapsed
        time between the acquisitions in that row.

Notes:
    Same-overpass adjacent-frame acquisitions (seconds apart) are collapsed to
    one per track before pairs/triplets are formed -- see
    _drop_same_pass_duplicates.
    The default 3x3 window (9 samples) positively biases the coherence
    estimator, worst at low true coherence (Touzi et al. 1999); increase
    COHERENCE_WINDOW_HALF to trade spatial resolution for a less-biased
    estimate.
"""

import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from util.nisar_station_timeseries import ADJACENT_PASS_GAP_SECONDS


COHERENCE_WINDOW_HALF = 1  # 3x3 window
MASK_INVALID_VALUES = (0, 255)
POLARIZATIONS = ("HH", "HV")

FILENAME_RE = re.compile(
    r"^(?P<station>.+?)_(?P<dt>\d{8}T\d{6})_track(?P<track>\d+)_(?P<dir>ASC|DESC)_(?P<native_id>.+)\.nc$"
)


def _center_index(size: int) -> int:
    return size // 2


def _window_slice(arr: np.ndarray, half_y: int, half_x: int) -> np.ndarray:
    # y and x centered independently, NOT a single shared index -- a saved window
    # isn't necessarily square (extract_and_save_gslc_window sizes window_size_x/y
    # independently from each granule's own anisotropic ~10m range x ~5m azimuth
    # spacing), so using arr.shape[0] for both axes can silently mis-center on x or
    # raise IndexError once x and y differ enough.
    cy = _center_index(arr.shape[0])
    cx = _center_index(arr.shape[1])
    return arr[cy - half_y: cy + half_y + 1, cx - half_x: cx + half_x + 1]


def _effective_half(arrays: list[np.ndarray], half: int) -> tuple[int, int]:
    """
    Largest half-window (half_y, half_x), capped at `half`, that fits inside
    every array in `arrays` without running past an edge.

    Args:
        arrays (list[np.ndarray]): Arrays being combined (2 for a pair, 3 for
            a triplet).
        half (int): Requested half-window size.

    Returns:
        tuple[int, int]: (half_y, half_x), the largest usable half-window <= `half`.

    Notes:
        Different acquisitions of the same station can have different saved
        window sizes (GSLC's per-granule pixel spacing is anisotropic), so
        the smallest half across every array being combined keeps them all
        the same shape. Assumes each array's true center is at
        `_center_index(size)` -- true for a window that wasn't edge-clipped
        near its granule's grid boundary.
    """
    half_y, half_x = half, half
    for arr in arrays:
        cy = _center_index(arr.shape[0])
        cx = _center_index(arr.shape[1])
        half_y = min(half_y, cy, arr.shape[0] - 1 - cy)
        half_x = min(half_x, cx, arr.shape[1] - 1 - cx)
    return half_y, half_x


def _complex_pol(ds: xr.Dataset, pol: str) -> np.ndarray | None:
    real_name, imag_name = f"{pol}_real", f"{pol}_imag"
    if real_name not in ds or imag_name not in ds:
        return None
    return ds[real_name].values + 1j * ds[imag_name].values


def load_station_acquisitions(gslc_dir: Path, station_name: str) -> dict[int, list[dict]]:
    """
    Load every saved GSLC window .nc file for `station_name` in `gslc_dir`,
    grouped by relative orbit/track and sorted by acquisition time within
    each group.

    Args:
        gslc_dir (Path): Directory of
            {station}_{datetime}_track{NNN}_{dir}_{native_id}.nc files (see
            util/nisar_gslc_windows.py).
        station_name (str): Station to load files for (matched against the
            filename prefix).

    Returns:
        dict[int, list[dict]]: {track_number: [acquisition, ...]}, sorted by
            time within each track. Each acquisition dict holds "datetime"
            (pd.Timestamp), "orb_dir", "native_id", "track", "mask" (2D array
            or None), "path", one entry per POLARIZATIONS (2D complex array
            or None), and "local_incidence_angle_center" (float).
    """
    groups: dict[int, list[dict]] = {}
    for path in sorted(gslc_dir.glob(f"{station_name}_*.nc")):
        m = FILENAME_RE.match(path.name)
        if not m or m.group("station") != station_name:
            continue
        with xr.open_dataset(path) as ds:
            mask = ds["mask"].values if "mask" in ds else None
            local_incidence = ds["local_incidence_angle_deg"].values
            # y and x centered independently -- see _window_slice's comment above;
            # local_incidence isn't necessarily square either.
            cy = _center_index(local_incidence.shape[0])
            cx = _center_index(local_incidence.shape[1])
            acq = {
                "datetime": pd.to_datetime(str(ds.attrs["acquisition_datetime_utc"])),
                "orb_dir": str(ds.attrs["orbit_pass_direction"]),
                "native_id": str(ds.attrs["native_id"]),
                "track": int(ds.attrs["track_number"]),
                "mask": mask,
                "local_incidence_angle_center": float(local_incidence[cy, cx]),
                "path": path,
            }
            for pol in POLARIZATIONS:
                acq[pol] = _complex_pol(ds, pol)
        groups.setdefault(acq["track"], []).append(acq)

    for track in groups:
        groups[track].sort(key=lambda a: a["datetime"])
        groups[track] = _drop_same_pass_duplicates(groups[track])
    return groups


def _drop_same_pass_duplicates(acqs: list[dict]) -> list[dict]:
    """
    Collapse same-overpass adjacent-frame acquisitions (gap <=
    ADJACENT_PASS_GAP_SECONDS) down to the first one seen, dropping the rest.

    Args:
        acqs (list[dict]): One track's acquisitions, already sorted by datetime.

    Returns:
        list[dict]: Acquisitions with same-pass duplicates dropped.

    Notes:
        A station near a frame boundary can get a GSLC window from each of
        two adjacent frames of the same physical pass, seconds apart; left
        in, a pair/triplet built from two such near-identical looks has a
        trivially high coherence and an in-situ-SM correlation that's
        structurally biased (found via GraniteCrk tracks 87/92, 0-8s apart).
    """
    if not acqs:
        return acqs
    kept = [acqs[0]]
    for acq in acqs[1:]:
        gap_seconds = (acq["datetime"] - kept[-1]["datetime"]).total_seconds()
        if gap_seconds > ADJACENT_PASS_GAP_SECONDS:
            kept.append(acq)
    return kept


def _valid_mask(masks: list[np.ndarray | None], half_y: int, half_x: int) -> np.ndarray:
    """
    Intersection of validity across the given acquisitions' mask windows; an
    acquisition with no saved mask contributes no restriction.

    Args:
        masks (list[np.ndarray | None]): Per-acquisition mask windows (2 for
            a pair, 3 for a triplet).
        half_y (int): Window half-height.
        half_x (int): Window half-width.

    Returns:
        np.ndarray: Boolean validity mask, shape (2*half_y+1, 2*half_x+1).
    """
    valid = np.ones((2 * half_y + 1, 2 * half_x + 1), dtype=bool)
    for mask in masks:
        if mask is not None:
            valid &= ~np.isin(_window_slice(mask, half_y, half_x), MASK_INVALID_VALUES)
    return valid


def pair_coherence(acq1: dict, acq2: dict, pol: str, half: int = COHERENCE_WINDOW_HALF) -> float:
    """
    Sample coherence |sum(s1 * conj(s2))| / sqrt(sum(|s1|^2) * sum(|s2|^2))
    over a window centered on the station pixel, excluding pixels flagged
    invalid in either acquisition's mask.

    Args:
        acq1 (dict): First acquisition (see load_station_acquisitions).
        acq2 (dict): Second acquisition.
        pol (str): Polarization, e.g. "HH".
        half (int): Window half-size (see _effective_half for the actual
            size used -- smaller than `half` if either saved window doesn't
            reach that far).

    Returns:
        float: Coherence magnitude, or NaN if `pol` wasn't saved for either
            acquisition or no valid pixel remains.
    """
    s1, s2 = acq1.get(pol), acq2.get(pol)
    if s1 is None or s2 is None:
        return np.nan
    half_y, half_x = _effective_half([s1, s2], half)
    valid = _valid_mask([acq1["mask"], acq2["mask"]], half_y, half_x)
    v1 = _window_slice(s1, half_y, half_x)[valid]
    v2 = _window_slice(s2, half_y, half_x)[valid]
    if v1.size == 0:
        return np.nan
    denom = np.sqrt(np.sum(np.abs(v1) ** 2) * np.sum(np.abs(v2) ** 2))
    if denom == 0:
        return np.nan
    return float(np.abs(np.sum(v1 * np.conj(v2))) / denom)


def pair_coherence_and_phase(
    acq1: dict, acq2: dict, pol: str, half: int = COHERENCE_WINDOW_HALF
) -> tuple[float, complex]:
    """
    Like pair_coherence, but also returns the complex interferometric value
    `coh * exp(i * angle(sum(s1 * conj(s2))))`, for De Zan-style
    phase-closure SM inversion (calc_SM_coherences/invert_sm), which needs
    per-pair phase, not just coherence magnitude.

    Args:
        acq1 (dict): First acquisition.
        acq2 (dict): Second acquisition.
        pol (str): Polarization, e.g. "HH".
        half (int): Window half-size.

    Returns:
        tuple[float, complex]: (coherence, coh * exp(i * phase)); (nan,
            nan+0j) under the same conditions as pair_coherence.

    Notes:
        Shares pair_coherence's window/mask handling exactly, so the two
        never disagree on which pixels a pair's coherence was computed over.
    """
    s1, s2 = acq1.get(pol), acq2.get(pol)
    if s1 is None or s2 is None:
        return np.nan, np.nan + 0j
    half_y, half_x = _effective_half([s1, s2], half)
    valid = _valid_mask([acq1["mask"], acq2["mask"]], half_y, half_x)
    v1 = _window_slice(s1, half_y, half_x)[valid]
    v2 = _window_slice(s2, half_y, half_x)[valid]
    if v1.size == 0:
        return np.nan, np.nan + 0j
    denom = np.sqrt(np.sum(np.abs(v1) ** 2) * np.sum(np.abs(v2) ** 2))
    if denom == 0:
        return np.nan, np.nan + 0j
    z = np.sum(v1 * np.conj(v2))
    coh = float(np.abs(z) / denom)
    ph = coh * np.exp(1j * np.angle(z))
    return coh, ph


def triplet_closure_phase(
    acq1: dict, acq2: dict, acq3: dict, pol: str, half: int = COHERENCE_WINDOW_HALF
) -> float:
    """
    Closure phase angle(z12 * z23 * conj(z13)) in radians, where
    z_ij = sum(s_i * conj(s_j)) over the window.

    Args:
        acq1 (dict): First acquisition.
        acq2 (dict): Second acquisition.
        acq3 (dict): Third acquisition.
        pol (str): Polarization, e.g. "HH".
        half (int): Window half-size (see pair_coherence).

    Returns:
        float: Closure phase in (-pi, pi], or NaN under the same conditions
            as pair_coherence. Uses pixels valid across all three
            acquisitions (mask intersection), so the same pixel set
            contributes to every leg of the loop.
    """
    s1, s2, s3 = acq1.get(pol), acq2.get(pol), acq3.get(pol)
    if s1 is None or s2 is None or s3 is None:
        return np.nan
    half_y, half_x = _effective_half([s1, s2, s3], half)
    valid = _valid_mask([acq1["mask"], acq2["mask"], acq3["mask"]], half_y, half_x)
    v1 = _window_slice(s1, half_y, half_x)[valid]
    v2 = _window_slice(s2, half_y, half_x)[valid]
    v3 = _window_slice(s3, half_y, half_x)[valid]
    if v1.size == 0:
        return np.nan
    z12 = np.sum(v1 * np.conj(v2))
    z23 = np.sum(v2 * np.conj(v3))
    z13 = np.sum(v1 * np.conj(v3))
    triple = z12 * z23 * np.conj(z13)
    if triple == 0:
        return np.nan
    return float(np.angle(triple))


def _safe_nanmin(values: list[float]) -> float:
    finite = [v for v in values if np.isfinite(v)]
    return float(np.min(finite)) if finite else np.nan


def _finalize(rows: list[dict], columns: list[str], date_columns: list[str], path: Path) -> Path:
    df = pd.DataFrame(rows, columns=columns)
    for col in date_columns:
        df[col] = pd.to_datetime(df[col]).dt.strftime("%Y-%m-%dT%H:%M:%S")
    df.to_csv(path, index=False)
    return path


def _baseline_days(acq_a: dict, acq_b: dict) -> float:
    """
    Elapsed time in days between two acquisitions.

    Args:
        acq_a (dict): Earlier acquisition.
        acq_b (dict): Later acquisition.

    Returns:
        float: Signed, exact baseline in days (not an assumed nominal repeat
            interval -- a missed acquisition between them makes this larger
            than the nominal cycle).
    """
    return (acq_b["datetime"] - acq_a["datetime"]).total_seconds() / 86400.0


def process_station(
    station_name: str, gslc_dir: Path, output_dir: Path, window_half: int = COHERENCE_WINDOW_HALF
) -> tuple[Path, Path]:
    """
    Build and save the full-pairwise-matrix coherence CSV and full-triplet
    closure-phase CSV for one station.

    Args:
        station_name (str): Station to process (matches {station_name}_*.nc
            in gslc_dir).
        gslc_dir (Path): Directory of saved GSLC window .nc files (see
            util/nisar_gslc_windows.py).
        output_dir (Path): Directory to write the two CSVs into.
        window_half (int): Half-size of the window centered on the station's
            pixel (1 = 3x3, 2 = 5x5, ...); see COHERENCE_WINDOW_HALF.

    Returns:
        tuple[Path, Path]: (coherence_csv_path, closure_phase_csv_path).

    Notes:
        Every pair (N*(N-1)/2) and triplet (C(N,3)) in a track gets a row --
        not a sequential-only subset -- matching
        03a_estimate_sm_nisar.ipynb's full-matrix inversion input.
        `coherence_<pol> * exp(i * phase_<pol>)` reconstructs the complex
        value calc_SM_coherences/invert_sm need directly from this CSV.
        `temporal_baseline_days` (pairs) / `_12_days`/`_23_days`/`_13_days`
        (triplets) records the exact elapsed time between acquisitions, so an
        irregular gap is visible downstream rather than assumed uniform.
        Both CSVs are fully rebuilt from scratch each run (local/in-memory,
        no network cost).
    """
    groups = load_station_acquisitions(gslc_dir, station_name)

    coherence_rows = []
    closure_rows = []

    for track, acqs in sorted(groups.items()):
        n = len(acqs)

        for i, j in combinations(range(n), 2):
            a, b = acqs[i], acqs[j]
            row = {
                "date_1": a["datetime"],
                "date_2": b["datetime"],
                "temporal_baseline_days": _baseline_days(a, b),
                "orbit_dir": a["orb_dir"],
                "orbit_num": track,
                "local_incidence_angle_deg": np.mean(
                    [a["local_incidence_angle_center"], b["local_incidence_angle_center"]]
                ),
            }
            for pol in POLARIZATIONS:
                coh_val, ph_val = pair_coherence_and_phase(a, b, pol, half=window_half)
                row[f"coherence_{pol}"] = coh_val
                row[f"phase_{pol}"] = float(np.angle(ph_val)) if not np.isnan(coh_val) else np.nan
            coherence_rows.append(row)

        for i, j, k in combinations(range(n), 3):
            a, b, c = acqs[i], acqs[j], acqs[k]
            row = {
                "date_1": a["datetime"],
                "date_2": b["datetime"],
                "date_3": c["datetime"],
                "temporal_baseline_12_days": _baseline_days(a, b),
                "temporal_baseline_23_days": _baseline_days(b, c),
                "temporal_baseline_13_days": _baseline_days(a, c),
                "orbit_dir": a["orb_dir"],
                "orbit_num": track,
                "local_incidence_angle_deg": np.mean([
                    a["local_incidence_angle_center"],
                    b["local_incidence_angle_center"],
                    c["local_incidence_angle_center"],
                ]),
            }
            for pol in POLARIZATIONS:
                row[f"closure_phase_{pol}"] = triplet_closure_phase(a, b, c, pol, half=window_half)
                row[f"min_coherence_{pol}"] = _safe_nanmin([
                    pair_coherence(a, b, pol, half=window_half),
                    pair_coherence(b, c, pol, half=window_half),
                    pair_coherence(a, c, pol, half=window_half),
                ])
            closure_rows.append(row)

    output_dir.mkdir(parents=True, exist_ok=True)

    coherence_columns = [
        "date_1", "date_2", "temporal_baseline_days", "orbit_dir", "orbit_num", "local_incidence_angle_deg",
    ]
    for pol in POLARIZATIONS:
        coherence_columns += [f"coherence_{pol}", f"phase_{pol}"]
    closure_columns = (
        ["date_1", "date_2", "date_3",
         "temporal_baseline_12_days", "temporal_baseline_23_days", "temporal_baseline_13_days",
         "orbit_dir", "orbit_num", "local_incidence_angle_deg"]
        + [f"closure_phase_{pol}" for pol in POLARIZATIONS]
        + [f"min_coherence_{pol}" for pol in POLARIZATIONS]
    )

    coherence_path = _finalize(
        coherence_rows, coherence_columns, ["date_1", "date_2"],
        output_dir / f"{station_name}_NISAR_GSLC_coherence.csv",
    )
    closure_path = _finalize(
        closure_rows, closure_columns, ["date_1", "date_2", "date_3"],
        output_dir / f"{station_name}_NISAR_GSLC_closure_phase.csv",
    )
    return coherence_path, closure_path
