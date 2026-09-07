# Earsel26_NISAR-main

Pipeline for joining ISMN in-situ soil moisture stations with NISAR L-band NISAR GCOV
backscatter, GSLC image chips, and derived interferometric coherence / closure phase.
Originally built for the Spain (Catalonia) pilot site — adapted from the
[ASFOpenSARlab/NISAR_Cookbook](https://github.com/ASFOpenSARlab/NISAR_Cookbook) — the
same notebooks (`02`–`05`, filenames still say "spain") now also drive the working
Alaska pilot; see the parent [`NISAR4SM`](../README.md) repo's own README for which
pilot is currently wired up and the full notebook execution order across both repos.

**This directory is a subfolder of the parent `NISAR4SM` repo, not a standalone
project** — it started as its own repo (hence its own `.gitignore`/`pixi.toml`
below) and was absorbed in directly. In practice, everything here now runs under
the parent repo's own `nisar4sm-pipeline` conda environment
(`../environment.yml`), not the `pixi`/`isce3` env described below — that setup
is kept for reference / in case this is ever split back out, but isn't what's
actually in use.

## Prerequisites

1. **NASA Earthdata Login** account, with the `Alaska Satellite Facility Data Access`
   EULA accepted. Register at https://urs.earthdata.nasa.gov/, then accept the EULA
   under the `EULAs` tab in your profile. `earthaccess.login()` (used in notebooks
   02/03) will prompt for these credentials the first time and cache them in
   `~/.netrc`.
2. **ISMN account** at https://ismn.earth, used for the manual data download step in
   `01_download_ismn_data.ipynb`.
3. [**pixi**](https://pixi.sh) installed (`curl -fsSL https://pixi.sh/install.sh | bash`)
   — only if running this directory standalone, outside the parent repo (see above;
   not the current setup).

## Setup

Use the parent repo's environment — see [`../README.md`](../README.md) §0. The
`pixi`-based standalone setup below is legacy, kept for reference:

```bash
git clone <this-repo-url>
cd Earsel26_NISAR-main
pixi install
pixi run -e isce3 jupyter lab
```

This creates a `isce3` pixi environment (Python 3.13) with everything the notebooks
and `util`/`ISMN_funcs` packages need: `earthaccess`, `rasterio`, `pyproj`,
`geopandas`, `xarray`/`h5py`/`h5netcdf`, `pygeogrids`, `ipyleaflet`, etc. — see
`pixi.toml`. Notebooks already carry an `isce3` kernel selection, so it should be
preselected when you open them in Jupyter Lab.

## Layout

```
Earsel26_NISAR-main/            # NISAR4SM/Earsel26_NISAR-main/ from the repo root
├── notebooks/            # run these, in order (see below)
│   ├── ismn_data/         # (gitignored) created by 01_download_ismn_data.ipynb
│   └── nisar_data/         # (gitignored) created by notebooks 02-04
├── util/                 # NISAR GCOV/GSLC/coherence-closure helper functions
├── ISMN_funcs/            # vendored ISMN reader (github.com/TUW-GEO/ismn)
├── pixi.toml / pixi.lock  # legacy standalone environment definition (see above)
└── .gitignore
```

`notebooks/` and the `util`/`ISMN_funcs` packages import each other via
`REPO_ROOT = os.path.dirname(os.getcwd())` (i.e. this directory's own root), so
keep that layout intact — don't move the notebooks out of `notebooks/`. (This
used to be `Earsel26_NISAR-main/Earsel26_NISAR-main/` — a doubled "zip extracted
into a same-named folder" artifact from however this was pulled in — collapsed
to the single level shown above; every notebook that referenced the old doubled
path has been updated to match.)

## Notebooks, in order

Written for the Spain pilot originally (hence the `spain`-named output
subfolders below, and the still-Spain-flavored prose in a few of the
docstring/markdown cells); notebooks `02`–`04` were renamed to drop "spain"
from the filename once Alaska became the working pilot, since they now drive
both. Currently pointed at Alaska (`ismn_data/alaska/`,
`nisar_data/alaska_GSLC_coherence_closure/`, etc.) — which AOI a given run
produces depends on the ISMN data downloaded in step 1 and each notebook's
own configuration, not the filenames.

1. **`01_download_ismn_data.ipynb`** — manual step: download an ISMN data bundle
   (zip) for the pilot site from https://ismn.earth, then parse it into
   `ismn_data/<site>/stations.geojson` plus one `sm_<station>_<instrument>_<depth>.csv`
   per station/instrument/depth.
2. **`01_prepare_ismn_data.ipynb`** — reshapes the raw per-station CSVs into
   `ismn_data/<site>/station_data/{station_ID}_metadata_{depth}.csv` and
   `{station_ID}_ts_{depth}.csv`, ready to join against NISAR time series.
3. **`02_get_NISAR_data.ipynb`** — for each station, builds a NISAR L-band GCOV
   (HH/HV) backscatter time series over a 3x3 window, written to
   `nisar_data/<site>/{station}_NISAR_raw.csv`.
4. **`03_get_NISAR_GSLC_windows.ipynb`** — saves a small NISAR GSLC
   (complex, geocoded SLC) pixel window per acquisition per station to
   `nisar_data/<site>_GSLC_windows/`, for interferometric analysis.
5. **`04_calc_coh_phase_closure_NISAR.ipynb`** — computes pairwise coherence and
   triplet closure phase from the GSLC windows, written to
   `nisar_data/<site>_GSLC_coherence_closure/`.

All data these notebooks produce is gitignored — rerun them locally to regenerate it.

## Methodology

### Data sources

- **ISMN** (`ismn_data/`): manually downloaded per-station soil moisture time
  series + metadata (depth, instrument, soil texture, land cover, climate class)
  for the Spain network, via the vendored `ISMN_funcs` reader
  (github.com/TUW-GEO/ismn).
- **NISAR L2 GCOV** (`NISAR_L2_GCOV_PROVISIONAL_V1`): geocoded polarimetric
  covariance (backscatter) product, frequency A (finer-resolution swath) only —
  frequency B (coarser, ~80 m) and S-band are not requested at all, since this
  product line carries no `/science/SSAR` data.
- **NISAR L2 GSLC** (`NISAR_L2_GSLC_PROVISIONAL_V1`, short name unconfirmed
  against a live CMR collection search — see `util/nisar_gslc_windows.py`):
  geocoded complex single-look-complex chips, used only for the interferometric
  coherence/closure-phase step.
- **Copernicus GLO-30 DEM**: streamed on demand per station (`/vsicurl/` +
  `rasterio`, one COG tile per 1°x1° cell, cached in-process) to compute local
  terrain slope/aspect — no DEM data is downloaded in bulk.

Both NISAR product lines are still `PROVISIONAL_V1` (pre-operational/early-mission
data) — this pipeline makes no assumption that calibration or the mask/algorithm
conventions below are final.

### Station-to-pixel matching

Each station's (lon, lat) is transformed into the granule's own projected CRS
(read from the granule's `projection` dataset, via `pyproj.Transformer`), then
matched to the nearest pixel *arithmetically* — `round((x - x_origin) / x_spacing)`
— rather than by reading and searching the full coordinate array (frequencyA's
`xCoordinates`/`yCoordinates` are read only for their origin, spacing, and length).
A granule where the station falls outside the frequency A grid extent entirely is
skipped for that station (no row emitted). GCOV reads a fixed 3x3 window centered
on that pixel, clipped to whatever part falls inside the grid for a station near
an edge. GSLC instead saves a configurable `window_size x window_size` chip
(default 15x15, `DEFAULT_WINDOW_SIZE`); for an even window size the extra pixel is
placed after the center, since a window can't straddle a single pixel evenly.

### Backscatter extraction (GCOV, notebook 02)

Each granule is opened **remotely** via `fsspec` + `h5py` (no full-granule
download — each is ~1.7 GB) with retry-with-backoff on transient `5xx` responses
from the presigned S3/CloudFront URLs (observed mid-stream under load, not just on
open). Per station per granule:

- `NISAR_L_HH_db_mean3x3`/`NISAR_L_HV_db_mean3x3` are computed by averaging the
  window's raw HHHH/HVHV values — **linear power**, per the GCOV spec — and only
  then converting the mean to dB, since dB values cannot be arithmetically
  averaged directly (`10*log10(mean(power))`, not `mean(10*log10(power))`).
- `..._std3x3` is the *sample* standard deviation (ddof=1) computed directly on
  the per-pixel dB values, as a texture/spread descriptor — there's no clean
  power-domain equivalent of a spread statistic the way there is for a mean.
- Pixels are excluded per the GCOV product spec's own `mask` layer convention
  (JPL D-102274 Rev E, §4.3.3, verified against the primary source): `0` =
  "invalid or partially focused" (the averaging ensemble contains at least one
  bad radar sample), `255` = fill (pixel geographically outside the RSLC
  acquisition bounds), `1`–`5` = valid, indicating the dominant sub-swath — only
  `1`–`5` are kept, along with a non-positive-power check per pixel. A
  non-HH/HV-polarized granule (e.g. VV/VH mode) yields an all-NaN row rather than
  being skipped, so orbit/angle metadata for that acquisition isn't lost.
- **Version dedup**: CMR can return multiple CRID (Composite Release Identifier)
  reprocessing versions of the same acquisition; only the highest-CRID one is
  kept, per unique acquisition parsed from `native-id`. The comparison is
  numeric-aware (`_crid_sort_key`, splitting the letter prefix from the digit
  run and comparing the digits as an int) rather than a plain string comparison,
  which would otherwise misorder differently-zero-padded CRIDs (e.g. `"D100"` <
  `"D99"` lexicographically, despite 100 > 99).
- **Track filtering**: a relative orbit/track covering the station fewer than
  `MIN_ACQUISITIONS_PER_TRACK` (default 3) times in the search window is dropped
  entirely, before any granule is opened — parsed straight from `native-id`'s
  track number and start time, so this costs no network I/O.
- **Same-pass merging**: a station near a frame boundary can be covered by two
  granules from the same physical overpass, seconds apart (same track). These
  are collapsed into one row: the mean columns are re-averaged in linear power;
  the std columns are averaged directly (a documented approximation — a
  statistically rigorous pooled-variance combination wasn't judged worth the
  complexity for this rare edge case).
- Re-running a station is incremental: a granule already represented in that
  station's CSV (same track, start time within the same-overpass tolerance) is
  never re-opened.

### Local (terrain-corrected) incidence angle

GCOV's own `incidenceAngle`/LOS fields are referenced to the smooth WGS84
**ellipsoid** normal, not actual terrain — the same ellipsoidal angle can
correspond to very different real illumination geometry depending on whether a
slope faces toward or away from the radar (a slope facing the sensor is
illuminated more directly — smaller true local incidence angle, typically
brighter — than the same slope facing away, at *identical* ellipsoidal incidence
angle). This was motivated by exploratory analysis elsewhere finding a
per-orbit-track effect on raw backscatter that the ellipsoidal incidence angle
(and elevation angle) didn't explain.

- **Terrain slope/aspect**, computed once per station (not per acquisition):
  a 3x3 window of the Copernicus GLO-30 DEM centered on the station is read, and
  slope/aspect are derived via Horn's method (the standard 3x3 Sobel-like
  gradient kernel), with the DEM's geographic (degrees) pixel spacing converted
  to meters (longitude spacing scaled by `cos(latitude)`) so the gradient comes
  out in real slope units rather than degrees/pixel. Aspect is the downhill
  compass bearing; care is taken to negate *both* gradient components when
  converting to a bearing (a common one-component-only bug). This becomes a unit
  surface normal in local East-North-Up (ENU) coordinates. A station within one
  pixel (~30 m) of a DEM tile boundary has its window clamped to stay inside the
  single opened tile rather than fetching the neighbor, which very slightly
  off-centers the terrain sample in that edge case.
- **Per-acquisition local incidence angle**: GCOV's `radarGrid` metadata cube
  provides the ellipsoidal incidence angle and the LOS unit vector's
  East/North components at the (coarser-resolution) nearest radar-grid cell,
  linearly interpolated over height only between the two bracketing height
  levels (reading only those two levels rather than the full ~43 MB/field cube —
  verified to match a full interpolation to within ~0.01°). The LOS vector's
  vertical ("Up") component is reconstructed as `cos(ellipsoidal incidence
  angle)` (exact by the definition of that angle against the ellipsoid normal,
  which *is* "Up" in a local ENU frame), the 3-vector is renormalized (since
  independently height-interpolating 3 fields can leave small norm error), and
  the true local incidence angle is `arccos(LOS_vector · station_surface_normal)`
  — the angle between the actual look direction and the real terrain surface,
  written out as `NISAR_local_incidence_angle_deg` in place of the ellipsoidal
  angle.

### GSLC window extraction (notebook 03)

Saves one NetCDF file per station per acquisition — a `window_size x window_size`
complex chip (`<pol>_real`/`<pol>_imag`, since NetCDF has no native complex type)
for every polarization present, plus a per-pixel geometry cube
(`local_incidence_angle_deg` computed the same way as above, but per-pixel across
the window rather than once at the station's center pixel; terrain slope/aspect
itself is *not* recomputed per pixel, since it varies negligibly across a chip
this size) and the `mask` layer if present. This module's GSLC HDF5 layout
assumptions (complex dtype readable directly by h5py vs. an `r`/`i` compound
type, `mask` presence, `radarGrid` structure, the `projection` EPSG convention)
are inferred by analogy to the verified GCOV pipeline and are **not confirmed
against a real granule** — `inspect_gslc_structure()` (notebook 03, Section 2)
prints one granule's actual layout and should be run once before trusting batch
output.

### Interferometric coherence & closure phase (GSLC, notebook 04)

For each station/track, saved acquisitions are sorted by time and grouped into
sequential pairs `(i, i+1)` and sequential, sliding/overlapping triplets
`(i, i+1, i+2)` — e.g. for acquisitions A/B/C/D: pairs AB/BC/CD and triplets
ABC/BCD. Both are computed over a window centered on the station pixel
(`COHERENCE_WINDOW_HALF` in `util/nisar_coherence_closure.py`, exposed as an
editable variable in notebook 04; `1` = 3x3 default, `2` = 5x5, etc.), using only
pixels valid across *every* acquisition in the pair/triplet (mask intersection),
so the same pixel set contributes to every leg of a triplet's phase loop.

- **Coherence**: the standard sample estimator
  `γ = |Σ(s₁·conj(s₂))| / sqrt(Σ|s₁|² · Σ|s₂|²)`.
- **Closure phase**: `Ξ = angle(z₁₂ · z₂₃ · conj(z₁₃))` where
  `z_ij = Σ(s_i·conj(s_j))` — algebraically equivalent to wrapping
  `angle(z₁₂) + angle(z₂₃) - angle(z₁₃)` to `(-π, π]`, the standard InSAR
  closure-phase quantity used to detect nonlinear (non-geometric) phase
  behavior.
- Every triplet is kept — none dropped for low coherence — and carries
  `min_coherence_<pol>`, the minimum of its three legs' pairwise coherence, so
  low-coherence triplets can be screened out downstream (e.g.
  `min_coherence_HH >= 0.3`) without recomputing anything.
- Because "sequential" means consecutive by *saved acquisition index*, not a
  fixed repeat cycle, a missing acquisition (failed download/processing, or a
  granule that didn't pass track filtering) silently lengthens the real gap
  between two "sequential" acquisitions. `temporal_baseline_days` (pairs) and
  `temporal_baseline_12_days`/`_23_days`/`_13_days` (triplets) record the
  *actual* elapsed time between each pair of dates in a row, so such gaps are
  visible/filterable downstream instead of being silently assumed uniform.

### Known limitations

- **Small-sample coherence bias.** The default 3x3 window is only 9 samples;
  the sample coherence magnitude estimator is known to be positively biased at
  this sample size, worst at low true coherence (Touzi et al., 1999; Lee et al.,
  1994) — apparent coherence over decorrelated (e.g. densely vegetated or
  actively-changing soil-moisture) pixels will read higher than the true value.
  Increasing `COHERENCE_WINDOW_HALF` trades spatial resolution around the
  station for a less-biased estimate; no de-biasing correction is applied.
- **No baseline-consistency enforcement.** The pipeline exposes the actual
  temporal baseline per row (see above) but doesn't itself flag or filter
  irregular gaps — a downstream consumer comparing coherence/closure phase
  across stations or time should check these columns rather than assume a
  uniform repeat cycle.
- **GSLC layout unverified.** See the GSLC section above — `inspect_gslc_structure()`
  should be run against a real granule before trusting batch GSLC/coherence/
  closure-phase output.
- **DEM tile-boundary clamping.** A station within one DEM pixel of a Copernicus
  GLO-30 tile boundary gets a slightly off-center terrain sample rather than a
  DEM read spanning two tiles (see terrain section above) — a small, edge-case-only
  effect on local incidence angle.
- **Provisional product line.** Both GCOV and GSLC are `PROVISIONAL_V1` — mask
  semantics, calibration, and algorithm conventions referenced here reflect the
  product spec as of this writing and may change as the mission moves to
  operational products.
- **No join/modeling logic here.** Joining these per-station CSVs against ISMN
  soil moisture, and any statistical modeling (correlation, regression, land-cover
  stratification) against backscatter/coherence/closure phase, happens in a
  separate downstream repo — none of that logic lives in this one.
