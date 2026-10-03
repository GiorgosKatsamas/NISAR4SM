# NISAR4SM

Pipeline for surface soil moisture (SSM) retrieval from NISAR (L-band) and Sentinel-1 (C-band) InSAR data, using InSAR4SM functionalities.

The repo ships a fully worked example over an **Alaska** SNOTEL/USCRN
network AOI (see [§2]) -- every notebook
reads/writes AOI-labeled paths like `ismn_data/alaska/` or
`nisar_data/<AOI>_GSLC_windows/`; re-running the pipeline for a different AOI
means pointing those same notebooks at your own ISMN download and letting
them write under your own AOI's label instead.

## 1. Execution order

| Notebook | Stage | What it does |
|---|---|---|
| [`01a_download_ismn_data.ipynb`](01a_download_ismn_data.ipynb) | acquisition | Manual step: parses a downloaded ISMN data bundle into `ismn_data/<AOI>/stations.geojson` + raw per-station/instrument/depth CSVs |
| [`01b_get_nisar_gslc_windows.ipynb`](01b_get_nisar_gslc_windows.ipynb) | acquisition | Saves a NISAR GSLC pixel-window NetCDF per acquisition per station |
| [`01c_get_hls_ndvi.ipynb`](01c_get_hls_ndvi.ipynb) | acquisition | Pulls HLS (Landsat/Sentinel-2) NDVI at each station, at NISAR/S1 acquisition dates |
| [`01d_get_s1_burst_stations.ipynb`](01d_get_s1_burst_stations.ipynb) | acquisition | Finds and downloads Sentinel-1 burst SLCs matched to each station's NISAR acquisition dates |
| [`01e_get_sme2.ipynb`](01e_get_sme2.ipynb) | acquisition (optional) | Downloads NISAR L3 SME2 reference SM granules per station -- only feeds the SME2 comparison in `04` |
| [`02a_prepare_ismn_data.ipynb`](02a_prepare_ismn_data.ipynb) | pre-processing | Reshapes `01a`'s raw per-station CSVs into `station_data/{station}_ts_{depth}.csv` + metadata, ready to join against NISAR/S1 time series |
| [`02b_calc_nisar_coherence_closure.ipynb`](02b_calc_nisar_coherence_closure.ipynb) | pre-processing | Computes NISAR pairwise coherence + triplet closure phase from `01b`'s saved GSLC windows |
| [`02c_calc_s1_coherence_closure.ipynb`](02c_calc_s1_coherence_closure.ipynb) | pre-processing | Computes the Sentinel-1 equivalent via HyP3 InSAR-pair jobs on `01d`'s downloaded bursts |
| [`03a_estimate_sm_nisar.ipynb`](03a_estimate_sm_nisar.ipynb) | estimation | De Zan SM inversion from NISAR coherence/closure phase, exported for validation in `04` |
| [`03b_estimate_sm_s1.ipynb`](03b_estimate_sm_s1.ipynb) | estimation | Same inversion from Sentinel-1 coherence/closure phase |
| [`04_station_overview.ipynb`](04_station_overview.ipynb) | checks | Per-station dashboard: NDVI/in-situ SM, coherence/closure-phase-vs-SM-difference screening, NISAR+S1 estimate-vs-in-situ timeseries, track-quality bar plots and selection thresholds, pooled scatter summaries |

If you grabbed the example data package (§3), skip straight to
`03a`/`03b`/`04`. Otherwise, `01a` needs an [ISMN](https://ismn.earth)
account (manual download step); `01b`-`01e` need a NASA Earthdata Login with
the *Alaska Satellite Facility Data Access* EULA accepted
(`earthaccess.login()` prompts for these and caches them in `~/.netrc`).

## 2. Example data package (Alaska)

**Not committed to this repo** (108 MB zipped -- over git's own 100 MB
per-file push limit, and not something a git history should carry anyway).
Download `refined_data_alaska.zip` from this repo's
[Releases](../../releases) page instead, and unzip it at the repo root
(it already mirrors the paths in the table below).

Run `python scripts/package_alaska_data.py` to (re)assemble `refined_data_alaska/`
-- a copy of everything `03a`/`03b`/`04` need in one folder, ready to zip and
attach to a release. It's a copy, not a move: the notebooks keep reading from
their original locations below, so this is safe to rerun any time the
underlying data changes.

| Path | Contents | Needed by |
|---|---|---|
| `ismn_data/alaska/` | `stations.geojson` + `station_data/{station}_ts_*.csv`/`_metadata_*.csv` (raw `sm_*.csv` only needed to re-run `02a`) | `03a`, `03b`, `04` |
| `nisar_data/alaska_GSLC_coherence_closure/` | Per-station NISAR coherence + closure-phase CSVs | `03a`, `04` |
| `nisar_data/alaska_HLS_ndvi/` | Per-station HLS NDVI CSVs | `04` |
| `S1_coherence_closure_alaska/` | Per-(station, orbit) Sentinel-1 coherence + closure-phase CSVs | `03b`, `04` |
| `soilgrids_data/alaska/` | SoilGrids sand/clay rasters, cropped to a buffered bbox around the stations (**manual download** -- see `03a` §3 for the exact files/names needed) | `03a`, `03b` |
| `SME2_data_alaska/station_manifest.json` + `sme2_cache.csv` | Per-station nearest-pixel SME2 SM% already extracted from the raw granules | `04`'s optional SME2 comparison |
| `nisar_data/alaska_NISAR_SM_estimates.csv` / `S1_SM_estimates_alaska.csv` | Exported inversion results (station, track, depth, date, sm_pct, is_dry_anchor) from `03a`/`03b` | `04` |

**Not needed for `03a`/`03b`/`04` at all** (raw inputs to the acquisition/
pre-processing chain -- regenerable from Earthdata/ASF):
`S1_burst_data_alaska/` (large raw burst SLCs), `nisar_data/alaska_GSLC_windows/`
(intermediate between `01b`/`02b`).

## 3. Methodology

### Data sources

- **ISMN** (`ismn_data/`): per-station soil moisture time series + metadata
  (depth, instrument, soil texture, land cover, climate class), via the
  vendored `ISMN_funcs` reader ([github.com/TUW-GEO/ismn](https://github.com/TUW-GEO/ismn)).
- **NISAR L2 GSLC** (`NISAR_L2_GSLC_PROVISIONAL_V1`): geocoded complex
  single-look-complex chips, used for the interferometric coherence/closure-phase
  step (see `util/nisar_gslc_windows.py`).
- **Sentinel-1**: burst SLCs matched to each station's NISAR acquisition
  dates, processed through HyP3 InSAR-pair jobs.
- **Copernicus GLO-30 DEM**: streamed on demand per station (`/vsicurl/` +
  `rasterio`) to compute local terrain slope/aspect -- no DEM data is
  downloaded in bulk.
- **HLS NDVI**, **NISAR L3 SME2** (reference product): each pulled per
  station, see the corresponding acquisition notebook.
- **SoilGrids** texture: **manual download**, no acquisition notebook --
  see `03a_estimate_sm_nisar.ipynb` §3 for where to get it and the exact
  filenames it must be saved under.

Both NISAR product lines are `PROVISIONAL_V1` (pre-operational/early-mission
data) -- this pipeline makes no assumption that calibration or mask/algorithm
conventions are final.

### Station-to-pixel matching

Each station's (lon, lat) is transformed into a granule's own projected CRS
(read from the granule's `projection` dataset via `pyproj.Transformer`), then
matched to the nearest pixel *arithmetically* -- `round((x - x_origin) / x_spacing)`
-- rather than by searching the full coordinate array. GSLC saves a
configurable `window_size x window_size` chip (default 15x15,
`DEFAULT_WINDOW_SIZE` in `util/nisar_gslc_windows.py`).

### Local (terrain-corrected) incidence angle

NISAR's own `incidenceAngle`/LOS fields are referenced to the smooth WGS84
**ellipsoid** normal, not actual terrain -- the same ellipsoidal angle can
correspond to very different real illumination geometry depending on whether
a slope faces toward or away from the radar. `util/nisar_station_timeseries.py`'s
`compute_station_terrain` derives slope/aspect from a 3x3 Copernicus GLO-30
window per station (Horn's method), builds a unit surface normal in local
East-North-Up coordinates, and reports the true local incidence angle --
the angle between the actual look direction and the real terrain surface --
in place of the ellipsoidal one.

### Interferometric coherence & closure phase

For each station/track, saved acquisitions are sorted by time; coherence and
closure phase are computed over a window centered on the station pixel
(`COHERENCE_WINDOW_HALF` in `util/nisar_coherence_closure.py`), using only
pixels valid across every acquisition in the pair/triplet (mask
intersection):

- **Coherence**: the standard sample estimator
  `gamma = |sum(s1 * conj(s2))| / sqrt(sum(|s1|^2) * sum(|s2|^2))`.
- **Closure phase**: `Xi = angle(z12 * z23 * conj(z13))` where
  `z_ij = sum(s_i * conj(s_j))` -- the standard InSAR closure-phase quantity
  used to detect nonlinear (non-geometric) phase behavior.
- Every triplet carries `min_coherence_<pol>`, the minimum of its three legs'
  pairwise coherence, so low-coherence triplets can be screened out
  downstream without recomputing anything.
- `temporal_baseline_days` (pairs) / `temporal_baseline_12_days` etc.
  (triplets) record the *actual* elapsed time between acquisitions, since a
  missing acquisition silently lengthens the gap between two
  "sequential" ones.

### De Zan SM inversion (`03a`/`03b`)

Per (station, track), one joint SLSQP minimization over the full n x n
coherence/closure-phase matrix (`INSAR4SM-main`'s `calc_SM_coherences`/
`invert_sm`) -- not just sequential-pair coherence -- with an in-situ-driven
per-station dry-date selection and a dry-anchor equality constraint. Run at 
both 5cm and 20cm depth where in-situ data supports it.

### Known limitations

- **Small-sample coherence bias**: the sample coherence magnitude estimator
  is positively biased at small window sizes, worst at low true coherence.
  Increasing `COHERENCE_WINDOW_HALF` trades spatial resolution around the
  station for a less-biased estimate; no de-biasing correction is applied.
- **No baseline-consistency enforcement**: the pipeline exposes the actual
  temporal baseline per row but doesn't itself flag or filter irregular
  gaps.
- **Provisional product line**: mask semantics, calibration, and algorithm
  conventions reflect the product spec as of this writing and may change as
  the mission moves to operational products.

## References

- Karamvasis, K., & Karathanassi, V. (2023). Soil moisture estimation from
  Sentinel-1 interferometric observations over arid regions. *Computers &
  Geosciences*, 178, 105410. https://doi.org/10.1016/j.cageo.2023.105410
- Fattahi, H., Agram, P., & Simons, M. (2017). A Network-Based Enhanced
  Spectral Diversity Approach for TOPS Time-Series Analysis. *IEEE TGRS*,
  55(2), 777-786.
- Dorigo, W., et al. (2021). The International Soil Moisture Network:
  serving Earth system science for over a decade. *HESS*, 25(11), 5749-5804.
