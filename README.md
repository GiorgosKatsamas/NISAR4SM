# NISAR4SM

Pipeline for surface soil moisture (SSM) retrieval and closure-phase signature
analysis from NISAR (L-band) and Sentinel-1 (C-band) InSAR data, built on top
of two vendored tools in this repo (trimmed to what's actually used — see
§3's cleanup checklist for the full audit):

- [`INSAR4SM-main/`](INSAR4SM-main/) — validated Sentinel-1 SSM inversion
  pipeline (De Zan coherence-loss model): `coh_funcs.calc_SM_coherences`,
  `inversion.invert_sm`.
- [`Earsel26_NISAR-main/ISMN_funcs/`](Earsel26_NISAR-main/ISMN_funcs/) —
  a vendored copy of the [`ismn`](https://github.com/TUW-GEO/ismn) package
  (ISMN station data reader), fixed up for this repo's pandas version.

(A third, `nisar_pytools/`, used to be vendored here too, for NISAR HDF5 I/O
— but the only thing this project ever actually called from it was
`find_nisar`, one ASF/CMR search function used by `SME2_alaska.ipynb`. That
function is now inlined directly in that notebook instead, and the ~5,000-line
package was dropped. If the Catalonia pilot (§2) resumes and needs GSLC/GUNW
reading, baseline, or local-incidence-angle primitives, `nisar_pytools` would
need to be re-vendored for that — it's not a dependency of anything currently
working.)

This targets the EARSeL 45th Symposium (Athens, 2026) abstract "*Toward
Improved Soil Moisture Retrieval Using NISAR InSAR Closure Phase Signatures*"
(Karamvasis & Karathanassi) — see
[`manuscript/EARSeL-template_abstract-submission_45th-Symp-Athens-2026.pdf`](manuscript/EARSeL-template_abstract-submission_45th-Symp-Athens-2026.pdf).
Two pilot sites: **Alaska** (SNOTEL/USCRN
network, working end-to-end with real NISAR + Sentinel-1 data — see
[§1](#1-alaska-pilot-working-end-to-end)) and **Catalonia** (XMS-CAT network,
planning-stage only, no driver notebook written yet — see
[§2](#2-catalonia-pilot-exploratory)).

## 0. Setup

**Environment.** `mamba env create -f environment.yml` (or `conda`), then
select the `nisar4sm-pipeline` env as your Jupyter kernel. `INSAR4SM-main/`
and `ISMN_funcs/` are vendored directly in this repo — every
notebook puts them on `sys.path` itself, so they don't need a separate
install.

**Data.** None of the raw/large datasets are committed to the repo (see
`.gitignore`) — NISAR/ERA5/ISMN/SME2 data carry their own usage terms, and
several are tens of GB. What you need depends on which pilot and which
notebooks you're running — see the data tables in §1 and §2 below.

## 1. Alaska pilot (working end-to-end)

### 1.1 Main notebooks

These four are the actual deliverables — everything else in this section
feeds them:

| Notebook | What it does |
|---|---|
| [`NISAR_SM_alaska.ipynb`](NISAR_SM_alaska.ipynb) | De Zan SSM inversion from NISAR GSLC coherence/closure phase, HH, SoilGrids texture, 5cm & 20cm, validated against ISMN + SME2 |
| [`S1_SM_alaska.ipynb`](S1_SM_alaska.ipynb) | Same inversion from Sentinel-1 (HyP3 burst) coherence/closure phase, VV, SoilGrids texture, 5cm & 20cm |
| [`Checks_alaska.ipynb`](Checks_alaska.ipynb) | De Zan model sanity check: does NISAR coherence loss / closure phase actually track real ISMN SM differences? Reads the same precomputed coherence/closure CSVs as `NISAR_SM_alaska.ipynb`, no inversion |
| [`Checks_alaska_S1.ipynb`](Checks_alaska_S1.ipynb) | Same sanity check for Sentinel-1 |

Supplementary: [`Station_overview_alaska.ipynb`](Station_overview_alaska.ipynb) — a
per-station dashboard pulling together everything above (NDVI/in-situ SM,
per-track coherence/closure-phase-vs-SM-diff regressions, NISAR+SME2/S1
estimate timeseries, local incidence angle, terrain slope), one page per
station, with axis ranges and track/orbit colors held fixed across every
station and panel so they're directly comparable. Needs
`NISAR_SM_alaska.ipynb`/`S1_SM_alaska.ipynb` to have been run at least once
first (it reads their exported SM-estimate CSVs, not recomputing the
inversion itself).

### 1.2 Auxiliary notebooks (data preparation — run once, in order)

You only need these if you're regenerating the data from scratch. **If you
received the Alaska data package (§1.3), skip straight to the main
notebooks above.**

```
01_download_ismn_data.ipynb  ──▶  01_prepare_ismn_data.ipynb  ──▶  (station_data/, stations.geojson)
        (01b_inspect_ismn_raw.ipynb — optional raw-data QC, not part of the chain)

S1_burst_stations.ipynb  ──▶  S1_coh_phase_closure.ipynb  ──▶  (S1_coherence_closure_alaska/)

02_get_NISAR_data.ipynb  ──▶  03_get_NISAR_GSLC_windows.ipynb  ──▶  04_calc_coh_phase_closure_NISAR.ipynb  ──▶  (alaska_GSLC_coherence_closure/)

05_get_HLS_ndvi.ipynb  ──▶  (alaska_HLS_ndvi/)

SME2_alaska.ipynb  ──▶  (SME2_data_alaska/)   — optional, only feeds the SME2 reference-product comparison
```

(Notebooks `02`–`04` originally had "spain" in their filenames — renamed since
they're parameterized and are what actually produced the Alaska GSLC
coherence/closure data above; see [`Earsel26_NISAR-main/README.md`](Earsel26_NISAR-main/README.md)
for their per-notebook docs.)

Requires Earthdata/ASF credentials for the NISAR and Sentinel-1 downloads,
and an ISMN account for the station data — none of that is needed if you're
using the provided data package instead.

### 1.3 Alaska data package (provided separately, not in git)

*[Karolos: fill in where colleagues get this — cloud link, filesize once
finalized, etc.]*

Run `python scripts/package_alaska_data.py` to assemble `refined_data_alaska/`
— a copy of everything below in one folder, ready to zip and send. It's a
copy, not a move: the notebooks keep reading from their original locations
below, so this is safe to rerun any time the underlying data changes.

| Path | Contents | Needed by | Size |
|---|---|---|---|
| `Earsel26_NISAR-main/notebooks/ismn_data/alaska/` | `stations.geojson` + `station_data/{station}_ts_*.csv`/`_metadata_*.csv` (also the raw `sm_*.csv`, only needed if re-running `01_prepare_ismn_data.ipynb`) | all 4 main notebooks | ~3.6 MB |
| `Earsel26_NISAR-main/notebooks/nisar_data/alaska_GSLC_coherence_closure/` | Per-station NISAR coherence + closure-phase CSVs | `NISAR_SM_alaska.ipynb`, `Checks_alaska.ipynb` | ~0.25 MB |
| `Earsel26_NISAR-main/notebooks/nisar_data/alaska_HLS_ndvi/` | Per-station HLS NDVI CSVs | `Checks_alaska.ipynb`, `Checks_alaska_S1.ipynb` | ~0.03 MB |
| `S1_coherence_closure_alaska/` | Per-(station, orbit) Sentinel-1 coherence + closure-phase CSVs | `S1_SM_alaska.ipynb`, `Checks_alaska_S1.ipynb` | ~0.1 MB |
| `soilgrids_data/alaska/` | SoilGrids sand/clay rasters, cropped to a 0.25°-buffered bbox around the 20 Alaska stations (was 343.5 MB full-resolution; sampled values verified pixel-identical for every station) | `NISAR_SM_alaska.ipynb`, `S1_SM_alaska.ipynb` | 107.7 MB |
| `SME2_data_alaska/station_manifest.json` + `sme2_cache.csv` | Per-station nearest-pixel SME2 SM% already extracted from the raw granules (both notebooks now read this cache first if present, falling back to the raw HDF5s otherwise) — **do not send the raw `.h5` files**, only these two small files | `NISAR_SM_alaska.ipynb` §8b, `S1_SM_alaska.ipynb` §8b (optional comparison only) | ~10 KB |
| `Earsel26_NISAR-main/notebooks/nisar_data/alaska_NISAR_SM_estimates.csv` | Exported `sm_map` from `NISAR_SM_alaska.ipynb`'s §12 (station, track, depth, date, sm_pct, is_dry_anchor) | `Station_overview_alaska.ipynb` | ~0.05 MB |
| `S1_SM_estimates_alaska.csv` | Same, exported from `S1_SM_alaska.ipynb`'s §12 | `Station_overview_alaska.ipynb` | ~0.02 MB |

Everything the 4 main notebooks need, including the two rows above, is now
about 111 MB total — down from ~26.6 GB. Both changes were verified to
produce byte-identical downstream results (inversion output counts and
values unchanged before/after).

**Not needed for the main notebooks at all** (raw inputs to the auxiliary
chain — regenerable from Earthdata/ASF, don't send these): `S1_burst_data_alaska/`
(~24.5 GB raw burst SLCs), `alaska_GSLC_windows/` (~11 MB, intermediate
between notebooks 03/04).

## 2. Catalonia pilot (exploratory)

**Planning-stage only — no notebook exists yet.** The NISAR HDF5 → De Zan
inversion bridge isn't built, and there is no `NISAR_SM.ipynb` (or any other
driver notebook) in this repo for Catalonia — what follows is prep work and
an intended approach, not working code. `ERA5_Land_data/`, `SME2_data/`
(root), and their planned `soilgrids_data/` tiles were downloaded toward this
but aren't read by anything currently in the repo (verified: no notebook
references any of these three paths). If this pilot isn't being picked back
up, those three directories (47 MB + 858 MB, no AOI-cropped soilgrids
tiles were ever pulled) are safe to delete; if it is, keep them — they took
Earthdata/Copernicus credentials and time to fetch.

**Data**, none committed — download and place under `BASE` once a driver
notebook exists:

| Folder | Contents | Source |
|---|---|---|
| `NISAR_data/` | NISAR L2 GSLC `.h5` products | [NASA Earthdata](https://search.earthdata.nasa.gov/) / ASF |
| `ERA5_Land_data/` | ERA5-Land GRIB files (skin temperature, precipitation), one file per month | [Copernicus CDS](https://cds.climate.copernicus.eu/) |
| `soilgrids_data/` | `sand.tif`, `clay_cont.tif` for the AOI | [soilgrids.org](https://soilgrids.org/) |
| `Data_separate_files_header_.../` | In-situ ISMN `.stm` files (XMS-CAT network) | [ISMN](https://ismn.earth/) |
| `SME2_data/` | NISAR L3 SME2 products | Same source as `NISAR_data/` |

**Status / open items:**
- AOI: 0.802–3.138°E, 41.702–42.606°N (Catalonia), 13 Nov 2025 – 14 Jan 2026.
- Bridge strategy decided: reuse `INSAR4SM-main`'s numpy functions
  (`Covar_2_phase_closures`, `calc_SM_coherences`, `invert_sm`, all still
  vendored here) directly on arrays extracted from NISAR HDF5, rather than
  reimplementing in dask/xarray — the science is peer-reviewed, only the I/O
  layer (ISCE2 vs. NISAR HDF5) differs. The I/O side itself isn't decided:
  `nisar_pytools` was the assumed tool for this, but it's no longer vendored
  here (see the intro's note on why) — re-vendor it, or pick something else,
  whenever this pilot actually resumes.
- Still open: land-cover source for bare/cropland/forest stratification
  (Thread A), ISMN station availability inside the Catalonia AOI, and
  whether to rerun once the Feb 2026 NISAR release defines the final study
  area.

## 3. Before uploading — cleanup checklist

- [x] Shrink `soilgrids_data/alaska/` (343.5 MB → 107.7 MB, cropped to the
      station bbox, verified pixel-identical) and `SME2_data_alaska/`
      (26.3 GB → ~10 KB cache, verified identical `sme2_series_by_station`
      output) — done, see §1.3.
- [x] Deleted `NISAR_data_summer/` (~113 GB), `SME2_data/` root (~858 MB),
      `ERA5_Land_data/` (~47 MB) — confirmed nothing in the repo reads any of
      them (see §2: the Catalonia pilot these were downloaded for has no
      driver notebook). `.kilo/` also removed (another tool's local state).
- [x] Collapsed the doubled `Earsel26_NISAR-main/Earsel26_NISAR-main/` nesting
      (a "zip extracted into a same-named folder" artifact) to a single
      `Earsel26_NISAR-main/` — every path constant across all 9 affected
      notebooks (4 main + `01_download_ismn_data`, `01b_inspect_ismn_raw`,
      `S1_coh_phase_closure`, `S1_burst_stations`, `SME2_alaska`) updated and
      each notebook re-run to confirm identical results.
- [x] `refined_data_alaska/` — a single folder with everything the 4 main
      notebooks need, assembled by `scripts/package_alaska_data.py` (copies
      from the working locations, never moves them — rerun any time the data
      changes). See §1.3.
- [x] Found and fixed a real data bug while re-verifying the above: 7
      stations' raw `sm_*.csv` files (RockyPoint, AtigunPass, Coldfoot,
      EagleSummit, ImnaviatCreek, McGrath, NukaGlacier) had their soil
      moisture values comma-decimal formatted (`0,447` instead of `0.447`,
      likely from an Excel round-trip at some point) instead of period. Under
      the pandas version now in the env this silently produced a `str`-dtype
      column instead of `object`, which crashed every notebook that touched
      one of those stations' dry-date selection. Fixed by normalizing the
      raw CSVs (comma → period, values otherwise untouched) and regenerating
      `station_data/` — all 4 main notebooks + `01_prepare_ismn_data`
      re-verified afterward with identical (non-crashing) results.
- [x] `git init` + confirm `.gitignore` covers every data directory now
      present (several were added since the last time it was edited),
      `.venv/`, `.kilo/`, `hyp3-docs-develop/` (a full clone of a public
      repo — link to it instead of vendoring), and the ISMN zip files at
      the repo root. Done in this pass — worth a final review before `git add`.
      (`git init` itself still needs to be run — that's yours to do.)
- [x] Consolidated the vendored-library dead code (audit was already done;
      this pass executed it, with every affected notebook re-run afterward to
      confirm byte-identical results):
      - `INSAR4SM-main/insar4sm/`: trimmed to the 5 modules actually used
        (`coh_funcs.py`, `inversion.py`, `forward_modelling_funcs.py`,
        `phase_closure_funcs.py`, `dielectic_models.py` — verified their own
        cross-imports don't reach anything else first). Dropped `classes.py`,
        `ERA5/`, `Soilgrids/`, `prep_*`, `sorting_funcs.py`, `DS_funcs.py`,
        `calc_SM_index.py`, `calc_int_data.py`, `amp_funcs.py`, `gridding.py`,
        `meteo_funcs.py`, `joblib_progress_bar.py`, `sm_estimation_funcs.py`,
        plus `INSAR4SM_app.ipynb` (the upstream demo notebook, since it
        depended on the modules just removed — its own README already links
        the canonical copy on GitHub) and `.readthedocs.yaml`/`INSAR4SM_env.yml`.
        `NISAR_SM_alaska.ipynb`/`S1_SM_alaska.ipynb` re-run afterward: 158 and
        44 inversion results respectively, identical to before trimming.
      - `nisar_pytools/` (~5,000+ lines): the only thing ever used from it —
        `find_nisar` plus its two validators, one cell in `SME2_alaska.ipynb`
        — was inlined directly into that notebook (~230 lines, self-contained,
        no more vendored-package dependency), then the whole package was
        deleted. Verified the inlined version against the *live* ASF/CMR
        catalog (found 25 SME2 URLs for RockyPoint — no auth needed for
        search), not just a syntax check. Added `asf_search` to
        `environment.yml` (previously only pulled in transitively).
      - `Earsel26_NISAR-main/ISMN_funcs/`: dropped `custom.py` (unused
        extension point) and `cli.py` (broken — imported the wrong package
        name). **Caught a mistake here**: I first also deleted `citations.txt`
        as "unused package metadata" per the original audit, but it's actually
        read at import time by `const.py` — broke `ISMN_Interface` entirely.
        Recovered it from the upstream `ismn` PyPI package (byte-for-byte the
        same file) and re-verified both `01_download_ismn_data.ipynb` and
        `01b_inspect_ismn_raw.ipynb` end-to-end (44 sensors, same station list
        as before) before moving on.
      - `hyp3-docs-develop/` (532 files, zero code coupling): already gone —
        removed directly at some point, same as the other large unused
        directories earlier in this cleanup.
- [x] Renamed `02_get_NISAR_spain_data.ipynb` → `02_get_NISAR_data.ipynb`,
      `03_get_NISAR_spain_GSLC_windows.ipynb` → `03_get_NISAR_GSLC_windows.ipynb`,
      `04_calc_coh_phase_closure_NISAR_spain.ipynb` → `04_calc_coh_phase_closure_NISAR.ipynb`
      — they drive both pilots now. Every reference updated (both READMEs,
      `Checks_alaska.ipynb`, `S1_burst_stations.ipynb`, `S1_coh_phase_closure.ipynb`,
      `05_get_HLS_ndvi.ipynb`, and notebooks 03/04's own cross-links to each
      other) — confirmed these are markdown/doc links only, nothing opens
      another notebook by filename programmatically, so no execution risk.
- [ ] Fill in the Alaska data package distribution link — placeholder still
      needed at the top of §1.3 (removed by accident in an earlier pass).
- [ ] **License — deliberately deferred.** `INSAR4SM-main/LICENSE` is **GPL**
      (Karamvasis, 2022–2025) and both main notebooks import from it directly
      (`calc_SM_coherences`, `invert_sm`), so this project already incorporates
      GPL-licensed code — worth keeping in mind whenever licensing does get
      decided, since it constrains the options. Not blocking this upload: the
      repo isn't going public yet. Karolos is handling this one directly.
- [x] Credential hygiene checked — `S1_burst_stations.ipynb`, `S1_coh_phase_closure.ipynb`,
      `SME2_alaska.ipynb` all use interactive/`.netrc`-based auth
      (`getpass`, `hyp3_sdk`'s `prompt="password"`, `earthaccess.login()`);
      no hardcoded credentials or tokens found anywhere in the repo.
- [x] Absolute local path (`D:\Docs\RS_2023\NISAR4SM\Data_separate_files_header_...zip`,
      the only one anywhere in the repo) hardcoded in `01_download_ismn_data.ipynb`
      and `01b_inspect_ismn_raw.ipynb` — replaced with `PROJECT_ROOT = Path.cwd().parents[1]`,
      matching this repo's existing relative-path convention. Both notebooks
      re-run to confirm identical results.
- [x] Notebook outputs: decided to keep them committed inline (not stripping
      with `nbstripout`) — colleagues can see results without re-running, which
      matches how these notebooks have been used as a record of findings, not
      just code. Diff noise on re-run is an accepted tradeoff.

## References

- Karamvasis, K., & Karathanassi, V. (2023). Soil moisture estimation from
  Sentinel-1 interferometric observations over arid regions. *Computers &
  Geosciences*, 178, 105410. https://doi.org/10.1016/j.cageo.2023.105410
- Fattahi, H., Agram, P., & Simons, M. (2017). A Network-Based Enhanced
  Spectral Diversity Approach for TOPS Time-Series Analysis. *IEEE TGRS*,
  55(2), 777–786.
- Dorigo, W., et al. (2021). The International Soil Moisture Network:
  serving Earth system science for over a decade. *HESS*, 25(11), 5749–5804.
