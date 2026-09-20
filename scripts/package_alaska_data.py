"""Assemble refined_data_alaska/: a single folder with everything the 4 main
Alaska notebooks (NISAR_SM_alaska.ipynb, S1_SM_alaska.ipynb, Checks_alaska.ipynb,
Checks_alaska_S1.ipynb) and Station_overview_alaska.ipynb need, for handing
off to colleagues.

Copies (never moves) from the working locations those notebooks actually read
from -- rerunning this script is always safe and picks up any data changes.
See README.md section 1.3 for what's in each subfolder and why.

Usage: python scripts/package_alaska_data.py
"""
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "refined_data_alaska"

# (destination subfolder, source directory, optional filename filter)
SOURCES = [
    ("ismn", REPO_ROOT / "Earsel26_NISAR-main" / "notebooks" / "ismn_data" / "alaska", None),
    ("nisar_coherence_closure",
     REPO_ROOT / "Earsel26_NISAR-main" / "notebooks" / "nisar_data" / "alaska_GSLC_coherence_closure", None),
    ("nisar_ndvi",
     REPO_ROOT / "Earsel26_NISAR-main" / "notebooks" / "nisar_data" / "alaska_HLS_ndvi", None),
    ("s1_coherence_closure", REPO_ROOT / "S1_coherence_closure_alaska", None),
    ("soilgrids", REPO_ROOT / "soilgrids_data" / "alaska", None),
    # SME2: only the small cache + manifest, never the ~26 GB of raw .h5 granules.
    ("sme2", REPO_ROOT / "SME2_data_alaska", {"station_manifest.json", "sme2_cache.csv"}),
]

# (destination path, source file) -- individual files rather than whole directories.
FILES = [
    ("nisar_sm_estimates.csv",
     REPO_ROOT / "Earsel26_NISAR-main" / "notebooks" / "nisar_data" / "alaska_NISAR_SM_estimates.csv"),
    ("s1_sm_estimates.csv", REPO_ROOT / "S1_SM_estimates_alaska.csv"),
]


def copy_dir(src: Path, dst: Path, name_filter):
    if not src.exists():
        print(f"  [skip] {src} does not exist")
        return 0, 0
    n_files, n_bytes = 0, 0
    for item in src.rglob("*"):
        if not item.is_file():
            continue
        if name_filter is not None and item.name not in name_filter:
            continue
        rel = item.relative_to(src)
        out_path = dst / rel
        out_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(item, out_path)
        n_files += 1
        n_bytes += item.stat().st_size
    return n_files, n_bytes


def main():
    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)
    OUT_DIR.mkdir(parents=True)

    total_files, total_bytes = 0, 0
    for name, src, name_filter in SOURCES:
        dst = OUT_DIR / name
        n_files, n_bytes = copy_dir(src, dst, name_filter)
        total_files += n_files
        total_bytes += n_bytes
        print(f"  {name}: {n_files} file(s), {n_bytes / 1e6:.1f} MB  <- {src}")

    for name, src in FILES:
        if not src.exists():
            print(f"  [skip] {src} does not exist")
            continue
        dst = OUT_DIR / name
        shutil.copy2(src, dst)
        n_bytes = src.stat().st_size
        total_files += 1
        total_bytes += n_bytes
        print(f"  {name}: {n_bytes / 1e6:.1f} MB  <- {src}")

    print(f"\n{OUT_DIR}: {total_files} file(s) total, {total_bytes / 1e6:.1f} MB")
    print("Zip this folder to hand off to colleagues.")


if __name__ == "__main__":
    main()
