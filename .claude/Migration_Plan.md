# Migration Plan: Rebinning + ICA Code into HST Paper Repo

**Created:** 2026-04-29
**Author:** Alexandros Pratsos + Claude
**Source:** `Richards/Trevor Code/`
**Destination:** `Richards/HST Paper/`

---

## Overview

Migrate the validated rebinning pipeline (2026AP) and ICA analysis code from
`Trevor Code/` into the `HST Paper/` GitHub repo. The end goal is a single
reproducible repo that can download data, rebin/coadd it, run ICA analysis,
and produce publication figures.

### Guiding Principles

1. **No functional changes** — algorithm logic stays identical; only structure,
   imports, and paths change.
2. **Rename for clarity** — files get clean Python names, but each file keeps a
   header comment recording its original `Trevor Code/` filename.
3. **Phased approach** — each phase is a separate commit with logged changes.
4. **Traceability** — all changes recorded in `Migration_Log.md`.
5. **Portable / collaborator-friendly** — the repo must work on any machine,
   not just AP's. All paths must be relative to the repo root (no hardcoded
   `/Users/...` or `C:\Users\...` paths). Any external data locations (e.g.,
   raw FITS files not stored in the repo) are specified via a single
   configuration variable or CLI argument, documented in README and CLAUDE.md,
   so a collaborator only needs to clone the repo, `pip install -e .`, point
   to their data directory, and run.

---

## Target Repo Structure (after all phases)

```
HST Paper/
├── .claude/
│   ├── Migration_Plan.md        # This file
│   └── Migration_Log.md         # Running log of all changes
├── CLAUDE.md                    # Codebase instructions for Claude
├── pyproject.toml               # Installable package (pip install -e .)
├── requirements.txt             # (updated)
│
├── rebinning/                   # Phase 1 — rebinning pipeline
│   ├── __init__.py
│   ├── coadd.py                 # from HST_coadd2_2026AP.py
│   ├── read_spec_data.py        # from Read_spec_data_2026AP.py
│   ├── cut_edge_pix.py          # from Cut_Edge_Pix_TVM_NoDQ.py
│   ├── lower_res_rebin.py       # from LowerResHSTRebin_TVM.py
│   ├── small_pix_filter.py      # from Small_Pix_Filter_TVM.py
│   ├── spec_cuts_fos.py         # from SpecCuts_FOS.py
│   ├── spec_cuts_hsla.py        # from SpecCuts_HSLA.py
│   └── run_rebin.py             # unified runner (from run_rebin_2026AP_*.py)
│
├── ica/                         # Phase 2 — ICA analysis pipeline
│   ├── __init__.py
│   ├── run_ica.py               # from ICA Scripts/run_ICA_r20_components.py
│   ├── fit_composites.py        # from ICA Scripts/fit_composites.py
│   ├── plot_ica.py              # from ICA Scripts/plot_ICA.py
│   ├── spec_morph.py            # from ICA Scripts/spec_morph.py
│   ├── bal_regions.py           # from ICA Scripts/BAL_regions.py
│   ├── civ_bal_regions.py       # from CIV_BAL_regions.py
│   ├── morphing_edges.py        # merged from MorphingEdges*.py (5 files)
│   ├── manual_fix.py            # from ica_manual_fix_refactored.py
│   ├── manual_fix_config.py     # from manual_fix_config.py
│   ├── quality_classification.json  # generated from ICA_full/ folder scan
│   └── run_all_objects.py       # from Run_ICA_All_Objects.py
│
├── Data/                        # Data files (existing + new)
│   ├── components/              # ICA component files (copied from Trevor Code/components/)
│   │   ├── amy_12603000_10c_180421.comp
│   │   ├── amy_12653000_hew_hsn_7c_190302.comp
│   │   ├── amy_12753000_lowew_10c_181101_v1.comp
│   │   ├── qsomod_z250_ng_continuumSDSSrebin.npy
│   │   ├── wav_12603000.dat
│   │   ├── wav_12653000.dat
│   │   └── wav_12753000.dat
│   ├── HST_CIV_Sulentic2007_HSLA2018_finalprops.csv   # (copied; eventually from pipeline)
│   └── ... (existing data files)
│
├── Figures/                     # (existing)
├── Notebooks/                   # (existing)
└── Plotting Code/               # (existing — will be reconnected later)
```

---

## Phase 0 — Project Scaffolding

**Goal:** Set up the repo infrastructure before copying any science code.

### Steps

0.1. Create `CLAUDE.md` for `HST Paper/` with repo conventions, package
     structure, and key principles.

0.2. Create `pyproject.toml` at `HST Paper/` root defining the installable
     package. Includes:
     - Package name (e.g., `hst_civ_paper`)
     - Python >=3.12
     - Dependencies: numpy, pandas, scipy, matplotlib, astropy, lmfit,
       weightedstats, palettable, richardsplot
     - Editable install: `pip install -e .`

0.3. Create empty `rebinning/__init__.py` and `ica/__init__.py`.

0.4. Create `Data/components/` directory.

0.5. Create `.claude/Migration_Log.md` with header (empty log body).

0.6. Update `.gitignore` if needed (e.g., add `__pycache__/`, `*.pyc`,
     `RebinnedSpec_*/`, `ICA_Plots_*/`).

0.7. **Commit:** "Phase 0: project scaffolding for code migration"

---

## Phase 1 — Rebinning Pipeline Migration

**Goal:** Copy the validated 2026AP rebinning code into `rebinning/`.

### Source → Destination Map

| Source (`Trevor Code/`)           | Destination (`HST Paper/rebinning/`) |
|-----------------------------------|--------------------------------------|
| `HST_coadd2_2026AP.py`           | `coadd.py`                           |
| `Read_spec_data_2026AP.py`       | `read_spec_data.py`                  |
| `Cut_Edge_Pix_TVM_NoDQ.py`       | `cut_edge_pix.py`                    |
| `LowerResHSTRebin_TVM.py`        | `lower_res_rebin.py`                 |
| `Small_Pix_Filter_TVM.py`        | `small_pix_filter.py`                |
| `SpecCuts_FOS.py`                 | `spec_cuts_fos.py`                   |
| `SpecCuts_HSLA.py`                | `spec_cuts_hsla.py`                  |
| `run_rebin_2026AP_FOS_STIS.py` + `run_rebin_2026AP_HSLA.py` | `run_rebin.py` (merged) |

### Steps

1.1. Copy each source file to its destination with:
     - A header comment recording the original filename and date of copy
     - No changes to algorithm logic

1.2. Update imports within the copied files:
     - Replace `sys.path.append` hacks with relative package imports
       (e.g., `from rebinning.cut_edge_pix import ...`)
     - Replace hardcoded Trevor paths with paths relative to a
       configurable `DATA_DIR`

1.3. Merge the two runner scripts (`run_rebin_2026AP_FOS_STIS.py` and
     `run_rebin_2026AP_HSLA.py`) into a single `run_rebin.py` that
     handles all instruments. Keep the same loop logic.

1.4. Update `rebinning/__init__.py` with key public imports.

1.5. Copy `HST_CIV_Sulentic2007_HSLA2018_finalprops.csv` to `Data/`
     (if not already there — check for the existing copy first).

1.6. Log all changes to `Migration_Log.md`.

1.7. **Commit:** "Phase 1: migrate rebinning pipeline (2026AP)"

---

## Phase 2 — ICA Analysis Pipeline Migration

**Goal:** Copy the ICA analysis code into `ica/`.

### Source → Destination Map

| Source (`Trevor Code/`)                        | Destination (`HST Paper/ica/`) |
|------------------------------------------------|--------------------------------|
| `ICA Scripts/run_ICA_r20_components.py`        | `run_ica.py`                   |
| `ICA Scripts/fit_composites.py`                | `fit_composites.py`            |
| `ICA Scripts/plot_ICA.py`                      | `plot_ica.py`                  |
| `ICA Scripts/spec_morph.py`                    | `spec_morph.py`                |
| `ICA Scripts/BAL_regions.py`                   | `bal_regions.py`               |
| `CIV_BAL_regions.py`                           | `civ_bal_regions.py`           |
| `MorphingEdges.py` + 4 instrument variants     | `morphing_edges.py` (merged)   |
| `ica_manual_fix_refactored.py`                 | `manual_fix.py`                |
| `manual_fix_config.py`                         | `manual_fix_config.py`         |
| `Run_ICA_All_Objects.py`                       | `run_all_objects.py`           |
| `components/*.comp`, `*.dat`, `*.npy`          | `Data/components/`             |

### Steps

2.1. Copy each source file to its destination with:
     - Header comment recording origin
     - No algorithm changes

2.2. Update imports:
     - Replace `sys.path.append` hacks with package imports
       (e.g., `from ica.fit_composites import ...`,
        `from ica.spec_morph import ...`)
     - The ICA code imports `Small_Pix_Filter_TVM` — update to
       `from rebinning.small_pix_filter import SDSS_pixel_filter`

2.3. Merge all five `MorphingEdges*.py` files into a single
     `morphing_edges.py`. Each instrument's edges become a function
     or dict within the merged file. No logic changes.

2.4. Generate `quality_classification.json` by scanning
     `Trevor Code/ICA_full/{Good,Bad,Fixable,Probably Good}/`
     and mapping each object name to its quality folder. This replaces
     the runtime `os.listdir` calls in `ica_manual_fix_refactored.py`.

2.5. Update `manual_fix.py` to load quality folder from
     `quality_classification.json` instead of scanning `ICA_full/`.

2.6. Update component file paths to point to `Data/components/`.

2.7. Copy component files to `Data/components/`.

2.8. Update `ica/__init__.py` with key public imports.

2.9. Log all changes to `Migration_Log.md`.

2.10. **Commit:** "Phase 2: migrate ICA analysis pipeline"

---

## Phase 3 — Integration Testing

**Goal:** Verify the migrated code produces identical results.

### Steps

3.1. Install the package in editable mode: `pip install -e .`

3.2. Run the rebinning pipeline on a small set of test objects
     (e.g., 3-5 objects spanning FOS, STIS, HSLA) using raw data
     from `Trevor Code/SulenticAllData/` and `Trevor Code/HSLA_coadds_wCIV/`.
     - Compare output to `Trevor Code/RebinnedSpec_2026AP/`
     - Use tolerance-based comparison (not pixel-exact, per known finding)

3.3. Run ICA on the rebinned test output.
     - Compare results to existing ICA output from `Trevor Code/`

3.4. Fix any import/path issues discovered. Log fixes.

3.5. **Commit:** "Phase 3: integration testing fixes" (if any fixes needed)

---

## Phase 4 — CLAUDE.md and Documentation

**Goal:** Create the CLAUDE.md for the HST Paper repo and update README.

### Steps

4.1. Write `CLAUDE.md` covering:
     - Repo overview and purpose
     - Package structure (`rebinning/`, `ica/`, `Data/`, `Plotting Code/`)
     - How to install and run
     - Key principles (reproducibility, no functional changes to algorithms,
       versioned files)
     - Data flow: raw data → rebinning → ICA → plots
     - Critical dependencies

4.2. Update `README.md` to reflect the new repo structure (rebinning + ICA
     sections in addition to existing plotting section).

4.3. Update `requirements.txt` to match `pyproject.toml` dependencies.

4.4. **Commit:** "Phase 4: CLAUDE.md and documentation updates"

---

## Phase 5 — Connect to Master Pipeline

**Goal:** Download calibrated HST spectra from MAST for every object in
`master_catalog_v21.csv` where `has_civ=True`, then wire the downloaded
files into the existing rebinning pipeline.

### Context

`HST_SDSS_Master_Pipeline_v21_script.py` (Richards) outputs
`pipeline_output/master_catalog_v21.csv` — a metadata catalog of MAST
observations (coordinates, instruments, redshifts, CIV coverage flag).
It queries MAST via `astroquery.mast.Observations.query_criteria()` but
never downloads actual FITS files.  The individual `obs_id` values
present in the per-instrument cache CSVs (`v21_cache_mast_cos.csv`,
`v21_cache_mast_stis.csv`, `v21_cache_mast_fos.csv`) are aggregated away
in the master catalog, so the download step must re-query MAST by sky
coordinates to recover them.

We want **individual calibrated exposures** (`calib_level` 2 or 3), not
HASP/HLSP co-added products (`calib_level` 4).  Richards' v21 script
does not filter by `calib_level` in its queries (it retains the column
but applies no cut), so the level choice is entirely ours in the download
step.

**Environment:** run under the `hasp-env` conda environment, which has
`astroquery` installed.

---

### New files introduced in this phase

```
HST Paper/
├── pipeline/
│   ├── __init__.py
│   ├── download_spectra.py    # main download script (Step 5.2)
│   └── catalog_bridge.py     # catalog → rebinning adapter (Step 5.3)
└── rebinning/
    └── run_rebin.py           # updated with --catalog CLI arg (Step 5.4)
```

---

### Downloaded data folder structure

```
RAW_DATA_DIR/              # configurable; default HST_PAPER_DATA_DIR env var
└── {common_name}/         # e.g. "3C273", "Fairall9"
    ├── COS/
    │   └── {obs_id}_x1d.fits
    ├── STIS/
    │   └── {obs_id}_x1d.fits
    └── FOS/
        └── {obs_id}.c1h.fits
```

`common_name` is taken directly from the `common_name` column of
`master_catalog_v21.csv` (already cleaned by Richards in Cell 13).
The per-instrument sub-folder (`COS/`, `STIS/`, `FOS/`) matches each
observation's `inst_family` tag used throughout the v21 pipeline.

---

### Steps

**5.1. Create `pipeline/` package**

Create `HST Paper/pipeline/__init__.py` (empty).

---

**5.2. Write `pipeline/download_spectra.py`**

This script is modelled directly on Richards' v21 script:
- Imports and constants copied verbatim from Cell 1 (paths, `MATCH_SEP`,
  `DEDUP_SEP`, etc.) and adapted for the `HST Paper/` package layout.
- Obs-id recovery uses Richards' cache files directly (Cells 2–3 output)
  rather than re-querying MAST, guaranteeing we use exactly the same
  observations Richards identified.

**Script structure (cells / sections):**

```
Section 1 — Imports and configuration
    astroquery.mast.Observations
    CATALOG_PATH    = pipeline_output/master_catalog_v21.csv
    RAW_DATA_DIR    = env var HST_PAPER_DATA_DIR  (CLI override: --data-dir)
    CACHE_COS       = pipeline_output/v21_cache_mast_cos.csv   (Richards Cell 2)
    CACHE_STIS      = pipeline_output/v21_cache_mast_stis.csv  (Richards Cell 3)
    CACHE_FOS       = pipeline_output/v21_cache_mast_fos.csv   (Richards Cell 3)
    DEDUP_SEP       = 2.0 arcsec  (same value Richards uses in Cell 4)
    CALIB_LEVELS    = [2, 3]
    PRODUCT_FILTERS per instrument:
        COS  → productSubGroupDescription in ['X1D', 'X1DSUM']
        STIS → productSubGroupDescription in ['X1D', 'SX1']
        FOS  → productType == 'SCIENCE', calib_level in [2, 3]
               (FOS uses non-standard subgroup names; filter by type+level)

Section 2 — Load catalog
    Read master_catalog_v21.csv.
    Filter to has_civ == True.
    Report N objects.

Section 3 — Recover obs_ids from Richards' cache files
    Load CACHE_COS, CACHE_STIS, CACHE_FOS.
    Stack into a single obs_df; tag each row with inst_family ('COS'/'STIS'/'FOS').
    For each has_civ object (ra_deg, dec_deg, common_name):
        Spatially cross-match against obs_df using DEDUP_SEP (same as Cell 4).
        Collect all obs_df rows within DEDUP_SEP → these are Richards' obs_ids
        for this object.
        Tag matched rows with common_name.
    Filter matched rows to calib_level in CALIB_LEVELS.
    Warn for any has_civ object with zero matched obs_ids.
    Result: obs_df_matched — one row per (common_name, obs_id).

    *** Verification A — object-level completeness (pre-download) ***
    catalog_names = set(catalog[has_civ].common_name)
    matched_names = set(obs_df_matched.common_name)
    missing = catalog_names - matched_names   # has_civ objects with no obs_ids
    extra   = matched_names - catalog_names   # should always be empty
    Print: N matched, N missing, N extra.
    Raise a warning (not an error) if missing is non-empty so the user
    can investigate before committing to the download.
    Script continues only when missing is empty or user passes --force.

Section 4 — Get product lists (with caching)
    PROD_CACHE = pipeline_output/v21_cache_download_products.csv
    If PROD_CACHE exists: load it.
    Else: for each unique obs_id in obs_df_matched:
        products = Observations.get_product_list(obs_id)
        tag with common_name, inst_family
        append to prod_rows list
    Filter product list:
        productType == 'SCIENCE'
        calib_level in CALIB_LEVELS
        productSubGroupDescription matches per-instrument filter above
    Save filtered products to PROD_CACHE.

Section 5 — Download
    for each common_name:
        sub = filtered products for this object
        for inst in ['COS', 'STIS', 'FOS']:
            inst_sub = sub rows for this instrument
            if empty: skip
            dest = RAW_DATA_DIR / common_name / inst
            dest.mkdir(parents=True, exist_ok=True)
            Observations.download_products(
                inst_sub,
                download_dir=str(dest),
                flat=True,          # no extra nested subdirs
            )
    Print per-object download summary (N files, total MB).

Section 6 — Download manifest and verification
    Write pipeline_output/download_manifest.csv:
        common_name | inst_family | obs_id | filename | file_size_mb | status

    *** Verification B — object-level completeness (post-download) ***
    downloaded_names = set(manifest[manifest.status == 'ok'].common_name)
    missing_after    = catalog_names - downloaded_names
    Print: N successfully downloaded, N missing.
    If missing_after is non-empty: list the object names explicitly so
    they can be investigated or re-run individually.

    Used by catalog_bridge.py to build the input list for run_rebin.py.
```

---

**5.3. Write `pipeline/catalog_bridge.py`** ✓ *Complete — 2026-05-02*

Reads `master_catalog_v21.csv` and `download_manifest.csv`, then
produces the per-object input dict that `rebinning/run_rebin.py` expects:

```python
{
    'name':       common_name,
    'redshift':   best_z,
    'instrument': inst_family,          # 'COS', 'STIS', 'FOS'
    'data_path':  RAW_DATA_DIR / common_name / inst_family,
}
```

Key logic:
- One dict per (common_name, inst_family) pair (multi-instrument objects
  produce multiple entries, matching how the existing pipeline handles
  e.g. FOS+STIS objects).
- Validates that `data_path` exists and contains at least one FITS file;
  warns and skips if not.
- Returns a list of dicts; callers can iterate directly or pass to
  `run_rebin.py`.

Implementation notes (2026-05-02):
- `master_catalog_v21.csv` carries no SDSS plate/MJD/fiber columns, so
  `fn_sdss=None` for all objects when `coadd.rebin()` is called. SDSS
  coadding is not supported for master-catalog objects at this stage.
- `download_spectra.py` is NOT imported by `catalog_bridge.py` (to avoid
  pulling in the `astroquery` dependency at rebin time); the two files
  define their own `CATALOG_PATH` / `MANIFEST_PATH` constants independently.
- No quality classification (`Good`/`Bad`/`Probably Good`) is used.
  The old `ICA_full/` folder categories apply only to the Sulentic-2007
  sample; new objects from the master catalog have no such category, and
  the ICA pipeline must run on all objects without pre-filtering by quality.
- The old Sulentic catalog (`HST_CIV_Sulentic2007_HSLA2018_finalprops.csv`)
  is obsolete for all new pipeline runs; `master_catalog_v21.csv` is the
  sole source of object metadata from Phase 5 onward.

---

**5.4. Update `rebinning/run_rebin.py` with `--catalog` CLI argument**

Add a `--catalog` argument (path to `master_catalog_v21.csv`) alongside
`--data-dir`.  When `--catalog` is supplied, the runner calls
`catalog_bridge.build_object_list()` instead of reading from the
hardcoded Trevor path layout.  The existing hardcoded path logic remains
as a fallback so the Phase 3 integration tests are not broken.

No changes to any algorithm logic — only the object-list construction
at the top of the runner is affected.

---

**5.5. Log all changes to `Migration_Log.md`.**

---

**5.6. End-to-end smoke test**

Run on 3–5 representative objects (one FOS-only, one STIS-only, one
COS-only, one multi-instrument) and confirm:
- Files download to the correct folder structure.
- `catalog_bridge.py` produces valid input dicts for each.
- `run_rebin.py --catalog ... --data-dir ...` completes without error.
- Output FITS files in `RebinnedSpec_*/` are plausible (non-zero flux,
  correct wavelength range for the object's redshift).

---

**5.7. Commit:** `"Phase 5: download pipeline and catalog bridge"`

---

## File-Level Dependency Graph (for reference)

### Rebinning

```
run_rebin.py
  └── coadd.py (HST_coadd2_2026AP)
        ├── read_spec_data.py (Read_spec_data_2026AP)
        │     ├── cut_edge_pix.py (Cut_Edge_Pix_TVM_NoDQ)
        │     ├── lower_res_rebin.py (LowerResHSTRebin_TVM)
        │     ├── spec_cuts_fos.py (SpecCuts_FOS)
        │     └── spec_cuts_hsla.py (SpecCuts_HSLA)
        └── small_pix_filter.py (Small_Pix_Filter_TVM)  [SDSS coadd path]
```

### ICA

```
run_all_objects.py (Run_ICA_All_Objects)
  └── manual_fix.py (ica_manual_fix_refactored)
        ├── manual_fix_config.py
        ├── quality_classification.json
        ├── run_ica.py (run_ICA_r20_components)
        │     ├── fit_composites.py
        │     ├── plot_ica.py
        │     └── spec_morph.py
        ├── bal_regions.py (BAL_regions)
        ├── civ_bal_regions.py (CIV_BAL_regions)
        ├── morphing_edges.py (MorphingEdges_*)
        └── rebinning.small_pix_filter (cross-package import)
```
