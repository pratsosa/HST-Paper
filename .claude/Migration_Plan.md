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

## Phase 5 — Connect to Master Pipeline (future, loosely defined)

**Goal:** Wire the v21 master pipeline output to the rebinning entry point.

### What We Know

- `HST_SDSS_Master_Pipeline_v21_script.py` outputs
  `pipeline_output/master_catalog_v21.csv` — a catalog of observations
  queried from MAST, with coordinates, instruments, redshifts, etc.
- It does **not** download the actual FITS files.
- A data download step is needed (to be added to v21 or as a bridge script).
- Once data is downloaded into an organized folder structure, `run_rebin.py`
  needs to be updated to read from it instead of from Trevor's
  `SulenticAllData/` layout.

### Steps (to be detailed later)

5.1. Define the download folder structure with Prof. Richards.

5.2. Write or extend the pipeline to download FITS files from MAST
     into the defined structure.

5.3. Write a bridge module that reads `master_catalog_v21.csv` and
     produces the inputs that `run_rebin.py` expects (object name,
     redshift, instrument, SDSS spec path, data path).

5.4. End-to-end test: pipeline → download → rebin → ICA → plots.

5.5. **Commit:** "Phase 5: connect to master pipeline"

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
