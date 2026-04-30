# Migration Log

**Project:** Rebinning + ICA Code Migration into HST Paper Repo
**Started:** 2026-04-29

All changes made during each phase are logged here with dates, file names,
and descriptions of what was done.

---

## Known Issues / Future Work

### Two diverged versions of `spec_morph.py`

**Logged:** 2026-04-30

There are two distinct versions of `spec_morph.py` in `Trevor Code/`, and both have been
migrated into this repo as separate files:

| Location in repo | Source in Trevor Code | Used by |
|------------------|-----------------------|---------|
| `rebinning/spec_morph.py` | `Trevor Code/spec_morph.py` | Rebinning pipeline (`coadd.py`) |
| `ica/spec_morph.py` *(Phase 2, not yet migrated)* | `Trevor Code/ICA Scripts/spec_morph.py` | ICA pipeline (`run_ICA_r20_components.py`) |

**Key differences between the two versions:**
- The ICA version adds `sys.path.append` (removed on migration) and loads
  `qsomod_z250_ng_continuumSDSSrebin.npy` from `components/` (rebinning version loads
  it from the same relative path, now `Data/components/`).
- The rebinning version has AP's defensive fixes to `cont_filtered()` (bounds checking
  on emission line interpolation) and a standalone `cont_filtered_AP()` function.
  These are absent from the ICA version.
- `continuum_fit()` signatures differ: the ICA version takes extra `value` and
  `markers_arr` parameters; the rebinning version uses the original fixed signature.
- The ICA version uses `np.polynomial.polynomial.Polynomial.fit` for the continuum
  fit; the rebinning version uses `np.polyfit`.

**Why this matters:** If a bug is fixed or an improvement is made to one version, it
may need to be manually ported to the other. Before doing so, confirm which version
was actually used to produce the validated `RebinnedSpec_2026AP/` and `ICA_full/`
results (the split is intentional — the rebinning and ICA runs each used their own
`spec_morph.py`).

**Future fix options:**
1. Consolidate into a single `spec_morph.py` (requires careful review of both versions
   to ensure no scientific differences are lost).
2. Keep them separate but add a note in each file cross-referencing the other.
3. Factor out common functions into a shared `utils/spec_morph_common.py`.

Option 2 (keep separate, document clearly) is the safest path until the two versions
are fully understood.

---

## Phase 0 — Project Scaffolding

**Date:** 2026-04-29

### Files Created

| File | Description |
|------|-------------|
| `CLAUDE.md` | Codebase guide: repo structure, principles, data flow, migration status table |
| `pyproject.toml` | Installable package `hst_civ_paper`; Python >=3.12; all dependencies listed |
| `rebinning/__init__.py` | Empty package init (to be populated in Phase 1) |
| `ica/__init__.py` | Empty package init (to be populated in Phase 2) |
| `Data/components/.gitkeep` | Placeholder so git tracks the empty components directory |

### Files Modified

| File | Change |
|------|--------|
| `.gitignore` | Added `RebinnedSpec_*/`, `ICA_Plots_*/`, `ICA_full/` to ignore large pipeline output directories |

### Notes
- `Migration_Log.md` header was already present (created alongside `Migration_Plan.md`)
- `.gitignore` already had `*.egg-info/`, `dist/`, `build/` from previous scaffold; no duplication needed

---

## Phase 1 — Rebinning Pipeline Migration

**Date:** 2026-04-30

### New files created in `rebinning/`

| Destination | Source (`Trevor Code/`) | Changes from source |
|-------------|------------------------|---------------------|
| `spec_cuts_fos.py` | `SpecCuts_FOS.py` | Header added only |
| `spec_cuts_stis.py` | `SpecCuts_STIS.py` | Header added only. **Note:** not in original Phase 1 plan; added because `cut_edge_pix.py` and `coadd.py` have STIS code paths and the FOS+STIS runner processes STIS objects. |
| `spec_cuts_cos.py` | `SpecCuts_COS.py` | Header added only. Included so the COS import path in `cut_edge_pix.py` resolves (no COS objects in 2026AP run). |
| `spec_cuts_hsla.py` | `SpecCuts_HSLA.py` | Header added only |
| `cut_edge_pix.py` | `Cut_Edge_Pix_TVM_NoDQ.py` | Instrument SpecCuts imports updated to `from rebinning import spec_cuts_* as SpecCuts` |
| `small_pix_filter.py` | `Small_Pix_Filter_TVM.py` | Header added only |
| `lower_res_rebin.py` | `LowerResHSTRebin_TVM.py` | `fits.open('spec-0266-51630-0080.fits')` → `fits.open(str(_SDSS_REF_SPEC))` where `_SDSS_REF_SPEC` resolves to `Data/spec-0266-51630-0080.fits` via `Path(__file__)`. Added `from pathlib import Path`. |
| `bal_regions.py` | `BAL_regions.py` | Header added only. BAL_regions.py is identical between `Trevor Code/` and `Trevor Code/ICA Scripts/`. Placed in `rebinning/` because `spec_morph.py` (rebinning version) imports it. Phase 2 ICA can import via `from rebinning import bal_regions`. |
| `spec_morph.py` | `spec_morph.py` (Trevor Code root) | This is the **rebinning version**; differs from `ICA Scripts/spec_morph.py` (ICA version, for Phase 2). Import changes: `import Small_Pix_Filter_TVM` → `from rebinning import small_pix_filter as Small_Pix_Filter_TVM`; `import BAL_regions` → `from rebinning import bal_regions as BAL_regions`. Path changes: `np.load("qsomod_...")` → `np.load(str(_QSOMOD_FILE))` (Data/components/); `np.load("/Users/Trevor1/...")` → `np.load(str(_VANDENBERK_FILE))` (Data/components/ — only needed if `morph()` is called; `morph()` is NOT called by the rebinning pipeline). |
| `read_spec_data.py` | `Read_spec_data_2026AP.py` | Import changes: `import Cut_Edge_Pix_TVM_NoDQ as Cut_Edge_Pix_TVM` → `from rebinning import cut_edge_pix as Cut_Edge_Pix_TVM`; `import SpecCuts_HSLA` → `from rebinning import spec_cuts_hsla as SpecCuts_HSLA`. Path change in `read_hsla()`: `path+"original\\%s"%Identifier` → `os.path.join(path, "original", Identifier)` (cross-platform). Added `import os`. |
| `coadd.py` | `HST_coadd2_2026AP.py` | Import changes: `import Small_Pix_Filter_TVM/LowerResHSTRebin_TVM/Read_spec_data_2026AP/spec_morph` → `from rebinning import ...` equivalents. Unused imports commented out: `Large_pix_filter_TVM`, `Plot_HST`, `Cut_Edge_Pix_TVM_NoDQ` (not directly called in this file). `MorphingEdges_*` imports inside `rebin()` commented out (vestigial — imported but never called; Phase 2 will add). Path changes: `fits.open('spec-0266...')` → `fits.open(str(_SDSS_REF_SPEC))`; `fits.open("SDSS_spec/lite/%s")` → `fits.open(os.path.join(sdss_spec_dir, fn_sdss))`; `t.write("RebinnedSpec_2026AP/...")` → `t.write(os.path.join(output_dir, ...))`. Signature: added `output_dir` and `sdss_spec_dir` keyword arguments. |
| `run_rebin.py` | `run_rebin_2026AP_FOS_STIS.py` + `run_rebin_2026AP_HSLA.py` | Merged into single script. `import HST_coadd2_2026AP` → `from rebinning import coadd`. All hardcoded paths replaced with `DATA_DIR`-relative paths via `argparse` (`--data-dir`, `--output-dir`, `--sdss-spec-dir`, `--catalog`, `--instruments`). `DATA_DIR` defaults to `HST_PAPER_DATA_DIR` env var. Default catalog: `Data/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv` in repo root. |
| `__init__.py` | (new — was empty placeholder) | Added docstring and public imports for all rebinning submodules. |

### Data files copied to `Data/`

| File | Source |
|------|--------|
| `Data/spec-0266-51630-0080.fits` | `Trevor Code/spec-0266-51630-0080.fits` |
| `Data/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv` | `Trevor Code/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv` |
| `Data/components/qsomod_z250_ng_continuumSDSSrebin.npy` | `Trevor Code/components/qsomod_z250_ng_continuumSDSSrebin.npy` |

### Notes
- `vandenberk01_continuum.npy` (used only by `spec_morph.morph()`, which is NOT called
  by the rebinning pipeline) was not copied. The file lives in `Richards/HST/ICA/` (archival).
  If `morph()` is ever needed, copy it to `Data/components/`.
- `Trevor Code/` source files are **untouched**. All Phase 1 files are new creations.
- The STIS and COS SpecCuts files were added beyond the original plan because they are
  required import targets; this is noted in the plan as a deviation.

---

## Phase 2 — ICA Analysis Pipeline Migration

*(not yet started)*

---

## Phase 3 — Integration Testing

*(not yet started)*

---

## Phase 4 — CLAUDE.md and Documentation

*(not yet started)*

---

## Phase 5 — Connect to Master Pipeline

*(not yet started)*
