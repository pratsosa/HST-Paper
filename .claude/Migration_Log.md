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

**Date:** 2026-04-29

### New files created in `ica/`

| Destination | Source (`Trevor Code/`) | Changes from source |
|-------------|------------------------|---------------------|
| `spec_morph.py` | `ICA Scripts/spec_morph.py` | Removed `sys.path.append`; `import Small_Pix_Filter_TVM` → `from rebinning import small_pix_filter as Small_Pix_Filter_TVM`; `import BAL_regions` → `from rebinning import bal_regions as BAL_regions`; both `np.load('components/qsomod_...')` calls → `np.load(str(_QSOMOD_FILE))` via `Path(__file__)`. NOTE: this is the ICA version — differs from `rebinning/spec_morph.py` (see Known Issues). |
| `fit_composites.py` | `ICA Scripts/fit_composites.py` | Removed `sys.path.append`; dead imports `plot_ICA`, `spec_morph` commented out; `path_priors` default changed from `"./c3_priors/"` to `str(_C3_PRIORS_DIR) + "/"` where `_C3_PRIORS_DIR = Path(__file__).parent.parent / "Data" / "c3_priors"`. |
| `plot_ica.py` | `ICA Scripts/plot_ICA.py` | Header added only. |
| `run_ica.py` | `ICA Scripts/run_ICA_r20_components.py` | Removed `sys.path.append`; `import fit_composites/plot_ICA/spec_morph` → `from ica import ...`; `load_ICA()` and all functions with `ica_path` now default to `str(_COMPONENTS_DIR) + "/"` via `Path(__file__)` instead of `"./"`. `sep="\s+"` → `sep=r"\s+"` (silence FutureWarning). |
| `civ_bal_regions.py` | `CIV_BAL_regions.py` | Header added only. |
| `morphing_edges.py` | `MorphingEdges.py` + `MorphingEdges_FOS/STIS/COS/HSLA.py` | All five files merged. `MorphingEdges` dict kept as `MorphingEdges`. Instrument-specific `MorphMarkers` dicts renamed: `FOS_MorphMarkers`, `STIS_MorphMarkers`, `COS_MorphMarkers`, `HSLA_MorphMarkers`. No algorithm changes. |
| `quality_classification.json` | `ICA_full/{Good,Bad,Fixable,Probably Good}/` | Generated by scanning folder filenames; maps `Spec_Name` prefix → quality folder. 207 entries (all objects). |
| `manual_fix.py` | `ica_manual_fix_refactored.py` | Removed `sys.path.append`; `from Small_Pix_Filter_TVM import` → `from rebinning.small_pix_filter import`; `import MorphingEdges_*` → `from ica import morphing_edges`; `import CIV_BAL_regions` → `from ica import civ_bal_regions as CIV_BAL_regions`; `import run_ICA_r20_components` → `from ica import run_ica as run_ICA_r20_components`; `import spec_morph` → `from ica import spec_morph`; dead `import components` removed; `from manual_fix_config import` → `from ica.manual_fix_config import`; `__init__` catalog path → `_DATA_DIR`-relative default; `__init__` loads `quality_classification.json` into `self.quality_map`; `get_folder()` replaced `os.listdir` scan with JSON lookup. |
| `manual_fix_config.py` | `manual_fix_config.py` | Header added only. |
| `run_all_objects.py` | `Run_ICA_All_Objects.py` | Removed `sys.path.append`; `from ica_manual_fix_refactored import` → `from ica.manual_fix import`; catalog CSV path → `Path(__file__)`-relative. |
| `__init__.py` | (was empty placeholder) | Added docstring and public imports for all ica submodules. |

### Data files copied to `Data/`

| File | Source |
|------|--------|
| `Data/components/amy_12603000_10c_180421.comp` | `Trevor Code/components/` |
| `Data/components/amy_12653000_hew_hsn_7c_190302.comp` | `Trevor Code/components/` |
| `Data/components/amy_12753000_lowew_10c_181101_v1.comp` | `Trevor Code/components/` |
| `Data/components/wav_12603000.dat` | `Trevor Code/components/` |
| `Data/components/wav_12653000.dat` | `Trevor Code/components/` |
| `Data/components/wav_12753000.dat` | `Trevor Code/components/` |
| `Data/c3_priors/` (full directory) | `Trevor Code/c3_priors/` |

### Notes
- `Trevor Code/` source files are **untouched**. All Phase 2 files are new copies.
- `morphing_edges.py` was built by concatenating all 5 source files; the `MorphMarkers` dicts were renamed to avoid name collision.
- `quality_classification.json` was generated programmatically; the JSON key is the `Spec_Name` prefix (text before `_{inst}` in the ICA_full filename).
- The `get_folder()` method in `manual_fix.py` no longer requires `ICA_full/` to be present at runtime — it reads from the committed JSON instead.
- `c3_priors/` added to `Data/` beyond the original Phase 2 plan (AP confirmed this is correct).

---

## Phase 3 — Integration Testing

*(not yet started)*

---

## Phase 4 — CLAUDE.md and Documentation

*(not yet started)*

---

## Phase 5 — Connect to Master Pipeline

*(not yet started)*
