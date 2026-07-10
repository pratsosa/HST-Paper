# Migration Log

**Project:** Rebinning + ICA Code Migration into HST Paper Repo
**Started:** 2026-04-29

All changes made during each phase are logged here with dates, file names,
and descriptions of what was done.

---

## Known Issues / Future Work

### HSLA-only objects missed by Phase 5 download (19 objects) — RESOLVED

**Logged:** 2026-05-02
**Resolved:** 2026-05-02 (re-run with new recovery logic recovered all 19; 0 missing)

**Symptom:** First end-to-end run of `pipeline/download_spectra.py` succeeded for 429 of
448 has_civ=True objects in `master_catalog_v21.csv`. The 19 missing objects were listed
in `missing_obj.txt` at the repo root.

**Root cause:**
1. Each of the 19 objects matches *exactly one* row in `v21_cache_mast_cos.csv` at
   <0.001 arcsec separation, so Verification A passed and they entered `obs_df_matched`.
2. **Their only matched obs_id is an HSLA HLSP co-add (`hst_hsla_*`)**. No underlying
   individual COS exposure rows are present in the cache for these objects.
3. Why the cache lacks individual exposures: Richards' Cell 2 COS query filters by
   `target_classification` against AGN keywords (`*Active galaxy*`, `*Quasar*`,
   `*Seyfert*`, ...). The HSLA HLSP entry's classification matches (HSLA tags rich
   classifications like `'UAT=BL Lacertae objects; TIER1=Galaxy; TIER2=AGN'`), but the
   underlying *individual* exposures often have a sparser classification (e.g. `Galaxy`
   or empty) and so were never returned by Cell 2's classification-filtered query.
4. When `download_spectra.py` calls `Observations.get_product_list('hst_hsla_*')` the
   call fails (HSLA HLSP obs_ids do not behave like standard exposure obs_ids in the
   MAST product API). The fallback `query_criteria(obs_id=...)` also fails. The
   exception handler prints a warning and `continue`s, silently dropping the object.
5. End result: 19 objects entered cross-match → 0 products in
   `v21_cache_download_products.csv` → 0 rows in `download_manifest.csv` → no folder
   under `RAW_DATA_DIR/{common_name}/`.

**Confirmed pattern:** Across the full has_civ=True catalog (448 objects), exactly 19
objects have all-HSLA matches. Every other object has at least one non-HSLA match,
so the silent skip does not affect them.

**Fix implemented (2026-05-02):** option (b) from triage — re-query MAST by sky
position for the underlying individual COS exposures of HSLA-only objects. See the
Phase 5 entry below for the new `recover_hsla_only_exposures()` function in
`pipeline/download_spectra.py` and the new `v21_cache_hsla_recovery.csv`.
After the re-run, all 19 objects were recovered and Verification B showed 0 missing.

**Open question (still pending):** AP is checking with Richards whether HSLA HLSP
co-adds should also be downloaded as a parallel data source (mirroring Trevor's
existing `Trevor Code/HSLA_coadds_wCIV/`). If yes, a separate HSLA download path
will be added. For now we proceed with individual-exposure rebinning only.

---

### COS NUV IndexError in `read_cos_flat` — RESOLVED

**Logged:** 2026-06-01
**Resolved:** 2026-06-01

**Symptom:** `python -m rebinning.run_rebin --master-catalog` failed for some COS objects
(e.g. `2MASS J10115298+5442063`) with:
```
IndexError: boolean index did not match indexed array along dimension 0;
dimension is 16384 but corresponding boolean dimension is 1274
```
at `rebinning/read_spec_data.py:426` (`masks[i,:][(err_wmask == 0.0)] = 1`).

**Root cause:** The dimension mismatch `16384 vs 1274` is COS detector geometry:
- COS FUV (G130M/G140L/G160M): 16384 px per segment
- COS NUV (G185M/G225M/G230L/G285M): 1274 px per stripe

When an object has a mix of FUV and NUV exposures, `read_cos_flat()` sets
`array_len = max(array_sizes) = 16384` and allocates `waves/fluxes/flux_errs/masks`
at that width. For the NUV file, `wavelength`/`fluxerr` are length 1274, but the
write target `masks[i,:]` is length 16384 — the boolean index has the wrong shape.

The FOS reader (lines 350-359) already handles the same situation correctly with
`[:len(err_wmask)]` slicing; the COS flat reader was missing that pattern. The
padding-with-zeros convention is intended — `coadd.py:144` filters with
`waves[waves!=0]`.

**Why this didn't surface in the 2026AP rebin run:** Trevor's original
`read_cos()` in `Read_spec_data_2026E.py` has the same latent bug (line 167), but
none of the 207 objects in `RebinnedSpec_2022Aug11/` combined COS FUV and NUV
exposures, so it never triggered. The new master-catalog flow surfaces it because
the broader object set includes mixed-detector COS data.

**Fix implemented:** In `rebinning/read_spec_data.py` `read_cos_flat()`, replaced
all `masks[i,:]` / `waves[i,:]` / `fluxes[i,:]` / `flux_errs[i,:]` write targets
in the per-file loop with `[:len(err_wmask)]` slicing, matching the FOS reader
pattern. Inline comment added at the assignment site.

**Open follow-up (not fixed):** The `else` branch at line 418
(`data['WAVELENGTH'][0]`) reads only the first stripe/segment for any non-E140M
grating. For COS NUV x1d files this drops NUVB and NUVC (only NUVA is read); for
COS FUV it would drop segment B if both segments lived in one row. This matches
Trevor's original behavior and was left unchanged. To capture all stripes, the
E140M multi-row concat pattern (lines 408-417) would need to be extended to any
`data.size > 1`. Confirm with AP / Richards whether NUV stripe coverage matters
for the paper's science scope before changing.

---

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

**Date:** 2026-04-30

### Bug fixed: dead `wave_empty` allocation in `lower_res_rebin.py`

**File:** `rebinning/lower_res_rebin.py` (line 42)

**Error:** `numpy.core._exceptions._ArrayMemoryError: Unable to allocate 687. MiB for an array with shape (90000001,) and data type float64`

**Root cause:** `HSTLowResRebin()` contained the line:
```python
wave_empty = np.arange(1000., 10000.+sdss_c1, sdss_c1)
```
`sdss_c1 ≈ 0.0001` (SDSS log-lambda spacing), making the array 90,000,001 elements (687 MiB). **This variable is never read again** — it is dead code in both the original `LowerResHSTRebin_TVM.py` and the migrated file. The 2026AP run presumably succeeded because more RAM was free at that time. STIS objects (smaller data footprint) happened to leave enough memory; HSLA coadd files (larger) did not.

**Fix:** Commented out the dead line:
```python
# wave_empty = np.arange(...)  # AP 2026-04-30: dead code — never read; removed to avoid 687 MiB allocation
```

**Scope:** No algorithmic change. `wave_empty` was created and discarded immediately in the original code. Removing it does not affect any output.

**Note:** The same dead line exists in the original `Trevor Code/LowerResHSTRebin_TVM.py`. It is NOT fixed there (per the policy of keeping Trevor Code untouched). The migrated version is now strictly better than the original on this point.

### Test objects run successfully

| Instrument | Object | Status |
|------------|--------|--------|
| FOS | J04232-0120 | OK |
| STIS | J07086-4933 | OK |
| HSLA | 1H1613-097 | OK (after fix above) |

---

## Phase 4 — CLAUDE.md and Documentation

**Date:** 2026-05-01

### Files Modified

| File | Change |
|------|--------|
| `CLAUDE.md` | Migration status table updated: phases 1, 2, 3 marked Complete |
| `README.md` | Full rewrite: now covers complete pipeline (rebinning → ICA → plotting), installation, data setup, and all CLI options. Previous version covered plotting figures only. |
| `requirements.txt` | Added `lmfit>=1.0.0` and `weightedstats>=0.4.0` to match `pyproject.toml` dependencies |

### Notes
- `QUICKSTART.md` retained as-is (still accurate for the plotting-only quick start path)
- `CLAUDE.md` required no structural changes beyond the status table

---

## Phase 5 — Connect to Master Pipeline

**Date:** 2026-05-01

### Steps completed: 5.1, 5.2 (download stage complete 2026-05-02 after HSLA-only recovery patch)

### Data copied into repo

| Source (Richards/) | Destination (HST Paper/) |
|--------------------|--------------------------|
| `pipeline_output/v21_cache_mast_cos.csv` | `pipeline_output/v21_cache_mast_cos.csv` |
| `pipeline_output/v21_cache_mast_stis.csv` | `pipeline_output/v21_cache_mast_stis.csv` |
| `pipeline_output/v21_cache_mast_fos.csv` | `pipeline_output/v21_cache_mast_fos.csv` |
| `master_catalog_v21.csv` | `pipeline_output/master_catalog_v21.csv` |

### Files created

| File | Description |
|------|-------------|
| `pipeline/__init__.py` | Empty package init (Step 5.1) |
| `pipeline/download_spectra.py` | MAST download script (Step 5.2) |

### Design decisions and filtering rationale

**Object-level filtering**: `master_catalog_v21.csv` with `has_civ=True` is the
direct output of Richards' Cells 1–8 (coordinate validation, redshift priority,
z-conflict resolution, UV ceiling cap, CIV coverage criterion). No replication
needed — we start from it directly.

**Observation-level coordinate filtering**: Mirrors Cell 4 exactly:
```python
valid = (ra.notna() & dec.notna() &
         (ra >= 0) & (ra <= 360) &
         (dec >= -90) & (dec <= 90))
```
This drops calibration-lamp rows (MAST target_name='WAVE', no sky coordinates)
and the MAST sentinel value -1.0. Cell 4 does NOT filter by calib_level at the
observation stage, and neither do we.

**Spatial cross-match**: The catalog RA/Dec values are the mean cluster positions
from Cell 4's greedy friends-of-friends. A radius search with the same
`DEDUP_SEP=2.0"` recovers exactly the obs_ids Richards grouped into each cluster.

**Product-level filtering**: `productType='SCIENCE'`, `calib_level in [2, 3]`
(excludes HASP/HLSP calib_level=4). Per-instrument subgroup filter:
COS → X1D/X1DSUM; STIS → X1D/SX1; FOS → type+level only (non-standard names).

**Environment**: must run under `hasp-env` conda environment (astroquery).

### Steps remaining in Phase 5

- 5.3: `pipeline/catalog_bridge.py`
- 5.4: `rebinning/run_rebin.py` `--catalog` CLI argument
- 5.5: Log (this entry)
- 5.6: End-to-end smoke test
- 5.7: Commit

### 2026-05-02 — HSLA-only recovery added to `download_spectra.py`

Triggered by the post-run finding (see Known Issues) that 19 has_civ objects had been
silently skipped because their only matched observations in the cache were HSLA HLSP
co-add obs_ids (`hst_hsla_*`).

**Code changes in `pipeline/download_spectra.py`:**

| Item | Description |
|------|-------------|
| New constant `HSLA_RECOVERY_CACHE` | `pipeline_output/v21_cache_hsla_recovery.csv` |
| New constant `HSLA_RECOVERY_RADIUS_ARCSEC = 30.0` | Cone-search radius (matches Cell 9b) |
| New function `recover_hsla_only_exposures(catalog, obs_df_matched)` | Detects HSLA-only objects, cone-searches MAST (`obs_collection='HST'`, `instrument_name='COS*'` and `'STIS*'`, `dataproduct_type='spectrum'`, no classification filter), drops HLSP/HSLA rows from results, returns recovered individual exposures. Cached to `HSLA_RECOVERY_CACHE`. |
| Refactored `_filter_products(...)` and `_query_products_for_obs(...)` | Helper functions extracted from `get_filtered_products()` so the bulk and incremental paths share logic. |
| `get_filtered_products()` is now incremental | If `PROD_CACHE` exists, only obs_ids in `obs_df_matched` that are NOT yet in the cache are queried. New filtered products are appended and saved. This is what lets the recovered HSLA-only obs_ids actually get queried on a re-run without re-querying everything. |
| New Section 3.5 in `main()` | Calls `recover_hsla_only_exposures(...)`, drops the original HSLA-only rows for objects covered by recovery, and concats the recovered individual-exposure rows into `obs_df_matched`. Warns if any object remains HSLA-only after recovery. |

**No algorithm change to existing successful path.** Objects with non-HSLA matches
behave identically to the prior run; the cache is reused as-is.

**Reproducibility:** All MAST queries during recovery are cached to
`pipeline_output/v21_cache_hsla_recovery.csv`. Subsequent runs read the cache and
skip the MAST round-trip.

**Verification:** Re-run on 2026-05-02 recovered all 19 HSLA-only objects;
Verification B reported 0 missing post-download. Stage 5.2 (download script) is
now considered complete.

**Open question (deferred to Richards):** Whether HSLA HLSP coadds should be
downloaded as a parallel data product (mirroring `Trevor Code/HSLA_coadds_wCIV/`).
If yes, a separate HSLA download path will be added — this recovery does not
address that question; it only ensures the 19 HSLA-only objects can be
processed via the individual-exposure rebinning path.

---

### 2026-06-01 — NaN guard added to `HSTLowResRebin` in `lower_res_rebin.py`

**Symptom:** Every object (FOS, STIS, COS) failed with:
```
Intel MKL ERROR: Parameter 6 was incorrect on entry to DGELSD.
LinAlgError: SVD did not converge in Linear Least Squares
```
Traceback pointed to `lower_res_rebin.py` line 100:
```python
m_flux, b_flux = np.polyfit([wave1, wave2], [flux1, flux2], 1)
```

**Root cause:** `Cut_Edge_Pix` sets the flux array to NaN beyond `end_stop_index`
(the trailing edge of the good-data region), but it does **not** zero out the
corresponding wavelengths — those remain non-zero beyond `end_stop_index`.
`HSTLowResRebin` builds its rebinned grid using `max(wavelength[wavelength!=0])`,
which extends all the way to the last non-zero wavelength, past `end_stop_index`.
When the inner loop reaches those trailing positions, `flux[arg1]` or `flux[arg2]`
is NaN. `np.polyfit([wave1, wave2], [NaN, flux2], 1)` passes NaN directly to
LAPACK's DGELSD, which rejects the input with the MKL parameter error rather than
returning NaN cleanly.

**Why this did not fail with Trevor's original data:** Trevor's curated
`SulenticAllData/` exposures had good S/N throughout nearly the full wavelength
array, so `end_stop_index` was close to the actual array end, leaving little or no
NaN tail. MAST-downloaded spectra can have larger degraded edge regions
(more bad pixels toward the red/blue limits), pushing `end_stop_index` well before
the end of the wavelength array and creating a significant NaN tail that the grid
walks into.

**Fix implemented (Option A) — `rebinning/lower_res_rebin.py`:**
Added a NaN guard at the top of the inner loop in `HSTLowResRebin`. If either
lookup pixel has non-finite flux, the rebinned pixel is set to NaN flux, 0 error,
and mask=1 (bad), then skipped:
```python
if not (np.isfinite(flux1) and np.isfinite(flux2)):
    flux_HST_rebin[i]    = np.nan
    fluxerr_HST_rebin[i] = 0.
    masks_HST_rebin[i]   = 1
    continue
```
Pixels with error=0 get weight=1/0²=inf, which the co-add in `coadd.py` immediately
zeros out — so these pixels are silently excluded from the variance-weighted median
and do not affect the final co-added spectrum.

**Scope:** The NaN guard fires only for the trailing NaN tail (and any sparse
interior NaN pixels, which are rare for properly calibrated STIS/FOS data). For
Trevor's data it would never fire. No science output is changed for clean data.

**Option B (not implemented — for future reference):**
A complementary fix would restrict the rebinned grid extent to the valid flux
region by changing the two `npix_from_*` lines in `HSTLowResRebin`:
```python
# Current code:
npix_from_startHST = (min(np.log10(wavelength[wavelength!=0])) - sdss_c0) // sdss_c1
npix_from_endHST   = (max(np.log10(wavelength[wavelength!=0])) - sdss_c0) // sdss_c1

# Option B replacement:
valid = (wavelength != 0) & np.isfinite(flux)
npix_from_startHST = (min(np.log10(wavelength[valid])) - sdss_c0) // sdss_c1
npix_from_endHST   = (max(np.log10(wavelength[valid])) - sdss_c0) // sdss_c1
```
This stops the grid from being built over the NaN tail entirely, producing a
slightly shorter but cleaner rebinned spectrum. Option B alone is not sufficient
(interior NaN pixels would still reach polyfit) but combined with Option A it
would be the most semantically correct solution. Option B is an algorithmic change
in the sense that output spectra would be shorter by the NaN tail; for Trevor's
data the change would be negligible. Implement Option B if the extra masked tail
pixels in the output FITS files cause issues downstream (e.g., in the ICA pipeline).

---

## Hardening of flat-layout COS/STIS readers against MAST-archived empty x1d files

**Date:** 2026-06-01
**Files modified:** `rebinning/read_spec_data.py` (`read_cos_flat`, `read_stis_flat`)

### Symptom

Master-catalog runs (`python -m rebinning.run_rebin --master-catalog --data-dir raw_data`)
crashed for some COS objects with:
```
IndexError: index 0 is out of bounds for axis 0 with size 0
  at read_spec_data.py:388 → array_sizes.append(len(data['Wavelength'][0]))
```
Reproducible failures on `2MASS J00435499+4234302` and `2MASS J00505073+3536429`.

### Root cause

Two distinct problems, both originating in the bulk MAST download path (not in
hand-curated `SulenticAllData/` folders, which is why Trevor's pipeline never
hit either):

**(1) Empty CalCOS extraction products.** When CalCOS at STScI fails to extract
a spectrum (guide-star loss, lamp not firing, count rate too low to extract,
aperture issues, etc.), it still writes out a structurally valid `_x1d.fits`
with all the expected header metadata but a **zero-row BinTable** in extension 1.
MAST archives these failed products and bulk downloads grab them. The classic
pattern: the entire failed visit gets re-observed under a new ASN_ID and both
visits are delivered. Example: `2MASS J00505073+3536429` has visit `lcxv03`
(20 empty files) followed by re-observation `lcxv26` (good data).

**(2) Same photons counted twice.** The `read_cos_flat` glob picked up both
`*_x1d.fits` (per-exposure) and `*_x1dsum.fits` (CalCOS per-visit weighted
coadd of the same exposures in one association). Including both means the
downstream `coadd.py` weights the same exposures in once via the individual
`_x1d` files and again via the `_x1dsum` coadd. Trevor's legacy readers never
referenced `_x1dsum.fits` — hand-curated SulenticAllData folders only ever
contained `_x1d.fits` per exposure.

### STIS investigation — _sx1.fits is NOT analogous to _x1dsum.fits

Before changing the STIS reader, checked the STIS raw data: 190 object folders,
424 `_sx1.fits` files across 67 objects. Two key findings refuted the initial
"same problem as COS" assumption:

- **Zero rootname collisions.** No single STIS dataset is delivered as both
  `_x1d.fits` and `_sx1.fits`. They are mutually exclusive.
- **Detector-specific, not cross-exposure sum.** `_x1d.fits` is produced for
  MAMA modes (G140L/M, G230L/M, E140M, E230M, PRISM); `_sx1.fits` is produced
  for CCD modes with CR-SPLIT (G430L/M, G750L/M, G230LB, G230MB) — the "s" is
  cosmic-ray-summed across CR-SPLIT sub-exposures, still per-exposure. They are
  complementary products covering different gratings, not duplicates.
- **67 objects have only `_sx1.fits`** (all-CCD-mode observations). Dropping
  `_sx1.fits` from the glob would silently lose all their data.

Conclusion: keep both globs in `read_stis_flat`; only add the empty-file filter
as a defensive hardening (no empty STIS file is currently observed in the data,
but the structural vulnerability is identical to COS's).

### Fix

**`read_cos_flat`:**
1. Removed `*_x1dsum.fits` from the glob (matches Trevor's legacy convention).
2. Added a filter step that drops any file where extension 1's `data` is `None`
   or has length 0, raising `FileNotFoundError` if no non-empty files remain.

**`read_stis_flat`:**
1. Kept the dual glob (`*_x1d.fits` + `*_sx1.fits`) — these are complementary
   detector-specific products, not duplicates. Added a comment explaining the
   distinction.
2. Added the same empty-file filter as a defensive precaution.

**FOS reader (`read_fos_flat`):** intentionally left alone. FOS uses a
multi-file format (`_c0f` / `_c1f` / `_c2f` / `_cqf`) and has no per-visit
summed product analog. Checked 182 FOS object folders: 0 empty `_c0f.fits`
files. No structural vulnerability triggered by current data.

### Verification

Smoke test of the two originally failing objects after the fix:
- `2MASS J00505073+3536429` — reads cleanly with shape `(16, 16384)`.
- `2MASS J00435499+4234302` — clears the `IndexError`, but progresses to a
  *different* downstream error in `Cut_Edge_Pix` (`boolean index did not match
  indexed array along axis 0; size of axis is 16384 but size of corresponding
  boolean axis is 1274`). This is a separate, pre-existing issue unrelated to
  the empty-file IndexError and is NOT addressed by this change. Logged here
  for future investigation.

---

## TEMP: COS exposure cap for Mrk 817 & NGC 5548 — REVISIT

**Date:** 2026-06-17
**Files modified:** `rebinning/read_spec_data.py` (`read_cos_flat` + new
`_COS_EXPOSURE_CAP` constant and `_cos_median_snr` helper)
**Status:** ⚠️ **STOPGAP — this decision must be revisited before final results.**

### Symptom

`python -m rebinning.run_rebin --master-catalog` appeared to "hang" immediately
after the `Mrk 771` objects. It was not hung: the next object, **Mrk 817 (COS)**,
has **1,598 individual `_x1d.fits` exposures**, and **NGC 5548 (COS)** has **795**.
Every other COS object in the sample has ≤ 73. Both are intensively monitored
reverberation-mapping AGN (e.g. AGN STORM / STORM 2), so MAST holds hundreds–
thousands of individual COS exposures for them.

Loading all of them stacks a `(Nexp, Npix)` array and runs the pure-Python
per-exposure morph loop and the per-pixel weighted-median double loop in
`coadd.py` (lines ~161 and ~209-211) over ~800–1,600 exposures — minutes-to-hours
of work that, combined with un-flushed `print()` output in a non-TTY console,
looks like a freeze.

### Stopgap applied

For **Mrk 817** and **NGC 5548 only** (hardcoded by `common_name` in
`_COS_EXPOSURE_CAP = {"Mrk 817": 100, "NGC 5548": 100}`), `read_cos_flat` now
keeps only the **100 highest-S/N `_x1d` exposures** before the heavy per-file
work. Ranking metric: median per-pixel `flux/error` across each exposure
(`_cos_median_snr`, flattens all segments/stripes). The cap block sits right
after the existing empty-file filter, so empties never count toward the 100.
All other objects hit `name not in _COS_EXPOSURE_CAP` → **zero behavior change
and zero overhead** for the rest of the sample.

### Why this is NOT the real fix — REVISIT

- **Arbitrary threshold.** "100 highest-S/N exposures" is a runtime convenience,
  not a science-justified selection. We are discarding the majority of the
  available photons for these two objects. The final coadd S/N and any CIV
  measurements for Mrk 817 / NGC 5548 depend on this arbitrary cut and must not
  be treated as final until the decision is revisited (ideally with G. Richards).
- **Hardcoded object names** will silently fail to protect any *future*
  monitoring target added to the catalog.

### Proper fixes to evaluate when revisiting

1. **Switch COS to `_x1dsum`-only** (per-visit CalCOS coadds). Investigated
   2026-06-17: all **244/244** COS objects have `_x1dsum` files, so no object
   would be dropped; Mrk 817 drops 1,598 → 613, NGC 5548 similarly. The original
   reason `_x1dsum` was excluded was **double-counting when globbed *alongside*
   `_x1d`** (see "Same photons counted twice", this log, ~line 503) — using
   `_x1dsum` *exclusively* resolves that. This is a sample-wide methodology
   change from Trevor's per-exposure convention and needs sign-off.
2. **Vectorize the `coadd.py` per-pixel weighted-median loop** so exposure count
   stops being the bottleneck, removing the need for any cap.

---

## Phase 5.8 — ICA `--master` mode (run ICA over all rebinned spectra)

**Date:** 2026-06-18
**Files modified:** `ica/run_all_objects.py`, `ica/manual_fix.py`

### Motivation

`rebinning/run_rebin.py --master-catalog` produces `raw_data/RebinnedSpec_master/`
(576 FITS files: 229 COS, 167 FOS, 180 STIS), named `{common_name}_{inst}.fits`.
The ICA batch runner, however, was **catalog-driven**: `run_all_objects.py` iterated
the 207 rows of `HST_CIV_Sulentic2007_HSLA2018_finalprops.csv` and
`manual_fix.setup_object` looked each up by `Spec_Name`/`Inst_final`
(`glob("{rebin}/{Spec_Name}*{Inst}.fits")`). Pointing `--rebin-dir` at the master
folder therefore only matched **28 of 576** spectra (the few FOS/STIS objects whose
old Sulentic Spec_Name/Inst coincide with the new filenames) and processed **zero
COS** objects, because COS is absent from the Sulentic catalog entirely.

Chosen approach: **Option A** — add a `master_mode` flag to `ICAManualFixProcessor`
that switches only the object-selection / file-resolution / plot-folder layer to be
folder-driven, leaving the ICA science path (`run_ica.main_ICA` and everything it
calls) byte-for-byte identical. (Option B, a standalone runner, was rejected to keep
a single code path and mirror how `run_rebin.py` got its `--master-catalog` flag.)

### Changes — `ica/manual_fix.py`

| Item | Description |
|------|-------------|
| `ICAManualFixProcessor.__init__` | Added `master_mode=False` and `output_path=None` params. `self.master_mode` stored; `self.output_path` defaults to module `OUTPUT_PATH`. Sulentic catalog + `quality_classification.json` are still loaded in master mode (cheap, unused) to avoid branching the constructor. |
| New class attr `_MASTER_FOLDER = "All"` | Single plot subfolder used for every master object (quality classification only exists for the 207 Sulentic objects). |
| New method `list_master_objects()` | Returns sorted FITS filename stems (`{common_name}_{inst}`) of every `*.fits` in `rebin_path`. One entry per spectrum → multi-instrument objects fit independently. |
| New method `_ensure_output_dirs()` | Creates the plot dir + subfolders for the active mode. Replaces the two duplicated `makedirs` blocks (previously in `process_object` and `batch_process`), which hardcoded the four quality folders. |
| `setup_object()` | Refactored so **only file resolution differs** by mode: master mode sets `spec_name = name` (the stem) and `fn = {rebin_path}/{name}.fits` (direct path — no glob, so spaces in COS names are safe); Sulentic mode keeps the `Spec_Name`/`Inst_final` lookup + glob. Redshift is read from the FITS in both modes (`spec[1].data["Redshift"][0]`), so master mode needs no redshift catalog. The shared body (flux/error normalisation, morphing, BAL masking, return tuple) is unchanged. |
| `get_folder()` | Master mode returns `_MASTER_FOLDER` instead of the JSON lookup (master objects are unclassified, so the JSON would raise). |
| `process_object` / `batch_process` / `create_diagnostic_plot` | `OUTPUT_PATH` references → `self.output_path`; duplicated `makedirs` → `self._ensure_output_dirs()`. |

### Changes — `ica/run_all_objects.py`

| Item | Description |
|------|-------------|
| New `--master` flag | When set: rebin dir defaults to `raw_data/RebinnedSpec_master` (overridable with `--rebin-dir`); constructs the processor with `master_mode=True, output_path="ICA_Plots_Rebin_master"`; iterates `processor.list_master_objects()`; writes `ICA_Results_Rebin_master.csv`. Separate output filenames so a master run never clobbers the Sulentic run's `ICA_Results_Rebin_AP.csv` / `ICA_Plots_Rebin_AP/`. |
| Default (no `--master`) path | Unchanged. |
| Docstring / header | Documented both modes; header bumped to "Updated (Phase 5.8)". |

### Behaviour of Sulentic-keyed overrides in master mode

`MorphingEdges`, `CIV_BAL_regions`, and `MANUAL_FIX_CONFIG` are all keyed on old
Sulentic `Spec_Name`s / `Final_Name`s. In master mode the object key is the new
filename stem, so these dicts simply **do not match** new objects and the overrides
are skipped (no morphing, no manual BAL windows, no forced components, no custom
pixel masks). This is the intended first-pass behaviour; per-object overrides can be
re-added later keyed on the master stems. For the ~28 FOS/STIS objects whose stem
happens to equal an old Spec_Name, the override **would** apply — acceptable and
harmless, but worth noting if exact parity with a Sulentic run is ever needed.

### Pre-existing bug fixed in passing — `ica_path=""`

`process_object` called `main_ICA(..., ica_path="", ...)`. The Phase 2 migration
changed `run_ica.load_ICA` (and the other `ica_path` functions) to fall back to
`Data/components/` **only when `ica_path is None`**; an empty string slips past that
guard and makes `pd.read_csv("wav_12603000.dat")` look in the cwd, raising
`FileNotFoundError` unless run from a directory containing the component files. This
broke **both** modes when run from the repo root; it surfaced here because the master
smoke test was the first run from repo root. Fixed by passing `ica_path=None` so the
canonical `Data/components/` default applies. No algorithmic change; strictly makes
the documented `python -m ica.run_all_objects` invocation work from the repo root.

### Verification

- `list_master_objects()` → 576 spectra.
- End-to-end `process_object(..., save_plot=True)` on one object per instrument
  (`1RXS J185800.9+485020_COS`, `0115+027_FOS`, `0107-0019_STIS`) all succeed; plots
  written to `ICA_Plots_Rebin_master/All/`. (Test plots removed afterward.)
- Default Sulentic-mode construction unaffected: `master_mode=False`,
  `output_path="ICA_Plots_Rebin_AP"`, `get_folder("J04232-0120") → "Good"`.

### Usage

```bash
# Run ICA over ALL rebinned spectra in raw_data/RebinnedSpec_master:
python -m ica.run_all_objects --master

# Master mode with an explicit rebin directory:
python -m ica.run_all_objects --master --rebin-dir /path/to/RebinnedSpec_master
```
