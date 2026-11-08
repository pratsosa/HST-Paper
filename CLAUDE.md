# HST CIV Paper — Codebase Guide for Claude

## Overview

This repository contains the **clean, reproducible analysis pipeline** for the paper on the
relationship between quasar luminosity and CIV emission line properties. It is the successor
to `Richards/Trevor Code/` and follows a phased migration plan documented in
`.claude/Migration_Plan.md`.

**Authors:** Alexandros Pratsos (AP), with foundational code from Trevor McCaffrey (TVM)
**Advisor:** Prof. Gordon Richards, Drexel University

---

## Repository Structure

```
HST Paper/
├── .claude/
│   ├── Migration_Plan.md        # Phase-by-phase migration plan
│   └── Migration_Log.md         # Running log of all changes
├── CLAUDE.md                    # This file
├── pyproject.toml               # Installable package (pip install -e .)
├── requirements.txt
│
├── rebinning/                   # Phase 1 — rebinning + coadding pipeline
│   ├── __init__.py
│   ├── coadd.py                 # Main coadd entry point
│   ├── read_spec_data.py        # Instrument-specific readers (FOS, STIS, HSLA)
│   ├── cut_edge_pix.py          # Edge pixel cutting
│   ├── lower_res_rebin.py       # Rebin to lower resolution
│   ├── small_pix_filter.py      # SDSS small pixel filter
│   ├── spec_cuts_fos.py         # FOS spectral cuts
│   ├── spec_cuts_hsla.py        # HSLA spectral cuts
│   └── run_rebin.py             # Unified runner (all instruments)
│
├── ica/                         # Phase 2 — ICA analysis pipeline
│   ├── __init__.py
│   ├── run_ica.py               # Main ICA fitting pipeline
│   ├── fit_composites.py        # Component fitting
│   ├── plot_ica.py              # ICA diagnostic plots
│   ├── spec_morph.py            # Spectral morphing
│   ├── bal_regions.py           # BAL region definitions
│   ├── civ_bal_regions.py       # CIV BAL regions
│   ├── morphing_edges.py        # Merged morphing edges (all instruments)
│   ├── manual_fix.py            # Manual override logic
│   ├── manual_fix_config.py     # Manual override configuration
│   ├── quality_classification.json
│   └── run_all_objects.py       # Batch runner
│
├── Data/
│   ├── components/              # ICA component files (.comp, .dat, .npy)
│   └── ...                      # Catalog CSVs and other data
│
├── diagnostics/                 # Pipeline audits (not part of the pipeline itself)
│   ├── audit_civ_coverage.py    # Census: raw vs rebinned C IV coverage, all objects
│   ├── trace_exposures.py       # Per-exposure trace for one object
│   ├── rebin_patched.py         # Candidate fixes applied by monkey-patching
│   ├── compare_rebins.py        # Compare two rebin output directories
│   └── CIV_Data_Loss_Investigation.md
│
├── Figures/                     # Output PDFs/PNGs (not committed by default)
├── Notebooks/                   # Exploratory notebooks
└── Plotting Code/               # Publication figure scripts
```

---

## Key Principles

1. **No functional changes to algorithms** — all science logic is identical to
   `Trevor Code/`. Only structure, imports, and paths change during migration.
2. **File provenance** — every migrated file starts with a header comment naming
   the original `Trevor Code/` source file and copy date.
3. **Portable paths** — no hardcoded `/Users/Trevor1/` or `C:\Users\...` paths
   anywhere. External data location is specified via a single `DATA_DIR` variable
   or CLI argument.
4. **Versioned files** — new versions get a suffix (e.g., `_2026E`, `_2026AP`)
   rather than overwriting working files.
5. **HST/ folder is archival** — never modify files in `Richards/HST/`.
6. **Reproducibility** — every step must be reproducible from raw data alone.

---

## Installation

```bash
pip install -e .
```

This installs the package in editable mode. The package is named `hst_civ_paper`.

---

## Data Flow

```
Raw FITS (SulenticAllData/ or MAST download)
    └── rebinning/run_rebin.py
            └── RebinnedSpec output (one FITS per object)
                    └── ica/run_all_objects.py
                            └── ICA fits, plots, CIV measurements
                                    └── Plotting Code/ → publication figures
```

---

## Configuring the Data Directory

All pipeline scripts accept a `DATA_DIR` argument pointing to the root directory
containing the raw spectral data. This defaults to the environment variable
`HST_PAPER_DATA_DIR` if set.

```bash
export HST_PAPER_DATA_DIR=/path/to/your/data
python -m rebinning.run_rebin --data-dir $HST_PAPER_DATA_DIR
```

---

## Known Issues

> **Start here:** `PIPELINE_CHANGES_2026-09.md` is the self-contained summary of the
> September 2026 fixes — which output files to use, what changed and why, what is still
> outstanding, and the caveats that matter for downstream figures. The sections below
> and the `diagnostics/*.md` logs are the supporting detail.


### C IV data loss — fixed 2026-09-10, spectra in `RebinnedSpec_master_v2/`

Four bugs in `rebinning/` cost 57 spectra their C IV coverage, three of them silently.
All four are now fixed in `rebinning/read_spec_data.py` (flat readers) and
`rebinning/spec_morph.py` (`cont_filtered`). Diagnosis and evidence:
`diagnostics/CIV_Data_Loss_Investigation.md`; change log: `.claude/Migration_Log.md`.

**`RebinnedSpec_master/` and `ICA_Results_Rebin_master.csv` are the OLD, affected
results** — kept for comparison. Corrected spectra are in `RebinnedSpec_master_v2/`.
The biggest fix (reading all COS detector segments, not just row 0) changes 91% of COS
spectra, so all COS ICA fits must be re-run against the v2 directory.

Still open: the `_COS_EXPOSURE_CAP` stopgap for Mrk 817 and NGC 5548 (AP to raise with
G. Richards), and re-review of the 36 COS manual ICA overrides.

### Two measurement defects in the ICA stage — FIXED 2026-09-23

Both predated this work and affected Trevor's published values as well as ours.

1. **C IV equivalent widths were biased low.** `get_CIV_parameters` read the continuum
   ~240 A too blue (`continuum[i]` while the flux used `flux_r[EW][i]`). Fixed by
   masking once into `cont_EW = continuum[EW]`, in `get_CIV_parameters` and
   `plot_CIV_analysis` (`ica/manual_fix.py`) and `get_CIV` (`ica/run_ica.py`).
   EW rises by a median factor of 1.495; blueshift changes for ~66% of objects
   (median 137 km/s). See `diagnostics/CIV_EW_Indexing_Investigation.md`.

2. **47 fits collapsed to an all-zero reconstruction** (`CIV_blueshift = 9574.1`).
   Cause: pixels with near-zero error give `ivar = 1/err**2` ~ 1e29 and pin the fit at
   zero weights. Fixed by a guard in `ICA_fit` (`ica/run_ica.py`) dropping non-finite
   ivar and ivar above `_IVAR_MAX_FACTOR` (1e4) x median.
   See `diagnostics/Collapsed_Fits_Investigation.md`.

**Current best results: `ICA_Results_Rebin_master_v2_fixed.csv` /
`ICA_Plots_Rebin_master_v2_fixed/`** (rebin v2 + both measurement fixes).
`ICA_Results_Rebin_master.csv` (original) and `..._v2.csv` (rebin fixed only) are kept
for comparison. Degenerate rows went 47 -> 6 with no regressions.

### REGRESSION found 2026-09-29 — 36 v2 spectra have an un-normalised block

`cont_filtered()` in `rebinning/spec_morph.py` fell back to a continuum of **1.0** when
an exposure's continuum was NaN everywhere (the `else 1.0` added 2026-09-10). Those
exposures were never divided by a continuum, so a contiguous block of the co-add sits in
raw physical units (~1e-14) while the rest is normalised (~1) — a factor ~1e13.

- **36 of 594** spectra in `RebinnedSpec_master_v2/` are affected (0 of 576 before).
- It is also the origin of most of the collapsed ICA fits: 30 of the 47 trace to it;
  only 17 are genuine.
- **Fixed** — `cont_filtered()` now uses `_flat_continuum_level()`, the exposure's own
  flux scale, never 1.0.
- **The 36 cannot be regenerated until the raw data on F: is reattached.** Until then,
  exclude the objects in `diagnostics/affected_unnormalised.csv` from
  `ICA_Results_Rebin_master_v2_fixed.csv`.

Also still to do: a **quality cut independent of fit success** — some objects return a
number from a spectrum with no data in the C IV window. All 6 remaining `EW = 0` rows
(incl. `PG0953+415_STIS`, `LEDA 138607_COS`, `Ton 488_COS`) have **zero** good pixels in
1500-1600 A and should simply be dropped.

### Two diverged versions of `spec_morph.py`

There are **two separate `spec_morph.py` files** in this repo, one per pipeline:

| File | Used by | Source |
|------|---------|--------|
| `rebinning/spec_morph.py` | Rebinning pipeline (`coadd.py`) | `Trevor Code/spec_morph.py` |
| `ica/spec_morph.py` | ICA pipeline (`run_ica.py`) | `Trevor Code/ICA Scripts/spec_morph.py` |

They diverged from a common base. Key differences: the rebinning version has AP's
defensive `cont_filtered()` fixes; the ICA version uses a different `continuum_fit()`
signature and `np.polynomial` fitting. **If you fix a bug in one, check whether the
other needs the same fix.** See `Migration_Log.md` → "Known Issues" for full details
and future fix options.

---

## Critical Dependencies

```python
numpy, pandas, scipy, matplotlib
astropy          # FITS I/O, cosmology
lmfit            # Robust spectral fitting (Rankine+2020)
weightedstats    # Weighted statistics utilities
palettable       # Color palettes
richardsplot     # Custom plot styles (Richards group)
```

---

## Migration Status

See `.claude/Migration_Plan.md` for the full plan and `.claude/Migration_Log.md`
for the running log of completed work.

| Phase | Description              | Status        |
|-------|--------------------------|---------------|
| 0     | Project scaffolding      | Complete      |
| 1     | Rebinning pipeline       | Complete      |
| 2     | ICA pipeline             | Complete      |
| 3     | Integration testing      | Complete      |
| 4     | Documentation            | Complete      |
| 5     | Connect to master pipeline | Future      |

---

## Scientific Context

- **Instrument coverage:** FOS, STIS, HSLA (HST); SDSS (optical complement)
- **Rebinning target:** 69 km/s log-space resolution (FWHM ≈ 162 km/s)
- **ICA components:** Three sets (mod/low/high EW), wavelength range 1260–3000 Å
- **Reference rebin run:** `Trevor Code/RebinnedSpec_2022Aug11/` (207 objects, Aug 2022)
- **Validated rebin run:** `Trevor Code/RebinnedSpec_2026AP/` (identical algorithm, AP 2026)
