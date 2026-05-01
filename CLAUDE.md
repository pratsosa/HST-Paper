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
