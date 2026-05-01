# HST Quasar CIV Paper — Analysis Pipeline

This repository contains the full, reproducible analysis pipeline for the paper on the
relationship between quasar luminosity and CIV emission line properties. It covers
everything from raw HST spectra through rebinning, ICA fitting, and publication figures.

**Authors:** Alexandros Pratsos, with foundational code from Trevor McCaffrey
**Advisor:** Prof. Gordon Richards, Drexel University

---

## Repository Structure

```
HST Paper/
├── rebinning/          # Rebinning and coadding pipeline (HST + SDSS)
├── ica/                # ICA spectral fitting pipeline
├── Plotting Code/      # Publication figure scripts
├── Data/               # Catalogs, component files, reference spectra
├── Figures/            # Output PDFs/PNGs
├── Notebooks/          # Exploratory analysis notebooks
├── pyproject.toml      # Installable package
└── requirements.txt
```

---

## Installation

Clone the repo and install in editable mode:

```bash
git clone <repository-url>
cd "HST Paper"
pip install -e .
```

This installs the `hst_civ_paper` package, making `rebinning` and `ica` importable
from anywhere in your environment.

---

## Data Setup

The pipeline reads raw spectral data from an external directory (not stored in this
repo). Set the location with the `HST_PAPER_DATA_DIR` environment variable, or pass
`--data-dir` on the command line:

```bash
export HST_PAPER_DATA_DIR=/path/to/your/data
```

The data directory must contain:
- `SulenticAllData/FOS/` — raw FOS FITS files
- `SulenticAllData/STIS/` — raw STIS FITS files
- `HSLA_coadds_wCIV/original/` — HSLA coadd FITS files
- `SDSS_spec/lite/` — SDSS spectra (for SDSS–HST coadds)

The file `Data/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv` (committed to the repo)
is the master object catalog and is used by default.

---

## Pipeline

### Step 1 — Rebinning

Rebin and coadd raw HST (and SDSS) spectra to a standard 69 km/s log-space grid.
Outputs one FITS file per object into a `RebinnedSpec/` directory.

```bash
# All instruments (FOS, STIS, HSLA):
python -m rebinning.run_rebin --data-dir $HST_PAPER_DATA_DIR --output-dir RebinnedSpec/

# FOS + STIS only:
python -m rebinning.run_rebin --instruments fos_stis --data-dir $HST_PAPER_DATA_DIR

# HSLA only:
python -m rebinning.run_rebin --instruments hsla --data-dir $HST_PAPER_DATA_DIR
```

Key options:
- `--data-dir` — root directory containing raw data (default: `$HST_PAPER_DATA_DIR`)
- `--output-dir` — where to write rebinned FITS files (default: `RebinnedSpec/`)
- `--sdss-spec-dir` — path to `SDSS_spec/lite/` if not under `--data-dir`
- `--catalog` — path to object catalog CSV (default: `Data/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv`)
- `--instruments` — `fos_stis`, `hsla`, or `all` (default: `all`)

### Step 2 — ICA Fitting

Run ICA spectral decomposition on the rebinned spectra. Outputs component fits and
diagnostic plots per object.

```bash
python -m ica.run_all_objects --rebinned-dir RebinnedSpec/ --output-dir ICA_output/
```

To process a single object:

```bash
python -m ica.run_all_objects --object 3C273 --rebinned-dir RebinnedSpec/ --output-dir ICA_output/
```

### Step 3 — Publication Figures

All figure scripts live in `Plotting Code/` and are run from that directory.

```bash
cd "Plotting Code"
```

**Figure 1: Luminosity vs. Redshift**

```bash
python make_fig1.py
# Output: ../Figures/Fig1_L_versus_z.pdf

python make_fig1.py --data-version temple
# Output: ../Figures/Fig1_L_versus_z_Temple.pdf
```

Options: `--data-version rankine|temple`, `--output-dir`

**Figure 2: Optical EV1 Diagram (R(FeII) vs. H-beta FWHM)**

```bash
python make_fig2.py
# Output: ../Figures/Fig2_EV1_optical.pdf
```

Options: `--output-dir`

Run any script with `--help` for the full option list.

---

## Key Dependencies

| Package | Purpose |
|---------|---------|
| `numpy`, `scipy` | Numerical computing |
| `pandas` | Catalog handling |
| `matplotlib` | Plotting |
| `astropy` | FITS I/O, cosmology |
| `lmfit` | Robust spectral fitting (Rankine+2020) |
| `weightedstats` | Weighted statistics |
| `palettable` | Color palettes |
| `richardsplot` | Richards group plot styles |

Install everything at once:

```bash
pip install -r requirements.txt
```

---

## Scientific Context

- **Sample:** ~207 HST quasar spectra (FOS, STIS, HSLA instruments)
- **Rebinning resolution:** 69 km/s log-space (FWHM ≈ 162 km/s)
- **ICA components:** Three component sets covering 1260–3000 Å (mod / low-EW / high-EW)
- **Reference rebin run:** `Trevor Code/RebinnedSpec_2022Aug11/` (Aug 2022, T. McCaffrey)
- **Validated rebin run:** `Trevor Code/RebinnedSpec_2026AP/` (AP 2026, identical algorithm)
- **Emission line focus:** CIV λ1549 Å properties vs. luminosity at 2500 Å

---

## Contributing

When adding a new figure:
1. Create `make_figX.py` in `Plotting Code/`
2. Import shared functions from `plotting_utils.py` where applicable
3. Add `argparse` CLI arguments following the pattern of existing scripts
4. Update this README with a usage block

When modifying pipeline code (`rebinning/` or `ica/`):
- No functional changes to algorithm logic without explicit discussion
- Log all changes in `.claude/Migration_Log.md`
- Keep paths relative and `DATA_DIR`-configurable — no hardcoded absolute paths
- See `CLAUDE.md` for the full set of conventions

---

## Project Documentation

- [CLAUDE.md](CLAUDE.md) — codebase guide for Claude (structure, principles, data flow)
- [.claude/Migration_Plan.md](.claude/Migration_Plan.md) — phased migration plan
- [.claude/Migration_Log.md](.claude/Migration_Log.md) — running log of all changes
