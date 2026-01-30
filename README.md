# HST Quasar Analysis Paper - Figure Generation

This repository contains scripts to generate publication-quality figures for the HST quasar spectroscopy paper.

## Repository Structure

```
HST Paper/
├── Data/                    # Data files required for plotting
├── Figures/                 # Output directory for generated figures
├── Notebooks/              # Original analysis notebooks
├── Plotting Code/          # Scripts to generate figures
│   └── make_fig1.py        # Figure 1: L vs z
└── README.md               # This file
```

## Requirements

### Python Dependencies

Install all required packages using the provided requirements file:

```bash
pip install -r requirements.txt
```

Or install packages individually:

```bash
pip install numpy pandas scipy matplotlib astropy palettable richardsplot
```

Or using conda:

```bash
conda install numpy pandas scipy matplotlib astropy
pip install palettable richardsplot
```

### Required Packages:
- `numpy` - Numerical computing
- `pandas` - Data manipulation
- `scipy` - Scientific computing (statistics, interpolation)
- `matplotlib` - Plotting
- `astropy` - Cosmology calculations and FITS file handling
- `palettable` - Color palettes
- `richardsplot` - Custom plotting configurations

## Usage

### Figure 1: Luminosity versus Redshift

Generate Figure 1 showing the relationship between L_2500 (luminosity at 2500 Angstroms) and redshift for HST, SDSS, and GNIRS-DQS quasar samples.

**Basic usage:**
```bash
cd "Plotting Code"
python make_fig1.py
```

This generates the figure using the Rankine+20 SDSS data and saves it to `../Figures/Fig1_L_versus_z.pdf`.

**Using Temple+21 SDSS data:**
```bash
cd "Plotting Code"
python make_fig1.py --data-version temple
```

This generates the figure using Temple+21 SDSS data and saves it to `../Figures/Fig1_L_versus_z_Temple.pdf`.

**Custom output directory:**
```bash
cd "Plotting Code"
python make_fig1.py --output-dir ./my_figures/
```

**Command-line options:**
- `--data-version`: Choose SDSS data source
  - `rankine` (default): Use Rankine+20 data
  - `temple`: Use Temple+21 data
- `--output-dir`: Specify output directory (default: `../Figures/`)

**Help:**
```bash
python make_fig1.py --help
```

## Output

Each script generates a PDF figure in the specified output directory with publication-quality resolution (300 DPI). The figures are created with:
- Proper axis labels and font sizes
- Legend with sample identification
- Secondary y-axis showing bolometric luminosity
- Tight bounding box for easy inclusion in papers

## Contributing

When adding new figures:
1. Create a new script `make_figX.py` in the `Plotting Code/` folder following the template of `make_fig1.py`
2. Include proper documentation and command-line arguments
3. Update this README with usage instructions
4. Test the script to ensure it runs correctly

## Contact

For questions or issues, please contact the paper authors or open an issue in the repository.

## Citation

If you use this code or data, please cite the associated paper:

[Citation information to be added upon publication]
