# Quick Start Guide

## First Time Setup

1. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd "HST Paper"
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

## Generate Figures

All commands should be run from the `Plotting Code/` directory:

```bash
cd "Plotting Code"
```

### Figure 1: L vs z (Rankine+20 SDSS data)
```bash
python make_fig1.py
```
Output: `../Figures/Fig1_L_versus_z.pdf`

### Figure 1: L vs z (Temple+21 SDSS data)
```bash
python make_fig1.py --data-version temple
```
Output: `../Figures/Fig1_L_versus_z_Temple.pdf`

### Figure 2: EV1 Optical Diagram
```bash
python make_fig2.py
```
Output: `../Figures/Fig2_EV1_optical.pdf`

## Need Help?

Run any script with `--help` for usage information:
```bash
python make_fig1.py --help
python make_fig2.py --help
```
