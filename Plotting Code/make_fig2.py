"""
Figure 2: Optical Emission Line Properties (EV1 Diagram)
=========================================================

This script generates Figure 2 for the paper, showing the relationship between
R(FeII) (ratio of FeII to H-beta equivalent widths) and H-beta FWHM for various 
quasar samples. This is the Eigenvector 1 (EV1) optical diagram.

Usage:
    python make_fig2.py [--output-dir DIR]

Arguments:
    --output-dir: Directory to save the figure (default: ../Figures/)

Examples:
    python make_fig2.py
    python make_fig2.py --output-dir ./output/
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import palettable
import richardsplot as rplot
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

# Import shared plotting utilities
from plotting_utils import plot_contour

# Initialize cosmology
cosmo = FlatLambdaCDM(Om0=0.3, H0=70)

# Configure matplotlib parameters
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True


def load_hst_data(data_dir):
    """
    Load HST quasar data with H-beta FWHM and FeII measurements.
    
    Returns:
    --------
    dict : Dictionary containing HST data arrays
    """
    hst_file = os.path.join(data_dir, "HST_CIV_Sulentic2007_HSLA2018_1Jan2023.csv")
    hst = pd.read_csv(hst_file)
    
    # Filter for objects with H-beta FWHM measurements
    hst = hst[~np.isnan(hst["Hbeta_FWHM_sul07"])]
    
    Hb_fwhm_hst = hst["Hbeta_FWHM_sul07"].values
    RFeII_hst = hst["RFeII_sul07"].values
    
    # Create mask for objects with BH mass measurements
    maskMBH = (~np.isnan(hst["log M_BH"].values))
    
    return {
        'RFeII': RFeII_hst,
        'Hb_fwhm': Hb_fwhm_hst,
        'mask_mbh': maskMBH,
        'mask_no_mbh': ~maskMBH
    }


def load_sdss_data(data_dir):
    """
    Load SDSS data from Rakshit et al. 2020.
    
    Returns:
    --------
    dict : Dictionary containing SDSS data arrays
    """
    sdss_file = os.path.join(data_dir, "Rakshit2020.fit")
    sdss = Table.read(sdss_file).to_pandas()
    
    # Filter for objects with complete FeII and H-beta measurements
    mask = ((~np.isnan(sdss["EWFe-OP"])) & 
            (~np.isnan(sdss["EWHb-BR"])) & 
            (~np.isnan(sdss["FWHMHb-BR"])))
    sdss = sdss[mask]
    
    EW_feii_sdss = sdss["EWFe-OP"].values
    EW_Hb_sdss = sdss["EWHb-BR"].values
    RFeII_sdss = EW_feii_sdss / EW_Hb_sdss
    Hb_fwhm_sdss = sdss["FWHMHb-BR"].values
    
    return {
        'RFeII': RFeII_sdss,
        'Hb_fwhm': Hb_fwhm_sdss
    }


def load_gnirs_data(data_dir):
    """
    Load GNIRS-DQS data.
    
    Returns:
    --------
    dict : Dictionary containing GNIRS data arrays
    """
    gnirs_file = os.path.join(data_dir, "gnirsdqs_wRankine_wdr16Lum.csv")
    gnirs = pd.read_csv(gnirs_file)
    
    EW_feii_gnirs = gnirs["EW-Fe2"].values
    EW_Hb_gnirs = gnirs["EW-Hb"].values
    RFeII_gnirs = EW_feii_gnirs / EW_Hb_gnirs
    Hb_fwhm_gnirs = gnirs["FWHM-Hb"].values
    
    return {
        'RFeII': RFeII_gnirs,
        'Hb_fwhm': Hb_fwhm_gnirs
    }


def create_figure(hst_data, sdss_data, gnirs_data):
    """
    Create the EV1 optical diagram figure.
    
    Parameters:
    -----------
    hst_data : dict
        HST data dictionary
    sdss_data : dict
        SDSS data dictionary
    gnirs_data : dict
        GNIRS data dictionary
    
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    # Set up color palette
    cs = palettable.colorbrewer.qualitative.Dark2_5.mpl_colors
    
    fig, ax = plt.subplots(1, 1, figsize=(9, 9))
    
    # Plot HST sample
    ax.plot(hst_data['RFeII'][hst_data['mask_no_mbh']], 
            hst_data['Hb_fwhm'][hst_data['mask_no_mbh']], 
            marker="^", markersize=12, markerfacecolor=cs[4], 
            markeredgecolor="k", label="HST, no $M_{BH}$", 
            linestyle="", zorder=3)
    
    ax.plot(hst_data['RFeII'][hst_data['mask_mbh']], 
            hst_data['Hb_fwhm'][hst_data['mask_mbh']], 
            marker="o", markersize=15, markerfacecolor=cs[3], 
            markeredgecolor="k", label="HST, yes $M_{BH}$", 
            linestyle="", zorder=3)
    
    # Plot SDSS sample with density contours
    print("Plotting SDSS data...")
    plot_contour(sdss_data['RFeII'], sdss_data['Hb_fwhm'], 
                 c=cs[2], mark="o", nlevels=[0.05, 0.25, 0.50, 0.75, 0.95], 
                 ax=ax, linewidths=1.5, s=3, alpha=0.4, label="SDSS")
    
    # Plot GNIRS-DQS
    ax.plot(gnirs_data['RFeII'], gnirs_data['Hb_fwhm'], 
            marker="v", markersize=11, markerfacecolor=cs[2], 
            markeredgecolor="k", label="GNIRS-DQS", 
            linestyle="", zorder=2)
    
    # Configure plot aesthetics
    ax.set_xlabel('R(FeII)', fontsize=35)
    ax.set_ylabel(r'FWHM\;H\,$\beta$', fontsize=35)
    ax.tick_params(axis='both', which='major', labelsize=27.5)
    ax.legend(frameon=True, labelspacing=0.1, loc="upper right", 
              prop={"size": 25})
    
    # Set axis limits
    ax.set_ylim(0, 20000)
    ax.set_xlim(-0.3, 5)
    
    plt.tight_layout()
    
    return fig, ax


def main():
    """Main function to generate Figure 2."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Generate Figure 2: Optical Emission Line Properties',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--output-dir',
                        default='../Figures/',
                        help='Directory to save the figure')
    
    args = parser.parse_args()
    
    # Set up paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)  # Go up to HST Paper directory
    data_dir = os.path.join(parent_dir, 'Data')
    
    # Handle relative output directory path
    if args.output_dir.startswith('..'):
        output_dir = os.path.join(parent_dir, args.output_dir)
    else:
        output_dir = os.path.join(parent_dir, args.output_dir)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 70)
    print("Figure 2: Optical Emission Line Properties (EV1 Diagram)")
    print("=" * 70)
    print(f"Output directory: {output_dir}")
    print()
    
    # Load data
    print("Loading data...")
    hst_data = load_hst_data(data_dir)
    print(f"  HST: {len(hst_data['RFeII'])} objects")
    
    sdss_data = load_sdss_data(data_dir)
    print(f"  SDSS: {len(sdss_data['RFeII'])} objects")
    
    gnirs_data = load_gnirs_data(data_dir)
    print(f"  GNIRS-DQS: {len(gnirs_data['RFeII'])} objects")
    print()
    
    # Create figure
    print("Creating figure...")
    fig, ax = create_figure(hst_data, sdss_data, gnirs_data)
    
    # Save figure
    output_filename = "Fig2_EV1_optical.pdf"
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved: {output_path}")
    print()
    print("=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
