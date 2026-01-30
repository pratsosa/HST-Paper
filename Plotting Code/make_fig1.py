"""
Figure 1: Luminosity versus Redshift Plot
==========================================

This script generates Figure 1 for the paper, showing the relationship between
L_2500 (luminosity at 2500 Angstroms) and redshift for various quasar samples.

Usage:
    python make_fig1.py [--data-version VERSION] [--output-dir DIR]

Arguments:
    --data-version: Choose 'rankine' or 'temple' for SDSS data source (default: rankine)
    --output-dir: Directory to save the figure (default: ./Figures/)

Examples:
    python make_fig1.py
    python make_fig1.py --data-version temple
    python make_fig1.py --data-version rankine --output-dir ./output/
"""

import os
import argparse
import pandas as pd
import numpy as np
from scipy import stats
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import palettable
import richardsplot as rplot
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

# Initialize cosmology
cosmo = FlatLambdaCDM(Om0=0.3, H0=70)

# Configure matplotlib parameters
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True


def plot_contour_fast2(xdata, ydata, c="k", mark=".", nlevels=3, ax=None,
                       linewidths=0.5, s=3, alpha=1, label=None):
    """
    Plot density contours with scatter points for low-density regions.
    
    Parameters:
    -----------
    xdata, ydata : array-like
        Data coordinates
    c : str or color
        Color for contours and points
    mark : str
        Marker style for scatter points
    nlevels : int or list
        Number or list of density levels for contours
    ax : matplotlib axis
        Axis to plot on (default: current axis)
    linewidths : float
        Width of contour lines
    s : float
        Size of scatter points
    alpha : float
        Transparency of scatter points
    label : str
        Label for legend
    """
    if ax is None:
        ax = plt.gca()

    dx = 0.1 * (xdata.max() - xdata.min())
    dy = 0.1 * (ydata.max() - ydata.min())

    xmin, xmax = xdata.min() - dx, xdata.max() + dx
    ymin, ymax = ydata.min() - dy, ydata.max() + dy

    X, Y = np.mgrid[xmin:xmax:60j, ymin:ymax:60j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([xdata, ydata])

    kernel = stats.gaussian_kde(values, bw_method=0.2)
    Z = kernel(positions).reshape(X.shape)

    Z = gaussian_filter(Z, sigma=1.0)  # enforce smooth topology
    Z /= Z.max()                       # normalize safely

    cs = ax.contour(X, Y, Z, levels=nlevels,
                    linewidths=linewidths, colors=[c])
    levels = cs.levels

    interp = RegularGridInterpolator(
        (X[:, 0], Y[0, :]), Z,
        bounds_error=False, fill_value=Z.min()
    )

    z = interp(np.column_stack((xdata, ydata)))
    mask = z > levels[0]

    ax.scatter(xdata[~mask], ydata[~mask],
               s=s, color=c, marker=mark, alpha=alpha, label=label)


def L2500_to_Lbol(L2500):
    """Convert log L_2500 to log L_bol using Richards+06 relation."""
    return 0.9869 * (np.log10(2.998e18/1549.48) + L2500) + 1.051


def Lbol_to_L2500(Lbol):
    """Convert log L_bol to log L_2500 using Richards+06 relation."""
    return (Lbol - 1.051) / 0.9869 - np.log10(2.998e18/1549.48)


def load_hst_data(data_dir):
    """
    Load HST quasar data.
    
    Returns:
    --------
    dict : Dictionary containing HST data arrays
    """
    hst_file = os.path.join(data_dir, "HST_CIV_Sulentic2007_HSLA2018_finalprops_f2500AP.csv")
    hst = pd.read_csv(hst_file)
    
    z_hst = hst["z"].values
    DLcm_hst = cosmo.luminosity_distance(z_hst).value * 3.086e24  # Mpc to cm
    
    # Compute luminosities in erg/s/Hz (convert f_lambda to f_nu)
    fnu2500_hst = hst["F2500_AP"].values * 2500. / (2.998e18/2500)
    logLnu2500_hst = (np.log10(4.*np.pi) + np.log10(fnu2500_hst) + 
                      2*np.log10(DLcm_hst) + np.log10(1.+z_hst))
    
    # Create masks for different BH mass categories
    maskMBH = (~np.isnan(hst["log M_BH"].values))
    maskMBH_rm = maskMBH & (hst["MBH_source_qual"].values < 2)
    maskMBH_sr = maskMBH & (hst["MBH_source_qual"].values >= 2)
    
    return {
        'z': z_hst,
        'logL': logLnu2500_hst,
        'mask_rm': maskMBH_rm,
        'mask_sr': maskMBH_sr,
        'mask_no_mbh': ~maskMBH
    }


def load_sdss_rankine_data(data_dir):
    """
    Load SDSS data from Rankine+20.
    
    Returns:
    --------
    dict : Dictionary containing SDSS data arrays
    """
    sdss_file = os.path.join(data_dir, "Rankine20_CIV_HeII_wDelGI_wdr16luminosity.fits")
    sdss_tbl = Table.read(sdss_file)
    names = [name for name in sdss_tbl.colnames if len(sdss_tbl[name].shape) <= 1]
    sdss = sdss_tbl[names].to_pandas()
    
    # Filter objects with valid luminosity measurements and good quality
    sdss = sdss[(sdss["LOGL2500_dr16"] > 0) & (sdss["good"] == True)]
    
    lognuLnu2500_sdss = sdss["LOGL2500_dr16"].values
    logLnu2500_sdss = lognuLnu2500_sdss - np.log10(2.9918e18/2500)
    z_sdss = sdss["z_paul"].values
    
    return {
        'z': z_sdss,
        'logL': logLnu2500_sdss
    }


def load_sdss_temple_data(data_dir):
    """
    Load SDSS data with Temple+21 luminosities.
    
    Returns:
    --------
    dict : Dictionary containing SDSS data arrays
    """
    sdss_file = os.path.join(data_dir, "Fig1Data_sdss_withTempleL3000.csv")
    sdss_temple = pd.read_csv(sdss_file)
    
    z_sdss_temple = sdss_temple["Z_SDSS"].values
    # Convert from L3000 to L2500
    logLnu2500_sdss_temple = sdss_temple["L3000_THB21"].values - 14.9775
    
    # Filter out NaN values
    mask_sdss_temple = ~np.isnan(logLnu2500_sdss_temple)
    z_sdss_temple = z_sdss_temple[mask_sdss_temple]
    logLnu2500_sdss_temple = logLnu2500_sdss_temple[mask_sdss_temple]
    
    return {
        'z': z_sdss_temple,
        'logL': logLnu2500_sdss_temple
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
    gnirs = gnirs[(gnirs["good"] == True) & (gnirs["LOGL2500_dr16"].values > 0)]
    
    lognuLnu2500_gnirs = gnirs["LOGL2500_dr16"].values
    logLnu2500_gnirs = lognuLnu2500_gnirs - np.log10(2.998e18/2500)
    z_gnirs = gnirs["z_paul"].values
    
    return {
        'z': z_gnirs,
        'logL': logLnu2500_gnirs
    }


def create_figure(hst_data, sdss_data, gnirs_data, data_version='rankine'):
    """
    Create the L vs z figure.
    
    Parameters:
    -----------
    hst_data : dict
        HST data dictionary
    sdss_data : dict
        SDSS data dictionary
    gnirs_data : dict
        GNIRS data dictionary
    data_version : str
        Version identifier ('rankine' or 'temple')
    
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    # Set up color palette
    cs = palettable.colorbrewer.qualitative.Dark2_5.mpl_colors
    
    # Configure matplotlib for this figure
    plt.rcParams['ytick.right'] = False
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 9))
    
    # Plot HST sample
    ax.plot(hst_data['z'][hst_data['mask_rm']], 
            hst_data['logL'][hst_data['mask_rm']], 
            marker="o", markersize=15, markerfacecolor=cs[3], 
            markeredgecolor="k", label="HST, RM $M_{BH}$", 
            linestyle="", zorder=3)
    
    ax.plot(hst_data['z'][hst_data['mask_sr']], 
            hst_data['logL'][hst_data['mask_sr']], 
            marker="o", markersize=13.5, markerfacecolor="w", 
            markeredgecolor=cs[3], mew=1.8, 
            label="HST, Scaling Relation $M_{BH}$", 
            linestyle="", zorder=2)
    
    ax.plot(hst_data['z'][hst_data['mask_no_mbh']], 
            hst_data['logL'][hst_data['mask_no_mbh']], 
            marker="^", markersize=16.5, markerfacecolor=cs[4], 
            markeredgecolor="k", label="HST, no $M_{BH}$", 
            linestyle="", zorder=2)
    
    # Plot SDSS sample with density contours
    print(f"Plotting SDSS data ({data_version} version)...")
    plot_contour_fast2(sdss_data['z'], sdss_data['logL'], 
                       c=cs[2], mark="o", 
                       nlevels=[0.05, 0.25, 0.50, 0.75, 0.95], 
                       ax=ax, linewidths=1.5, s=3, alpha=0.4, 
                       label="SDSS")
    
    # Plot GNIRS-DQS
    ax.plot(gnirs_data['z'], gnirs_data['logL'], 
            marker="v", markersize=11, markerfacecolor=cs[2], 
            markeredgecolor="k", label="GNIRS-DQS", 
            linestyle="", zorder=3)
    
    # Configure plot aesthetics
    ax.set_xlabel('Redshift', fontsize=35)
    ax.set_ylabel(r'log\,$L_{2500\AA}$ [erg\,s$^{-1}$\,Hz$^{-1}$]', 
                  fontsize=35)
    ax.tick_params(axis='both', which='major', labelsize=27.5)
    ax.legend(frameon=True, labelspacing=0.1, loc="lower right", 
              prop={"size": 25})
    
    # Add secondary y-axis for bolometric luminosity
    secax = ax.secondary_yaxis('right', functions=(L2500_to_Lbol, Lbol_to_L2500))
    secax.set_ylabel(r'log\,$L_\mathrm{Bol}$ [erg\,s$^{-1}$]', fontsize=35)
    secax.tick_params(labelsize=27.5)
    
    # Add threshold line at L_bol = 3e45 erg/s
    L2500_threshold = Lbol_to_L2500(np.log10(3e45))
    ax.plot([-1e5, 1e5], [L2500_threshold, L2500_threshold], 
            "-b", zorder=1)
    
    # Set axis limits
    ax.set_ylim(27.1, 32.5)
    ax.set_xlim(-0.16575, 3.25)
    
    plt.tight_layout()
    
    return fig, ax


def main():
    """Main function to generate Figure 1."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Generate Figure 1: Luminosity vs Redshift',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--data-version', 
                        choices=['rankine', 'temple'],
                        default='rankine',
                        help='SDSS data source: rankine (Rankine+20) or temple (Temple+21)')
    parser.add_argument('--output-dir',
                        default='./Figures/',
                        help='Directory to save the figure')
    
    args = parser.parse_args()
    
    # Set up paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)  # Go up to HST Paper directory
    data_dir = os.path.join(parent_dir, 'Data')
    output_dir = os.path.join(parent_dir, args.output_dir)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 70)
    print("Figure 1: Luminosity versus Redshift")
    print("=" * 70)
    print(f"Data version: {args.data_version}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Load data
    print("Loading data...")
    hst_data = load_hst_data(data_dir)
    print(f"  HST: {len(hst_data['z'])} objects")
    
    if args.data_version == 'rankine':
        sdss_data = load_sdss_rankine_data(data_dir)
        print(f"  SDSS (Rankine+20): {len(sdss_data['z'])} objects")
    else:
        sdss_data = load_sdss_temple_data(data_dir)
        print(f"  SDSS (Temple+21): {len(sdss_data['z'])} objects")
    
    gnirs_data = load_gnirs_data(data_dir)
    print(f"  GNIRS-DQS: {len(gnirs_data['z'])} objects")
    print()
    
    # Create figure
    print("Creating figure...")
    fig, ax = create_figure(hst_data, sdss_data, gnirs_data, args.data_version)
    
    # Save figure
    if args.data_version == 'rankine':
        output_filename = "Fig1_L_versus_z.pdf"
    else:
        output_filename = "Fig1_L_versus_z_Temple.pdf"
    
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved: {output_path}")
    print()
    print("=" * 70)
    print("Done!")
    print("=" * 70)
    
    # Optionally display the figure
    # plt.show()


if __name__ == "__main__":
    main()
