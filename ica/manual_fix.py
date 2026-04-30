# Originally: Trevor Code/ica_manual_fix_refactored.py
# Copied: 2026-04-29 (Phase 2 migration)
# Changes from source:
#   - Removed sys.path.append hacks
#   - from Small_Pix_Filter_TVM import -> from rebinning.small_pix_filter import
#   - import MorphingEdges_* -> from ica import morphing_edges (all 5 files merged)
#   - import CIV_BAL_regions -> from ica import civ_bal_regions as CIV_BAL_regions
#   - import run_ICA_r20_components -> from ica import run_ica as run_ICA_r20_components
#   - import spec_morph -> from ica import spec_morph
#   - import components -> removed (dead import, never used in function bodies)
#   - from manual_fix_config import -> from ica.manual_fix_config import
#   - get_folder(): replaced os.listdir scan with JSON lookup (quality_classification.json)
#   - __init__: added json load for quality_map; catalog path -> Data/ relative
#   - REBIN_PATH: kept as-is (caller must set or pass via env)
"""
Refactored Manual ICA Fix System

This module provides a systematic approach to applying manual fixes to ICA analysis
based on configuration settings. It replaces the ad-hoc approach in Manual Object Fix.ipynb.
"""

from pathlib import Path
import json
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib
from astropy.io import fits

# Import existing modules
from rebinning.small_pix_filter import SDSS_pixel_filter
from ica import morphing_edges
from ica import morphing_edges as MorphingEdges         # MorphingEdges.MorphingEdges dict
from ica import civ_bal_regions as CIV_BAL_regions
from lmfit import minimize, Parameters
from ica import run_ica as run_ICA_r20_components
from ica import spec_morph
import warnings

# Import our manual fix configuration
from ica.manual_fix_config import MANUAL_FIX_CONFIG

_DATA_DIR = Path(__file__).resolve().parent.parent / "Data"
_QUALITY_JSON = Path(__file__).resolve().parent / "quality_classification.json"

warnings.filterwarnings('ignore')

OUTPUT_PATH = "ICA_Plots_Rebin_AP"
# REBIN_PATH = "RebinnedSpec_2022Aug11"
REBIN_PATH = "RebinnedSpec_2026AP"

class ICAManualFixProcessor:
    """
    Main class for processing objects with manual ICA fixes
    """
    
    def __init__(self, catalog_file=None):
        """
        Initialize the processor with catalog data

        Parameters:
        -----------
        catalog_file : str or None
            Path to the CSV file containing object information.
            Defaults to Data/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv in repo root.
        """
        if catalog_file is None:
            catalog_file = str(_DATA_DIR / "HST_CIV_Sulentic2007_HSLA2018_finalprops.csv")
        self.dat = pd.read_csv(catalog_file)
        self.names = self.dat["Final_Name"].values
        self.names_spec = self.dat["Spec_Name"].values
        self.inst_final = self.dat["Inst_final"].values

        # Clean up names (remove spaces)
        for i in range(len(self.names)):
            self.names[i] = "".join(self.names[i].split(" "))

        # Load quality classification map (replaces runtime os.listdir of ICA_full/)
        with open(_QUALITY_JSON) as f:
            self.quality_map = json.load(f)
    
    def get_folder(self, name):
        """
        Get the quality folder for a specific object from quality_classification.json.
        Replaces the original os.listdir scan of ICA_full/ subdirectories.

        Parameters:
        -----------
        name : str
            Object name (Final_Name / Spec_Name prefix)

        Returns:
        --------
        str
            Folder name ('Good', 'Bad', 'Fixable', 'Probably Good')
        """
        folder = self.quality_map.get(name)
        if folder is None:
            raise Exception(f'Object "{name}" not found in quality_classification.json')
        return folder
    
    def plot_HST(self, wave, flux, mask, ax):
        iplot_start = 0
        #iplot_end   = 0
        mask_curr   = mask[0]
        for i in range(len(wave)):
            #go until mask is a different value
            if mask[i]!=mask_curr or i==len(wave)-1:
                if mask_curr==3:
                    #plot BAL region
                    ax.plot(wave[iplot_start:i], flux[iplot_start:i], "-m", zorder=2)
                elif mask_curr==2:
                    #plot NAL region
                    ax.plot(wave[iplot_start:i], flux[iplot_start:i], "-y", zorder=2)
                elif mask_curr==1:
                    #plot badpix region
                    ax.plot(wave[iplot_start:i], flux[iplot_start:i], "-m", zorder=2)
                else:
                    #plot BAL region
                    ax.plot(wave[iplot_start:i], flux[iplot_start:i], "-k", zorder=1)

                mask_curr = mask[i]
                iplot_start = i-1
        ax.plot([1e7,1e7], [1e7,1e8], "-m", zorder=1, alpha=0.1, label="Input Bad Pixels")
        ax.legend()                                     

    def custom_pix_filter(self, wave_orig, mask_orig, filter_wavelengths):
        """
        Apply custom pixel filtering based on wavelength array
        
        Parameters:
        -----------
        wave_orig : array
            Original wavelength array
        mask_orig : array
            Original mask array
        filter_wavelengths : array
            Wavelengths to mask
            
        Returns:
        --------
        mask_orig : array
            Updated mask with custom pixels flagged
        indices : array
            Indices of masked pixels
        """
        indices = []
        for wavelength in filter_wavelengths:
            index = np.argmin(np.abs(wave_orig - wavelength))
            indices.append(index)
        
        indices = np.unique(np.array(indices))
        mask_orig[indices] = 1
        return mask_orig, indices
    
    def setup_object(self, name):
        """
        Set up an object for analysis - loads data and applies morphing/BAL regions
        
        Parameters:
        -----------
        name : str
            Object name
            
        Returns:
        --------
        wave_orig, flux_orig, z, errs_orig, mask_orig, spec_name : tuple
            Processed spectral data
        """
        name_final = name
        index = np.where(self.names == name_final)[0][0]
        spec_name = self.names_spec[index]
        inst = self.inst_final[index]

        fn_list = glob.glob("%s/%s*%s.fits"%(REBIN_PATH, spec_name, inst))
        fn = fn_list[0]
        
        spec = fits.open(fn)
        # Print the column names for debugging
        print(f'File column names: {spec[1].data.columns}') # Debugging line
        
        wave_orig = spec[1].data["Rest-Frame Wavelength"]
        flux_orig = spec[1].data["Coadded Flux (Arbitrary Units)"]
        z         = spec[1].data["Redshift"][0]
        flux_orig /= np.nanmedian(flux_orig)
        errs_orig = spec[1].data["Coadded Flux Errors"] / np.nanmedian(flux_orig)
        
        # Apply morphing if specified
        if spec_name in MorphingEdges.MorphingEdges:
            try: 
                print(f"Applying morphing for {spec_name}: {MorphingEdges.MorphingEdges[spec_name]}")
                flux_orig, morph_coeff, unique_markers = spec_morph.morph(
                    wave_orig*(1+z), flux_orig, z, np.nanmedian(flux_orig), 
                    unique_markers=MorphingEdges.MorphingEdges[spec_name], 
                    return_markers=True)
                errs_orig *= morph_coeff
            except:
                print(f"Morphing failed for: {spec_name}, final name: {name_final}")
                pass
        
        mask_orig = spec[1].data["Bad Pixel Mask"]
        
        # Apply manual BAL regions if specified
        if spec_name in CIV_BAL_regions.BAL_regions:
            print(f"{spec_name} is in CIV_BAL_regions")
            for window in CIV_BAL_regions.BAL_regions[spec_name]:
                mask_orig[((wave_orig>=window[0])&(wave_orig<=window[1]))] = 1
        
        return wave_orig, flux_orig, z, errs_orig, mask_orig, spec_name
    
    def apply_manual_fix(self, name, wave_orig, mask_orig):
        """
        Apply manual fix based on configuration
        
        Parameters:
        -----------
        name : str
            Object name
        wave_orig : array
            Wavelength array
        mask_orig : array
            Original mask
            
        Returns:
        --------
        mask_orig : array
            Updated mask
        comps_use : str or None
            Component set to force
        """
        comps_use = None
        
        if name in MANUAL_FIX_CONFIG:
            config = MANUAL_FIX_CONFIG[name]
            
            # Apply custom pixel masking
            if config['custom_mask_pixels'] is not None:
                print(f"Applying custom pixel masking for {name}")
                mask_orig, indices = self.custom_pix_filter(
                    wave_orig, mask_orig, config['custom_mask_pixels'])
                print(f"Masked {len(indices)} pixels")
            
            # Set forced components
            if config['forced_components'] is not None:
                comps_use = config['forced_components']
                print(f"Forcing components to: {comps_use}")
        else:
            print(f"No manual fix configuration found for {name}")
        
        return mask_orig, comps_use
    
    def get_CIV_parameters(self, wave, flux, wave_r, flux_r, name, 
                          EW_region=[1500,1600], cont_region=[[1445,1465],[1700,1705]]):
        """
        Extract CIV parameters from ICA fit
        
        Parameters:
        -----------
        wave : array
            Original wavelength
        flux : array
            Original flux
        wave_r : array
            Reconstructed wavelength
        flux_r : array
            Reconstructed flux
        name : str
            Object name
        EW_region : list
            Wavelength region for EW calculation
        cont_region : list
            Continuum fitting regions
            
        Returns:
        --------
        CIV_blue : float
            CIV blueshift in km/s
        CIV_EW : float
            CIV equivalent width in Angstroms
        """
        # Fit continuum
        cont1 = ((wave_r>=cont_region[0][0])&(wave_r<=cont_region[0][1]))
        cont2 = ((wave_r>=cont_region[1][0])&(wave_r<=cont_region[1][1]))
        
        m, b = np.polyfit(np.concatenate((wave_r[cont1], wave_r[cont2])), 
                         np.concatenate((flux_r[cont1], flux_r[cont2])), 1)
        continuum = wave_r * m + b
        
        # Calculate EW
        EW = ((wave_r>=EW_region[0])&(wave_r<=EW_region[1]))
        
        CIV_EW = 0
        ew_list = [0.]
        for i in range(len(wave_r[EW])):
            try:
                CIV_EW += max(((flux_r[EW][i] - continuum[i]) / continuum[i]) * 
                             (wave_r[EW][i+1] - wave_r[EW][i]), 0)
            except IndexError:
                CIV_EW += max(((flux_r[EW][i] - continuum[i]) / continuum[i]) * 
                             (wave_r[EW][i] - wave_r[EW][i-1]), 0)
            ew_list.append(CIV_EW)
        
        # Calculate blueshift
        ind_half_flux = abs((CIV_EW / 2) - np.array(ew_list)).argmin()
        CIV_blue = ((1549.48 - wave_r[EW][ind_half_flux]) / 1549.48) * 3e5
        
        return CIV_blue, CIV_EW
    
    
    def create_diagnostic_plot(self, wave_arb, flux_arb, errs_arb, mask_arb, 
                              wave_ica, flux_ica, name, spec_name, CIV_blue, CIV_EW, show=False, save=False):
        """
        Create diagnostic plots for the ICA fit
        """
        c4mask = wave_arb > 1400
        ylow, yhigh = max(0, np.nanpercentile(flux_arb, 1)), \
                      np.nanpercentile(flux_arb[c4mask], 99) + np.nanmedian(flux_arb[c4mask])
        ylow_err, yhigh_err = max(0, np.nanpercentile(errs_arb, 1)), np.nanpercentile(errs_arb, 99)
        
        xlow, xhigh = 1500, 1600
        xlow0, xhigh0 = max(min(wave_arb), 1250), min(max(wave_arb), max(wave_ica))
        
        fig = plt.figure(figsize=(18, 10.5), constrained_layout=True)
        gs = GridSpec(7, 12, figure=fig)
        
        ax0 = fig.add_subplot(gs[:3, :8])  # Full spectrum
        ax1 = fig.add_subplot(gs[3:6, :8])  # CIV region
        ax2 = fig.add_subplot(gs[3:, 8:])   # CIV fit
        ax3 = fig.add_subplot(gs[6, :8], sharex=ax1)  # Errors
        plt.subplots_adjust(hspace=0)

        # Plot full spectrum
        ax0.plot(wave_arb, flux_arb, "-k", alpha=0.6)#, label="Data")
        ax0.plot(wave_ica, flux_ica, "-r")#, label="ICA fit")
        ax0.set_xlim(xlow0, xhigh0)
        ax0.set_ylim(ylow, yhigh + 5)
        ax0.set_ylabel("Flux Density (Arb. Units)", fontsize=20)
        ax0.set_title(f"{name} - ICA Fit", fontsize=18)
        ax0.legend()

        if spec_name in CIV_BAL_regions.BAL_regions:
            for window in CIV_BAL_regions.BAL_regions[spec_name]:
                ax1.axvspan(window[0], window[1], color="r", alpha=0.3)

        self.plot_HST(wave_arb, flux_arb, mask_arb, ax0)
        

        # Plot CIV region
        ax1.plot(wave_arb, flux_arb, "-k", alpha=0.6)
        ax1.plot(wave_ica, flux_ica, "-r")
        ax1.set_xlim(xlow, xhigh)
        ax1.set_ylim(ylow, yhigh + 5)
        ax1.set_ylabel("Flux Density (Arb. Units)", fontsize=20)
        self.plot_HST(wave_arb, flux_arb, mask_arb, ax1)

        # Plot errors
        ax3.plot(wave_arb, errs_arb, "-k", alpha=0.6)
        ax3.set_xlim(xlow, xhigh)
        ax3.set_ylim(ylow_err, yhigh_err)
        ax3.set_xlabel("Wavelength (Å)", fontsize=20)
        ax3.set_ylabel("Error", fontsize=17.5)
        
        # Plot CIV analysis
        self.plot_CIV_analysis(ax2, wave_arb, flux_arb, wave_ica, flux_ica, name, CIV_blue, CIV_EW)
        
        plt.tight_layout()
        
        if save:
            config_type = "nomask"
            if name in MANUAL_FIX_CONFIG:
                if MANUAL_FIX_CONFIG[name]['custom_mask_pixels'] is not None:
                    config_type = "mask"
                if MANUAL_FIX_CONFIG[name]['forced_components'] is not None:
                    config_type += f"_comps{MANUAL_FIX_CONFIG[name]['forced_components']}"
            
            # plt.savefig(f'ICA_Plots_Refactored/{spec_name}_{config_type}.png', dpi=150)
            print(f'Spec_name is {spec_name}')
            plt.savefig(f'{OUTPUT_PATH}/{self.get_folder(spec_name)}/{spec_name}_{config_type}.png', dpi=150)
        if show:
            plt.show()
    
    def plot_CIV_analysis(self, ax, wave, flux, wave_r, flux_r, name, CIV_blue, CIV_EW,
                        EW_region=[1500,1600], cont_region=[[1445,1465],[1700,1705]]):
        """
        Plot CIV analysis on given axis
        """
        # ylow, yup = max(0, np.percentile(flux_r, 1) - np.nanmedian(flux_r)/5), \
                    # np.percentile(flux_r, 99) + np.nanmedian(flux_r)
        ylow, yup = max(0, np.percentile(flux_r, 1)-np.nanmedian(flux_r)/5), \
                        np.percentile(flux_r, 99)+np.nanpercentile(flux_r, 99)

        ax.plot(wave, flux, "-k", alpha=0.5)#, label="Data")
        ax.plot(wave_r, flux_r, "-r", lw=1.8)#, label="ICA fit")
        ax.plot([1549.48, 1549.48], [ylow, yup], "--k", label="CIV Rest")
        
        # Fit and plot continuum
        cont1 = ((wave_r>=cont_region[0][0])&(wave_r<=cont_region[0][1]))
        cont2 = ((wave_r>=cont_region[1][0])&(wave_r<=cont_region[1][1]))
        ax.axvspan(cont_region[0][0], cont_region[0][1], alpha=0.5, color='grey')
        ax.axvspan(cont_region[1][0], cont_region[1][1], alpha=0.5, color='grey')
        
        m, b = np.polyfit(np.concatenate((wave_r[cont1], wave_r[cont2])), 
                         np.concatenate((flux_r[cont1], flux_r[cont2])), 1)
        continuum = wave_r * m + b
        ax.plot(wave_r, continuum, "-m")#, label="Continuum")
        
        # Show EW region and blueshift
        EW = ((wave_r>=EW_region[0])&(wave_r<=EW_region[1]))
        CIV_EW = 0
        ew_list = [0.]
        for i in range(len(wave_r[EW])):
            try:
                CIV_EW += max(( (flux_r[EW][i] - continuum[i]) / continuum[i] ) * ( wave_r[EW][i+1] - wave_r[EW][i] ), 0) #no absorption
            except IndexError:
                CIV_EW += max(( (flux_r[EW][i] - continuum[i]) / continuum[i] ) * ( wave_r[EW][i] - wave_r[EW][i-1] ), 0)

            ew_list.append(CIV_EW)

        ind_half_flux = abs((CIV_EW / 2) - np.array(ew_list)).argmin()

        ax.fill_between(wave_r[EW], continuum[EW], flux_r[EW], color="blue", alpha=0.2, 
                       label=f"EW = {CIV_EW:.1f} Å")
        ax.plot([EW_region[0], EW_region[0]], [ylow, yup], "-b", alpha=0.7)
        ax.plot([EW_region[1], EW_region[1]], [ylow, yup], "-b", alpha=0.7)
        ax.plot([wave_r[EW][ind_half_flux],wave_r[EW][ind_half_flux]], [ylow,yup], "--b", label="%.1f km/s" % (CIV_blue))

        ax.set_xlim(1435, 1710)
        ax.set_ylim(ylow, yup)
        ax.set_xlabel("Wavelength (Å)", fontsize=20)
        ax.set_title(name, fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=17.5)
        ax.tick_params(axis='both', which='minor', labelsize=17.5)
        ax.tick_params(labelleft=False)
        ax.legend(loc="best", prop={"size":12.5})

    def process_object(self, name, plot=False, save_plot=False):
        """
        Process a single object with manual fixes
        
        Parameters:
        -----------
        name : str
            Object name
        plot : bool
            Whether to create diagnostic plots
        save_plot : bool
            Whether to save plots
            
        Returns:
        --------
        results : dict
            Dictionary containing all fit results
        """
        print(f"\n=== Processing {name} ===")
        
        # Setup object data
        wave_orig, flux_orig, z, errs_orig, mask_orig, spec_name = self.setup_object(name)
        print(f'Median of flux_orig: {np.nanmedian(flux_orig)}')  # Debugging line
        # Apply manual fixes
        mask_orig, comps_use = self.apply_manual_fix(name, wave_orig, mask_orig)
        
        # Run ICA analysis
        wave_arb, flux_arb, errs_arb, mask_arb, wave_ica, flux_ica, f2500_ica = \
            run_ICA_r20_components.main_ICA(wave_orig, flux_orig, errs_orig, mask_orig, z, 
                                           name="", ica_path="", comps_use=comps_use)
        
        # Extract CIV parameters
        CIV_blue, CIV_EW = self.get_CIV_parameters(wave_arb, flux_arb, wave_ica, flux_ica, name)
        
        # Compile results
        results = {
            'object_name': name,
            'spec_name': spec_name,
            'redshift': z,
            'CIV_blueshift': CIV_blue,
            'CIV_EW': CIV_EW,
            'f2500': f2500_ica,
            'components_used': comps_use,
            'manual_masking_applied': name in MANUAL_FIX_CONFIG and MANUAL_FIX_CONFIG[name]['custom_mask_pixels'] is not None,
            'forced_components': comps_use is not None
        }
        
        print(f"Results: CIV_blue={CIV_blue:.1f} km/s, CIV_EW={CIV_EW:.2f} Å, F2500={f2500_ica:.3e}")
        
        # Create diagnostic plots if requested
        if save_plot:
            os.makedirs(OUTPUT_PATH, exist_ok=True)
            os.makedirs(f'{OUTPUT_PATH}/Probably Good', exist_ok=True)
            os.makedirs(f'{OUTPUT_PATH}/Good', exist_ok=True)
            os.makedirs(f'{OUTPUT_PATH}/Fixable', exist_ok=True)
            os.makedirs(f'{OUTPUT_PATH}/Bad', exist_ok=True)
        
        if plot or save_plot:
            self.create_diagnostic_plot(wave_arb, flux_arb, errs_arb, mask_arb, 
                                       wave_ica, flux_ica, name, spec_name, CIV_blue, CIV_EW, 
                                       show=plot, save=save_plot)
        
        return results
    
    def batch_process(self, object_list=None, output_file="manual_fix_results.csv", 
                     plot=False, save_plots=False):
        """
        Process multiple objects and save results
        
        Parameters:
        -----------
        object_list : list or None
            List of object names to process. If None, process all in config
        output_file : str
            Output CSV filename
        plot : bool
            Whether to show plots
        save_plots : bool
            Whether to save plots
            
        Returns:
        --------
        results_df : pandas.DataFrame
            DataFrame with all results
        """
        if object_list is None:
            object_list = list(MANUAL_FIX_CONFIG.keys())
        
        # Create output directory for plots if needed
        if save_plots:
            os.makedirs(OUTPUT_PATH, exist_ok=True)
            os.makedirs(f'{OUTPUT_PATH}/Probably Good', exist_ok=True)
            os.makedirs(f'{OUTPUT_PATH}/Good', exist_ok=True)
            os.makedirs(f'{OUTPUT_PATH}/Fixable', exist_ok=True)
            os.makedirs(f'{OUTPUT_PATH}/Bad', exist_ok=True)
        all_results = []
        
        for obj_name in object_list:
            try:
                result = self.process_object(obj_name, plot=plot, save_plot=save_plots)
                all_results.append(result)
            except Exception as e:
                print(f"ERROR processing {obj_name}: {str(e)}")
                continue
        
        # Convert to DataFrame and save
        results_df = pd.DataFrame(all_results)
        results_df.to_csv(output_file, index=False)
        print(f"\nResults saved to {output_file}")
        print(f"Successfully processed {len(all_results)}/{len(object_list)} objects")
        
        return results_df

# Example usage
if __name__ == "__main__":
    # Initialize processor
    processor = ICAManualFixProcessor()
    
    # Process a single object for testing
    result = processor.process_object('RBS1763', plot=False, save_plot=False)
    
    # Batch process all objects
    # results = processor.batch_process(plot=False, save_plots=True)


