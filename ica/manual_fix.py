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
import os
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
# GUI-authored per-object overrides (manual_fix_gui.py Save writes this). Keyed on
# master-mode FITS stems (e.g. "Mrk 1148_COS"); the going-forward override source
# for master / overrides-only runs. MANUAL_FIX_CONFIG remains only for the legacy
# Sulentic default run.
_OVERRIDES_JSON = Path(__file__).resolve().parent / "manual_fix_overrides.json"

warnings.filterwarnings('ignore')


def _override_kwargs(entry):
    """Map a manual_fix_overrides.json entry to fit_with_overrides() kwargs.

    Returns an empty dict for a missing/empty entry, so the fit reproduces the
    plain pipeline default (== the batch baseline fit)."""
    if not entry:
        return {}
    return dict(
        mask_ranges=entry.get("mask_ranges") or [],
        mask_pixels=entry.get("mask_pixels") or [],
        unmask_ranges=entry.get("unmask_ranges") or [],
        unmask_pixels=entry.get("unmask_pixels") or [],
        comps_use=entry.get("forced_components"),
    )


def _override_has_mask(entry):
    return bool(entry and (entry.get("mask_ranges") or entry.get("mask_pixels")))


def _override_has_unmask(entry):
    return bool(entry and (entry.get("unmask_ranges") or entry.get("unmask_pixels")))


def _override_tag(entry):
    """Filename tag describing which override types an entry carries.

    "baseline"           -> no override (or empty entry)
    "override_mask"      -> only pixels/ranges masked
    "override_unmask"    -> only pixels/ranges force-included
    "override_maskunmask_high" -> masked + unmasked + forced 'high' components
    Mask/unmask are distinguished (per user choice); components append as a token.
    """
    if not entry:
        return "baseline"
    mask_tok = ""
    if _override_has_mask(entry):
        mask_tok += "mask"
    if _override_has_unmask(entry):
        mask_tok += "unmask"
    toks = []
    if mask_tok:
        toks.append(mask_tok)
    comps = entry.get("forced_components")
    if comps in ("mod", "low", "high"):
        toks.append(comps)
    if not toks:
        return "baseline"  # entry present but no actionable override
    return "override_" + "_".join(toks)

OUTPUT_PATH = "ICA_Plots_Rebin_AP"
# REBIN_PATH = "RebinnedSpec_2022Aug11"
# REBIN_PATH = "RebinnedSpec_2026AP"   ← old default; now a constructor arg
_DEFAULT_REBIN_PATH = os.environ.get("HST_PAPER_REBIN_DIR", "RebinnedSpec_2026AP")

class ICAManualFixProcessor:
    """
    Main class for processing objects with manual ICA fixes
    """
    
    # Plot subfolder used for every object in master mode (the Good/Bad/Fixable/
    # "Probably Good" quality classification only exists for the 207 Sulentic objects).
    _MASTER_FOLDER = "All"

    def __init__(self, catalog_file=None, rebin_path=None, master_mode=False,
                 output_path=None, overrides_path=None):
        """
        Initialize the processor with catalog data

        Parameters:
        -----------
        catalog_file : str or None
            Path to the CSV file containing object information.
            Defaults to Data/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv in repo root.
        overrides_path : str or None
            Path to the GUI-authored override store (manual_fix_overrides.json).
            Defaults to the file next to this module. Loaded into self.overrides;
            used as the override source for master / overrides-only runs.
        rebin_path : str or None
            Directory containing rebinned FITS files to run ICA on.
            Defaults to HST_PAPER_REBIN_DIR env var, then "RebinnedSpec_2026AP".
            Use "RebinnedSpec_master" for spectra rebinned from master_catalog_v21.csv.
        master_mode : bool
            If True, drive object selection from the FITS files in rebin_path
            (named {common_name}_{inst}.fits, e.g. from master_catalog_v21.csv)
            instead of the Sulentic catalog. The object name is the FITS filename
            stem; redshift is read from each FITS file; and morphing-edges /
            CIV-BAL / manual-fix overrides (all keyed on old Sulentic Spec_Names)
            are skipped for objects not present in those dicts. Diagnostic plots
            are saved into a single "All" subfolder rather than quality folders.
        output_path : str or None
            Directory to write diagnostic plots into. Defaults to the module-level
            OUTPUT_PATH ("ICA_Plots_Rebin_AP"). run_all_objects passes a separate
            directory in master mode to avoid clobbering Sulentic-run plots.
        """
        if catalog_file is None:
            catalog_file = str(_DATA_DIR / "HST_CIV_Sulentic2007_HSLA2018_finalprops.csv")
        self.dat = pd.read_csv(catalog_file)
        self.names = self.dat["Final_Name"].values
        self.names_spec = self.dat["Spec_Name"].values
        self.inst_final = self.dat["Inst_final"].values

        self.rebin_path = rebin_path if rebin_path is not None else _DEFAULT_REBIN_PATH
        self.master_mode = master_mode
        self.output_path = output_path if output_path is not None else OUTPUT_PATH

        # Clean up names (remove spaces)
        for i in range(len(self.names)):
            self.names[i] = "".join(self.names[i].split(" "))

        # Load quality classification map (replaces runtime os.listdir of ICA_full/).
        # Unused in master mode (objects are not classified), but cheap to load.
        with open(_QUALITY_JSON) as f:
            self.quality_map = json.load(f)

        # Load GUI-authored overrides (keyed on master-mode FITS stems).
        self.overrides_path = overrides_path if overrides_path is not None else str(_OVERRIDES_JSON)
        self.overrides = {}
        if os.path.exists(self.overrides_path):
            with open(self.overrides_path) as f:
                self.overrides = json.load(f)

    def list_master_objects(self):
        """
        Return the list of object identifiers to process in master mode: the
        FITS filename stems ({common_name}_{inst}) of every .fits file in
        rebin_path. One entry per spectrum, so multi-instrument objects yield
        one identifier per instrument (each fit independently).
        """
        files = sorted(glob.glob(os.path.join(self.rebin_path, "*.fits")))
        return [os.path.splitext(os.path.basename(f))[0] for f in files]

    def _ensure_output_dirs(self):
        """Create the plot output directory and its subfolders for the active mode."""
        os.makedirs(self.output_path, exist_ok=True)
        if self.master_mode:
            os.makedirs(os.path.join(self.output_path, self._MASTER_FOLDER), exist_ok=True)
        else:
            for sub in ("Probably Good", "Good", "Fixable", "Bad"):
                os.makedirs(os.path.join(self.output_path, sub), exist_ok=True)

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
            Folder name ('Good', 'Bad', 'Fixable', 'Probably Good'), or the
            single "All" folder in master mode (objects are unclassified).
        """
        if self.master_mode:
            return self._MASTER_FOLDER
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
        if self.master_mode:
            # Master mode: the identifier IS the FITS filename stem
            # ({common_name}_{inst}). No Sulentic catalog lookup; redshift is
            # read from the FITS below. spec_name = stem, so the morphing-edges
            # and CIV-BAL dicts (keyed on old Sulentic Spec_Names) simply do not
            # match new master objects and are skipped.
            spec_name = name
            fn = os.path.join(self.rebin_path, name + ".fits")
            if not os.path.exists(fn):
                raise FileNotFoundError(
                    f'No rebinned FITS for "{name}" in {self.rebin_path}')
        else:
            index = np.where(self.names == name_final)[0][0]
            spec_name = self.names_spec[index]
            inst = self.inst_final[index]

            fn_list = glob.glob("%s/%s*%s.fits"%(self.rebin_path, spec_name, inst))
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

    def fit_with_overrides(self, name, mask_ranges=None, mask_pixels=None,
                           unmask_ranges=None, unmask_pixels=None, comps_use=None):
        """Run the ICA pipeline from *explicit* overrides and return fit arrays.

        This is the single shared fit engine used by both the interactive GUI
        (ica.manual_fix_gui.GuiFixProcessor.fit_for_gui delegates here) and the
        batch runner, so a GUI fit and a batch fit for the same override are
        numerically identical. It reproduces process_object()'s numerical path
        (setup_object -> main_ICA -> get_CIV_parameters); only the source of the
        mask/comps overrides differs. No matplotlib is touched, so it is safe to
        call off a UI thread.

        Parameters
        ----------
        name : str
            Object identifier (master-mode FITS stem, or Sulentic Final_Name).
        mask_ranges : list of [lo, hi] or None
            Rest-frame wavelength ranges to mask (bad-pixel flag).
        mask_pixels : list of float or None
            Wavelengths whose nearest pixel is masked.
        unmask_ranges, unmask_pixels :
            Ranges / pixels to force-*include* even if the pipeline masked them
            (NAL/BAL/input bad). Honouring an unmask requires re-fitting with the
            pipeline's recomputed mask cleared at those positions (see below).
        comps_use : {'mod', 'low', 'high', None, '', 'auto'}
            Force a component set; None/''/'auto' means auto-select by chi-squared.

        Returns
        -------
        dict with name, spec_name, z, the diagnostic arrays (wave_arb, flux_arb,
        errs_arb, mask_arb, wave_ica, flux_ica), f2500, civ_blue, civ_ew.
        """
        mask_ranges = mask_ranges or []
        mask_pixels = mask_pixels or []
        unmask_ranges = unmask_ranges or []
        unmask_pixels = unmask_pixels or []

        # Reload from disk every fit -> guaranteed-clean mask state, no cross-fit
        # mutation bugs. Cheap next to the ~5-10 s fit itself.
        wave, flux, z, errs, mask, spec_name = self.setup_object(name)
        mask = np.asarray(mask, dtype=float).copy()

        for lo, hi in mask_ranges:
            mask[(wave >= lo) & (wave <= hi)] = 1
        for wl in mask_pixels:
            mask[np.argmin(np.abs(wave - wl))] = 1

        cu = None if comps_use in (None, "", "auto") else comps_use

        (wave_arb, flux_arb, errs_arb, mask_arb,
         wave_ica, flux_ica, f2500) = run_ICA_r20_components.main_ICA(
            wave, flux, errs, mask, z, name="", ica_path=None, comps_use=cu)

        # Force-unmask pixels the pipeline masked (NAL/BAL/input bad). Those masks
        # are recomputed inside main_ICA, so honouring an unmask means re-fitting
        # with them cleared. We reuse main_ICA's own outputs: flux_arb carries the
        # *real* (not median-replaced) flux at NAL pixels, which is exactly what we
        # want at an unmasked pixel; still-masked pixels are ignored either way.
        # When nothing is unmasked we skip this entirely -> identical to the batch.
        if unmask_ranges or unmask_pixels:
            mask_ref = np.asarray(mask_arb, dtype=float).copy()
            for lo, hi in unmask_ranges:
                mask_ref[(wave_arb >= lo) & (wave_arb <= hi)] = 0
            for wl in unmask_pixels:
                mask_ref[np.argmin(np.abs(wave_arb - wl))] = 0

            wave_ica2, flux_ica2 = run_ICA_r20_components.get_ICA(
                wave_arb, flux_arb, errs_arb, mask_ref, z,
                ica_path=None, use_priors=False, comps_use=cu)

            # Rescale f2500 by the original fit's real-unit factor at 2500 A
            # (norm/morph scaling is fixed by the data, not the fit).
            i0 = np.argmin(np.abs(wave_ica - 2500.))
            i2 = np.argmin(np.abs(wave_ica2 - 2500.))
            if flux_ica[i0]:
                f2500 = float(flux_ica2[i2] * (f2500 / flux_ica[i0]))
            wave_ica, flux_ica = wave_ica2, flux_ica2
            # Show the mask actually used for the fit so unmasked pixels no longer
            # render as NAL/BAL (they were set to 0 in mask_ref).
            mask_arb = mask_ref

        civ_blue, civ_ew = self.get_CIV_parameters(
            wave_arb, flux_arb, wave_ica, flux_ica, name)

        return dict(
            name=name, spec_name=spec_name, z=float(z),
            wave_arb=wave_arb, flux_arb=flux_arb, errs_arb=errs_arb, mask_arb=mask_arb,
            wave_ica=wave_ica, flux_ica=flux_ica,
            f2500=float(f2500), civ_blue=float(civ_blue), civ_ew=float(civ_ew),
        )

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

        # AP 2026-09-23: continuum must be indexed with the SAME mask as the flux.
        # The loop below runs over the EW-selected subset, so `continuum[i]` read the
        # full-array element i -- i.e. the continuum ~240 A too blue (the ICA grid
        # starts at 1260 A, the window at 1500 A). Because quasar continua fall with
        # wavelength the divisor was too large, biasing every EW low by a median
        # factor of 1.46. Masking once here also avoids re-masking inside the loop.
        # See diagnostics/CIV_EW_Indexing_Investigation.md
        cont_EW = continuum[EW]
        
        CIV_EW = 0
        ew_list = [0.]
        for i in range(len(wave_r[EW])):
            try:
                CIV_EW += max(((flux_r[EW][i] - cont_EW[i]) / cont_EW[i]) * 
                             (wave_r[EW][i+1] - wave_r[EW][i]), 0)
            except IndexError:
                CIV_EW += max(((flux_r[EW][i] - cont_EW[i]) / cont_EW[i]) * 
                             (wave_r[EW][i] - wave_r[EW][i-1]), 0)
            ew_list.append(CIV_EW)
        
        # Calculate blueshift
        ind_half_flux = abs((CIV_EW / 2) - np.array(ew_list)).argmin()
        CIV_blue = ((1549.48 - wave_r[EW][ind_half_flux]) / 1549.48) * 3e5
        
        return CIV_blue, CIV_EW
    
    
    def create_diagnostic_plot(self, wave_arb, flux_arb, errs_arb, mask_arb,
                              wave_ica, flux_ica, name, spec_name, CIV_blue, CIV_EW, show=False, save=False,
                              file_tag=None):
        """
        Create diagnostic plots for the ICA fit

        file_tag : str or None
            If given, the saved filename is "{spec_name}_{file_tag}.png" (used by
            the override-aware path, e.g. "baseline" / "override_maskunmask_high").
            If None, the legacy MANUAL_FIX_CONFIG-derived tag is used (default mode).
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
            if file_tag is not None:
                config_type = file_tag
            else:
                config_type = "nomask"
                if name in MANUAL_FIX_CONFIG:
                    if MANUAL_FIX_CONFIG[name]['custom_mask_pixels'] is not None:
                        config_type = "mask"
                    if MANUAL_FIX_CONFIG[name]['forced_components'] is not None:
                        config_type += f"_comps{MANUAL_FIX_CONFIG[name]['forced_components']}"

            # plt.savefig(f'ICA_Plots_Refactored/{spec_name}_{config_type}.png', dpi=150)
            print(f'Spec_name is {spec_name}')
            plt.savefig(f'{self.output_path}/{self.get_folder(spec_name)}/{spec_name}_{config_type}.png', dpi=150)
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
        cont_EW = continuum[EW]   # AP 2026-09-23: see get_CIV_parameters above
        CIV_EW = 0
        ew_list = [0.]
        for i in range(len(wave_r[EW])):
            try:
                CIV_EW += max(( (flux_r[EW][i] - cont_EW[i]) / cont_EW[i] ) * ( wave_r[EW][i+1] - wave_r[EW][i] ), 0) #no absorption
            except IndexError:
                CIV_EW += max(( (flux_r[EW][i] - cont_EW[i]) / cont_EW[i] ) * ( wave_r[EW][i] - wave_r[EW][i-1] ), 0)

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
        # Master / overrides-driven path: source the override from the JSON store
        # (going-forward override source) and fit through the shared engine, which
        # also handles mask ranges and force-unmasking. Legacy default (Sulentic)
        # mode keeps the MANUAL_FIX_CONFIG path below.
        if self.master_mode:
            return self.process_object_override(
                name, self.overrides.get(name), plot=plot, save_plot=save_plot)

        print(f"\n=== Processing {name} ===")

        # Setup object data
        wave_orig, flux_orig, z, errs_orig, mask_orig, spec_name = self.setup_object(name)
        print(f'Median of flux_orig: {np.nanmedian(flux_orig)}')  # Debugging line
        # Apply manual fixes
        mask_orig, comps_use = self.apply_manual_fix(name, wave_orig, mask_orig)
        
        # Run ICA analysis
        wave_arb, flux_arb, errs_arb, mask_arb, wave_ica, flux_ica, f2500_ica = \
            run_ICA_r20_components.main_ICA(wave_orig, flux_orig, errs_orig, mask_orig, z,
                                           name="", ica_path=None, comps_use=comps_use)
        
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
            self._ensure_output_dirs()

        if plot or save_plot:
            self.create_diagnostic_plot(wave_arb, flux_arb, errs_arb, mask_arb,
                                       wave_ica, flux_ica, name, spec_name, CIV_blue, CIV_EW,
                                       show=plot, save=save_plot)

        return results

    def process_object_override(self, name, override, plot=False, save_plot=False,
                                tag=None, is_override=None):
        """
        Process a single object through the shared override-aware fit engine.

        Parameters
        ----------
        name : str
            Object identifier (master-mode FITS stem).
        override : dict or None
            A manual_fix_overrides.json entry, or None for the plain baseline fit.
        tag : str or None
            Filename tag for the saved plot. Defaults to _override_tag(override)
            (i.e. "baseline" when override is None, else "override_<types>").
        is_override : bool or None
            Value recorded in the CSV 'override' column. Defaults to
            (override is not None); pass explicitly to force a baseline row.

        Returns
        -------
        results : dict
            One CSV row for this fit.
        """
        if tag is None:
            tag = _override_tag(override)
        if is_override is None:
            is_override = override is not None

        print(f"\n=== Processing {name} [{tag}] ===")

        res = self.fit_with_overrides(name, **_override_kwargs(override))
        spec_name = res["spec_name"]
        comps_use = (_override_kwargs(override) or {}).get("comps_use")

        results = {
            'object_name': name,
            'spec_name': spec_name,
            'redshift': res["z"],
            'CIV_blueshift': res["civ_blue"],
            'CIV_EW': res["civ_ew"],
            'f2500': res["f2500"],
            'components_used': comps_use,
            'override': bool(is_override),
            'override_type': tag,
            'manual_masking_applied': _override_has_mask(override) or _override_has_unmask(override),
            'forced_components': comps_use is not None,
        }

        print("Results [%s]: CIV_blue=%.1f km/s, CIV_EW=%.2f A, F2500=%.3e"
              % (tag, res["civ_blue"], res["civ_ew"], res["f2500"]))

        if save_plot:
            self._ensure_output_dirs()
        if plot or save_plot:
            self.create_diagnostic_plot(
                res["wave_arb"], res["flux_arb"], res["errs_arb"], res["mask_arb"],
                res["wave_ica"], res["flux_ica"], name, spec_name,
                res["civ_blue"], res["civ_ew"], show=plot, save=save_plot, file_tag=tag)

        return results

    def batch_process_overrides(self, object_list=None,
                                output_file="ICA_Results_Rebin_overrides.csv",
                                plot=False, save_plots=False):
        """
        Run only objects that have a manual override, fitting each TWICE: once
        with no override (baseline) and once with the override applied. Both fits
        are plotted (distinct filename tags) and both are written to the results
        CSV (one row each, distinguished by the 'override' column), so the two
        can be compared directly.

        Parameters
        ----------
        object_list : list or None
            Object identifiers to run. If None, uses every key in self.overrides.
        output_file : str
            Output CSV filename (two rows per object).
        plot, save_plots : bool
            Show / save the diagnostic plots.

        Returns
        -------
        results_df : pandas.DataFrame
        """
        if object_list is None:
            object_list = sorted(self.overrides.keys())

        if save_plots:
            self._ensure_output_dirs()

        all_results = []
        for name in object_list:
            entry = self.overrides.get(name)
            # Baseline (no override) then the overridden fit.
            for is_override in (False, True):
                ov = entry if is_override else None
                tag = _override_tag(entry) if is_override else "baseline"
                try:
                    result = self.process_object_override(
                        name, ov, plot=plot, save_plot=save_plots,
                        tag=tag, is_override=is_override)
                    all_results.append(result)
                except Exception as e:
                    print(f"ERROR processing {name} [{tag}]: {str(e)}")
                    continue

        results_df = pd.DataFrame(all_results)
        results_df.to_csv(output_file, index=False)
        print(f"\nResults saved to {output_file}")
        n_obj = len(object_list)
        print(f"Processed {len(all_results)} fits for {n_obj} override object(s) "
              f"(baseline + override each).")

        return results_df

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
            self._ensure_output_dirs()
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


