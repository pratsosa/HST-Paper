# Originally: Trevor Code/ICA Scripts/fit_composites.py
# Copied: 2026-04-29 (Phase 2 migration)
# Changes from source:
#   - Removed sys.path.append hacks
#   - import plot_ICA, import spec_morph -> commented out (dead imports; neither is called in this file)
#   - path_priors default: "./" -> str(_C3_PRIORS_DIR) + "/" (points to Data/c3_priors/)
"""
6 February 2023
Philly
Trevor McCaffrey

Fit Paul Hewett's composite spectra to an input spectrum.
After the best-fit spectrum is determined, return the weights
used to fit that composite to use as priors for the spectrum
that is actually being fit.  Note that this code assumes
the wavelength resolution of SDSS (shouldn't be a problem
for the ICA module since the re-binning is a part of that code)
"""

#necessary imports - should only need these
#see Rankine+2020 for lmfit installation
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from lmfit import minimize, Parameters
from astropy.io import fits
import glob

# import plot_ICA   # dead import in source; not used in any function body
# import spec_morph # dead import in source; morph2 call is commented out

_C3_PRIORS_DIR = Path(__file__).resolve().parent.parent / "Data" / "c3_priors"


def residual(params, composite, data, eps_data):
    '''
    compute residuals between data and model,
    weighted by the variance array;
    used by minimize() function in ICA_fit()
    '''
    model = sum(params["A"].value*composite)
    return (data-model)**2 * eps_data

def fit_comp(composite, wave, spectrum, ivar, flags):
    # FIXME: explain this
    #this is called ICA_fit, but we get the component weights from minimize()
    '''
    Instansiate Parameters and add empty weight values;
    then compute the ICA weights using minimize()
    '''
    params = Parameters()
    params.add('A', value=0)
    ivar[flags>0] = 0.
    return minimize(residual, params, args=(composite, spectrum, ivar), method="leastsq", nan_policy='omit')

def get_fit(wave, flux, errs, mask, wave_comp, flux_comp):
    #Apply mask and get fits
    wavecoverage_mask  = (wave>=1260)&(wave<=max(wave_comp))
    comp_mask  = (wave_comp>min(wave[wavecoverage_mask]))&(wave_comp<max(wave[wavecoverage_mask]))

    #Adjust array sizes
    wave_fit        = wave_comp[comp_mask]
    composite_fit       = flux_comp[comp_mask]

    if len(wave[wavecoverage_mask])==len(wave_fit):
        fit_y  = fit_comp(composite_fit, wave[wavecoverage_mask], flux[wavecoverage_mask], 1/(errs[wavecoverage_mask]**2), mask[wavecoverage_mask])
    else:
        #input wave is greater by one pixel, cut components byb one
        if len(wave[wavecoverage_mask])-len(wave_fit) == -1:
            fit_y  = fit_comp(composite_fit[:,:-1], wave[wavecoverage_mask], flux[wavecoverage_mask], 1/(errs[wavecoverage_mask]**2), mask[wavecoverage_mask])
        elif len(wave[wavecoverage_mask])-len(wave_fit) == 1:
            fit_y  = fit_comp(composite_fit, wave[wavecoverage_mask][:-1], flux[wavecoverage_mask][:-1], 1/(errs[wavecoverage_mask][:-1]**2), mask[wavecoverage_mask][:-1])

    #compute weights
    weight  = fit_y.params["A"].value
    return fit_y.redchi, weight


def main(waveSpec, fluxSpec, errsSpec, maskSpec, z, path_priors=None):
    #main() should return the three sets of prior weights for the best-fit composite
    if path_priors is None:
        path_priors = str(_C3_PRIORS_DIR) + "/"

    wave = waveSpec.copy()
    flux = fluxSpec.copy()
    errs = errsSpec.copy()
    mask = maskSpec.copy() ; mask[mask>0] = 1

    #normalize spectrum
    #norm_coeff = np.nanmedian(flux)
    #flux /= norm_coeff
    #errs /= norm_coeff

    fn_comp_list = []
    redchi_list  = []
    for fn in glob.glob(path_priors+"MattS_CIV_tplte/*.fits"):
        #print(fn) #---- AP
        comp = fits.open(fn)
        c0 = comp[0].header["COEFF0"]
        c1 = comp[0].header["COEFF1"]
        npix = comp[0].header["NGOOD"]
        wave_composite = 10.**np.linspace(c0, c0+c1*(npix-1), npix)
        flux_composite = comp[0].data[0,:]
        #flux_composite_morph, _ = spec_morph.morph2(wave_composite*(1+z), flux_composite, flux_composite.copy(), z, "")
        #redchi_i, weight_i = get_fit(wave, flux, errs, mask, wave_composite, flux_composite_morph)
        redchi_i, weight_i = get_fit(wave, flux, errs, mask, wave_composite, flux_composite)
        fn_comp_list.append(fn.split("/")[-1])
        redchi_list.append(redchi_i)
    #print(fn_comp_list)
    #print(redchi_list)
    fn_comp_list = np.array(fn_comp_list)
    redchi_list  = np.array(redchi_list)
    fn_comp_best = fn_comp_list[np.argmin(redchi_list)].split('\\')[-1]
    #print(fn_comp_best) ----AP

    df_mod = pd.read_csv(path_priors+"prior_weights_fromcomposites/modEW_priors.csv")
    df_low = pd.read_csv(path_priors+"prior_weights_fromcomposites/lowEW_priors.csv")
    df_high= pd.read_csv(path_priors+"prior_weights_fromcomposites/highEW_priors.csv")

    # print(df_mod['Composite_Name'])
    # print(df_low['Composite_Name'])
    # print(df_high['Composite_Name'])

    #print(df_mod["Composite_Name"].values==fn_comp_best) #AP
    #print(df_low["Composite_Name"].values==fn_comp_best) #AP
    #print(df_high["Composite_Name"].values==fn_comp_best) #AP
    weights_mod = df_mod.iloc[df_mod["Composite_Name"].values==fn_comp_best,1:].values
    weights_low = df_low.iloc[df_low["Composite_Name"].values==fn_comp_best,1:].values
    weights_high= df_high.iloc[df_high["Composite_Name"].values==fn_comp_best,1:].values

    # print(weights_low)  #AP
    # print(weights_mod)  #AP
    # print(weights_high) #AP
    return weights_low.flatten(), weights_mod.flatten(), weights_high.flatten()
