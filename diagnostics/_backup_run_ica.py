# Originally: Trevor Code/ICA Scripts/run_ICA_r20_components.py
# Copied: 2026-04-29 (Phase 2 migration)
# Changes from source:
#   - Removed sys.path.append hacks (Trevor's Dropbox paths)
#   - import fit_composites -> from ica import fit_composites
#   - import plot_ICA -> from ica import plot_ica as plot_ICA
#   - import spec_morph -> from ica import spec_morph
#   - load_ICA() default ica_path: "./" -> str(_COMPONENTS_DIR) + "/"
#     (points to Data/components/ via Path(__file__))
#   - All function signatures that have ica_path="./": unchanged (callers can override)
#   - No algorithm changes
"""
13 January 2023
Seattle en route to Philly, Post-AAS241
Trevor McCaffrey

Running ICA on different data sets can be annoying because
just loading everything up in a notebook and the functions
take up loads of space, making it kind of a mess to work with.

Beyond that, it'd be good to share a script like this with
people like Keri and Gordon so they can quickly run ICA
without waiting on me.  Hopefully this works!
"""

#necessary imports - should only need these
#see Rankine+2020 for lmfit installation
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib
import numpy as np
import pandas as pd
from lmfit import minimize, Parameters
from ica import fit_composites
from ica import plot_ica as plot_ICA
from ica import spec_morph

import inspect

_COMPONENTS_DIR = Path(__file__).resolve().parent.parent / "Data" / "components"

def arg2min(arr):
    #return index of second closest point
    return np.argsort(arr)[1]

def rebin(wave, flux, errs, mask, orig_width_kms):
    #create re-binned wavelength array - this is the easy part!
    #resolution is constant in log-space for SDSS, so easier to stick with that
    error_multiple = np.sqrt(orig_width_kms/69.)

    logwave = np.log10(wave)
    loglam_rebin = np.arange(logwave[0], logwave[-1], 0.0001)
    loglam_rebin_edges = loglam_rebin - 0.0001/2
    loglam_rebin_edges = np.append(loglam_rebin_edges, loglam_rebin[-1]+0.0001/2)

    #start empty
    flux_rebin = np.nan*np.zeros(len(loglam_rebin))
    errs_rebin = np.nan*np.zeros(len(loglam_rebin))
    mask_rebin = np.zeros(len(loglam_rebin))

    #and fill them in
    for i in range(len(loglam_rebin_edges)-1):
        ledge_rebin, redge_rebin = loglam_rebin_edges[i], loglam_rebin_edges[i+1]
        #orig edges
        #oldwidth =

        maskOrig = (logwave>ledge_rebin)&(logwave<redge_rebin)
        if not maskOrig.sum()==0:
            #flux_rebin[i] = ws.numpy_weighted_median(flux[maskOrig], weights=1./(errs[maskOrig]**2))
            flux_rebin[i] = np.nanmedian(flux[maskOrig])
            errs_rebin[i] = np.nanmedian(errs[maskOrig])
            mask_rebin[i] = np.nanmax(mask[maskOrig])

    if (np.isnan(flux_rebin)).any():
        #initial resolution could be low enough such that some re-binned pixels weren't assigned a value
        #for those, interpolate with nearest neighbors
        for j in range(len(flux_rebin)):
            if np.isnan(flux_rebin[j]): #should make sure here it's not a known "gap", instead just ones missed from above
                #arg1 = np.argmin(loglam_rebin[~np.isnan(flux_rebin)]-logwave[j])
                arg1 = np.argmin(np.abs(loglam_rebin[~np.isnan(flux_rebin)]-loglam_rebin[j]))
                arg2 = arg2min(np.abs(loglam_rebin[~np.isnan(flux_rebin)]-loglam_rebin[j]))
                #print(10.**loglam_rebin[arg1], 10.**loglam_rebin[arg2])
                loglam1, flux1, errs1 = loglam_rebin[~np.isnan(flux_rebin)][arg1], \
                                        flux_rebin[~np.isnan(flux_rebin)][arg1], \
                                        errs_rebin[~np.isnan(flux_rebin)][arg1]
                loglam2, flux2, errs2 = loglam_rebin[~np.isnan(flux_rebin)][arg2], \
                                        flux_rebin[~np.isnan(flux_rebin)][arg2], \
                                        errs_rebin[~np.isnan(flux_rebin)][arg2]
                oldwidth = abs(loglam2-loglam1)
                #do linear fit on nearest non-nan neighbors
                m_flux, b_flux = np.polyfit([loglam1, loglam2], [flux1, flux2], 1)
                #m_flux, b_flux = np.polynomial.polynomial.Polynomial.fit([loglam1, loglam2], [flux1, flux2], 1)
                # Changed to np.polynomial.polynomial.Polynomial.fit

                flux_rebin[j] = m_flux*loglam_rebin[j] + b_flux
                m_errs, b_errs = np.polyfit([loglam1, loglam2], [errs1, errs2], 1)
                # m_errs, b_errs = np.polynomial.polynomial.Polynomial.fit([loglam1, loglam2], [errs1, errs2], 1)
                # Changed to np.polynomial.polynomial.Polynomial.fit
                errs_rebin[j] = (m_errs*loglam_rebin[j] + b_errs) #* np.sqrt(oldwidth/0.0001)
                #leave mask[j]==0 for now

    wave_rebin = 10.**loglam_rebin
    errs_rebin *= error_multiple
    return wave_rebin, flux_rebin, errs_rebin, mask_rebin

def blueshift(wave_half_flux):
    '''
    compute blueshift given wavelength where half the EW
    of CIV has been accumulated; units Angstroms
    '''
    return ((1549.48 - wave_half_flux) / 1549.48) * 3e5

def get_CIV(wave, flux, wave_r, flux_r, name, ax, EW_region=[1500,1600], cont_region=[[1445,1465],[1700,1705]]):
    #print("Getting CIV")
    #Going to save plots, but not show with magic above
    #"""
    ylow, yup = max(0, np.percentile(flux_r, 1)-np.nanmedian(flux_r)/5), np.percentile(flux_r, 99)+np.nanmedian(flux_r)
    #fig, ax = plt.subplots(figsize=(9,9))
    #fig = plt.figure(figsize=(9,9))
    ax.plot(wave, flux, "-k", alpha=0.5)#, label="HST co-add")
    ax.plot(wave_r, flux_r, "-b", lw=1.8)#, label="ICA fit")
    ax.plot([1549.48,1549.48], [ylow,yup], "--k", label="CIV Rest Wavelength")
    #"""

    #Fit continuum -- ADD MANUAL HERE
    cont1 = ((wave_r>=cont_region[0][0])&(wave_r<=cont_region[0][1]))
    cont2 = ((wave_r>=cont_region[1][0])&(wave_r<=cont_region[1][1]))
    ax.axvspan(cont_region[0][0], cont_region[0][1], alpha=0.5, color='grey')
    ax.axvspan(cont_region[1][0], cont_region[1][1], alpha=0.5, color='grey')

    m,b = np.polyfit(np.concatenate((wave_r[cont1], wave_r[cont2])), np.concatenate((flux_r[cont1], flux_r[cont2])), 1)
    continuum = wave_r*m + b
    ax.plot(wave_r, continuum, "-m")

    #Plot EW region under flux_r -- ADD MANUAL HERE
    EW = ((wave_r>=EW_region[0])&(wave_r<=EW_region[1]))

    #And compute the EW
    CIV_EW = 0
    ew_list = [0.]
    for i in range(len(wave_r[EW])):
        try:
            CIV_EW += max(( (flux_r[EW][i] - continuum[i]) / continuum[i] ) * ( wave_r[EW][i+1] - wave_r[EW][i] ), 0) #no absorption
        except IndexError:
            CIV_EW += max(( (flux_r[EW][i] - continuum[i]) / continuum[i] ) * ( wave_r[EW][i] - wave_r[EW][i-1] ), 0)

        ew_list.append(CIV_EW)

    ax.fill_between(wave_r[EW], continuum[EW], flux_r[EW], color="blue", alpha=0.2, label="Equivalent Width = %.1f Å"%CIV_EW)

    ind_half_flux = abs((CIV_EW / 2) - np.array(ew_list)).argmin()
    CIV_blue = blueshift(wave_r[EW][ind_half_flux])
    #"""
    #for w in np.arange(1450,1530,10): ax.plot([w,w],[ylow,yup],"--",label=str(w))
    #for w in np.arange(1610,1640,10): ax.plot([w,w],[ylow,yup],"--",label=str(w))
    ax.plot([EW_region[0],EW_region[0]], [ylow,yup], "-b", alpha=0.7)
    ax.plot([EW_region[1],EW_region[1]], [ylow,yup], "-b", alpha=0.7)
    ax.plot([wave_r[EW][ind_half_flux],wave_r[EW][ind_half_flux]], [ylow,yup], "--b", label="Blueshift = %.1f km/s" % (CIV_blue))
    ax.set_xlim(1435, 1710)
    ax.set_ylim(ylow, yup)
    ax.set_xlabel("Rest Wavelength (Å)", fontsize=30)
    #ax.set_ylabel("Flux (Arbitrary Units)", fontsize=20)
    #ax.set_title("Blue, EW = %.2f km/s, %.2f Å" % (CIV_blue,CIV_EW), fontsize=20)
    ax.set_title(name, fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=30)
    ax.tick_params(axis='both', which='minor', labelsize=30)
    ax.tick_params(labelleft=False)
    ax.legend(loc="best", prop={"size":20})
    return CIV_blue, CIV_EW

def pixel_filter(wavelength,data,npix=20):

    filterval=npix//2
    filterval2=npix//2+1

    #print(len(data), data.dtype, data)
    indices = ~np.isnan(data) #indices we'll look at

    end        = data[indices][-filterval2-1:]
    end_wave   = wavelength[indices][-filterval2-1:]
    start      = data[indices][:filterval+1]
    start_wave = wavelength[indices][:filterval+1]
    right_median = np.nanmedian(end)
    left_median  = np.nanmedian(start)

    new_left_wave   = np.zeros(filterval)
    new_right_wave  = np.zeros(filterval)
    new_left_flux   = np.zeros(filterval)
    new_right_flux  = np.zeros(filterval)
    wave_devs_left  = np.zeros(filterval)
    wave_devs_right = np.zeros(filterval)
    flux_devs_left  = np.zeros(filterval)
    flux_devs_right = np.zeros(filterval)

    #Now flipping the mirrored spectra over the left and right medians:
    for k in range(filterval):
        wave_devs_left[-k-1]  = start_wave[k+1] - min(start_wave) #mirror across 0th pixel
        new_left_wave[-k-1]   = -2*wave_devs_left[-k-1] + start_wave[k+1]
        flux_devs_left[-k-1]  = start[k+1] - start[0] #left_median - TVM changed to edge pixel since I think that's right?
        new_left_flux[-k-1]   = -2*flux_devs_left[-k-1]+start[k+1] #want to add to mirrored array starting at end and moving left

        wave_devs_right[k]    = end_wave[-k-2] - max(end_wave)
        new_right_wave[k]     = -2*wave_devs_right[k] + end_wave[-k-2]
        flux_devs_right[k]    = end[-k-2] - end[-1] #right_median
        new_right_flux[k]     = -2*flux_devs_right[k]+end[-k-2]

    extended_data_array=np.concatenate((new_left_flux,data,new_right_flux),axis=0)
    continuum = np.array([np.nanmedian(extended_data_array[l-filterval:l+filterval]) for l in np.arange(filterval,len(data)+filterval,1)])

    return continuum

def maskNAL(wave, flux, errs, mask):
    N = 3. #flag pixels N*sigma below median
    flux_NAL = flux.copy()
    mask_NAL = mask.copy()
    flux_median_61p = pixel_filter(wave, flux, 61)
    NALpix = (flux_median_61p-flux > N*errs)
    mask_NAL[NALpix] = 2.
    flux_NAL[NALpix] = flux_median_61p[NALpix]
    #Pad the NAL mask three pixels each side
    for i in range(len(wave)):
        if mask_NAL[i]==2.:
            mask_NAL[i-3:i+3+1] = 3. #this keeps from masking the entire spectrum in the for loop!
            flux_NAL[i-3:i+3+1] = flux_median_61p[i-3:i+3+1]
    mask_NAL[mask_NAL==3] = 2 #save them all as two
    return mask_NAL, flux_NAL, flux_median_61p

def get_f2500(wave, flux, errs, mask, z, norm_coeff, ica_path=None):
    if ica_path is None:
        ica_path = str(_COMPONENTS_DIR) + "/"
    #real units!  also using the final generated mask
    wave_norm_ica, flux_norm_ica = get_ICA(wave, flux, errs, mask, z,
                                 ica_path=ica_path, use_priors=False,
                                 plot_spectrum=False, comps_use=None)
    ind2500 = np.argmin(abs(wave_norm_ica-2500.))
    return flux_norm_ica[ind2500]*norm_coeff

def load_ICA(ica_path=None):
    """
    load in the ICA components; there are three sets:
    (1) amy_12603000_10c_180421.comp are the 10 "main" components
    (2) amy_12753000_lowew_10c_181101_v1.comp are 10 different components
        which typically work better on low EW spectra when the main components fail
    (3) amy_12653000_hew_hsn_7c_190302.comp are 7 components which typically work
        better on high EW spectra when the main components fail

    Note the difference in wavelength coverage for each set of components!

    First load components in as a pandas dataframe, then turn them into
    (Ncomponents)x(Npixels) arrays.
    """
    if ica_path is None:
        ica_path = str(_COMPONENTS_DIR) + "/"
    wave_modEW          = pd.read_csv(ica_path+"wav_12603000.dat", names=["wave"]).values.flatten()
    components_modEW_df = pd.read_csv(ica_path+"amy_12603000_10c_180421.comp", sep=r"\s+",
                          names=["c1","c2","c3","c4","c5","c6","c7","c8","c9","c10"])
    components_modEW    = np.array([components_modEW_df.iloc[:,i].values for i in range(components_modEW_df.shape[1])])

    wave_lowEW          = pd.read_csv(ica_path+"wav_12753000.dat", names=["wave"]).values.flatten()
    components_lowEW_df = pd.read_csv(ica_path+"amy_12753000_lowew_10c_181101_v1.comp", sep=r"\s+",
                          names=["c1","c2","c3","c4","c5","c6","c7","c8","c9","c10"])
    components_lowEW    = np.array([components_lowEW_df.iloc[:,i].values for i in range(components_lowEW_df.shape[1])])

    wave_highEW          = pd.read_csv(ica_path+"wav_12653000.dat", names=["wave"]).values.flatten()
    components_highEW_df = pd.read_csv(ica_path+"amy_12653000_hew_hsn_7c_190302.comp", sep=r"\s+",
                          names=["c1","c2","c3","c4","c5","c6","c7"])
    components_highEW    = np.array([components_highEW_df.iloc[:,i].values for i in range(components_highEW_df.shape[1])])

    return wave_modEW, components_modEW, wave_lowEW, components_lowEW, wave_highEW, components_highEW

def residual(params, comps, data, eps_data):
    '''
    compute residuals between data and model,
    weighted by the variance array;
    used by minimize() function in ICA_fit()
    '''
    model = sum(params[w].value*comps[n] for (n, w) in enumerate(params))
    return (data-model)**2 * eps_data

def ICA_fit(components, wave, spectrum, ivar, flags, priors):
    # FIXME: explain this
    #this is called ICA_fit, but we get the component weights from minimize()
    '''
    Instansiate Parameters and add empty weight values;
    then compute the ICA weights using minimize()
    '''
    params = Parameters()
    for i in range(len(components)):
        params.add('W%d' % (i+1), value=priors[i])
    #print(params)
    #badpix = masking(wave, spectrum, 1/np.sqrt(ivar))
    ivar[flags>0] = 0.
    return minimize(residual, params, args=(components, spectrum, ivar), method="lbfgsb", nan_policy='omit')

def get_ICA(wave, flux, errs, mask, z, ica_path=None, use_priors=True, plot_spectrum=True, comps_use=None):
    if ica_path is None:
        ica_path = str(_COMPONENTS_DIR) + "/"

    caller_frame = inspect.stack()[1]
    caller_name = caller_frame.function
    #print(f"Function '{caller_name}' called get_ICA()")
    #Debugging stuff I was doing - AP

    #most of this function just aligns the pixels from your spectrum and the ICA components
    wave_mod, components_mod, wave_low, components_low, wave_high, components_high = load_ICA(ica_path)

    #Apply mask and get fits
    mod_mask  = (wave>1260)&(wave<3000)
    low_mask  = (wave>1275)&(wave<3000)
    high_mask = (wave>1265)&(wave<3000)

    mod_comp_mask  = (wave_mod>min(wave[mod_mask]))&(wave_mod<max(wave[mod_mask]))
    low_comp_mask  = (wave_low>min(wave[low_mask]))&(wave_low<max(wave[low_mask]))
    high_comp_mask = (wave_high>min(wave[high_mask]))&(wave_high<max(wave[high_mask]))
    # I believe this is just finding the maximal region where both the inputted wavelengths and the ICA wavelengths are defined

    #print(mod_comp_mask.shape)
    ###########
    #Adjust array sizes
    wave_mod_fit        = wave_mod[mod_comp_mask]
    components_mod_fit  = components_mod[:,mod_comp_mask]
    wave_low_fit        = wave_low[low_comp_mask]
    components_low_fit  = components_low[:,low_comp_mask]
    wave_high_fit       = wave_high[high_comp_mask]
    components_high_fit = components_high[:,high_comp_mask]

    ###############
    #Now do fitting
    if use_priors:
        prior_low, prior_mod, prior_high = fit_composites.main(wave, flux, errs, mask, z)
    else:
        prior_low, prior_mod, prior_high = np.zeros(10), np.zeros(10), np.zeros(7)

    #print(prior_low)

    # mod EW
    #print(len(wave[mod_mask]), len(wave_mod_fit))
    if len(wave[mod_mask])==len(wave_mod_fit):
        fit_mod  = ICA_fit(components_mod_fit, wave[mod_mask], flux[mod_mask], 1/(errs[mod_mask]**2), mask[mod_mask], priors=prior_mod)
    else:
        #input wave is greater by one pixel, cut components byb one
        if len(wave[mod_mask])-len(wave_mod_fit) == -1:
            fit_mod  = ICA_fit(components_mod_fit[:,:-1], wave[mod_mask], flux[mod_mask], 1/(errs[mod_mask]**2), mask[mod_mask], priors=prior_mod)
        elif len(wave[mod_mask])-len(wave_mod_fit) == 1:
            fit_mod  = ICA_fit(components_mod_fit, wave[mod_mask][:-1], flux[mod_mask][:-1], 1/(errs[mod_mask][:-1]**2), mask[mod_mask][:-1], priors=prior_mod)
    #low EW
    if len(wave[low_mask])==len(wave_low_fit):
        fit_low  = ICA_fit(components_low_fit, wave[low_mask], flux[low_mask], 1/(errs[low_mask]**2), mask[low_mask], priors=prior_low)
    else:
        #input wave is greater by one pixel, cut components byb one
        if len(wave[low_mask])-len(wave_low_fit) == -1:
            fit_low  = ICA_fit(components_low_fit[:,:-1], wave[low_mask], flux[low_mask], 1/(errs[low_mask]**2), mask[low_mask], priors=prior_low)
        elif len(wave[low_mask])-len(wave_low_fit) == 1:
            fit_low  = ICA_fit(components_low_fit, wave[low_mask][:-1], flux[low_mask][:-1], 1/(errs[low_mask][:-1]**2), mask[low_mask][:-1], priors=prior_low)

    #high EW
    if len(wave[high_mask])==len(wave_high_fit):
        fit_high  = ICA_fit(components_high_fit, wave[high_mask], flux[high_mask], 1/(errs[high_mask]**2), mask[high_mask], priors=prior_high)
    else:
        #input wave is greater by one pixel, cut components byb one
        if len(wave[high_mask])-len(wave_high_fit) == -1:
            fit_high  = ICA_fit(components_high_fit[:,:-1], wave[high_mask], flux[high_mask], 1/(errs[high_mask]**2), mask[high_mask], priors=prior_high)
        elif len(wave[high_mask])-len(wave_high_fit) == 1:
            #print(wave_high_fit.shape, wave[high_mask].shape)
            fit_high  = ICA_fit(components_high_fit, wave[high_mask][:-1], flux[high_mask][:-1], 1/(errs[high_mask][:-1]**2), mask[high_mask][:-1], priors=prior_high)

    #this is what all branches above would look like w/o index problems
    #fit_mod   = ICA_fit(components_mod[:,mod_comp_mask], wave[mod_mask], flux[mod_mask], 1/(errs[mod_mask]**2), mask[mod_mask]) if len(wave[mod_mask])==len(wave_mod) else ICA_fit(components_mod[:,:-1], wave[mod_mask], flux[mod_mask], 1/(errs[mod_mask]**2), mask[mod_mask])
    #fit_low   = ICA_fit(components_low[:,mod_comp_mask], wave[low_mask], flux[low_mask], 1/(errs[low_mask]**2), mask[low_mask]) if len(wave[low_mask])==len(wave_low) else ICA_fit(components_low[:,:-1], wave[low_mask], flux[low_mask], 1/(errs[low_mask]**2), mask[low_mask])
    #fit_high  = ICA_fit(components_high[:,mod_comp_mask], wave[high_mask], flux[high_mask], 1/(errs[high_mask]**2), mask[high_mask]) if len(wave[high_mask])==len(wave_high) else ICA_fit(components_high[:,:-1], wave[high_mask], flux[high_mask], 1/(errs[high_mask]**2), mask[high_mask])

    #compute weights
    weights_mod  = [fit_mod.params[i].value for i in fit_mod.params]
    weights_low  = [fit_low.params[i].value for i in fit_low.params]
    weights_high = [fit_high.params[i].value for i in fit_high.params]

    if comps_use is None:
        if fit_mod.redchi < fit_high.redchi and fit_mod.redchi < fit_low.redchi:
            #print("using mod EW components")
            return wave_mod, np.dot(weights_mod, components_mod)
        elif fit_low.redchi < fit_high.redchi and fit_low.redchi < fit_mod.redchi:
            #print("using low EW components")
            return wave_low, np.dot(weights_low, components_low)
        else:
            #print("using high EW components")
            return wave_high, np.dot(weights_high, components_high)

    #enter this branch if you want to tell ICA which components to use
    else:
        if comps_use=="mod":
            return wave_mod, np.dot(weights_mod, components_mod)
        elif comps_use=="low":
            return wave_low, np.dot(weights_low, components_low)
        elif comps_use=="high":
            return wave_high, np.dot(weights_high, components_high)
        else:
            print("Didn't recognize comps_use")

def maskIterate(wave, flux, errs, mask, z, ica_path=None):
    if ica_path is None:
        ica_path = str(_COMPONENTS_DIR) + "/"
    #Iterative masking for broad absorption
    #It seems that, without this step, the reconstructions are in general underestimated
    #Start: if (reconstruction-flux)>N*sigma, mask this pixel
    N = 2.0
    #Break here if N > 4
    while N <= 4:
        #Mask *potential* broad absorption initially
        maskBAL = mask.copy()
        #maskBAL[(wave>=1295)&(wave<=1400)] = 99 #SiIV+OIV
        #maskBAL[(wave>=1430)&(wave<=1546)] = 99 #CIV
        #maskBAL[(wave>=1780)&(wave<=1880)] = 99 #SiIV+OIV
        wave_ica, flux_ica = get_ICA(wave, flux, errs, maskBAL, z,
                                      ica_path=ica_path, use_priors=True) #normally True - AP
        #-------------------------Change back to True
        """
        fig, ax = plt.subplots(1,1, figsize=(15,5))
        ax.plot(wave, flux)
        ax.plot(wave_ica, flux_ica)
        plt.show()
        """
        #Start from scratch with original wave, flux, etc
        niter = 0
        while niter < 10:
            nbad = 0
            #Work out the window from full spec we're working on
            wave_begin = max( min(wave_ica), min(wave) )
            wave_end   = min( max(wave_ica), max(wave) )
            ind_begin = np.argmin( abs(wave-wave_begin) )
            ind_end   = np.argmin( abs(wave-wave_end) )

            #Do the actual masking
            for i in range(ind_begin, ind_end+1):
                #ind_reconst = np.argmin( abs(wave[i]-wave_ica) )
                i_left  = max(ind_begin, i-30)
                i_right = min(ind_end, i+30)
                ir_left  = np.argmin( abs(wave[i_left]-wave_ica) )
                ir_right = np.argmin( abs(wave[i_right]-wave_ica) )
                #change here to only mask that pixel for now;
                #should only mask the 10 nearest after the initial bad pixels are established
                # print(ir_left)
                # print(ir_right)
                # print(i_left)
                # print(i_right)

                badpix = (flux_ica[ir_left:ir_right+1]-flux[i_left:i_right+1]) > N*errs[i_left:i_right+1]
                if badpix.sum() > (len(badpix)//2) and mask[i]==0 and maskBAL[i]==99:
                    #print("found a bad pix")
                    nbad += 1
                    #mask[max(ind_begin,i-10):min(ind_end,i+11)] = 1
                    mask[i] = 3

            #If no badpix, mask has converged and we can break
            if nbad == 0:
                break
            else:
                for i in range(len(mask)):
                    if mask[i]==3 and i<ind_end:
                        mask[i-10:i+11] = 4
                niter += 1
                #print("N = %.2f, Try Number %d" % (N,niter))
                wave_ica, flux_ica = get_ICA(wave, flux, errs, mask, z,
                                              ica_path=ica_path, use_priors=False)
        N += 0.25
        if nbad == 0:
            break
    mask[mask==4] = 3
    return mask

def f2500_ICA(waveSpec, fluxSpec, errsSpec, maskSpec, z, name="", ica_path=None, use_priors=False, comps_use=None, plot_process=False):
    if ica_path is None:
        ica_path = str(_COMPONENTS_DIR) + "/"
    #this one skips the morphing
    wave = waveSpec.copy()
    flux = fluxSpec.copy()
    errs = errsSpec.copy()
    mask = maskSpec.copy() ; mask[mask>0] = 1

    # check wavelength resolution, and rebin if it doesn't match components (~69 km/s)
    resolution = np.nanmedian([3.e5*((wave[i+1]-wave[i])/wave[i]) for i in range(len(wave)-1)])
    if round(resolution)!=69:
        #print("Input wavelength resolution %d km/s does not match 69 km/s.  Re-binning spectrum." % resolution)
        wave, flux, errs, mask = rebin(wave, flux, errs, mask, resolution)

    #normalize spectrum
    norm_coeff = np.nanmedian(flux)
    flux /= norm_coeff
    errs /= norm_coeff

    #Mask pixels >3sigma below pseudo continuum
    #initial mask is just from initial pipeline (e.g. SDSS, HST)
    mask_NAL, flux_NAL, flux_median_61p = maskNAL(wave, flux, errs, mask)

    #do initial fit without BAL mask
    #wave_init_ica, flux_init_ica = get_ICA(wave, flux_morph, errs_morph, mask_NAL, z, ica_path, use_priors=use_priors)

    #generate iterative BAL mask
    ## FIXME: for BAL quasars, can start with either manually defined mask covering troughs
    ##        or generic definition done in Rankine+20
    ##        don't want to mask generic BAL troughs for non-BALs though
    mask_witer = maskIterate(wave, flux_NAL, errs, mask_NAL, z, ica_path)

    #now do ICA fit with final mask applied
    wave_ica, flux_ica = get_ICA(wave, flux_NAL, errs, mask_witer, z,
                                 ica_path=ica_path, use_priors=use_priors, comps_use=comps_use)

    #get f2500 from reconstruction
    #de-morph and de-normalize the spectrum to get back to real units
    # FIXME: need to get indexing right here for edge cases - commenting out for now
    flux_ica_realunits = flux_ica * norm_coeff
    f2500_ica = flux_ica[np.argmin(abs(wave_ica-2500.))] * norm_coeff
    return wave, flux*norm_coeff, errs*norm_coeff, mask_witer, wave_ica, flux_ica*norm_coeff, f2500_ica

def main_ICA(waveSpec, fluxSpec, errsSpec, maskSpec, z, name="", ica_path=None, use_priors=False, comps_use=None, plot_process=False):
    if ica_path is None:
        ica_path = str(_COMPONENTS_DIR) + "/"
    wave = waveSpec.copy()
    flux = fluxSpec.copy()
    errs = errsSpec.copy()
    mask = maskSpec.copy() ; mask[mask>0] = 1

    # check wavelength resolution, and rebin if it doesn't match components (~69 km/s)
    resolution = np.nanmedian([3.e5*((wave[i+1]-wave[i])/wave[i]) for i in range(len(wave)-1)])
    if round(resolution)!=69:
        print("Input wavelength resolution %d km/s does not match 69 km/s.  Re-binning spectrum." % resolution)
        wave, flux, errs, mask = rebin(wave, flux, errs, mask, resolution)

    _, _, _, _, _, _, f2500_ica_debug = f2500_ICA(waveSpec, fluxSpec, errsSpec, maskSpec, z, name="", ica_path=ica_path, use_priors=False, comps_use=None, plot_process=False)
    print(f'Debug f2500_ica from f2500_ICA function: {f2500_ica_debug}')

    #normalize spectrum
    norm_coeff = np.nanmedian(flux)
    flux /= norm_coeff
    errs /= norm_coeff

    #Mask pixels >3sigma below pseudo continuum
    #initial mask is just from initial pipeline (e.g. SDSS, HST)
    mask_NAL, flux_NAL, flux_median_61p = maskNAL(wave, flux, errs, mask)

    #standard-ize spectrum shape by morphing
    flux_morph_NAL, morph_coeff = spec_morph.morph2(wave*(1+z), flux_NAL, errs, z, "") #note NAL pixels replaces in flux array
    flux_morph = flux * morph_coeff
    errs_morph = errs * morph_coeff

    #do initial fit without BAL mask
    #wave_init_ica, flux_init_ica = get_ICA(wave, flux_morph, errs_morph, mask_NAL, z, ica_path, use_priors=use_priors)

    #generate iterative BAL mask
    ## FIXME: for BAL quasars, can start with either manually defined mask covering troughs
    ##        or generic definition done in Rankine+20
    ##        don't want to mask generic BAL troughs for non-BALs though
    mask_witer = maskIterate(wave, flux_morph_NAL, errs_morph, mask_NAL, z)

    #now do ICA fit with final mask applied
    wave_ica, flux_ica = get_ICA(wave, flux_morph_NAL, errs_morph, mask_witer, z,
                                 ica_path=ica_path, use_priors=use_priors, comps_use=comps_use)

    #get f2500 from reconstruction
    #de-morph and de-normalize the spectrum to get back to real units
    # FIXME: need to get indexing right here for edge cases - commenting out for now
    try:
        iwave_start = np.argmin(abs(wave-min(wave_ica)))
        flux_ica_realunits = flux_ica / morph_coeff[iwave_start:iwave_start+len(wave_ica)] * norm_coeff
        f2500_ica = flux_ica_realunits[np.argmin(abs(wave_ica-2500.))]
    #Note the input units - normalized, but not morphed
    except ValueError:
        print("Using get_f2500 function to compute f2500 directly.") # AP - added this for debugging
        print("Median morph_coeff:", np.nanmedian(morph_coeff))
        f2500_ica = get_f2500(wave, flux, errs, mask_witer, z, norm_coeff, ica_path=ica_path)

    if plot_process:
        #fig, [axRebin, axNAL, axIter, axICA] = plt.subplots(4, 1, figsize=(16,20))
        """
        if round(resolution)!=69:
            fig = plt.figure(figsize=(32,14), constrained_layout=True)
            gs  = GridSpec(8, 24, figure=fig)
            #plt.subplots_adjust(hspace=0)
            axRebin = fig.add_subplot(gs[:4,:12]) #ica reconstruction
            axIter = fig.add_subplot(gs[4:,:12])   #CIV fitting
            axNAL = fig.add_subplot(gs[:4,12:])
            axICAflux = fig.add_subplot(gs[4:7,12:20])
            axICAerrs = fig.add_subplot(gs[7:,12:20], sharex=axICAflux)
            axCIV = fig.add_subplot(gs[4:,20:])
        """
        #else:
        fig = plt.figure(figsize=(24,24), constrained_layout=True)
        gs  = GridSpec(12, 18, figure=fig, hspace=2.5)
        #plt.subplots_adjust(hspace=0)
        axNAL = fig.add_subplot(gs[:4,:])   #CIV fitting
        axIter = fig.add_subplot(gs[4:8,:])
        #axICAflux = fig.add_subplot(gs[8:11,:12])
        #axICAerrs = fig.add_subplot(gs[11:,:12], sharex=axICAflux)
        axICAflux = fig.add_subplot(gs[8:,:12])
        axCIV = fig.add_subplot(gs[8:,12:])
        #plt.subplots_adjust(hspace=0)
        #plot initial rebin if necessary:
        """
        if round(resolution)!=69:
            axRebin.set_title("Rebin stage")
            axRebin.plot(waveSpec, fluxSpec, "-b", alpha=0.3)
            axRebin.plot(wave, flux, "-b")
            axRebin.plot(waveSpec, fluxSpec, "-k", alpha=0.3)
            axRebin.plot(wave, errs)
        """
        #plot NAL masking
        axNAL.plot(wave, flux, "-k", zorder=1, label="Normalized Spectrum")
        axNAL.scatter(wave[mask_NAL==1], flux[mask_NAL==1], s=20, color="r", zorder=2, label="Input bad pixels")
        axNAL.plot(wave, flux_median_61p, "-g", lw=3, zorder=1, label="Median Spectrum")
        axNAL.scatter(wave[mask_NAL==2], flux[mask_NAL==2], s=30, color="y", zorder=3, label="NAL bad pixels")
        axNAL.set_ylabel("Flux Density (Arbitrary units)", fontsize=32.5)
        axNAL.set_xlabel("Rest Wavelength (Å)", fontsize=32.5)
        axNAL.tick_params(axis='both', which='major', labelsize=30)
        axNAL.set_title("NAL Masking", fontsize=30)
        axNAL.legend(loc="best", prop={"size":25})
        ylims = np.nanpercentile(flux,5), np.nanpercentile(flux,99)
        axNAL.set_ylim(ylims[0],None)
        axNAL.set_xscale("log")
        axNAL.set_yscale("log")
        ticks = np.array([1200,1500,1750,2000,2200,2500,3000,3500,4000,4500,5000,6000,7000,8000])
        ticks_use = ticks[(ticks>min(wave))&(ticks<max(wave))]
        #ticks_use = ticks[(ticks>min(wave))&(ticks<red_expxlim)]
        axNAL.set_xticks(ticks_use)
        axNAL.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axNAL.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axNAL.get_xaxis().set_tick_params(which='minor', size=0)
        axNAL.get_xaxis().set_tick_params(which='minor', width=0)
        plt.minorticks_off()

        #plot iterative masking
        axIter.plot(wave, flux_morph, "-k", zorder=1, label="Normalized Spectrum")
        #axIter.plot(wave_init_ica, flux_init_ica, "-b", lw=3, zorder=1, label="Initial ICA Fit")
        axIter.scatter(wave[mask_witer==1], flux_morph[mask_witer==1], s=25, color="r", zorder=2, label="Input bad pixels")
        axIter.scatter(wave[mask_witer==2], flux_morph[mask_witer==2], s=25, color="y", zorder=2, label="NAL bad pixels")
        axIter.scatter(wave[mask_witer==3], flux_morph[mask_witer==3], s=25, color="m", zorder=2, label="Iteratively masked pixels")
        #axIter.set_ylim(0.8*min(flux_init_ica), 1.2*max(flux_init_ica))
        axIter.set_ylabel("Flux Density (Arbitrary units)", fontsize=32.5)
        axIter.set_xlabel("Rest Wavelength (Å)", fontsize=32.5)
        axIter.set_title("Post morphing + iterative masking", fontsize=30)
        axIter.tick_params(axis='both', which='major', labelsize=30)
        axIter.legend(loc="best", prop={"size":25})
        ylims = np.nanpercentile(flux,5), np.nanpercentile(flux,99)
        axIter.set_ylim(ylims[0],None)
        axIter.set_xscale("log")
        axIter.set_yscale("log")
        ticks = np.array([1200,1500,1750,2000,2200,2500,3000,3500,4000,4500,5000,6000,7000,8000])
        ticks_use = ticks[(ticks>min(wave))&(ticks<max(wave))]
        #ticks_use = ticks[(ticks>min(wave))&(ticks<red_expxlim)]
        axIter.set_xticks(ticks_use)
        axIter.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axIter.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axIter.get_xaxis().set_tick_params(which='minor', size=0)
        axIter.get_xaxis().set_tick_params(which='minor', width=0)
        plt.minorticks_off()

        #axICA
        axICAflux.plot(wave, flux_morph, "-k", alpha=0.7, label="Morphed spectrum")
        axICAflux.plot(wave_ica, flux_ica, "-b", label="Final ICA fit")
        axICAflux.plot(wave, errs_morph, color="g", alpha=0.7, label="Morphed errors")
        #axICAerrs.plot(wave, errs_morph, "-k", alpha=0.7)
        #axICAerrs.set_ylim(np.nanpercentile(errs_morph[((wave>=min(wave_ica))&(wave<=max(wave_ica)))], [1,99]))
        axICAflux_ylim_low = min(min(flux_ica)-np.nanmedian(flux_ica)/5, np.nanpercentile(errs_morph[((wave>=min(wave_ica))&(wave<=max(wave_ica)))], 1))
        #axICAflux.set_ylim(min(flux_ica)-np.nanmedian(flux_ica)/5, max(flux_ica)+np.nanmedian(flux_ica)/5)
        axICAflux.set_xlim(min(wave_ica), max(wave_ica))
        axICAflux.set_ylim(axICAflux_ylim_low, max(flux_ica)+np.nanmedian(flux_ica)/5)
        axICAflux.set_ylabel("Flux Density (Arbitrary units)", fontsize=30)
        axICAflux.set_xlabel("Rest Wavelength (Å)", fontsize=30)
        axICAflux.set_title("Final ICA with morphing", fontsize=30)
        axICAflux.legend(loc="best", prop={"size":22.5})
        axICAflux.tick_params(axis='both', which='major', labelsize=30)
        get_CIV(wave, flux_morph, wave_ica, flux_ica, name, ax=axCIV)

    return wave, flux_morph, errs_morph, mask_witer, wave_ica, flux_ica, f2500_ica #, morph_coeff #return flux with Arbitrary units as well for plotting
