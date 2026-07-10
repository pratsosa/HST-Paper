"""
rebinning/coadd.py
Migrated from: Trevor Code/HST_coadd2_2026AP.py
Migration date: 2026-04-30

Import/path changes:
  - `import Small_Pix_Filter_TVM`
      → `from rebinning import small_pix_filter as Small_Pix_Filter_TVM`
  - `import LowerResHSTRebin_TVM`
      → `from rebinning import lower_res_rebin as LowerResHSTRebin_TVM`
  - `import Read_spec_data_2026AP as Read_spec_data`
      → `from rebinning import read_spec_data`
  - `import spec_morph`
      → `from rebinning import spec_morph`
  - `import Large_pix_filter_TVM` — commented out (unused; not called anywhere in this file)
  - `import Plot_HST` — commented out (unused; not called anywhere in this file)
  - `import Cut_Edge_Pix_TVM_NoDQ as Cut_Edge_Pix_TVM` — commented out (not used
      directly in this file; used internally by read_spec_data)
  - Instrument-specific imports inside rebin() (`SpecCuts_*`, `MorphingEdges_*`)
      updated to use package-relative imports. MorphingEdges imports are commented
      out — they are vestigial (imported but never called in this file); they will
      be re-added when needed in Phase 2.

Path changes:
  - `fits.open('spec-0266-51630-0080.fits')` → path via _SDSS_REF_SPEC (Data/)
  - `fits.open("SDSS_spec/lite/%s" % fn_sdss)` → `fits.open(os.path.join(sdss_spec_dir, fn_sdss))`
  - `t.write("RebinnedSpec_2026AP/...")` → `t.write(os.path.join(output_dir, ...))`

Signature change:
  - `rebin(Identifier, z, data_origin, fn_sdss, data_path=None)` gains two new
    keyword arguments:
      output_dir  (str): directory to write output FITS files. Default "RebinnedSpec".
      sdss_spec_dir (str|None): root directory containing SDSS spec files
        (e.g. the folder that holds plate subdirectories). Required when fn_sdss
        is not None. Default None.

Original file header (HST_coadd2_2026AP.py):
  AP: 2026-03-08 — DQ==0 quality cut removed entirely (uses Cut_Edge_Pix_TVM_NoDQ
  via Read_spec_data_2026AP). All other settings identical to 2026E:
    - FOS multi-row error formula: sqrt(sum(e**2)) / (N-1)
    - SpecCuts wavelength edge cuts still applied
  Output goes to RebinnedSpec_2026AP/.
  Purpose: quantify the effect of DQ flagging on rebinned spectra.
"""

import os
import matplotlib
import weightedstats as ws
import decimal
import time
import math
import scipy
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from astropy.io import fits
from astropy.table import Table

from rebinning import small_pix_filter as Small_Pix_Filter_TVM
from rebinning import lower_res_rebin as LowerResHSTRebin_TVM
from rebinning import read_spec_data
from rebinning import spec_morph
# import Large_pix_filter_TVM  # unused in this file
# import Plot_HST               # unused in this file
# import Cut_Edge_Pix_TVM_NoDQ  # not used directly here; used inside read_spec_data

# Path to the reference SDSS spectrum (ships with the repo in Data/)
_SDSS_REF_SPEC = Path(__file__).parent.parent / "Data" / "spec-0266-51630-0080.fits"


def rebin(Identifier, z, data_origin, fn_sdss, data_path=None,
          output_dir="RebinnedSpec", sdss_spec_dir=None, flat=False):
    #Identifier: string; name of object; helps you find files
    #z: redshift; float; recorded by TVM (i.e. not updated by Paul)
    #data_origin: string; HST instrument, HSLA, Gordon/Angelica DR7 targets, SDSS-RM, GNIRS-DQS
    #fn_sdss: string; spec-PPPP-MMMMM-FFFF.fits plate/mjd/fiber SDSS name
    #output_dir: string; directory to write rebinned FITS output
    #sdss_spec_dir: string; root directory containing SDSS spec files (required if fn_sdss is not None)
    #flat: bool; if True, use flat-layout readers (MAST download structure, no NecessaryParams.csv)

    #Initial Instrument-dependent imports
    if data_origin == "FOS":
        from rebinning import spec_cuts_fos as SpecCuts
        # MorphingEdges_FOS — vestigial import; not called in this file (Phase 2)
    elif data_origin == "STIS":
        from rebinning import spec_cuts_stis as SpecCuts
        # MorphingEdges_STIS — vestigial import; not called in this file (Phase 2)
    elif data_origin == "COS":
        from rebinning import spec_cuts_cos as SpecCuts
        # MorphingEdges_COS — vestigial import; not called in this file (Phase 2)
    elif data_origin == "HSLA":
        from rebinning import spec_cuts_hsla as SpecCuts
        # MorphingEdges_HSLA — vestigial import; not called in this file (Phase 2)

    #read the data should be the same as before; only free parameter is the instrument/origin of data
    if flat:
        waves, fluxes, flux_errs, masks = read_spec_data.read_data_flat(Identifier, data_path, data_origin, z)
    else:
        waves, fluxes, flux_errs, masks = read_spec_data.read_data(Identifier, data_path, data_origin, z)

    """
    Step 1: Get a smoothed "median" flux and fluxerr spectrum.
            Determine noise spikes (and remove) by checking if
            a given pixel is N-sigma from the smoothed version.
    """
    for i in range(waves.shape[0]):
        median_flux = Small_Pix_Filter_TVM.pixel_filter(waves[i,:], fluxes[i,:], flux_errs[i,:], Identifier,
                                                  data_origin, "Fluxes")
        median_err  = Small_Pix_Filter_TVM.pixel_filter(waves[i,:], fluxes[i,:], flux_errs[i,:], Identifier,
                                                  data_origin, "Fluxerrs")

        em = spec_morph.emission_lines2(waves[i,:]/(1+z))
        for j in range(len(median_flux)):
            if ( (abs(fluxes[i,j]-median_flux[j]) > 5.*abs(median_err[j])) and ~em[j] ): #change condition as needed
                fluxes[i,j]    = median_flux[j]
                flux_errs[i,j] = median_err[j]
                masks[i,j]     = 2

    '''
    Step 2: Rebin the HST wave/flux/fluxerr/masks to the (often lower-)resolution
            SDSS scale.

    Start with an "empty" wave array covering pixels the spectrum *would* have
    if it were observed with HST resolution.  We want this to truly cover the full
    *potential* wavelength range so that co-adding from different observations
    later on is relatively seamless.  Then can insert fluxes at each pixel.

    c0: central value of first wavelength pixel
    c1: spacing from pixel to pixel; for SDSS c1 = 0.0001 loglam (constant in logarithmic space)
    '''
    #Generic SDSS spectrum to get the wavelength scale
    hdulist = fits.open(str(_SDSS_REF_SPEC))
    sdss_c0 = hdulist[0].header['coeff0'] #c0 is the central log10(wavelength) of the first pixel
    sdss_c1 = hdulist[0].header['coeff1'] #c1 is the log10 dispersion per pixel
    sdss_npix = hdulist[1].header['NAXIS2']
    sdss_loglam = sdss_c0 + sdss_c1 * np.arange(sdss_npix)
    sdss_wave = 10.**(sdss_loglam)

    #Number of SDSS pixels separating sdss_c0 and start/end of HST spectrum
    npix_from_startHST = (min(np.log10(waves[waves!=0])) - sdss_c0) // sdss_c1
    npix_from_endHST   = (max(np.log10(waves[waves!=0])) - sdss_c0) // sdss_c1
    #Now define where starting/ending HST wave lies in SDSS
    loglam_minHST_rebin = sdss_c0 + npix_from_startHST*sdss_c1
    loglam_maxHST_rebin = sdss_c0 + npix_from_endHST*sdss_c1
    loglam_HST_rebin    = np.arange(loglam_minHST_rebin, loglam_maxHST_rebin+sdss_c1, sdss_c1)
    hst_wave_new = 10.**(loglam_HST_rebin)

    #Now begin rebinning data - start by initializing
    old_binned_waves  = np.zeros((waves.shape[0], len(hst_wave_new)))
    old_binned_fluxes = np.zeros((waves.shape[0], len(hst_wave_new)))*np.nan
    old_binned_errs   = np.zeros((waves.shape[0], len(hst_wave_new)))
    old_binned_masks  = np.zeros((waves.shape[0], len(hst_wave_new)))

    f1450_list = np.array([])
    f1450_wts  = np.array([])
    f1450_s2n  = np.array([])
    for i in range(waves.shape[0]):
        #get re-binned arrays for individual exposure
        new_binned_waves, new_binned_fluxes, new_binned_errs, new_binned_masks = \
        LowerResHSTRebin_TVM.HSTLowResRebin(waves[i,:],fluxes[i,:],flux_errs[i,:],masks[i,:],Identifier,z)
        #and morph each individual exposure
        continuum = spec_morph.cont_filtered(new_binned_waves, new_binned_fluxes, z, Identifier)

        #if 1450Å in this exposure, add to list
        if min(new_binned_waves)/(1+z)<1450 and max(new_binned_waves)/(1+z)>1450:
            ind_f1450 = np.abs(new_binned_waves/(1+z) - 1450).argmin()
            f1450_list = np.append(f1450_list, np.nanmedian(continuum[new_binned_errs!=0][ind_f1450-5:ind_f1450+5]))
            f1450_wts  = np.append(f1450_wts, np.nanmedian(1/new_binned_errs[new_binned_errs!=0][ind_f1450-5:ind_f1450+5]))
            f1450_s2n  = np.append(f1450_s2n, np.nanmedian(continuum[new_binned_errs!=0][ind_f1450-5:ind_f1450+5]/new_binned_errs[new_binned_errs!=0][ind_f1450-5:ind_f1450+5]))

        new_binned_fluxes /= continuum
        new_binned_errs /= continuum

        #Save arrays with wavelength "in place"
        argstart = max(0, np.argmin(np.abs(new_binned_waves[0]-hst_wave_new))-1) #avoid ValueError below
        old_binned_waves[i,argstart:argstart+len(new_binned_waves)]   = new_binned_waves
        old_binned_fluxes[i,argstart:argstart+len(new_binned_fluxes)] = new_binned_fluxes
        old_binned_errs[i,argstart:argstart+len(new_binned_errs)]     = new_binned_errs
        old_binned_masks[i,argstart:argstart+len(new_binned_masks)]   = new_binned_masks

    """
    Step 3: Now that data are rebinned, co-add them
            to get the variance-weighted median spectrum.
    """
    if data_origin=="FOS" or data_origin=="STIS" or data_origin=="COS" or data_origin=="SDSS-RM":
        #Initialize - note we want weights for each exposure
        varweighted_flux = np.zeros(old_binned_waves.shape[1])*np.nan
        varweighted_errs = np.zeros(old_binned_waves.shape[1])
        varweighted_mask = np.zeros(old_binned_waves.shape[1])
        total_variance = np.zeros(old_binned_waves.shape[1])
        my_weights     = np.ones(old_binned_waves.shape)

        #compute weights for each pixel in each exposure
        for i in range(len(total_variance)):
            total_variance[i] = np.nansum(old_binned_errs[:,i]**2.)
            if total_variance[i] != 0:
                my_weights[:,i] = 1 / old_binned_errs[:,i]**2.  #should be inverse here - otherwise noisier spectra contribute more
                my_weights[np.isnan(my_weights[:,i]), i] = 0. #some errors are saved as nan - change weight to 0 to feed to weightedstats
                my_weights[np.isinf(my_weights[:,i]), i] = 0.
        #don't include masked points in the final co-add
        em = spec_morph.emission_lines2(old_binned_waves/(1+z))
        my_weights[(old_binned_masks>0)&(~em)] = 0.

        #the masking portion need not be co-added; just mask final pixel if that pixel from all exposures is masked
        for i in range(old_binned_waves.shape[1]):
            varweighted_flux[i] = ws.numpy_weighted_median(old_binned_fluxes[:,i], weights=my_weights[:,i])
            varweighted_errs[i] = ws.numpy_weighted_median(old_binned_errs[:,i], weights=my_weights[:,i])

            if np.isnan(varweighted_flux[i]):
                goodpix   = ( (~np.isnan(old_binned_fluxes[:,i])) & (old_binned_fluxes[:,i]!=0) )
                new_flux = old_binned_fluxes[goodpix,i]
                new_errs = old_binned_errs[goodpix,i]
                new_masks = old_binned_masks[goodpix,i]
                new_weights = 1 / (new_errs**2.)
                new_weights[((np.isinf(new_weights))|(np.isnan(new_weights)))] = 0.

                if np.nansum(new_errs) == 0.:
                    varweighted_flux[i]  = np.nan
                    varweighted_errs[i]  = 0
                    varweighted_mask[i] = 1
                else:
                    varweighted_flux[i]  = ws.numpy_weighted_median(new_flux,weights=new_weights)
                    # AP: Defensive check - weighted_median can return None if weights are all zero/nan/inf
                    med_err = ws.numpy_weighted_median(new_errs,weights=new_weights)
                    if med_err is None:
                        varweighted_errs[i] = 0
                        varweighted_mask[i] = 1
                    else:
                        varweighted_errs[i] = med_err / np.sqrt(len(new_errs))
                        if (new_masks==0).any():
                            varweighted_mask[i] = 0
                        else:
                            varweighted_mask[i] = 1
                    if np.isnan(varweighted_flux[i]):
                        pass

    else:
        #if we're not co-adding (e.g. HSLA), just take the rebinned spectrum
        varweighted_flux = old_binned_fluxes[0].copy()
        varweighted_errs = old_binned_errs[0].copy()
        varweighted_mask = old_binned_masks[0].copy()

    #data is messy blueward of ~1120Å, so cut it there
    lambdacut = 1120
    argstart = np.argmin( np.abs((hst_wave_new/(1+z))-lambdacut) )
    """
    Step 4: Load in SDSS if available.
    """
    if fn_sdss is not None:
        hdulist = fits.open(os.path.join(sdss_spec_dir, fn_sdss))
        sdss_flux    = hdulist[1].data["flux"]*10**-17
        sdss_fluxerr = (1/(np.sqrt(hdulist[1].data["ivar"])+1.e-28))*10**-17 #rest-frame flux uncertainty
        sdss_ANDMask = hdulist[1].data['and_mask'] #ANDmask
        sdss_ANDMask[sdss_fluxerr==0.0]=1

        #wavelength
        sdss_c0     = hdulist[0].header['COEFF0']
        sdss_c1     = hdulist[0].header['COEFF1']
        sdss_npix   = hdulist[1].header['naxis2']
        sdss_loglam = sdss_c0 + sdss_c1 * np.arange(sdss_npix)
        sdss_wave   = (10.**(sdss_loglam))

        #Normalize by continuum to stitch together better, THEN morph
        varweighted_continuum = spec_morph.cont_filtered(hst_wave_new, varweighted_flux, z, Identifier)
        varweighted_flux /= varweighted_continuum
        varweighted_errs /= varweighted_continuum
        sdss_continuum = spec_morph.continuum_fit2(sdss_wave, sdss_flux, z)
        sdss_flux /= sdss_continuum
        sdss_fluxerr /= sdss_continuum

        varweighted_flux, varweighted_morph_factor = spec_morph.morph2(hst_wave_new, varweighted_flux, varweighted_errs, z, Identifier)
        sdss_flux, sdss_morph_factor = spec_morph.morph2(sdss_wave, sdss_flux, sdss_fluxerr, z, Identifier)
        varweighted_errs *= varweighted_morph_factor
        sdss_fluxerr *= sdss_morph_factor
        """
        Step 5a: Now stitch together HST/SDSS and save.
        """
        hst_med_res = np.nanmedian([np.log10(hst_wave_new[i+1])-np.log10(hst_wave_new[i]) for i in range(len(hst_wave_new)-1)])
        sdss_med_res= np.nanmedian([np.log10(sdss_wave[i+1])-np.log10(sdss_wave[i]) for i in range(len(sdss_wave)-1)])
        final_res   = (hst_med_res+sdss_med_res) / 2.

        if min(sdss_wave) < max(hst_wave_new):
            #spectra overlap - cut HST spectrum
            argcut     = np.argmin( np.abs(hst_wave_new-min(sdss_wave)) )
            final_wave = np.concatenate((hst_wave_new[argstart:argcut], sdss_wave))
            final_flux = np.concatenate((varweighted_flux[argstart:argcut], sdss_flux))
            final_errs = np.concatenate((varweighted_errs[argstart:argcut], sdss_fluxerr))
            final_flux, final_morph_factor = spec_morph.morph2(final_wave, final_flux, final_errs, z, Identifier)
            final_errs *= final_morph_factor
            final_mask = np.concatenate((varweighted_mask[argstart:argcut], sdss_ANDMask))
        else:
            #there's a gap between HST and SDSS waves
            gap_wave    = 10.**np.arange(np.log10(max(hst_wave_new))+final_res, np.log10(min(sdss_wave)), final_res)
            gap_nans    = np.ones(len(gap_wave))*np.nan
            gap_zeros   = np.zeros(len(gap_wave))
            gap_ones    = np.ones(len(gap_wave))
            final_wave = np.concatenate((hst_wave_new[argstart:], gap_wave, sdss_wave))
            final_flux = np.concatenate((varweighted_flux[argstart:], gap_nans, sdss_flux))
            final_errs = np.concatenate((varweighted_errs[argstart:], gap_zeros, sdss_fluxerr))

            _, final_morph_factor = spec_morph.morph2(final_wave, final_flux, final_errs, z, Identifier)
            final_mask = np.concatenate((varweighted_mask[argstart:], gap_ones, sdss_ANDMask))
    else:
        """
        Step 5b: Morph variance-weighted-median-continuum-normalized HST spectrum.
        """
        final_wave = hst_wave_new[argstart:]
        varweighted_flux, varweighted_morph_factor = spec_morph.morph2(hst_wave_new, varweighted_flux, varweighted_errs, z, Identifier)
        final_flux = varweighted_flux[argstart:]
        final_errs = varweighted_errs[argstart:] * varweighted_morph_factor[argstart:]
        final_mask = varweighted_mask[argstart:]
    """
    Finally, save!
    """
    t = Table([final_wave/(1+z), final_flux, final_errs, final_mask, z*np.ones(len(final_wave))],
              names=('Rest-frame Wavelength','Coadded Flux (Arbitrary Units)','Coadded Flux Errors','Bad Pixel Mask','Redshift'))

    # AP: For HSLA data, Identifier already contains .fits extension, so strip it before appending _HSLA.fits
    if data_origin == "HSLA":
        Identifier_clean = Identifier.replace('.fits', '') if Identifier.endswith('.fits') else Identifier
        t.write(os.path.join(output_dir, "%s_%s.fits" % (Identifier_clean, data_origin)), overwrite=True)
    else:
        t.write(os.path.join(output_dir, "%s_%s.fits" % (Identifier, data_origin)), overwrite=True)
