"""
rebinning/small_pix_filter.py
Migrated from: Trevor Code/Small_Pix_Filter_TVM.py
Migration date: 2026-04-30
No algorithmic changes — only header added.

Original filename: 11_Pix_Filter.py
"""

import math
import numpy as np
import scipy
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt


#100 is necessary to fit the G130M images better...
def pixel_filter(wavelength, flux, fluxerr, Identifier, Instrument, Cont_Type, PLOT=False):
    if Instrument == "COS":
       filterval = 5
       filterval2 = 6
    elif Instrument == "HSLA":
       filterval = 20
       filterval2 = 21
    elif Instrument == "STIS":
       filterval = 5
       filterval2 = 6
    elif Instrument == "FOS":
       filterval = 30
       filterval2 = 31
    elif Instrument == "GHRS":
       filterval = 5
       filterval2 = 6
    elif Instrument == "SDSS-RM":
       filterval = 60
       filterval2 = 61

    indices = ~np.isnan(flux) #indices we'll look at

   # AP: Updated defensive check - mirror construction needs at least 2 valid pixels.
   # If n_valid <= 1, return NaNs so pathological exposures are skipped gracefully.
   # Otherwise shrink filter width when valid coverage is shorter than default window.
    n_valid = indices.sum()
    if n_valid <= 1:
       return np.full_like(flux, np.nan)
    if n_valid < filterval + 1:
       filterval = n_valid - 1
       filterval2 = filterval + 1

    if Cont_Type == "S/N":
       #Horizontally flip the fluxes:
       end        = np.flipud(flux[indices][-filterval2-1:])/np.flipud(fluxerr[indices][-filterval2-1:])
       end_wave   = wavelength[indices][-filterval2-1:]
       start      = (flux[indices][:filterval+1])/(fluxerr[indices][:filterval+1])
       start_wave = wavelength[indices][:filterval+1]

    elif Cont_Type == "Fluxerrs":
       end        = fluxerr[indices][-filterval2-1:]
       end_wave   = wavelength[indices][-filterval2-1:]
       start      = fluxerr[indices][:filterval+1]
       start_wave = wavelength[indices][:filterval+1]

    elif Cont_Type == "Fluxes":
       end        = flux[indices][-filterval2-1:]
       end_wave   = wavelength[indices][-filterval2-1:]
       start      = flux[indices][:filterval+1]
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

    if Cont_Type == "S/N":
       extended_flux_array = np.concatenate((new_left_flux, flux/fluxerr, new_right_flux), axis=0)
    elif Cont_Type == "Fluxerrs":
       extended_flux_array = np.concatenate((new_left_flux, fluxerr, new_right_flux), axis=0)
    elif Cont_Type == "Fluxes":
       extended_flux_array = np.concatenate((new_left_flux, flux, new_right_flux), axis=0)

    continuum = [np.nanmedian(extended_flux_array[l-filterval:l+filterval]) if 1 else 0 for l in np.arange(filterval, len(flux)+filterval, 1)]
    continuum = np.asarray(continuum)

    return continuum


def SDSS_pixel_filter(wavelength, data, npix=20):

    filterval  = npix//2
    filterval2 = npix//2+1

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

    extended_data_array = np.concatenate((new_left_flux, data, new_right_flux), axis=0)
    continuum = np.array([np.nanmedian(extended_data_array[l-filterval:l+filterval]) for l in np.arange(filterval, len(data)+filterval, 1)])

    return continuum
