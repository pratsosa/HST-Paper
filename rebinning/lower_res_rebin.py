"""
rebinning/lower_res_rebin.py
Migrated from: Trevor Code/LowerResHSTRebin_TVM.py
Migration date: 2026-04-30

Import/path changes:
  - `fits.open('spec-0266-51630-0080.fits')` updated to use a path relative to
    this file's location: Data/spec-0266-51630-0080.fits at the repo root.
No algorithmic changes.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from astropy.io import fits
from scipy import stats
from astropy.table import Table

# Path to the reference SDSS spectrum (ships with the repo in Data/)
_SDSS_REF_SPEC = Path(__file__).parent.parent / "Data" / "spec-0266-51630-0080.fits"


def arg2min(arr):
    #return index of second closest point
    return np.argsort(arr)[1]


def HSTLowResRebin(wavelength, flux, flux_error, masks, Identifier, z):
    HSTresolutionlower = 1

    ShortID = Identifier[0:9]

    hdulist = fits.open(str(_SDSS_REF_SPEC))
    sdss_c0 = hdulist[0].header['coeff0']   #central loglam of first pixel
    sdss_c1 = hdulist[0].header['coeff1']   #loglam spacing between pixels
    sdss_npix = hdulist[1].header['NAXIS2']
    sdss_loglam = sdss_c0 + sdss_c1 * np.arange(sdss_npix)
    sdss_wave = 10.**(sdss_loglam)

    #Create empty wavelength array with all possible wavelengths we might deal with
    wave_empty = np.arange(1000., 10000.+sdss_c1, sdss_c1)

    #Number of SDSS pixels separating sdss_c0 and start/end of HST spectrum
    if (wavelength==0).all():
        print(Identifier)
    npix_from_startHST = (min(np.log10(wavelength[wavelength!=0])) - sdss_c0) // sdss_c1
    npix_from_endHST   = (max(np.log10(wavelength[wavelength!=0])) - sdss_c0) // sdss_c1

    #Now define where starting/ending HST wave lies in SDSS
    #TVM changed these lines on 09/13/21
    loglam_minHST_rebin = sdss_c0 + npix_from_startHST*sdss_c1
    loglam_maxHST_rebin = sdss_c0 + npix_from_endHST*sdss_c1
    loglam_HST_rebin    = np.arange(loglam_minHST_rebin, loglam_maxHST_rebin+sdss_c1, sdss_c1)

    #Get "edges" of new loglam pixels
    loglam_HST_rebin_edges = loglam_HST_rebin-0.5*sdss_c1
    loglam_HST_rebin_edges = np.append(loglam_HST_rebin_edges, loglam_HST_rebin_edges[-1]+sdss_c1) #add final edge

    #Linear HST wavelength on SDSS scale
    wave_HST_rebin = 10.**(loglam_HST_rebin)
    wave_HST_rebin_edges = 10.**(loglam_HST_rebin_edges)

    #Now want to compute HST flux at each new rebinned wavelength pixel - initialize first
    flux_HST_rebin    = np.zeros(len(wave_HST_rebin))
    fluxerr_HST_rebin = np.zeros(len(wave_HST_rebin))
    masks_HST_rebin   = np.zeros(len(wave_HST_rebin))

    #Get original by-pixel spacing for HST spectrum
    lambdaSpacingSTIS = np.zeros(len(wavelength)-1)
    for i in range(len(wavelength)-1):
        lambdaSpacingSTIS[i] = wavelength[i+1] - wavelength[i]

    #Get the edges of pixel "bins" in HST spectrum
    lambdaEdgesSTIS = np.zeros(len(wavelength)+1)
    lambdaEdgesSTIS[0] = wavelength[0] - 0.5*lambdaSpacingSTIS[0]
    for i in range(1, len(wavelength)): lambdaEdgesSTIS[i] = wavelength[i-1] + 0.5*lambdaSpacingSTIS[i-1]
    lambdaEdgesSTIS[-1] = wavelength[-1] + 0.5*lambdaSpacingSTIS[-1]

    #######
    lambdaSpacingMedianSTIS = np.median(lambdaSpacingSTIS)

    #ABR:
    #Now we'll loop through each of the pixels in the new grid and see which original pixels overlap with it
    #(either completely or partially).  Then we'll compute the weighted sums of the flux densities and the
    #quadrature sum of the flux density errors.
    #TVM: no, do this instead; why would we use summed errors in quad??? No.
    for i in range(len(wave_HST_rebin)):
        arg1 = np.argmin(np.abs(wavelength-wave_HST_rebin[i]))
        arg2 = arg2min(np.abs(wavelength-wave_HST_rebin[i]))
        wave1, flux1, fluxerr1 = wavelength[arg1], flux[arg1], flux_error[arg1]
        wave2, flux2, fluxerr2 = wavelength[arg2], flux[arg2], flux_error[arg2]

        #Since there are effectively fewer photons captured in each
        #wavelength bin, the S/N should decrease by a factor of sqrt(old_width/new_width)
        old_width = lambdaEdgesSTIS[arg1+1] - lambdaEdgesSTIS[arg1] #old pixel width is that of nearest HST pixel
        new_width = wave_HST_rebin_edges[i+1] - wave_HST_rebin_edges[i] #new pixel width is according to the SDSS scale

        #Do linear fit on these two points
        m_flux, b_flux = np.polyfit([wave1, wave2], [flux1, flux2], 1)
        m_err, b_err   = np.polyfit([wave1, wave2], [fluxerr1, fluxerr2], 1)

        flux_HST_rebin[i]    = m_flux*wave_HST_rebin[i] + b_flux
        fluxerr_HST_rebin[i] = (m_err*wave_HST_rebin[i] + b_err) * np.sqrt(old_width/new_width)
        masks_HST_rebin[i]   = masks[arg1]

    return wave_HST_rebin, flux_HST_rebin, fluxerr_HST_rebin, masks_HST_rebin
