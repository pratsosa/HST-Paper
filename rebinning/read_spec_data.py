"""
rebinning/read_spec_data.py
Migrated from: Trevor Code/Read_spec_data_2026AP.py
Migration date: 2026-04-30

Import/path changes:
  - `import Cut_Edge_Pix_TVM_NoDQ as Cut_Edge_Pix_TVM`
      → `from rebinning import cut_edge_pix as Cut_Edge_Pix_TVM`
  - `import SpecCuts_HSLA`
      → `from rebinning import spec_cuts_hsla as SpecCuts_HSLA`
  - `spec = fits.open(path+"original\\%s"%Identifier)` in read_hsla()
      → `fits.open(os.path.join(path, "original", Identifier))`
      (backslash replaced with os.path.join for cross-platform portability)
No algorithmic changes.

Original file header (Read_spec_data_2026AP.py):
  AP: 2026-03-08 — 2026AP iteration
  This is Read_spec_data_2026E.py with one change:
    Imports Cut_Edge_Pix_TVM_NoDQ instead of Cut_Edge_Pix_TVM_DQOnly.
    The DQ==0 quality cut is removed entirely.
  All other logic is identical to 2026E:
    - FOS multi-row error formula: sqrt(sum(e**2)) / (N-1)  [hypothesized Aug 2022]
    - No SNR edge cut
    - SpecCuts wavelength edge cuts still applied
"""

import os
import pandas as pd
import numpy as np
from astropy.io import fits
import glob

from rebinning import cut_edge_pix as Cut_Edge_Pix_TVM
from rebinning import spec_cuts_hsla as SpecCuts_HSLA


def read_data(Identifier, path, data_origin, z):
    if data_origin == "FOS":
        return read_fos(Identifier, path, z)
    elif data_origin == "STIS":
        return read_stis(Identifier, path, z)
    elif data_origin == "COS":
        return read_cos(Identifier, path, z)
    elif data_origin == "HSLA":
        return read_hsla(Identifier, path, z)
    elif data_origin == "SDSS-RM":
        return read_sdssrm(Identifier, path, z)
    print("data_origin not recognized")


def read_sdssrm(Identifier, path, z):
    fn_list = glob.glob(path+"%s/*.fits"%Identifier)
    array_lens = []
    for i in range(len(fn_list)):
        wavelength = 10.**fits.open(fn_list[i])[1].data["LOGLAM"]
        array_lens.append(len(wavelength))
    array_len = max(array_lens)

    waves  = np.zeros((len(fn_list), array_len))
    fluxes = np.zeros((len(fn_list), array_len))
    errs   = np.zeros((len(fn_list), array_len))
    masks  = np.zeros((len(fn_list), array_len))

    for i in range(len(fn_list)):
        hdu = fits.open(fn_list[i])
        loglam = hdu[1].data["LOGLAM"]
        wave   = 10.**loglam
        flux   = hdu[1].data["FLUX"]
        err    = 1. / np.sqrt(hdu[1].data["IVAR"])
        mask   = hdu[1].data["AND_MASK"]
        waves[i,:len(wave)]  = wave
        fluxes[i,:len(wave)] = flux
        errs[i,:len(wave)]   = err
        masks[i,:len(wave)]  = mask

    return waves, fluxes, errs, masks


def read_hsla(Identifier, path, z):
    spec = fits.open(os.path.join(path, "original", Identifier))
    coadd_wave = spec[1].data["WAVE"]
    coadd_flux = spec[1].data["FLUX"]
    coadd_errs = spec[1].data["ERROR"]
    coadd_mask = np.zeros(len(coadd_wave))

    #Manually cut red end - big problem for HSLA co-adds
    obj_name = "_".join(Identifier.split("_")[:-1])
    if obj_name in SpecCuts_HSLA.RedEdges:
        lambdaend = SpecCuts_HSLA.RedEdges[obj_name]
        indstop  = np.argmin( np.abs((coadd_wave/(1+z))-lambdaend) )
        coadd_wave = coadd_wave[:indstop]
        coadd_flux = coadd_flux[:indstop]
        coadd_errs = coadd_errs[:indstop]
        coadd_mask = coadd_mask[:indstop]

    coadd_wave = Cut_Edge_Pix_TVM.Cut_Edge_Pix(np.zeros(len(coadd_wave)), coadd_wave, coadd_flux, coadd_errs, coadd_wave.copy(), len(coadd_wave), "wavelength", False, z, "%s"%(Identifier), "HSLA")
    coadd_flux = Cut_Edge_Pix_TVM.Cut_Edge_Pix(np.zeros(len(coadd_wave)), coadd_wave, coadd_flux, coadd_errs, coadd_flux.copy(), len(coadd_flux), "flux", False, z, "%s"%(Identifier), "HSLA")
    #HSLA saves empty values as zero; change to nan
    coadd_flux[coadd_flux==0] = np.nan
    coadd_errs = Cut_Edge_Pix_TVM.Cut_Edge_Pix(np.zeros(len(coadd_wave)), coadd_wave, coadd_flux, coadd_errs, coadd_errs.copy(), len(coadd_errs), "flux error", False, z, "%s"%(Identifier), "HSLA")
    coadd_mask = Cut_Edge_Pix_TVM.Cut_Edge_Pix(np.zeros(len(coadd_wave)), coadd_wave, coadd_flux, coadd_errs, coadd_mask.copy(), len(coadd_mask), "masks", False, z, "%s"%(Identifier), "HSLA")

    return np.array([coadd_wave]), np.array([coadd_flux]), np.array([coadd_errs]), np.array([coadd_mask])


def read_cos(Identifier, path, z):
    #get observation details
    try:
        obs_details = pd.read_csv(path+"%s/all_exposures.txt" %(Identifier), sep="\s+")
        gratings    = obs_details["Grating"].values
        spec_names  = obs_details["Rootname"].values
    except FileNotFoundError:
        obs_details = pd.read_csv(path+"%s/NecessaryParams.csv" %(Identifier))
        gratings    = obs_details["filters"].values
        spec_names  = obs_details["obs_id"].values

    #take max wave size for initializing below
    array_sizes = []
    for i in range(len(spec_names)):
        try:
            data = fits.open(path+'%s/%s/%s_x1d.fits' %(Identifier,spec_names[i],spec_names[i]))
        except FileNotFoundError:
            try:
                data = fits.open(path+'%s/%s.fits' %(Identifier,spec_names[i]))
            except FileNotFoundError:
                data = fits.open(path+'%s/Data/%s_x1d.fits' %(Identifier,spec_names[i]))
        data = data[1].data
        if gratings[i]=='E140M':
            array_sizes.append(data.size*1024)
        else:
            array_sizes.append(len(data['Wavelength'][0]))
    array_len = max(array_sizes)

    waves     = np.zeros((len(spec_names), array_len))
    fluxes    = np.zeros((len(spec_names), array_len))
    flux_errs = np.zeros((len(spec_names), array_len))
    masks     = np.zeros((len(spec_names), array_len))

    for i in range(len(spec_names)):
        Bad_list = []
        try:
            data = fits.open(path+'%s/%s/%s_x1d.fits' %(Identifier,spec_names[i],spec_names[i]))[1].data
        except FileNotFoundError:
            try:
                data = fits.open(path+'%s/%s.fits' %(Identifier,spec_names[i]))[1].data
            except FileNotFoundError:
                data = fits.open(path+'%s/Data/%s_x1d.fits' %(Identifier,spec_names[i]))[1].data

        if gratings[i]=="E140M":
            wavelength = []
            flux       = []
            fluxerr    = []
            DQ         = []
            for t in np.arange(0,data.size,1):
                wavelength = np.append(wavelength, data['WAVELENGTH'][t])
                flux       = np.append(flux, data['FLUX'][t])
                fluxerr    = np.append(fluxerr, data['ERROR'][t])
                DQ         = np.append(DQ, data['DQ'][t])
        else:
            wavelength = data['WAVELENGTH'][0]
            flux       = data['FLUX'][0]
            fluxerr    = data['ERROR'][0]
            DQ         = data['DQ'][0]

        flux_wmask = flux.copy()
        err_wmask  = fluxerr.copy()
        sel        = ((wavelength >= 1215.0) & (wavelength <= 1216.0))
        masks[i,:][(err_wmask == 0.0)] = 1
        waves[i,:]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                wavelength, array_len, "wavelength", False, z, "%s - %s"%(Identifier,spec_names[i]), "COS")
        fluxes[i,:]    = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                flux_wmask, array_len, "flux", False, z, "%s - %s"%(Identifier,spec_names[i]), "COS")
        flux_errs[i,:] = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                err_wmask, array_len, "flux error", False, z, "%s - %s"%(Identifier,spec_names[i]), "COS")
        masks[i,:]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                masks[i,:], array_len, "masks", False, z, "%s - %s"%(Identifier,spec_names[i]), "COS")

    return waves, fluxes, flux_errs, masks


def read_stis(Identifier, path, z):
    obs_details = pd.read_csv(path+"%s/NecessaryParams.csv" %(Identifier))
    gratings    = obs_details["filters"].values
    spec_names  = obs_details["obs_id"].values

    array_sizes = []
    for i in range(len(spec_names)):
        try:
            data = fits.open(path+'%s/%s/%s_x1d.fits' %(Identifier,spec_names[i],spec_names[i]))
        except FileNotFoundError:
            data = fits.open(path+'%s/%s/%s_sx1.fits' %(Identifier,spec_names[i],spec_names[i]))
        data = data[1].data
        if gratings[i]=='E140M':
            array_sizes.append(data.size*1024)
        else:
            array_sizes.append(len(data['Wavelength'][0]))
    array_len = max(array_sizes)

    waves     = np.zeros((len(spec_names), array_len))
    fluxes    = np.zeros((len(spec_names), array_len))
    flux_errs = np.zeros((len(spec_names), array_len))
    masks     = np.zeros((len(spec_names), array_len))

    for i in range(len(spec_names)):
        try:
            data = fits.open(path+'%s/%s/%s_x1d.fits' %(Identifier,spec_names[i],spec_names[i]))[1].data
        except FileNotFoundError:
            data = fits.open(path+'%s/%s/%s_sx1.fits' %(Identifier,spec_names[i],spec_names[i]))[1].data
        if gratings[i]=="E140M":
            wavelength = []
            flux       = []
            fluxerr    = []
            DQ         = []
            for t in np.arange(0,data.size,1):
                wavelength = np.append(wavelength, data['WAVELENGTH'][t])
                flux       = np.append(flux, data['FLUX'][t])
                fluxerr    = np.append(fluxerr, data['ERROR'][t])
                DQ         = np.append(DQ, data['DQ'][t])
        else:
            wavelength = data['WAVELENGTH'][0]
            flux       = data['FLUX'][0]
            fluxerr    = data['ERROR'][0]
            DQ         = data['DQ'][0]

        flux_wmask = flux.copy()
        err_wmask  = fluxerr.copy()
        sel        = ((wavelength >= 1215.0) & (wavelength <= 1216.0))
        masks[i,:len(err_wmask)][err_wmask==0.] = 1
        waves[i,:len(err_wmask)]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                wavelength, array_len, "wavelength", False, z, "%s - %s"%(Identifier,spec_names[i]), "STIS")
        fluxes[i,:len(err_wmask)]    = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                flux_wmask, array_len, "flux", False, z, "%s - %s"%(Identifier,spec_names[i]), "STIS")
        flux_errs[i,:len(err_wmask)] = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                err_wmask, array_len, "flux error", False, z, "%s - %s"%(Identifier,spec_names[i]), "STIS")
        masks[i,:len(err_wmask)]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                masks[i,:len(err_wmask)], array_len, "masks", False, z, "%s - %s"%(Identifier,spec_names[i]), "STIS")

    return waves, fluxes, flux_errs, masks


def read_fos(Identifier, path, z):
    # AP: HST original read_fos — no try-except, no bad_indices removal.

    obs_details = pd.read_csv(path+"%s/NecessaryParams.csv"%Identifier)
    gratings    = obs_details["filters"].values
    spec_names  = obs_details["obs_id"].values

    array_lens = np.array([], dtype=int)
    for spectrum in spec_names:
        wave = fits.open(path+"%s/%s/%s_c0f.fits" % (Identifier, spectrum, spectrum))[0].data
        array_lens = np.append(array_lens, wave.shape[-1])
    array_len = max(array_lens)

    waves     = np.zeros((len(spec_names), array_len))
    fluxes    = np.zeros((len(spec_names), array_len))
    flux_errs = np.zeros((len(spec_names), array_len))
    masks     = np.zeros((len(spec_names), array_len))

    def accum_flag(index):
        return np.mean(flux[index]) > np.mean(flux[index-1]) and np.mean(flux[index-1]) > np.mean(flux[index-2])

    for i in range(len(spec_names)):
        wave = fits.open(path+"%s/%s/%s_c0f.fits" % (Identifier, spec_names[i], spec_names[i]))[0].data

        ind = -1

        if len(wave.shape) > 1:
            nzero = (wave[ind]>0.)
            if (wave[ind][nzero][0] > wave[ind][nzero][-1]):
                wavelength = wave[ind][::-1]
            else:
                wavelength = wave[ind][:]
        else:
            nzero = (wave>0.)
            if (wave[nzero][0] > wave[nzero][-1]):
                wavelength = wave[::-1]
            else:
                wavelength = wave[:]

        flux    = np.zeros(len(wavelength))
        fluxerr = np.zeros(len(wavelength))
        DQ      = np.zeros(len(wavelength))

        obs_flux    = fits.open(path+'%s/%s/%s_c1f.fits' %(Identifier,spec_names[i],spec_names[i]))[0].data
        obs_fluxerr = fits.open(path+'%s/%s/%s_c2f.fits' %(Identifier,spec_names[i],spec_names[i]))[0].data
        obs_DQ      = fits.open(path+'%s/%s/%s_cqf.fits' %(Identifier,spec_names[i],spec_names[i]))[0].data

        if len(wave.shape) > 1:
            if accum_flag(ind):
                if (wave[ind][nzero][0] > wave[ind][nzero][-1]):
                   flux    = obs_flux[ind][::-1]
                   fluxerr = obs_fluxerr[ind][::-1]
                   DQ      = obs_DQ[ind][::-1]
                else:
                   flux    = obs_flux[ind]
                   fluxerr = obs_fluxerr[ind]
                   DQ      = obs_DQ[ind]
            else:
                for l in range(len(wavelength)):
                    if (wave[ind][nzero][0] > wave[ind][nzero][-1]):
                        flux[l]    = np.mean(obs_flux[:,-l-1])
                        fluxerr[l] = np.sqrt(sum(obs_fluxerr[:,-l-1]**2)) / (obs_fluxerr.shape[0]-1)  # AP: restored /N-1 (hypothesized Aug 2022)
                        DQ[l]      = sum(obs_DQ[:,-l-1])
                    else:
                        flux[l]    = np.mean(obs_flux[:,l])
                        fluxerr[l] = np.sqrt(sum(obs_fluxerr[:,l]**2)) / (obs_fluxerr.shape[0]-1)  # AP: restored /N-1 (hypothesized Aug 2022)
                        DQ[l]      = sum(obs_DQ[:,l])
        else:
            if (wave[nzero][0] > wave[nzero][-1]):
                flux    = obs_flux[::-1]
                fluxerr = obs_fluxerr[::-1]
                DQ      = obs_DQ[::-1]
            else:
                flux    = obs_flux
                fluxerr = obs_fluxerr
                DQ      = obs_DQ

        sel                       = (wavelength >= 1205) & (wavelength <= 1225)
        mask                      = DQ>0
        flux_wmask                = flux.copy()
        err_wmask                 = fluxerr.copy()
        masks[i,:len(err_wmask)][err_wmask==0.] = 1

        waves[i,:len(flux)]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                wavelength, array_len, "wavelength", False, z, "%s - %s"%(Identifier,spec_names[i]),"FOS")
        fluxes[i,:len(flux)]    = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                flux_wmask, array_len, "flux", False, z, "%s - %s"%(Identifier,spec_names[i]),"FOS")
        flux_errs[i,:len(flux)] = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                err_wmask, array_len, "flux error", False, z, "%s - %s"%(Identifier,spec_names[i]),"FOS")
        masks[i,:len(flux)]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask, \
                                                masks[i,:len(err_wmask)], array_len, "masks", False, z, "%s - %s"%(Identifier,spec_names[i]),"FOS")

    return waves, fluxes, flux_errs, masks
