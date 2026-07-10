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


# AP 2026-06-17: TEMPORARY exposure cap for COS reverberation-mapping monitoring
# targets. These objects have far more exposures than any other (Mrk 817 ~1600,
# NGC 5548 ~800), which makes the pure-Python coadd in coadd.py appear to hang.
# STOPGAP: keep only the N highest-S/N exposures. See Migration_Log.md
# "TEMP: COS exposure cap" -- REVISIT this decision.
_COS_EXPOSURE_CAP = {"Mrk 817": 100, "NGC 5548": 100}


def _cos_median_snr(fn):
    """Median per-pixel S/N of a COS x1d exposure, for ranking. Flattens all
    rows (segments/stripes) so it works for FUV (2 seg) and NUV (3 stripe) data."""
    data = fits.open(fn)[1].data
    flux = np.concatenate([np.asarray(r) for r in data['FLUX']])
    err  = np.concatenate([np.asarray(r) for r in data['ERROR']])
    good = np.isfinite(flux) & np.isfinite(err) & (err > 0)
    return float(np.nanmedian(flux[good] / err[good])) if good.any() else -np.inf


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


def read_data_flat(name, path, data_origin, z):
    """
    Read spectra from a flat directory layout produced by pipeline/download_spectra.py.

    Unlike read_data(), which expects NecessaryParams.csv and per-obs subdirectories,
    this function globs all relevant files directly from `path` and reads the grating
    from each file's FITS primary header (OPT_ELEM keyword).

    Parameters
    ----------
    name : str
        Object common name (used for logging in Cut_Edge_Pix calls).
    path : str
        Flat directory containing the FITS files (raw_data/{name}/{inst}/).
    data_origin : str
        Instrument: 'COS', 'STIS', or 'FOS'.
    z : float
        Redshift.
    """
    if data_origin == "COS":
        return read_cos_flat(name, path, z)
    elif data_origin == "STIS":
        return read_stis_flat(name, path, z)
    elif data_origin == "FOS":
        return read_fos_flat(name, path, z)
    print("data_origin not recognized for flat reader")


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


# ── Flat-layout readers (master catalog / MAST download) ─────────────────────
# Files land directly in the instrument folder with no object subdirectory
# and no NecessaryParams.csv.  Grating is read from the FITS primary header.

def read_cos_flat(name, path, z):
    """Read COS x1d files from a flat directory (no NecessaryParams.csv)."""
    # Only use per-exposure _x1d.fits products (matches Trevor's legacy convention).
    # _x1dsum.fits files are CalCOS per-visit coadds of the same exposures, so
    # including them would double-count the same photons in coadd.py.
    fn_list = sorted(glob.glob(os.path.join(path, '*_x1d.fits')))
    if not fn_list:
        raise FileNotFoundError("No COS x1d files found in %s" % path)

    # Drop files where CalCOS extraction failed and produced an empty BinTable
    # (zero rows). These are still archived by MAST and would crash the readers
    # downstream. Failed visits are typically re-observed under a new ASN_ID.
    fn_list = [fn for fn in fn_list if fits.open(fn)[1].data is not None
               and len(fits.open(fn)[1].data) > 0]
    if not fn_list:
        raise FileNotFoundError("No non-empty COS x1d files found in %s" % path)

    # AP 2026-06-17: TEMP cap for monitoring targets -- see _COS_EXPOSURE_CAP above.
    if name in _COS_EXPOSURE_CAP and len(fn_list) > _COS_EXPOSURE_CAP[name]:
        n_keep = _COS_EXPOSURE_CAP[name]
        ranked = sorted(fn_list, key=_cos_median_snr, reverse=True)
        fn_list = sorted(ranked[:n_keep])  # re-sort by name for deterministic order
        print("CAP COS %s: kept %d highest-S/N of %d exposures (TEMP -- see "
              "Migration_Log.md)" % (name, n_keep, len(ranked)), flush=True)

    gratings = []
    for fn in fn_list:
        hdr = fits.open(fn)[0].header
        gratings.append(hdr.get('OPT_ELEM', hdr.get('FILTER', 'UNKNOWN')))

    array_sizes = []
    for i, fn in enumerate(fn_list):
        data = fits.open(fn)[1].data
        if gratings[i] == 'E140M':
            array_sizes.append(data.size * 1024)
        else:
            array_sizes.append(len(data['Wavelength'][0]))
    array_len = max(array_sizes)

    waves     = np.zeros((len(fn_list), array_len))
    fluxes    = np.zeros((len(fn_list), array_len))
    flux_errs = np.zeros((len(fn_list), array_len))
    masks     = np.zeros((len(fn_list), array_len))

    for i, fn in enumerate(fn_list):
        spec_id = os.path.basename(fn).replace('_x1d.fits', '').replace('_x1dsum.fits', '')
        data = fits.open(fn)[1].data

        if gratings[i] == "E140M":
            wavelength = []
            flux       = []
            fluxerr    = []
            DQ         = []
            for t in range(data.size):
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
        # Slice to len(err_wmask) so that COS NUV stripes (1274 px) don't crash
        # when this object also has COS FUV segments (16384 px) and array_len is
        # the max. Matches the FOS/STIS reader pattern; padding stays as zeros
        # and coadd.py filters with waves[waves!=0].
        masks[i,:len(err_wmask)][(err_wmask == 0.0)] = 1
        waves[i,:len(err_wmask)]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                wavelength, array_len, "wavelength", False, z,
                                                "%s - %s" % (name, spec_id), "COS")
        fluxes[i,:len(err_wmask)]    = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                flux_wmask, array_len, "flux", False, z,
                                                "%s - %s" % (name, spec_id), "COS")
        flux_errs[i,:len(err_wmask)] = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                err_wmask, array_len, "flux error", False, z,
                                                "%s - %s" % (name, spec_id), "COS")
        masks[i,:len(err_wmask)]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                masks[i,:len(err_wmask)], array_len, "masks", False, z,
                                                "%s - %s" % (name, spec_id), "COS")

    return waves, fluxes, flux_errs, masks


def read_stis_flat(name, path, z):
    """Read STIS x1d/sx1 files from a flat directory (no NecessaryParams.csv)."""
    # _x1d.fits and _sx1.fits are detector-specific products, not duplicates:
    # MAMA modes (G140L/M, G230L/M, E140M, E230M, PRISM) → _x1d.fits;
    # CCD modes with CR-SPLIT (G430L/M, G750L/M, G230LB, G230MB) → _sx1.fits.
    # They are mutually exclusive per rootname, so include both.
    fn_list = sorted(
        glob.glob(os.path.join(path, '*_x1d.fits')) +
        glob.glob(os.path.join(path, '*_sx1.fits'))
    )
    if not fn_list:
        raise FileNotFoundError("No STIS x1d/sx1 files found in %s" % path)

    # Defensive: drop files where pipeline extraction produced an empty BinTable.
    # Not currently observed in STIS data, but matches the COS reader hardening.
    fn_list = [fn for fn in fn_list if fits.open(fn)[1].data is not None
               and len(fits.open(fn)[1].data) > 0]
    if not fn_list:
        raise FileNotFoundError("No non-empty STIS x1d/sx1 files found in %s" % path)

    gratings = []
    for fn in fn_list:
        hdr = fits.open(fn)[0].header
        gratings.append(hdr.get('OPT_ELEM', hdr.get('FILTER', 'UNKNOWN')))

    array_sizes = []
    for i, fn in enumerate(fn_list):
        data = fits.open(fn)[1].data
        if gratings[i] == 'E140M':
            array_sizes.append(data.size * 1024)
        else:
            array_sizes.append(len(data['Wavelength'][0]))
    array_len = max(array_sizes)

    waves     = np.zeros((len(fn_list), array_len))
    fluxes    = np.zeros((len(fn_list), array_len))
    flux_errs = np.zeros((len(fn_list), array_len))
    masks     = np.zeros((len(fn_list), array_len))

    for i, fn in enumerate(fn_list):
        spec_id = os.path.basename(fn).replace('_x1d.fits', '').replace('_sx1.fits', '')
        data = fits.open(fn)[1].data

        if gratings[i] == "E140M":
            wavelength = []
            flux       = []
            fluxerr    = []
            DQ         = []
            for t in range(data.size):
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
        masks[i,:len(err_wmask)][err_wmask == 0.] = 1
        waves[i,:len(err_wmask)]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                wavelength, array_len, "wavelength", False, z,
                                                "%s - %s" % (name, spec_id), "STIS")
        fluxes[i,:len(err_wmask)]    = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                flux_wmask, array_len, "flux", False, z,
                                                "%s - %s" % (name, spec_id), "STIS")
        flux_errs[i,:len(err_wmask)] = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                err_wmask, array_len, "flux error", False, z,
                                                "%s - %s" % (name, spec_id), "STIS")
        masks[i,:len(err_wmask)]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                masks[i,:len(err_wmask)], array_len, "masks", False, z,
                                                "%s - %s" % (name, spec_id), "STIS")

    return waves, fluxes, flux_errs, masks


def read_fos_flat(name, path, z):
    """
    Read FOS c0f/c1f/c2f/cqf files from a flat directory (no NecessaryParams.csv).
    Discovers exposures by globbing *_c0f.fits; algorithm is identical to read_fos().
    """
    wave_files = sorted(glob.glob(os.path.join(path, '*_c0f.fits')))
    if not wave_files:
        raise FileNotFoundError("No FOS c0f files found in %s" % path)

    spec_names = [os.path.basename(f).replace('_c0f.fits', '') for f in wave_files]

    array_lens = np.array([], dtype=int)
    for sn in spec_names:
        wave = fits.open(os.path.join(path, '%s_c0f.fits' % sn))[0].data
        array_lens = np.append(array_lens, wave.shape[-1])
    array_len = max(array_lens)

    waves     = np.zeros((len(spec_names), array_len))
    fluxes    = np.zeros((len(spec_names), array_len))
    flux_errs = np.zeros((len(spec_names), array_len))
    masks     = np.zeros((len(spec_names), array_len))

    def accum_flag(index, flux_arr):
        return (np.mean(flux_arr[index]) > np.mean(flux_arr[index-1]) and
                np.mean(flux_arr[index-1]) > np.mean(flux_arr[index-2]))

    for i, sn in enumerate(spec_names):
        wave        = fits.open(os.path.join(path, '%s_c0f.fits' % sn))[0].data
        obs_flux    = fits.open(os.path.join(path, '%s_c1f.fits' % sn))[0].data
        obs_fluxerr = fits.open(os.path.join(path, '%s_c2f.fits' % sn))[0].data
        cqf_path    = os.path.join(path, '%s_cqf.fits' % sn)
        if os.path.exists(cqf_path):
            obs_DQ = fits.open(cqf_path)[0].data
        else:
            obs_DQ = np.zeros_like(obs_flux)

        ind = -1

        if len(wave.shape) > 1:
            nzero = (wave[ind] > 0.)
            wavelength = wave[ind][::-1] if wave[ind][nzero][0] > wave[ind][nzero][-1] else wave[ind][:]
        else:
            nzero = (wave > 0.)
            wavelength = wave[::-1] if wave[nzero][0] > wave[nzero][-1] else wave[:]

        flux    = np.zeros(len(wavelength))
        fluxerr = np.zeros(len(wavelength))
        DQ      = np.zeros(len(wavelength))

        if len(wave.shape) > 1:
            if accum_flag(ind, obs_flux):
                if wave[ind][nzero][0] > wave[ind][nzero][-1]:
                    flux    = obs_flux[ind][::-1]
                    fluxerr = obs_fluxerr[ind][::-1]
                    DQ      = obs_DQ[ind][::-1]
                else:
                    flux    = obs_flux[ind]
                    fluxerr = obs_fluxerr[ind]
                    DQ      = obs_DQ[ind]
            else:
                for l in range(len(wavelength)):
                    if wave[ind][nzero][0] > wave[ind][nzero][-1]:
                        flux[l]    = np.mean(obs_flux[:, -l-1])
                        fluxerr[l] = np.sqrt(sum(obs_fluxerr[:, -l-1]**2)) / (obs_fluxerr.shape[0]-1)
                        DQ[l]      = sum(obs_DQ[:, -l-1])
                    else:
                        flux[l]    = np.mean(obs_flux[:, l])
                        fluxerr[l] = np.sqrt(sum(obs_fluxerr[:, l]**2)) / (obs_fluxerr.shape[0]-1)
                        DQ[l]      = sum(obs_DQ[:, l])
        else:
            if wave[nzero][0] > wave[nzero][-1]:
                flux    = obs_flux[::-1]
                fluxerr = obs_fluxerr[::-1]
                DQ      = obs_DQ[::-1]
            else:
                flux    = obs_flux
                fluxerr = obs_fluxerr
                DQ      = obs_DQ

        flux_wmask = flux.copy()
        err_wmask  = fluxerr.copy()
        masks[i,:len(err_wmask)][err_wmask == 0.] = 1

        waves[i,:len(flux)]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                wavelength, array_len, "wavelength", False, z,
                                                "%s - %s" % (name, sn), "FOS")
        fluxes[i,:len(flux)]    = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                flux_wmask, array_len, "flux", False, z,
                                                "%s - %s" % (name, sn), "FOS")
        flux_errs[i,:len(flux)] = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                err_wmask, array_len, "flux error", False, z,
                                                "%s - %s" % (name, sn), "FOS")
        masks[i,:len(flux)]     = Cut_Edge_Pix_TVM.Cut_Edge_Pix(DQ, wavelength, flux_wmask, err_wmask,
                                                masks[i,:len(err_wmask)], array_len, "masks", False, z,
                                                "%s - %s" % (name, sn), "FOS")

    return waves, fluxes, flux_errs, masks
