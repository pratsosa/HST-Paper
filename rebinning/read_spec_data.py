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

AP: 2026-09-10 -- the three flat-layout readers (read_cos_flat, read_stis_flat,
read_fos_flat) were rewritten to stop discarding data. See the comment block above
them, and diagnostics/CIV_Data_Loss_Investigation.md. The legacy readers
(read_cos, read_stis, read_fos, read_hsla, read_sdssrm) are untouched.
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

# AP 2026-09-10: floor below which a "wavelength" is not a wavelength solution.
# The bluest genuine value anywhere in this dataset is 703 A (COS G140L/800 FUVA);
# COS G130M/1055 FUVB reaches 892 A. Nothing real sits below 500 A, and 500 A
# observed is rest-frame <500 A for every object, far blueward of the 1120 A rest
# cut coadd.rebin() applies, so this can never remove data that reaches the output.
_MIN_PHYSICAL_WAVE = 500.0


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


# --- Flat-layout readers (master catalog / MAST download) -------------------
# Files land directly in the instrument folder with no object subdirectory
# and no NecessaryParams.csv.  Grating is read from the FITS primary header.
#
# AP 2026-09-10: rewritten to fix three ways spectra were being lost. The three
# public readers below keep their names and signatures; they now share one
# implementation because the bugs were duplicated across all three copies.
# See diagnostics/CIV_Data_Loss_Investigation.md for the evidence.
#
#   Fix A - every detector row is read, not just row 0.
#       The old code did `wavelength = data['WAVELENGTH'][0]`. A COS FUV x1d
#       table has two rows (FUVA, FUVB) and a COS NUV table has three
#       (NUVA, NUVB, NUVC); each row is a physically separate piece of the
#       spectrum, so row 0 alone discarded half to two-thirds of the coverage.
#       Every row now becomes its own entry in the exposure stack. They are kept
#       separate rather than concatenated because segments are neither
#       contiguous nor sorted by wavelength; coadd.py already handles
#       independent, overlapping exposures, so nothing downstream changes.
#       (E140M echelle orders ARE contiguous pieces of one spectrum and stay
#       concatenated, exactly as before.)
#       STIS first-order x1d tables have one row, so STIS output is unchanged.
#
#   Fix B - an exposure that fails Cut_Edge_Pix is skipped, not fatal.
#       Cut_Edge_Pix raises IndexError when an exposure has no usable pixel at
#       all. That exception used to propagate out of the reader and discard
#       every good exposure alongside the bad one, leaving the object with no
#       output file.
#
#   Fix C - FOS exposures with only two readout groups no longer crash.
#       accum_flag() tested mean(flux[-1]) > mean(flux[-2]) > mean(flux[-3]);
#       files with NAXIS2=2 have no third group and raised IndexError. With two
#       groups the test degrades to the single comparison it can make.
#
# The FOS multi-row error formula sqrt(sum(e**2))/(N-1) and the SpecCuts edge
# cuts are unchanged from the 2026AP lineage.


def _cos_stis_rows(fn, skipped=None):
    """
    Yield (spec_id, wave, flux, fluxerr, DQ) for each usable row of one
    COS/STIS x1d (or sx1) file. Empty extraction products yield nothing.
    Rows rejected for an unusable wavelength solution are appended to `skipped`
    (a list) so the caller can report them once rather than per file.
    """
    base = os.path.basename(fn)
    for suf in ("_x1d.fits", "_x1dsum.fits", "_sx1.fits"):
        base = base.replace(suf, "")
    with fits.open(fn) as h:
        grating = h[0].header.get("OPT_ELEM", h[0].header.get("FILTER", "UNKNOWN"))
        data = h[1].data
        if data is None or len(data) == 0:
            return
        if grating == "E140M":
            # echelle orders are contiguous pieces of one spectrum
            cat = lambda col: np.concatenate([data[col][t] for t in range(data.size)])
            yield (base, cat("WAVELENGTH"), cat("FLUX"), cat("ERROR"), cat("DQ"))
            return
        has_seg = "SEGMENT" in data.columns.names
        for t in range(len(data)):
            seg = data["SEGMENT"][t] if has_seg else str(t)
            wavelength = np.asarray(data["WAVELENGTH"][t])

            # AP 2026-09-10: two rejections, both added when Fix A started reading
            # rows that the old row-0-only reader never saw.
            #
            # COS G140L segment FUVB is not illuminated in first order at the
            # cenwaves used here (1230, 1280). CalCOS still writes the row, but the
            # wavelength solution is a placeholder running from ~0.02 A, and flux
            # IS recorded against those nonsense wavelengths. Including it stretches
            # the rebin grid in coadd.py over four decades of wavelength, which
            # overruns the Vanden Berk composite in spec_morph.morph2().
            if grating == "G140L" and seg == "FUVB":
                if skipped is not None:
                    skipped.append("%s %s (G140L FUVB, not first-order illuminated)"
                                   % (grating, seg))
                continue

            # General safety net for the same failure mode in a configuration we
            # have not seen: a row whose wavelengths start below _MIN_PHYSICAL_WAVE
            # does not carry a real wavelength solution.
            pos = wavelength > 0
            if pos.any() and wavelength[pos].min() < _MIN_PHYSICAL_WAVE:
                if skipped is not None:
                    skipped.append("%s %s (wavelengths start at %.3f A)"
                                   % (grating, seg, wavelength[pos].min()))
                continue

            yield ("%s-%s" % (base, seg), wavelength,
                   np.asarray(data["FLUX"][t]), np.asarray(data["ERROR"][t]),
                   np.asarray(data["DQ"][t]))


def _fos_rows(fn, path):
    """Yield (spec_id, wave, flux, fluxerr, DQ) for one FOS exposure."""
    sn = os.path.basename(fn).replace("_c0f.fits", "")
    part = lambda suf: fits.open(os.path.join(path, "%s_%s.fits" % (sn, suf)))[0].data
    wave, obs_flux, obs_fluxerr = part("c0f"), part("c1f"), part("c2f")
    cqf = os.path.join(path, "%s_cqf.fits" % sn)
    obs_DQ = fits.open(cqf)[0].data if os.path.exists(cqf) else np.zeros_like(obs_flux)

    if wave.ndim == 1:
        nzero = wave > 0.
        rev = wave[nzero][0] > wave[nzero][-1]
        sl = slice(None, None, -1) if rev else slice(None)
        yield sn, wave[sl], obs_flux[sl], obs_fluxerr[sl], obs_DQ[sl]
        return

    ind = -1
    nzero = wave[ind] > 0.
    rev = wave[ind][nzero][0] > wave[ind][nzero][-1]
    wavelength = wave[ind][::-1] if rev else wave[ind][:]

    # Fix C: the cumulative-group test needs three groups; with two, compare the
    # only pair there is. (Original: mean[-1] > mean[-2] > mean[-3].)
    means = [np.mean(obs_flux[k]) for k in range(obs_flux.shape[0])]
    if len(means) >= 3:
        accum = means[-1] > means[-2] and means[-2] > means[-3]
    else:
        accum = means[-1] > means[-2]

    sl = slice(None, None, -1) if rev else slice(None)
    if accum:
        yield sn, wavelength, obs_flux[ind][sl], obs_fluxerr[ind][sl], obs_DQ[ind][sl]
    else:
        n_rows = obs_fluxerr.shape[0]
        yield (sn, wavelength,
               np.mean(obs_flux, axis=0)[sl],
               (np.sqrt(np.sum(obs_fluxerr ** 2, axis=0)) / (n_rows - 1))[sl],
               np.sum(obs_DQ, axis=0)[sl])


def _report_skipped_rows(name, skipped):
    """One compact line per distinct reason, however many rows it applied to."""
    if not skipped:
        return
    counts = {}
    for r in skipped:
        counts[r] = counts.get(r, 0) + 1
    for reason, n in sorted(counts.items()):
        print("  %s: skipped %d rows -- %s" % (name, n, reason), flush=True)


def _stack_flat_rows(name, rows, z, Instrument):
    """
    Run Cut_Edge_Pix over each (spec_id, wave, flux, fluxerr, DQ) row and stack
    the survivors into the (n_rows, array_len) arrays coadd.rebin() expects.
    Rows that Cut_Edge_Pix rejects outright are skipped and reported (Fix B).
    """
    array_len = max(len(r[1]) for r in rows)
    waves     = np.zeros((len(rows), array_len))
    fluxes    = np.zeros((len(rows), array_len))
    flux_errs = np.zeros((len(rows), array_len))
    masks     = np.zeros((len(rows), array_len))
    keep, dropped = [], []

    for i, (spec_id, wavelength, flux, fluxerr, DQ) in enumerate(rows):
        n = len(fluxerr)
        flux_wmask = flux.copy()
        err_wmask  = fluxerr.copy()
        exp_id = "%s - %s" % (name, spec_id)
        try:
            masks[i, :n][err_wmask == 0.] = 1
            cut = lambda input_data, kind: Cut_Edge_Pix_TVM.Cut_Edge_Pix(
                DQ, wavelength, flux_wmask, err_wmask, input_data, array_len,
                kind, False, z, exp_id, Instrument)
            waves[i, :n]     = cut(wavelength, "wavelength")
            fluxes[i, :n]    = cut(flux_wmask, "flux")
            flux_errs[i, :n] = cut(err_wmask, "flux error")
            masks[i, :n]     = cut(masks[i, :n], "masks")
            keep.append(i)
        except IndexError as e:
            dropped.append("%s (%s)" % (spec_id, e))

    if dropped:
        print("  %s: skipped %d unusable exposure rows: %s"
              % (name, len(dropped), "; ".join(dropped[:3])), flush=True)
    if not keep:
        raise IndexError("All exposure rows failed the edge cut for %s" % name)
    return waves[keep], fluxes[keep], flux_errs[keep], masks[keep]


def read_cos_flat(name, path, z):
    """Read COS x1d files from a flat directory (no NecessaryParams.csv)."""
    # Only per-exposure _x1d.fits products (matches Trevor's legacy convention).
    # _x1dsum.fits are CalCOS per-visit coadds of the same exposures, so
    # including them would double-count the same photons in coadd.py.
    fn_list = sorted(glob.glob(os.path.join(path, '*_x1d.fits')))
    if not fn_list:
        raise FileNotFoundError("No COS x1d files found in %s" % path)

    # Drop files where CalCOS extraction failed and produced an empty BinTable
    # (zero rows). MAST archives these; failed visits are typically re-observed
    # under a new ASN_ID. _cos_stis_rows() skips them too, but the exposure cap
    # below ranks files by S/N and cannot read an empty one, so filter first.
    fn_list = [fn for fn in fn_list if fits.open(fn)[1].data is not None
               and len(fits.open(fn)[1].data) > 0]
    if not fn_list:
        raise FileNotFoundError("No non-empty COS x1d files found in %s" % path)

    # AP 2026-06-17: TEMPORARY exposure cap for monitoring targets. Applied to
    # FILES, before rows are expanded, so the cap still means "N exposures".
    # See Migration_Log.md -- still a stopgap, still needs revisiting.
    if name in _COS_EXPOSURE_CAP and len(fn_list) > _COS_EXPOSURE_CAP[name]:
        n_keep = _COS_EXPOSURE_CAP[name]
        ranked = sorted(fn_list, key=_cos_median_snr, reverse=True)
        fn_list = sorted(ranked[:n_keep])
        print("CAP COS %s: kept %d highest-S/N of %d exposures (TEMP -- see "
              "Migration_Log.md)" % (name, n_keep, len(ranked)), flush=True)

    skipped = []
    rows = [r for fn in fn_list for r in _cos_stis_rows(fn, skipped)]
    _report_skipped_rows(name, skipped)
    if not rows:
        raise FileNotFoundError("No usable COS x1d rows found in %s" % path)
    return _stack_flat_rows(name, rows, z, "COS")


def read_stis_flat(name, path, z):
    """Read STIS x1d/sx1 files from a flat directory (no NecessaryParams.csv)."""
    # _x1d.fits and _sx1.fits are detector-specific products, not duplicates:
    # MAMA modes (G140L/M, G230L/M, E140M, E230M, PRISM) -> _x1d.fits;
    # CCD modes with CR-SPLIT (G430L/M, G750L/M, G230LB, G230MB) -> _sx1.fits.
    # They are mutually exclusive per rootname, so include both.
    fn_list = sorted(
        glob.glob(os.path.join(path, '*_x1d.fits')) +
        glob.glob(os.path.join(path, '*_sx1.fits'))
    )
    if not fn_list:
        raise FileNotFoundError("No STIS x1d/sx1 files found in %s" % path)

    skipped = []
    rows = [r for fn in fn_list for r in _cos_stis_rows(fn, skipped)]
    _report_skipped_rows(name, skipped)
    if not rows:
        raise FileNotFoundError("No usable STIS x1d/sx1 rows found in %s" % path)
    return _stack_flat_rows(name, rows, z, "STIS")


def read_fos_flat(name, path, z):
    """
    Read FOS c0f/c1f/c2f/cqf files from a flat directory (no NecessaryParams.csv).
    Discovers exposures by globbing *_c0f.fits.
    """
    wave_files = sorted(glob.glob(os.path.join(path, '*_c0f.fits')))
    if not wave_files:
        raise FileNotFoundError("No FOS c0f files found in %s" % path)

    rows = [r for fn in wave_files for r in _fos_rows(fn, path)]
    if not rows:
        raise FileNotFoundError("No usable FOS exposures found in %s" % path)
    return _stack_flat_rows(name, rows, z, "FOS")
