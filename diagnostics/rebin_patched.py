"""
diagnostics/rebin_patched.py
AP: 2026-09-10

Candidate fixes for the three ways C IV data is lost between the raw MAST files
and RebinnedSpec_master/. Nothing in rebinning/ is modified: the fixes are applied
by monkey-patching the two functions coadd.rebin() looks up at call time
(`read_spec_data.read_data_flat` and `spec_morph.cont_filtered`), so the same
`coadd.rebin` runs with patched behaviour.

Fix A - all detector rows are read (read_data_flat)
    The COS/STIS readers take `data['WAVELENGTH'][0]`, i.e. only the FIRST row of
    the x1d table. A COS FUV exposure has two rows (FUVA, FUVB) and a COS NUV
    exposure has three (NUVA, NUVB, NUVC); each is a physically separate piece of
    the spectrum. Row 0 alone throws away half to two-thirds of the coverage.
    Here every row becomes its own entry in the exposure stack. They are kept
    separate rather than concatenated because segments are not contiguous and are
    not sorted by wavelength - coadd.py already handles independent, overlapping
    exposures, so this needs no other change.

Fix B - exposures that fail the edge cut are skipped, not fatal (read_data_flat)
    Cut_Edge_Pix raises IndexError("No good pixels found") when one exposure is
    entirely bad. That exception propagates out of the reader and kills the whole
    object, discarding every good exposure alongside it.

Fix C - FOS exposures with only two readout groups (read_data_flat)
    read_fos_flat's accum_flag() indexes flux[-1], flux[-2] and flux[-3] to test
    whether the groups are cumulative. FOS files with NAXIS2=2 have no third
    group, so it raises IndexError and the object produces no output at all.
    With two groups the test degrades to the single comparison it can make.

Fix D - continuum of zero inside an emission line (cont_filtered)
    cont_filtered() masks out emission-line windows, takes a running median of
    what is left, and interpolates back across the masked gaps. If an exposure
    lies ENTIRELY inside one window - e.g. a COS G160M exposure covering rest
    1475-1660 A, wholly inside the C IV window of 1460-1685 A - there is nothing
    left to take a median of, and the defensive branch added during the migration
    fills the continuum with 0.0. coadd.py then divides the exposure by that
    continuum, so every pixel becomes +/-inf and the exposure is silently dropped.
    (Before that defensive branch existed the same case raised IndexError, so the
    object failed loudly instead.) Here the continuum falls back to a constant -
    the exposure's own median flux - which preserves the data and puts it on the
    same normalised scale as the neighbouring exposures. The absolute scale does
    not matter: coadd.py morphs the co-added spectrum afterwards.

Usage:
    python -m diagnostics.rebin_patched --raw-data-dir "F:/.../raw_data" \
        --out-dir RebinnedSpec_master_patched --verdicts LOST_ROWS,LOST_COADD,NO_OUTPUT
"""

import argparse
import glob
import os
import sys
import traceback
import warnings

import numpy as np
from astropy.io import fits

warnings.simplefilter("ignore")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from rebinning import coadd, read_spec_data, spec_morph  # noqa: E402
from rebinning import cut_edge_pix as CEP  # noqa: E402


# ── Fix A/B: one entry per detector row, bad exposures skipped ───────────────

def _x1d_rows(fn, inst):
    """Yield (spec_id, wave, flux, err, DQ) for every usable row of one x1d file."""
    base = os.path.basename(fn)
    for suf in ("_x1d.fits", "_x1dsum.fits", "_sx1.fits"):
        base = base.replace(suf, "")
    with fits.open(fn) as h:
        grating = h[0].header.get("OPT_ELEM", h[0].header.get("FILTER", "UNKNOWN"))
        d = h[1].data
        if d is None or len(d) == 0:
            return
        if grating == "E140M":
            # echelle orders are contiguous pieces of one spectrum: keep concatenated
            yield (base,
                   np.concatenate([d["WAVELENGTH"][t] for t in range(d.size)]),
                   np.concatenate([d["FLUX"][t] for t in range(d.size)]),
                   np.concatenate([d["ERROR"][t] for t in range(d.size)]),
                   np.concatenate([d["DQ"][t] for t in range(d.size)]))
            return
        for t in range(len(d)):
            seg = d["SEGMENT"][t] if "SEGMENT" in d.columns.names else str(t)
            yield ("%s-%s" % (base, seg), np.asarray(d["WAVELENGTH"][t]),
                   np.asarray(d["FLUX"][t]), np.asarray(d["ERROR"][t]),
                   np.asarray(d["DQ"][t]))


def _fos_rows(fn, path):
    """Yield (spec_id, wave, flux, err, DQ) for one FOS exposure (Fix C)."""
    sn = os.path.basename(fn).replace("_c0f.fits", "")
    g = lambda suf: fits.getdata(os.path.join(path, "%s_%s.fits" % (sn, suf)))
    wave, obs_flux, obs_err = g("c0f"), g("c1f"), g("c2f")
    cqf = os.path.join(path, "%s_cqf.fits" % sn)
    obs_DQ = fits.getdata(cqf) if os.path.exists(cqf) else np.zeros_like(obs_flux)

    if wave.ndim == 1:
        rev = wave[wave > 0][0] > wave[wave > 0][-1]
        sl = slice(None, None, -1) if rev else slice(None)
        yield sn, wave[sl], obs_flux[sl], obs_err[sl], obs_DQ[sl]
        return

    nzero = wave[-1] > 0.
    rev = wave[-1][nzero][0] > wave[-1][nzero][-1]
    wavelength = wave[-1][::-1] if rev else wave[-1][:]

    # Fix C: the cumulative-group test needs three groups; with two, compare the
    # only pair available. (Original: mean[-1] > mean[-2] > mean[-3].)
    means = [np.mean(obs_flux[k]) for k in range(obs_flux.shape[0])]
    if len(means) >= 3:
        accum = means[-1] > means[-2] > means[-3]
    else:
        accum = means[-1] > means[-2]

    if accum:
        sl = slice(None, None, -1) if rev else slice(None)
        yield sn, wavelength, obs_flux[-1][sl], obs_err[-1][sl], obs_DQ[-1][sl]
    else:
        cols = slice(None, None, -1) if rev else slice(None)
        n = obs_err.shape[0]
        yield (sn, wavelength,
               np.mean(obs_flux, axis=0)[cols],
               (np.sqrt(np.sum(obs_err ** 2, axis=0)) / (n - 1))[cols],
               np.sum(obs_DQ, axis=0)[cols])


def read_data_flat_patched(name, path, data_origin, z):
    if data_origin == "FOS":
        files = sorted(glob.glob(os.path.join(path, "*_c0f.fits")))
        rowgen = lambda fn: _fos_rows(fn, path)
    else:
        pats = ["*_x1d.fits"] + (["*_sx1.fits"] if data_origin == "STIS" else [])
        files = sorted(f for p in pats for f in glob.glob(os.path.join(path, p)))
        rowgen = lambda fn: _x1d_rows(fn, data_origin)
    if not files:
        raise FileNotFoundError("No %s files found in %s" % (data_origin, path))

    cap = read_spec_data._COS_EXPOSURE_CAP.get(name) if data_origin == "COS" else None
    if cap and len(files) > cap:
        files = sorted(sorted(files, key=read_spec_data._cos_median_snr,
                              reverse=True)[:cap])

    rows = [r for fn in files for r in rowgen(fn)]
    if not rows:
        raise FileNotFoundError("No non-empty %s rows in %s" % (data_origin, path))
    array_len = max(len(r[1]) for r in rows)

    waves     = np.zeros((len(rows), array_len))
    fluxes    = np.zeros((len(rows), array_len))
    flux_errs = np.zeros((len(rows), array_len))
    masks     = np.zeros((len(rows), array_len))
    keep, dropped = [], []

    for i, (spec_id, wavelength, flux, fluxerr, DQ) in enumerate(rows):
        n = len(fluxerr)
        tag = "%s - %s" % (name, spec_id)
        args = (DQ, wavelength, flux.copy(), fluxerr.copy())
        try:
            masks[i, :n][fluxerr == 0.] = 1
            cut = lambda data, kind: CEP.Cut_Edge_Pix(
                *args, data, array_len, kind, False, z, tag, data_origin)
            waves[i, :n]     = cut(wavelength, "wavelength")
            fluxes[i, :n]    = cut(flux.copy(), "flux")
            flux_errs[i, :n] = cut(fluxerr.copy(), "flux error")
            masks[i, :n]     = cut(masks[i, :n], "masks")
            keep.append(i)
        except IndexError as e:      # Fix B: one dead exposure must not kill the object
            dropped.append("%s (%s)" % (spec_id, e))

    if dropped:
        print("  skipped %d unusable exposure rows: %s"
              % (len(dropped), "; ".join(dropped[:3])), flush=True)
    if not keep:
        raise IndexError("All exposure rows failed the edge cut for %s" % name)
    return waves[keep], fluxes[keep], flux_errs[keep], masks[keep]


# ── Fix D: never return a continuum of zero ─────────────────────────────────

_cont_filtered_orig = spec_morph.cont_filtered


def cont_filtered_patched(wave, flux, z, Identifier):
    em = spec_morph.emission_lines2(wave / (1 + z))
    if em.all():
        # whole exposure inside one emission window -> flat continuum at its own
        # median flux, so the exposure survives the division in coadd.py
        med = np.nanmedian(flux[np.isfinite(flux) & (flux != 0)]) if np.isfinite(flux).any() else np.nan
        if not np.isfinite(med) or med == 0:
            med = 1.0
        return np.full(len(wave), med)
    cont = _cont_filtered_orig(wave, flux, z, Identifier)
    bad = ~np.isfinite(cont) | (cont == 0)
    if bad.any():                       # leftover zeros from partial-coverage cases
        good = np.isfinite(cont) & (cont != 0)
        cont = cont.copy()
        cont[bad] = np.median(cont[good]) if good.any() else 1.0
    return cont


def apply_patches():
    read_spec_data.read_data_flat = read_data_flat_patched
    spec_morph.cont_filtered = cont_filtered_patched


# ── Runner ───────────────────────────────────────────────────────────────────

def main():
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-data-dir", required=True)
    ap.add_argument("--out-dir", default=os.path.join(REPO, "RebinnedSpec_master_patched"))
    ap.add_argument("--audit", default=os.path.join(REPO, "diagnostics", "civ_coverage_audit.csv"))
    ap.add_argument("--verdicts", default="LOST_ROWS,LOST_COADD,NO_OUTPUT")
    ap.add_argument("--only", default=None, help="comma-separated object names")
    a = ap.parse_args()

    apply_patches()
    os.makedirs(a.out_dir, exist_ok=True)

    aud = pd.read_csv(a.audit)
    sel = aud[aud.verdict.isin(a.verdicts.split(","))]
    if a.only:
        sel = sel[sel.name.isin(a.only.split(","))]
    want = {(r["name"], r.instrument) for _, r in sel.iterrows()}

    from pipeline.catalog_bridge import build_object_list
    ok = fail = 0
    for o in build_object_list(a.raw_data_dir, verbose=False):
        if (o["name"], o["instrument"]) not in want:
            continue
        print("%s %s" % (o["name"], o["instrument"]), flush=True)
        try:
            coadd.rebin(o["name"], float(o["redshift"]), o["instrument"], None,
                        data_path=str(o["data_path"]) + os.sep, output_dir=a.out_dir,
                        sdss_spec_dir=None, flat=True)
            ok += 1
        except Exception as e:
            tb = traceback.extract_tb(sys.exc_info()[2])[-1]
            print("  FAIL %s: %s (%s:%d)" % (type(e).__name__, str(e)[:120],
                                             os.path.basename(tb.filename), tb.lineno),
                  flush=True)
            fail += 1
    print("\npatched rebin: %d written, %d failed -> %s" % (ok, fail, a.out_dir))


if __name__ == "__main__":
    main()
