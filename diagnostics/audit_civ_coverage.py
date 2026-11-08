"""
diagnostics/audit_civ_coverage.py
AP: 2026-09-10

Audit: which objects have C IV (1549.48 A rest) on disk in the raw MAST data but
NOT in the rebinned spectrum written to RebinnedSpec_master/?

For every (object, instrument) pair the rebinning pipeline was asked to process,
this compares three wavelength coverages, all in the rest frame:

  raw_all   - every row of every exposure file the reader globs
  raw_read  - only the rows the reader actually keeps (row 0, except E140M)
  output    - pixels in the rebinned FITS with finite, non-zero flux

and classifies each pair:

  NO_OUTPUT     no rebinned file exists -> the rebin run crashed for this object
  LOST_ROWS     C IV signal in raw_all but not in raw_read -> reader dropped the
                detector segment / echelle order that holds it
  LOST_COADD    C IV signal in raw_read but not in output -> lost inside coadd.py
  NO_CIV_SIGNAL C IV inside the nominal wavelength array but the raw flux/error
                there is at or below the noise -> nothing worth recovering
  NO_CIV_ON_SKY C IV outside the nominal wavelength array
  OK            C IV present in the output

"Signal" means the median per-pixel flux/error over rest-frame 1520-1580 A exceeds
SNR_FLOOR. The measured value is written to the CSV (civ_raw_snr_all / _read) so the
threshold can be re-judged without re-running.

Reading every wavelength array would be slow (some objects have >1000 exposures),
so exposures are grouped by their instrument configuration (OPT_ELEM/CENWAVE from
the primary header, which is cheap to read) and one representative file per group
is opened. Different exposures in the same configuration have essentially the same
wavelength range, so this reproduces the coverage faithfully.

Usage:
    python -m diagnostics.audit_civ_coverage --raw-data-dir "F:/Richards-data/HST Paper/raw_data"
"""

import argparse
import glob
import os
import warnings

import numpy as np
import pandas as pd
from astropy.io import fits

warnings.simplefilter("ignore")

CIV = 1549.48
SNR_FLOOR = 0.3   # median per-pixel flux/error in the C IV window below which
                  # the raw data carries no usable signal there
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def reader_files(path, inst):
    """The file list each flat reader globs, after its empty-BinTable filter."""
    if inst == "FOS":
        return sorted(glob.glob(os.path.join(path, "*_c0f.fits")))
    if inst == "COS":
        fns = sorted(glob.glob(os.path.join(path, "*_x1d.fits")))
    else:  # STIS
        fns = sorted(glob.glob(os.path.join(path, "*_x1d.fits")) +
                     glob.glob(os.path.join(path, "*_sx1.fits")))
    keep = []
    for fn in fns:
        with fits.open(fn, memmap=True) as h:
            if h[1].data is not None and len(h[1].data) > 0:
                keep.append(fn)
    return keep


def config_key(fn):
    """Cheap instrument-configuration key used to group near-identical exposures."""
    with fits.open(fn, memmap=True) as h:
        hd = h[0].header
        return (hd.get("OPT_ELEM", hd.get("FILTER", "?")), hd.get("CENWAVE", "?"))


def _fos_rows(fn):
    """(wave, flux, err) rows for a FOS exposure; the reader uses the last accum row."""
    base = fn[:-len("_c0f.fits")]
    w = fits.getdata(fn)
    try:
        f = fits.getdata(base + "_c1f.fits")
        e = fits.getdata(base + "_c2f.fits")
    except FileNotFoundError:
        f = e = np.zeros_like(w)
    if w.ndim == 1:
        return [(w, f, e)]
    return [(w[-1], f[-1], e[-1])]


def file_rows(fn, inst):
    """(all_rows, read_rows) as lists of (wave, flux, err) arrays for one exposure."""
    if inst == "FOS":
        rows = _fos_rows(fn)
        return rows, rows
    with fits.open(fn, memmap=True) as h:
        grating = h[0].header.get("OPT_ELEM", h[0].header.get("FILTER", "?"))
        d = h[1].data
        allr = [(np.asarray(d["WAVELENGTH"][t]), np.asarray(d["FLUX"][t]),
                 np.asarray(d["ERROR"][t])) for t in range(len(d))]
    return allr, (allr if grating == "E140M" else allr[:1])


def spans(rows):
    return [(float(w[w > 0].min()), float(w[w > 0].max()))
            for w, _, _ in rows if (w > 0).any()]


def civ_snr(rows, z, lo=1520.0, hi=1580.0):
    """Median per-pixel flux/error in the rest-frame C IV window across `rows`."""
    vals = []
    for w, f, e in rows:
        s = (w / (1 + z) > lo) & (w / (1 + z) < hi) & (e > 0)
        if s.any():
            vals.append(f[s] / e[s])
    return float(np.median(np.concatenate(vals))) if vals else np.nan


def covers(intervals, x):
    return any(lo <= x <= hi for lo, hi in intervals)


def audit_pair(name, z, inst, path, rebin_dir):
    fns = reader_files(path, inst)
    row = dict(name=name, instrument=inst, z=z, n_exposures=len(fns),
               n_x1dsum_excluded=len(glob.glob(os.path.join(path, "*_x1dsum.fits")))
               if inst == "COS" else 0)
    if not fns:
        row.update(verdict="NO_RAW_FILES")
        return row

    # one representative exposure per instrument configuration
    reps = {}
    for fn in fns:
        reps.setdefault(config_key(fn), fn)
    all_rows, read_rows = [], []
    for fn in reps.values():
        a, r = file_rows(fn, inst)
        all_rows += a
        read_rows += r
    to_rest = lambda iv: [(lo / (1 + z), hi / (1 + z)) for lo, hi in iv]
    all_iv, read_iv = to_rest(spans(all_rows)), to_rest(spans(read_rows))

    row["n_configs"] = len(reps)
    row["raw_all_min"] = min(lo for lo, _ in all_iv)
    row["raw_all_max"] = max(hi for _, hi in all_iv)
    row["civ_in_raw_all"] = covers(all_iv, CIV)
    row["civ_in_raw_read"] = covers(read_iv, CIV)
    row["civ_raw_snr_all"] = civ_snr(all_rows, z)
    row["civ_raw_snr_read"] = civ_snr(read_rows, z)

    out = os.path.join(rebin_dir, "%s_%s.fits" % (name, inst))
    if not os.path.exists(out):
        row.update(verdict="NO_OUTPUT", civ_in_output=False)
        return row

    d = fits.getdata(out, 1)
    w, f = d["Rest-frame Wavelength"], d["Coadded Flux (Arbitrary Units)"]
    good = np.isfinite(f) & (f != 0)
    row["out_good_min"] = float(w[good].min()) if good.any() else np.nan
    row["out_good_max"] = float(w[good].max()) if good.any() else np.nan
    # require real data across the C IV line core, not just one stray pixel
    core = (w > 1500) & (w < 1600)
    row["civ_good_frac"] = float(good[core].mean()) if core.any() else 0.0
    row["civ_in_output"] = row["civ_good_frac"] > 0.5

    if row["civ_in_output"]:
        row["verdict"] = "OK"
    elif not row["civ_in_raw_all"]:
        row["verdict"] = "NO_CIV_ON_SKY"
    elif not (row["civ_raw_snr_all"] > SNR_FLOOR):
        row["verdict"] = "NO_CIV_SIGNAL"
    elif not row["civ_in_raw_read"] or not (row["civ_raw_snr_read"] > SNR_FLOOR):
        row["verdict"] = "LOST_ROWS"
    else:
        row["verdict"] = "LOST_COADD"
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-data-dir", required=True)
    ap.add_argument("--rebin-dir", default=os.path.join(REPO, "RebinnedSpec_master"))
    ap.add_argument("--out", default=os.path.join(REPO, "diagnostics", "civ_coverage_audit.csv"))
    args = ap.parse_args()

    import sys
    sys.path.insert(0, REPO)
    from pipeline.catalog_bridge import build_object_list

    objs = build_object_list(args.raw_data_dir, verbose=False)
    rows = []
    for k, o in enumerate(objs):
        try:
            rows.append(audit_pair(o["name"], float(o["redshift"]), o["instrument"],
                                   str(o["data_path"]), args.rebin_dir))
        except Exception as e:
            rows.append(dict(name=o["name"], instrument=o["instrument"],
                             z=o["redshift"], verdict="AUDIT_ERROR",
                             error="%s: %s" % (type(e).__name__, e)))
        if (k + 1) % 25 == 0:
            print("  %d/%d" % (k + 1, len(objs)), flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print("\nwrote %s (%d rows)\n" % (args.out, len(df)))
    print(df.groupby(["verdict", "instrument"]).size().to_string())


if __name__ == "__main__":
    main()
