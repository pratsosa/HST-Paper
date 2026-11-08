"""
diagnostics/v2_report.py
AP: 2026-09-10

Summary diff of RebinnedSpec_master_v2/ against RebinnedSpec_master/.

compare_rebins.py prints one line per object, which is unreadable across ~590
spectra. This aggregates instead: how many spectra are bit-identical, how many
changed and by how much, which appeared, which disappeared. Per-object detail
goes to a CSV.

Change metrics, over the wavelengths both versions cover:
  pix_ratio   good (finite, non-zero) pixels, new / old
  snr_ratio   median |flux/error|, new / old  -- depth, not coverage
  civ_aa_before / civ_aa_after
              Angstroms of the rest-frame C IV window (1500-1600 A) carrying
              good flux -- absolute, so it is comparable across runs whose grids
              span different wavelength ranges

Usage:
    python -m diagnostics.v2_report
    python -m diagnostics.v2_report --old RebinnedSpec_master --new RebinnedSpec_master_v2
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
from astropy.io import fits

warnings.simplefilter("ignore")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(p):
    d = fits.getdata(p, 1)
    return (np.asarray(d["Rest-frame Wavelength"], float),
            np.asarray(d["Coadded Flux (Arbitrary Units)"], float),
            np.asarray(d["Coadded Flux Errors"], float))


def civ_aa(w, g):
    """
    Angstroms of the rest-frame C IV window (1500-1600 A) that carry good flux.

    AP 2026-09-10: this replaces "fraction of grid pixels in the window that are
    good". That fraction is measured against whatever the grid happens to span, so
    it is NOT comparable between two runs with different coverage: a spectrum whose
    old grid stopped at 1530 A had only 77 pixels in the window and scored 0.88,
    while the same spectrum with a grid reaching 2300 A has 280 pixels in the window
    and scores 0.41 -- despite having MORE good data. Angstroms covered is absolute
    and comparable.
    """
    c = (w > 1500) & (w < 1600)
    gc = g & c
    if not gc.any():
        return 0.0
    dw = np.diff(w)
    dw = np.append(dw, dw[-1]) if len(dw) else np.zeros(len(w))
    return float(dw[gc].sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", default=os.path.join(REPO, "RebinnedSpec_master"))
    ap.add_argument("--new", default=os.path.join(REPO, "RebinnedSpec_master_v2"))
    ap.add_argument("--out", default=os.path.join(REPO, "diagnostics", "v2_vs_master.csv"))
    a = ap.parse_args()

    old_f = set(f for f in os.listdir(a.old) if f.endswith(".fits"))
    new_f = set(f for f in os.listdir(a.new) if f.endswith(".fits"))

    rows = []
    for fn in sorted(old_f & new_f):
        inst = fn[:-5].rsplit("_", 1)[1]
        wo, fo, eo = load(os.path.join(a.old, fn))
        wn, fnx, en = load(os.path.join(a.new, fn))
        go = np.isfinite(fo) & (fo != 0)
        gn = np.isfinite(fnx) & (fnx != 0)

        identical = (len(wo) == len(wn)
                     and np.array_equal(np.nan_to_num(fo), np.nan_to_num(fnx))
                     and np.array_equal(np.nan_to_num(eo), np.nan_to_num(en)))

        f_i = np.interp(wo, wn, np.where(gn, fnx, np.nan), left=np.nan, right=np.nan)
        e_i = np.interp(wo, wn, np.where(gn, en, np.nan), left=np.nan, right=np.nan)
        both = go & (eo > 0) & np.isfinite(f_i) & np.isfinite(e_i) & (e_i > 0)
        if both.sum() > 50:
            snr = (np.median(np.abs(f_i[both] / e_i[both]))
                   / np.median(np.abs(fo[both] / eo[both])))
        else:
            snr = np.nan

        rows.append(dict(
            name=fn[:-5], instrument=inst, status="identical" if identical else "changed",
            pix_old=int(go.sum()), pix_new=int(gn.sum()),
            pix_ratio=round(gn.sum() / max(go.sum(), 1), 3),
            snr_ratio=round(float(snr), 3) if np.isfinite(snr) else np.nan,
            civ_aa_before=round(civ_aa(wo, go), 1), civ_aa_after=round(civ_aa(wn, gn), 1)))

    for fn in sorted(new_f - old_f):
        wn, fnx, en = load(os.path.join(a.new, fn))
        gn = np.isfinite(fnx) & (fnx != 0)
        rows.append(dict(name=fn[:-5], instrument=fn[:-5].rsplit("_", 1)[1],
                         status="new", pix_old=0, pix_new=int(gn.sum()),
                         pix_ratio=np.nan, snr_ratio=np.nan,
                         civ_aa_before=0.0, civ_aa_after=round(civ_aa(wn, gn), 1)))
    for fn in sorted(old_f - new_f):
        wo, fo, eo = load(os.path.join(a.old, fn))
        go = np.isfinite(fo) & (fo != 0)
        rows.append(dict(name=fn[:-5], instrument=fn[:-5].rsplit("_", 1)[1],
                         status="LOST", pix_old=int(go.sum()), pix_new=0,
                         pix_ratio=np.nan, snr_ratio=np.nan,
                         civ_aa_before=round(civ_aa(wo, go), 1), civ_aa_after=0.0))

    d = pd.DataFrame(rows)
    d.to_csv(a.out, index=False)

    print("old: %d spectra   new: %d spectra\n" % (len(old_f), len(new_f)))
    print(d.pivot_table(index="status", columns="instrument", values="name",
                        aggfunc="count", fill_value=0).to_string())

    ch = d[d.status == "changed"]
    if len(ch):
        print("\nchanged spectra (%d):" % len(ch))
        print(ch[["pix_ratio", "snr_ratio"]].describe(
            percentiles=[.1, .5, .9]).round(3).to_string())

    # "usable C IV" = at least 30 A of the 100 A window carries good flux
    MIN_AA = 30.0
    gained = d[(d.civ_aa_before < MIN_AA) & (d.civ_aa_after >= MIN_AA)]
    lost   = d[(d.civ_aa_before >= MIN_AA) & (d.civ_aa_after < MIN_AA)]
    shrunk = d[d.civ_aa_after < d.civ_aa_before - 1.0]
    print("\nC IV window (rest 1500-1600 A), Angstroms carrying good flux:")
    print("  gained usable C IV (>=%.0f A, was less): %d" % (MIN_AA, len(gained)))
    print("  lost usable C IV:                       %d" % len(lost))
    print("  any decrease of more than 1 A:          %d" % len(shrunk))
    if len(shrunk):
        print(shrunk[["name", "instrument", "civ_aa_before", "civ_aa_after",
                      "pix_ratio", "snr_ratio"]].to_string(index=False))
    print("  median A covered: %.1f -> %.1f" % (d.civ_aa_before.median(),
                                                d.civ_aa_after.median()))
    print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
