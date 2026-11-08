"""
diagnostics/cos_fix_impact.py
AP: 2026-09-10

How much does Fix A (read every COS detector segment, not just row 0) change?

For every COS (object, instrument) pair, compares what the current reader sees
(`WAVELENGTH[0]` of each x1d table) with what all rows hold. Reports, in the rest
frame: wavelength coverage, coverage inside the ICA fitting range (1260-3000 A),
and the number of exposure rows that would feed the co-add - the latter being the
first-order driver of co-added S/N, which improves roughly as sqrt(N) wherever the
extra rows overlap what is already there.

Exposures are grouped by instrument configuration (OPT_ELEM/CENWAVE) and one
representative per group is opened, as in audit_civ_coverage.py; the row COUNT is
taken from every file, since that is cheap and exact.

Usage:
    python -m diagnostics.cos_fix_impact --raw-data-dir "F:/Richards-data/HST Paper/raw_data"
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
from astropy.io import fits

warnings.simplefilter("ignore")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from diagnostics.audit_civ_coverage import reader_files, config_key  # noqa: E402

ICA_LO, ICA_HI = 1260.0, 3000.0
GRID = np.arange(1000.0, 3200.0, 0.5)


def union_coverage(rows, z):
    """Boolean mask on GRID (rest frame) covered by any of `rows`."""
    cov = np.zeros(len(GRID), bool)
    for w in rows:
        w = np.asarray(w, float) / (1 + z)
        w = w[w > 0]
        if len(w):
            cov |= (GRID >= w.min()) & (GRID <= w.max())
    return cov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-data-dir", required=True)
    ap.add_argument("--out", default=os.path.join(REPO, "diagnostics", "cos_fix_impact.csv"))
    a = ap.parse_args()

    from pipeline.catalog_bridge import build_object_list
    objs = [o for o in build_object_list(a.raw_data_dir, verbose=False)
            if o["instrument"] == "COS"]

    out = []
    for k, o in enumerate(objs):
        name, z = o["name"], float(o["redshift"])
        fns = reader_files(str(o["data_path"]), "COS")
        if not fns:
            continue
        reps, nrows_total, nrows_row0, gratings = {}, 0, 0, set()
        for fn in fns:
            with fits.open(fn, memmap=True) as h:
                gratings.add(h[0].header.get("OPT_ELEM", "?"))
                nrows_total += len(h[1].data)
                nrows_row0 += 1
            reps.setdefault(config_key(fn), fn)

        all_w, row0_w = [], []
        for fn in reps.values():
            with fits.open(fn, memmap=True) as h:
                d = h[1].data
                all_w += [np.asarray(d["WAVELENGTH"][t]) for t in range(len(d))]
                row0_w.append(np.asarray(d["WAVELENGTH"][0]))

        c_all, c_row0 = union_coverage(all_w, z), union_coverage(row0_w, z)
        ica = (GRID >= ICA_LO) & (GRID <= ICA_HI)
        px = GRID[1] - GRID[0]
        out.append(dict(
            name=name, z=z, gratings="+".join(sorted(gratings)),
            n_files=len(fns), rows_now=nrows_row0, rows_fixed=nrows_total,
            cov_now_AA=round(c_row0.sum() * px, 1),
            cov_fixed_AA=round(c_all.sum() * px, 1),
            gain_AA=round((c_all.sum() - c_row0.sum()) * px, 1),
            cov_now_ica_AA=round((c_row0 & ica).sum() * px, 1),
            cov_fixed_ica_AA=round((c_all & ica).sum() * px, 1),
            gain_ica_AA=round(((c_all & ica).sum() - (c_row0 & ica).sum()) * px, 1),
        ))
        if (k + 1) % 40 == 0:
            print("  %d/%d" % (k + 1, len(objs)), flush=True)

    d = pd.DataFrame(out)
    d["affected"] = d.rows_fixed > d.rows_now
    d.to_csv(a.out, index=False)

    n = len(d)
    aff = d[d.affected]
    print("\nCOS pairs with raw files: %d" % n)
    print("  affected by Fix A (x1d has >1 row): %d (%.0f%%)" % (len(aff), 100 * len(aff) / n))
    print("  unaffected (single-row tables):     %d" % (n - len(aff)))
    print("\nExposure rows entering the co-add, over the affected pairs:")
    print("  now %d  ->  fixed %d  (x%.2f)" % (aff.rows_now.sum(), aff.rows_fixed.sum(),
                                               aff.rows_fixed.sum() / aff.rows_now.sum()))
    print("\nRest-frame coverage gain (A), affected pairs:")
    for col, lbl in [("gain_AA", "full spectrum"), ("gain_ica_AA", "inside ICA range 1260-3000")]:
        s = aff[col]
        print("  %-28s median %6.0f   mean %6.0f   max %6.0f   >0 for %d pairs"
              % (lbl, s.median(), s.mean(), s.max(), (s > 0).sum()))
    print("\nFractional coverage increase inside the ICA range:")
    r = (aff.cov_fixed_ica_AA / aff.cov_now_ica_AA.replace(0, np.nan))
    print("  median x%.2f   90th pct x%.2f" % (r.median(), r.quantile(0.9)))
    print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
