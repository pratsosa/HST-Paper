"""
diagnostics/compare_rebins.py
AP: 2026-09-10

Side-by-side comparison of two rebinned-spectrum directories. Reports, per object:
coverage of good (finite, non-zero) flux, the fraction of the rest-frame C IV core
(1500-1600 A) that is good, and - over the wavelengths both versions have in common -
the median ratio and scatter of the fluxes. A median ratio near 1 with small scatter
means the change did not disturb data that was already fine.

Usage:
    python -m diagnostics.compare_rebins RebinnedSpec_master RebinnedSpec_master_patched
"""

import argparse
import os

import numpy as np
from astropy.io import fits


def load(p):
    d = fits.getdata(p, 1)
    w = np.asarray(d["Rest-frame Wavelength"], float)
    f = np.asarray(d["Coadded Flux (Arbitrary Units)"], float)
    return w, f, np.isfinite(f) & (f != 0)


def describe(w, g):
    if not g.any():
        return "        all-NaN       "
    return "%7.0f-%-7.0f" % (w[g].min(), w[g].max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir_a")
    ap.add_argument("dir_b")
    a = ap.parse_args()

    names = sorted(set(os.listdir(a.dir_b)) & set(f for f in os.listdir(a.dir_a)))
    only_b = sorted(set(os.listdir(a.dir_b)) - set(os.listdir(a.dir_a)))

    print("%-30s %-16s %-16s %6s %6s %8s %7s"
          % ("object", "A good (rest)", "B good (rest)", "civ_A", "civ_B",
             "B/A med", "scatter"))
    for n in names:
        wa, fa, ga = load(os.path.join(a.dir_a, n))
        wb, fb, gb = load(os.path.join(a.dir_b, n))
        core = lambda w, g: float(g[(w > 1500) & (w < 1600)].mean()) \
            if ((w > 1500) & (w < 1600)).any() else 0.0
        # compare on the shared wavelength grid (both are the same log-lambda scale)
        fb_i = np.interp(wa, wb, np.where(gb, fb, np.nan), left=np.nan, right=np.nan)
        both = ga & np.isfinite(fb_i) & (fa != 0)
        if both.sum() > 20:
            r = fb_i[both] / fa[both]
            med, sc = np.nanmedian(r), np.nanpercentile(np.abs(r - 1), 68)
        else:
            med = sc = np.nan
        print("%-30s %-16s %-16s %6.2f %6.2f %8.3f %7.3f"
              % (n.replace(".fits", ""), describe(wa, ga), describe(wb, gb),
                 core(wa, ga), core(wb, gb), med, sc))

    for n in only_b:
        wb, fb, gb = load(os.path.join(a.dir_b, n))
        print("%-30s %-16s %-16s %6s %6.2f" % (n.replace(".fits", ""), "(absent)",
                                               describe(wb, gb), "-",
                                               float(gb[(wb > 1500) & (wb < 1600)].mean())))


if __name__ == "__main__":
    main()
