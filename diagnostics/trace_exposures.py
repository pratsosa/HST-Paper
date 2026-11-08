"""
diagnostics/trace_exposures.py
AP: 2026-09-10

Per-exposure trace of the rebinning pipeline for one (object, instrument) pair.

Reproduces the steps coadd.rebin() performs, but instead of writing a spectrum it
prints, for each exposure, where the data is and whether it survives:

  read      rest-frame range the reader returned for that exposure
  rebin     rest-frame range after HSTLowResRebin
  cont      the per-exposure continuum from spec_morph.cont_filtered, which the
            exposure is divided by. `cont<=0` counts pixels where that division
            destroys the flux (turns it into +/-inf or flips its sign).
  surv      pixels with finite non-zero flux AFTER the continuum division
  CIV       whether rest-frame 1549.48 A falls inside that exposure

Usage:
    python -m diagnostics.trace_exposures "LEDA 3095694" COS \
        --raw-data-dir "F:/Richards-data/HST Paper/raw_data"
"""

import argparse
import os
import sys
import warnings

import numpy as np
from astropy.io import fits

warnings.simplefilter("ignore")

CIV = 1549.48
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from rebinning import read_spec_data, small_pix_filter, lower_res_rebin, spec_morph  # noqa: E402


def trace(name, inst, z, path):
    reader = {"COS": read_spec_data.read_cos_flat,
              "STIS": read_spec_data.read_stis_flat,
              "FOS": read_spec_data.read_fos_flat}[inst]
    waves, fluxes, flux_errs, masks = reader(name, path + os.sep, z)
    print("%s %s  z=%.4f  %d exposures  (C IV rest 1549.48 = %.1f A observed)"
          % (name, inst, z, waves.shape[0], CIV * (1 + z)))

    # Step 1 of coadd.rebin(): noise-spike replacement
    for i in range(waves.shape[0]):
        mf = small_pix_filter.pixel_filter(waves[i], fluxes[i], flux_errs[i], name, inst, "Fluxes")
        me = small_pix_filter.pixel_filter(waves[i], fluxes[i], flux_errs[i], name, inst, "Fluxerrs")
        em = spec_morph.emission_lines2(waves[i] / (1 + z))
        bad = (np.abs(fluxes[i] - mf) > 5 * np.abs(me)) & (~em)
        fluxes[i][bad], flux_errs[i][bad], masks[i][bad] = mf[bad], me[bad], 2

    print("%3s %-19s %-19s %8s %8s %8s  %s"
          % ("i", "read (rest)", "rebin (rest)", "npix", "cont<=0", "surv", "CIV?"))
    for i in range(waves.shape[0]):
        nz = waves[i] > 0
        bw, bf, be, bm = lower_res_rebin.HSTLowResRebin(
            waves[i], fluxes[i], flux_errs[i], masks[i], name, z)
        cont = spec_morph.cont_filtered(bw, bf, z, name)
        norm = bf / cont
        r0, r1 = bw.min() / (1 + z), bw.max() / (1 + z)
        print("%3d %8.1f-%8.1f  %8.1f-%8.1f %8d %8d %8d  %s"
              % (i, waves[i][nz].min() / (1 + z), waves[i][nz].max() / (1 + z), r0, r1,
                 len(bw), int((cont <= 0).sum()), int((np.isfinite(norm) & (norm != 0)).sum()),
                 "yes" if r0 <= CIV <= r1 else ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("instrument", choices=["COS", "STIS", "FOS"])
    ap.add_argument("--raw-data-dir", required=True)
    ap.add_argument("--z", type=float, default=None,
                    help="redshift; looked up in master_catalog_v21.csv if omitted")
    a = ap.parse_args()

    z = a.z
    if z is None:
        from pipeline.catalog_bridge import build_object_list
        for o in build_object_list(a.raw_data_dir, verbose=False):
            if o["name"] == a.name and o["instrument"] == a.instrument:
                z = float(o["redshift"])
                break
        if z is None:
            raise SystemExit("redshift not found for %s/%s; pass --z" % (a.name, a.instrument))

    trace(a.name, a.instrument, z, os.path.join(a.raw_data_dir, a.name, a.instrument))


if __name__ == "__main__":
    main()
