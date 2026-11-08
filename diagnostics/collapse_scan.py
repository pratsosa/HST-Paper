"""
diagnostics/collapse_scan.py
AP: 2026-09-20

Why do only ~20 ICA fits collapse when 133 spectra carry negative errors?

For each object this records what actually reaches the ICA fit, rather than what
is in the rebinned file:

  file_max_absflux   largest |flux| in the rebinned spectrum (ICA range 1260-3000 A)
  file_n_negerr      negative-error pixels in the rebinned spectrum
  cont_spec_min      minimum of the continuum that morph2() divides by
  cont_spec_nonpos   how many of its pixels are <= 0
  morph_absmax       largest |morph_coeff| (= cont_morph / cont_spec)
  fit_absmax         largest |flux| in the array handed to ICA_fit
  fit_sum_weights    sum |component weights| the fit returned; 0 means collapse

A fit collapses when one pixel's residual dominates chi-squared, which needs a
large |flux| AND a not-correspondingly-large error. Negative errors are squared
into ivar = 1/err^2, so a huge negative error actually DOWN-weights its pixel --
which is why negative errors alone do not predict collapse.

Usage:
    python -m diagnostics.collapse_scan --rebin-dir RebinnedSpec_master_v2 \
        --names diagnostics/_collapse_names.txt --out diagnostics/collapse_scan.csv
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebin-dir", default="RebinnedSpec_master_v2")
    ap.add_argument("--names", default=None, help="file with one object stem per line")
    ap.add_argument("--out", default=os.path.join(REPO, "diagnostics", "collapse_scan.csv"))
    a = ap.parse_args()

    from ica import spec_morph as sm
    from ica import run_ica
    from ica.manual_fix import ICAManualFixProcessor as P

    cur = {}

    orig_morph = sm.morph2
    def morph_spy(wave, flux, errs, z, Identifier):
        out = orig_morph(wave, flux, errs, z, Identifier)
        try:
            cs = np.asarray(sm.cont_filtered(wave, flux, z, Identifier), float)
            mc = np.asarray(out[1], float)
            cur.setdefault("cont_spec_min", float(np.nanmin(cs)))
            cur.setdefault("cont_spec_absmin", float(np.nanmin(np.abs(cs))))
            cur.setdefault("cont_spec_nonpos", int((cs <= 0).sum()))
            cur.setdefault("morph_absmax", float(np.nanmax(np.abs(mc))))
        except Exception:
            pass
        return out
    sm.morph2 = morph_spy

    orig_fit = run_ica.ICA_fit
    def fit_spy(components, wave, spectrum, ivar, flags, priors):
        r = orig_fit(components, wave, spectrum, ivar, flags, priors)
        try:
            s = np.asarray(spectrum, float)
            w = np.abs([r.params[i].value for i in r.params]).sum()
            cur.setdefault("fit_absmax", float(np.nanmax(np.abs(s))))
            cur.setdefault("fit_redchi", float(r.redchi) if r.redchi is not None else np.nan)
            cur["fit_sum_weights"] = min(cur.get("fit_sum_weights", np.inf), float(w))
        except Exception:
            pass
        return r
    run_ica.ICA_fit = fit_spy

    proc = P(rebin_path=a.rebin_dir, master_mode=True, output_path="_collapse_plots")
    if a.names:
        names = [l.strip() for l in open(a.names, encoding="utf-8") if l.strip()]
    else:
        names = proc.list_master_objects()

    rows = []
    for k, nm in enumerate(names):
        cur.clear()
        p = os.path.join(a.rebin_dir, nm + ".fits")
        rec = dict(name=nm)
        try:
            d = fits.getdata(p, 1)
            w = np.asarray(d["Rest-frame Wavelength"], float)
            f = np.asarray(d["Coadded Flux (Arbitrary Units)"], float)
            e = np.asarray(d["Coadded Flux Errors"], float)
            ica = (w >= 1260) & (w <= 3000) & np.isfinite(f)
            rec["file_max_absflux"] = float(np.nanmax(np.abs(f[ica]))) if ica.any() else np.nan
            rec["file_n_negerr"] = int((e < 0).sum())
        except Exception as ex:
            rec["error"] = str(ex)[:80]
        try:
            res = proc.process_object(nm, plot=False, save_plot=False)
            rec["CIV_blueshift"] = res.get("CIV_blueshift")
            rec["CIV_EW"] = res.get("CIV_EW")
            rec["fit_error"] = ""
        except Exception as ex:
            rec["fit_error"] = "%s: %s" % (type(ex).__name__, str(ex)[:60])
        rec.update(cur)
        rows.append(rec)
        if (k + 1) % 25 == 0:
            print("  %d/%d" % (k + 1, len(names)), flush=True)

    d = pd.DataFrame(rows)
    b = pd.to_numeric(d.get("CIV_blueshift"), errors="coerce")
    ew = pd.to_numeric(d.get("CIV_EW"), errors="coerce")
    d["outcome"] = np.where(d.fit_error.astype(str).str.len() > 0, "no fit",
                    np.where((b.round(1) == 9574.1) | ew.isna() | (ew == 0),
                             "DEGENERATE", "fit OK"))
    d.to_csv(a.out, index=False)
    print("\nwrote %s (%d rows)" % (a.out, len(d)))
    print(d.outcome.value_counts().to_string())


if __name__ == "__main__":
    main()
