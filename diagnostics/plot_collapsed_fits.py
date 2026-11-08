"""
diagnostics/plot_collapsed_fits.py
AP: 2026-09-20

Visual check for the ICA fits that collapse to an all-zero reconstruction.

Question it answers: are those objects simply bad data we can drop, or is the
spectrum fine and the fit ruined by a few pathological pixels?

Mechanism under test: get_ICA passes ivar = 1/errs**2 to ICA_fit. A pixel whose
error is ~0 gets ivar ~1e29, so any model value differing from the data there
costs ~1e30 in chi-squared. The optimiser is pinned at all-zero weights, because
zero really is the least-squares optimum under that weighting. Collapsed objects
have median max(ivar) ~1.4e29; healthy controls ~8.4.

Each panel shows, per object:
  grey   the morphed spectrum handed to the ICA fit
  red    the pixels whose ivar exceeds CLIP_FACTOR x median(ivar) -- the culprits
  blue   the ICA reconstruction obtained with ONLY those pixels de-weighted

If the blue curve tracks the grey spectrum, the data is good and the failure is
the pixels. If the grey spectrum is itself junk, the object is genuinely bad and
dropping it is justified.

Usage:
    python -m diagnostics.plot_collapsed_fits
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

# Pre-validated steps from the dataviz reference palette (light surface).
SURFACE = "#fcfcfb"
INK = "#1a1a19"
MUTED = "#8a8a86"
SPEC = "#55555260"
FIT = "#2a78d6"      # categorical slot 1
BAD = "#d03b3b"      # status: critical -- these pixels are a fault, not a series
CLIP_FACTOR = 1e4    # ivar above this multiple of the median is not a real error


def main():
    from ica import run_ica
    from ica.manual_fix import ICAManualFixProcessor as P

    stash = {}

    orig_fit = run_ica.ICA_fit
    def fit_clipped(components, wave, spectrum, ivar, flags, priors):
        iv = np.asarray(ivar, float)
        good = np.isfinite(iv) & (iv > 0)
        if good.any():
            lim = np.nanmedian(iv[good]) * CLIP_FACTOR
            bad = good & (iv > lim)
            if bad.any():
                w = np.asarray(wave, float)
                stash.setdefault("bad_wave", w[bad])
                stash["n_bad"] = max(stash.get("n_bad", 0), int(bad.sum()))
                stash["max_ivar"] = max(stash.get("max_ivar", 0.0), float(np.nanmax(iv[bad])))
                iv = iv.copy()
                iv[bad] = 0.0
        return orig_fit(components, wave, spectrum, iv, flags, priors)
    run_ica.ICA_fit = fit_clipped

    orig_civ = P.get_CIV_parameters
    def civ_spy(self, wave, flux, wave_r, flux_r, name, *a, **k):
        stash["wave"] = np.asarray(wave, float)
        stash["flux"] = np.asarray(flux, float)
        stash["wave_r"] = np.asarray(wave_r, float)
        stash["flux_r"] = np.asarray(flux_r, float)
        return orig_civ(self, wave, flux, wave_r, flux_r, name, *a, **k)
    P.get_CIV_parameters = civ_spy

    d = pd.read_csv(os.path.join(REPO, "diagnostics", "degenerate_causes.csv"))
    d["nfit"] = pd.to_numeric(d.nfit, errors="coerce")
    d["sumw"] = pd.to_numeric(d.sumw, errors="coerce")
    targets = d[(d.outcome == "DEGENERATE") & (d.nfit >= 100) &
                ((d.sumw == 0) | d.sumw.isna())].name.head(6).tolist()

    proc = P(rebin_path=os.path.join(REPO, "RebinnedSpec_master_v2"),
             master_mode=True, output_path="_collapse_vis")

    fig, axes = plt.subplots(3, 2, figsize=(13, 10.5), facecolor=SURFACE)
    rows = []
    for ax, nm in zip(axes.ravel(), targets):
        stash.clear()
        try:
            proc.process_object(nm, plot=False, save_plot=False)
        except Exception as e:
            ax.text(.5, .5, "%s\nfit raised %s" % (nm, type(e).__name__),
                    ha="center", va="center", transform=ax.transAxes, color=MUTED)
            ax.set_facecolor(SURFACE)
            continue

        w, f = stash.get("wave"), stash.get("flux")
        wr, fr = stash.get("wave_r"), stash.get("flux_r")
        ax.set_facecolor(SURFACE)
        if w is None:
            ax.text(.5, .5, nm + "\nno spectrum captured", ha="center", va="center",
                    transform=ax.transAxes, color=MUTED)
            continue

        ax.plot(w, f, "-", color=SPEC, lw=0.9, label="spectrum fed to the fit")
        if fr is not None:
            ax.plot(wr, fr, "-", color=FIT, lw=2.0,
                    label="ICA fit, culprit pixels de-weighted")
        bw = stash.get("bad_wave")
        if bw is not None and len(bw):
            lo, hi = ax.get_ylim()
            ax.plot(bw, np.full(len(bw), lo + 0.02 * (hi - lo)), "|", color=BAD,
                    ms=9, mew=1.6, label="pixels with ivar > 1e4 x median")

        finite = np.isfinite(f)
        if finite.any():
            ax.set_ylim(np.nanpercentile(f[finite], 0.5),
                        np.nanpercentile(f[finite], 99.7) * 1.12)
        ax.set_xlim(1260, 3000)
        ax.set_title("%s   (%d culprit pixels, max ivar %.1e)"
                     % (nm, stash.get("n_bad", 0), stash.get("max_ivar", np.nan)),
                     fontsize=10, color=INK, loc="left")
        ax.grid(True, color="#00000010", lw=0.6)
        for s in ax.spines.values():
            s.set_color("#00000022")
        ax.tick_params(colors=MUTED, labelsize=8)
        rows.append(dict(name=nm, n_bad=stash.get("n_bad", 0),
                         max_ivar=stash.get("max_ivar", np.nan)))

    axes.ravel()[0].legend(loc="upper right", fontsize=8, frameon=False,
                           labelcolor=INK)
    fig.suptitle("Collapsed ICA fits: is the spectrum bad, or only a few pixels?",
                 fontsize=13, color=INK, x=0.01, ha="left")
    fig.supxlabel("Rest-frame wavelength (A)", fontsize=10, color=MUTED)
    fig.tight_layout(rect=[0, 0.01, 1, 0.97])
    out = os.path.join(REPO, "diagnostics", "collapsed_fits_visual_check.png")
    fig.savefig(out, dpi=140, facecolor=SURFACE)
    print("wrote %s" % out)
    if rows:
        print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
