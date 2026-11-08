"""
diagnostics/plot_guard_rescued.py
AP: 2026-09-23

Visual triage of the 29 fits rescued by the ivar guard in ICA_fit.

These objects previously returned an all-zero reconstruction (CIV_blueshift = 9574.1)
because a few pixels with near-zero error claimed infinite precision and pinned the
fit. The guard de-weights those pixels. This plots what the fit now looks like across
the C IV measurement region, so each can be accepted or rejected by eye.

Panels are sorted by median S/N, lowest first, and show:
  grey    the morphed spectrum the fit saw
  blue    the ICA reconstruction
  band    the 1500-1600 A window the equivalent width is integrated over
  dashed  C IV rest wavelength, 1549.48 A

Usage:
    python -m diagnostics.plot_guard_rescued
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

SURFACE = "#fcfcfb"
INK = "#1a1a19"
MUTED = "#8a8a86"
SPEC = "#55555255"
FIT = "#2a78d6"
BAND = "#2a78d61a"
CIV = 1549.48


def main():
    from ica.manual_fix import ICAManualFixProcessor as P

    stash = {}
    orig = P.get_CIV_parameters
    def spy(self, wave, flux, wave_r, flux_r, name, *a, **k):
        stash["w"] = np.asarray(wave, float)
        stash["f"] = np.asarray(flux, float)
        stash["wr"] = np.asarray(wave_r, float)
        stash["fr"] = np.asarray(flux_r, float)
        return orig(self, wave, flux, wave_r, flux_r, name, *a, **k)
    P.get_CIV_parameters = spy

    t = pd.read_csv(os.path.join(REPO, "diagnostics", "guard_rescued.csv"))
    proc = P(rebin_path=os.path.join(REPO, "RebinnedSpec_master_v2"),
             master_mode=True, output_path="_guard_vis")

    n = len(t)
    ncol, nrow = 5, int(np.ceil(n / 5))
    fig, axes = plt.subplots(nrow, ncol, figsize=(19, 3.0 * nrow), facecolor=SURFACE)
    axes = np.atleast_1d(axes).ravel()

    for ax, (_, r) in zip(axes, t.iterrows()):
        nm = r["name"]
        stash.clear()
        ax.set_facecolor(SURFACE)
        try:
            proc.process_object(nm, plot=False, save_plot=False)
        except Exception as e:
            ax.text(.5, .5, "%s\n%s" % (nm, type(e).__name__), ha="center",
                    va="center", transform=ax.transAxes, fontsize=7, color=MUTED)
            continue
        w, f, wr, fr = stash.get("w"), stash.get("f"), stash.get("wr"), stash.get("fr")
        if w is None:
            continue
        ax.axvspan(1500, 1600, color=BAND, lw=0)
        ax.plot(w, f, "-", color=SPEC, lw=0.8)
        ax.plot(wr, fr, "-", color=FIT, lw=1.7)
        ax.axvline(CIV, color=INK, ls="--", lw=0.8, alpha=0.5)
        ax.set_xlim(1440, 1710)

        sel = (w > 1440) & (w < 1710) & np.isfinite(f)
        selr = (wr > 1440) & (wr < 1710) & np.isfinite(fr)
        if sel.sum() > 5:
            lo = min(np.nanpercentile(f[sel], 2), np.nanmin(fr[selr]) if selr.any() else 0)
            hi = max(np.nanpercentile(f[sel], 98), np.nanmax(fr[selr]) if selr.any() else 1)
            pad = 0.12 * (hi - lo + 1e-9)
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_title("%s\nS/N %.1f  cov %d%%  EW %.1f  blue %.0f"
                     % (nm[:30], r["snr"], r["cov"], r["EW"], r["blue"]),
                     fontsize=7.5, color=INK, loc="left")
        ax.tick_params(colors=MUTED, labelsize=6.5)
        ax.grid(True, color="#00000010", lw=0.5)
        for s in ax.spines.values():
            s.set_color("#00000022")

    for ax in axes[n:]:
        ax.set_visible(False)

    fig.suptitle("Fits rescued by the ivar guard - C IV region, sorted by S/N "
                 "(grey = spectrum, blue = ICA fit, band = EW window)",
                 fontsize=12, color=INK, x=0.01, ha="left")
    fig.supxlabel("Rest-frame wavelength (A)", fontsize=9, color=MUTED)
    fig.tight_layout(rect=[0, 0.012, 1, 0.975])
    out = os.path.join(REPO, "diagnostics", "guard_rescued_visual_check.png")
    fig.savefig(out, dpi=130, facecolor=SURFACE)
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
