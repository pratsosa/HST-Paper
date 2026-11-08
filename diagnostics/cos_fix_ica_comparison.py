"""
diagnostics/cos_fix_ica_comparison.py
AP: 2026-09-30

Before/after ICA diagnostic plots for every COS spectrum changed by the COS reader
fix (fix 1 in PIPELINE_CHANGES_2026-09.md: read every detector row, not row 0).

Both sides are fit with the *current* ICA code (EW indexing fix + ivar guard,
2026-09-23), so the only difference between a before and an after plot is the
rebinned spectrum:

    before : RebinnedSpec_master     (old reader, row 0 only)
    after  : RebinnedSpec_master_v2  (all rows)

Objects that have a manual override in ica/manual_fix_overrides.json are fit
twice on each side (baseline, and with the override applied unchanged), so
four plots are written:

    {name}_1_before.png
    {name}_2_before_override.png     (override objects only)
    {name}_3_after.png
    {name}_4_after_override.png      (override objects only)

Sorting by filename puts each object's plots together in this order, for
flipping through in an image viewer.

Every fit gets a plot. Where the standard diagnostic plot cannot be drawn, a
simpler FALLBACK PLOT (titled as such) is written under the same filename:
  - "fit ran, standard plot failed": no data redward of 1400 A, so the standard
    plot's y-limits are NaN. The ICA fit ran and its EW/blueshift are recorded,
    but they are template extrapolation, not a measurement.
  - "ICA fit failed": the rebinned spectrum alone, no fit.
The CSV columns plot_kind, error and civ_good_px (good pixels in 1500-1600 A of
the rebinned spectrum) record which case applies.

Selection: COS spectra with status 'changed' in diagnostics/v2_vs_master.csv AND
affected=True in diagnostics/cos_fix_impact.csv (the reader fix added rows). The
20 COS spectra that are 'changed' but not affected differ only by floating-point
noise (max relative flux difference ~1e-9) and are left out. Spectra listed in
diagnostics/affected_unnormalised.csv are skipped: their after spectrum is
corrupted by the fix-2 regression (PIPELINE_CHANGES_2026-09.md section 5) and
must be regenerated first.

Outputs:
    ICA_Plots_COS_fix_comparison/All/*.png
    diagnostics/cos_fix_ica_comparison.csv    one row per fit
    diagnostics/cos_fix_ica_comparison_skipped.csv

Usage (from HST Paper/):
    python -m diagnostics.cos_fix_ica_comparison
    python -m diagnostics.cos_fix_ica_comparison --only "Mrk 1148_COS" "NGC 863_COS"
"""

import argparse
import json
import os
import sys
import time
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from astropy.io import fits  # noqa: E402

warnings.simplefilter("ignore")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from ica.manual_fix import ICAManualFixProcessor, _override_kwargs, _override_tag  # noqa: E402

DIAG = os.path.join(REPO, "diagnostics")
EW_LO, EW_HI = 1500.0, 1600.0
VERSIONS = [  # (label, rebin dir, plot-title description, file-tag prefix)
    ("before", "RebinnedSpec_master", "BEFORE COS fix (RebinnedSpec_master)", "1_before"),
    ("after", "RebinnedSpec_master_v2", "AFTER COS fix (RebinnedSpec_master_v2)", "3_after"),
]


def select_objects():
    """Return (selected names, skipped DataFrame)."""
    v = pd.read_csv(os.path.join(DIAG, "v2_vs_master.csv"))
    c = pd.read_csv(os.path.join(DIAG, "cos_fix_impact.csv"))
    affected = set(c.loc[c.affected, "name"] + "_COS")
    unnorm = set(pd.read_csv(os.path.join(DIAG, "affected_unnormalised.csv")).name)

    changed = v[(v.instrument == "COS") & (v.status == "changed")].name
    cand = sorted(n for n in changed if n in affected)
    selected = [n for n in cand if n not in unnorm]
    skipped = pd.DataFrame({"name": [n for n in cand if n in unnorm],
                            "reason": "in affected_unnormalised.csv (after spectrum corrupted)"})
    return selected, skipped


def read_rebinned(rebin_dir, name):
    """(wave, flux, err, mask) of a rebinned FITS; mask 0 = good pixel."""
    with fits.open(os.path.join(rebin_dir, name + ".fits")) as h:
        d = h[1].data
        return tuple(np.asarray(d[c], float) for c in
                     ("Rest-frame Wavelength", "Coadded Flux (Arbitrary Units)",
                      "Coadded Flux Errors", "Bad Pixel Mask"))


def civ_good_pixels(rebin_dir, name):
    """Good pixels of the rebinned spectrum inside the 1500-1600 A EW window."""
    w, f, _, m = read_rebinned(rebin_dir, name)
    return int(np.sum((w >= EW_LO) & (w <= EW_HI) & np.isfinite(f) & (f != 0) & (m == 0)))


def fallback_plot(path, title, reason, wave, flux, wave_ica=None, flux_ica=None):
    """Simple full-spectrum + C IV-region plot, used when the standard diagnostic
    plot cannot be drawn (no data > 1400 A, or the ICA fit itself failed)."""
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(18, 10.5))
    good = np.isfinite(flux) & (flux != 0)
    for ax, (lo, hi) in ((ax0, (1150, 3100)), (ax1, (1435, 1710))):
        ax.plot(wave[good], flux[good], "-k", alpha=0.6, lw=0.8, label="rebinned spectrum")
        if wave_ica is not None:
            ax.plot(wave_ica, flux_ica, "-r", label="ICA fit")
        ax.axvspan(EW_LO, EW_HI, color="b", alpha=0.08, label="C IV EW window")
        ax.set_xlim(lo, hi)
        sel = good & (wave >= lo) & (wave <= hi)
        if sel.any():
            y = flux[sel]
            ax.set_ylim(min(0, np.nanpercentile(y, 1)), 1.3 * np.nanpercentile(y, 99))
        else:
            ax.text(0.5, 0.5, "no data in %d-%d A" % (lo, hi), transform=ax.transAxes,
                    ha="center", va="center", fontsize=20, color="0.4")
        ax.set_ylabel("Flux Density (Arb. Units)", fontsize=16)
    ax0.legend(loc="upper right")
    ax1.set_xlabel("Rest-frame Wavelength (A)", fontsize=16)
    ax0.set_title("%s\nFALLBACK PLOT - %s" % (title, reason), fontsize=16)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close("all")


def fit_and_plot(proc, name, entry, title, tag):
    """One fit + diagnostic plot; returns a results row. Always writes a plot:
    the standard diagnostic if possible, otherwise a fallback plot."""
    path = os.path.join(proc.output_path, proc._MASTER_FOLDER, "%s_%s.png" % (name, tag))
    row = dict(z=np.nan, CIV_blueshift=np.nan, CIV_EW=np.nan, f2500=np.nan,
               civ_good_px=civ_good_pixels(proc.rebin_path, name), plot_kind="standard",
               error="")
    try:
        res = proc.fit_with_overrides(name, **_override_kwargs(entry))
    except Exception as e:
        plt.close("all")
        w, f, _, _ = read_rebinned(proc.rebin_path, name)
        fallback_plot(path, title, "ICA fit failed: %s" % e, w, f)
        row.update(plot_kind="spectrum_only (fit failed)", error=str(e))
        return row

    row.update(z=res["z"], CIV_blueshift=res["civ_blue"], CIV_EW=res["civ_ew"],
               f2500=res["f2500"])
    try:
        proc.create_diagnostic_plot(
            res["wave_arb"], res["flux_arb"], res["errs_arb"], res["mask_arb"],
            res["wave_ica"], res["flux_ica"], title, res["spec_name"],
            res["civ_blue"], res["civ_ew"], show=False, save=True, file_tag=tag)
        plt.close("all")
    except Exception as e:
        plt.close("all")
        fallback_plot(path, title, "fit ran, standard plot failed (%s); EW=%.1f A, "
                      "blueshift=%.1f km/s" % (e, res["civ_ew"], res["civ_blue"]),
                      np.asarray(res["wave_arb"]), np.asarray(res["flux_arb"]),
                      res["wave_ica"], res["flux_ica"])
        row.update(plot_kind="fallback (standard plot failed)", error=str(e))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="+", default=None, help="restrict to these names")
    ap.add_argument("--plots-dir", default=os.path.join(REPO, "ICA_Plots_COS_fix_comparison"))
    ap.add_argument("--out", default=os.path.join(DIAG, "cos_fix_ica_comparison.csv"))
    a = ap.parse_args()

    names, skipped = select_objects()
    if a.only:
        names = [n for n in names if n in a.only]
    with open(os.path.join(REPO, "ica", "manual_fix_overrides.json")) as f:
        overrides = json.load(f)
    n_ov = sum(n in overrides for n in names)
    print("%d objects (%d with a manual override); %d skipped as un-normalised"
          % (len(names), n_ov, len(skipped)))

    procs = {lab: ICAManualFixProcessor(rebin_path=os.path.join(REPO, d), master_mode=True,
                                        output_path=a.plots_dir)
             for lab, d, _, _ in VERSIONS}
    procs["before"]._ensure_output_dirs()

    rows, errors = [], []
    t0 = time.time()
    for k, name in enumerate(names):
        entry = overrides.get(name)
        for lab, _, desc, prefix in VERSIONS:
            fits_to_do = [(None, "baseline", prefix)]
            if entry:
                fits_to_do.append((entry, _override_tag(entry),
                                   prefix.replace("1_", "2_").replace("3_", "4_") + "_override"))
            for ov, ov_type, tag in fits_to_do:
                title = "%s  |  %s  |  %s" % (name, desc, ov_type)
                try:
                    r = fit_and_plot(procs[lab], name, ov, title, tag)
                except Exception as e:  # even the fallback failed; record and keep going
                    plt.close("all")
                    errors.append(dict(name=name, reason="%s %s: no plot: %s" % (lab, ov_type, e)))
                    print("ERROR %s [%s %s]: %s" % (name, lab, ov_type, e))
                    continue
                if r["error"]:
                    print("FALLBACK %s [%s %s]: %s" % (name, lab, ov_type, r["error"]))
                rows.append(dict(name=name, version=lab, override_type=ov_type,
                                 plot=name + "_" + tag + ".png", **r))
        print("[%d/%d] %s  (%.0f s elapsed)" % (k + 1, len(names), name, time.time() - t0),
              flush=True)

    pd.DataFrame(rows).to_csv(a.out, index=False)
    skipped = pd.concat([skipped, pd.DataFrame(errors)], ignore_index=True)
    skipped.to_csv(a.out.replace(".csv", "_skipped.csv"), index=False)
    print("Wrote %d fits to %s; %d skipped/failed entries" % (len(rows), a.out, len(skipped)))


if __name__ == "__main__":
    main()
