"""
diagnostics/ew_bias_scan.py
AP: 2026-09-20

Measures the C IV equivalent-width bias caused by the continuum indexing in
get_CIV_parameters(), for every object, without changing any pipeline code.

The defect (inherited verbatim from Trevor's get_CIV, see
diagnostics/CIV_EW_Indexing_Investigation.md):

    EW = ((wave_r >= 1500) & (wave_r <= 1600))
    for i in range(len(wave_r[EW])):
        CIV_EW += max(((flux_r[EW][i] - continuum[i]) / continuum[i]) * dlam, 0)

`flux_r[EW][i]` indexes the selected window but `continuum[i]` indexes the FULL
array, so the continuum is read from the blue end of the ICA grid (~1260 A)
instead of from 1500-1600 A.

This script patches get_CIV_parameters to compute the EW and blueshift BOTH ways
for each object and writes them side by side. Nothing else is altered: the value
the pipeline actually uses is still the as-coded one, so the fits are unaffected.

Usage:
    python -m diagnostics.ew_bias_scan --rebin-dir RebinnedSpec_master --out diagnostics/ew_bias_master.csv
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter("ignore")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)


def both_ways(wave_r, flux_r, EW_region, cont_region):
    """(EW, blueshift) computed as-coded and with the continuum indexed correctly."""
    wr = np.asarray(wave_r, float)
    fr = np.asarray(flux_r, float)
    c1 = (wr >= cont_region[0][0]) & (wr <= cont_region[0][1])
    c2 = (wr >= cont_region[1][0]) & (wr <= cont_region[1][1])
    m, b = np.polyfit(np.concatenate((wr[c1], wr[c2])),
                      np.concatenate((fr[c1], fr[c2])), 1)
    cont = wr * m + b
    EW = (wr >= EW_region[0]) & (wr <= EW_region[1])
    n = int(EW.sum())
    out = {}
    for tag, c in (("as_coded", cont), ("corrected", cont[EW])):
        ew, cum = 0.0, [0.0]
        for i in range(n):
            dlam = (wr[EW][i + 1] - wr[EW][i]) if i + 1 < n else (wr[EW][i] - wr[EW][i - 1])
            ew += max(((fr[EW][i] - c[i]) / c[i]) * dlam, 0)
            cum.append(ew)
        j = int(np.abs((ew / 2) - np.array(cum)).argmin())
        out[tag] = (ew, ((1549.48 - wr[EW][j]) / 1549.48) * 3e5)
    # context for interpreting the bias
    out["slope"] = float(m)
    out["grid_start"] = float(wr.min())
    out["cont_read_at"] = (float(wr[0]), float(wr[n - 1]))
    out["ew_window_at"] = (float(wr[EW][0]), float(wr[EW][-1]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebin-dir", default="RebinnedSpec_master")
    ap.add_argument("--out", default=os.path.join(REPO, "diagnostics", "ew_bias.csv"))
    a = ap.parse_args()

    from ica.manual_fix import ICAManualFixProcessor as P
    rows = []
    orig = P.get_CIV_parameters

    def spy(self, wave, flux, wave_r, flux_r, name,
            EW_region=[1500, 1600], cont_region=[[1445, 1465], [1700, 1705]]):
        try:
            r = both_ways(wave_r, flux_r, EW_region, cont_region)
            rows.append(dict(
                name=name,
                ew_as_coded=r["as_coded"][0], ew_corrected=r["corrected"][0],
                blue_as_coded=r["as_coded"][1], blue_corrected=r["corrected"][1],
                cont_slope=r["slope"], grid_start=r["grid_start"],
                cont_read_lo=r["cont_read_at"][0], cont_read_hi=r["cont_read_at"][1],
                ew_win_lo=r["ew_window_at"][0], ew_win_hi=r["ew_window_at"][1],
                ew_region_lo=EW_region[0], ew_region_hi=EW_region[1]))
        except Exception as e:
            rows.append(dict(name=name, error="%s: %s" % (type(e).__name__, e)))
        return orig(self, wave, flux, wave_r, flux_r, name, EW_region, cont_region)

    P.get_CIV_parameters = spy
    proc = P(rebin_path=a.rebin_dir, master_mode=True, output_path="_ew_scan_plots")
    names = proc.list_master_objects()
    print("scanning %d spectra in %s" % (len(names), a.rebin_dir), flush=True)
    for k, nm in enumerate(names):
        try:
            proc.process_object(nm, plot=False, save_plot=False)
        except Exception:
            pass
        if (k + 1) % 50 == 0:
            print("  %d/%d  (%d measured)" % (k + 1, len(names), len(rows)), flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(a.out, index=False)
    print("\nwrote %s (%d rows)" % (a.out, len(d)))


if __name__ == "__main__":
    main()
