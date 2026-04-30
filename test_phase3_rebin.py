"""
Phase 3 Integration Test — Rebinning
Run 3 test objects (1 FOS, 1 STIS, 1 HSLA), compare output to RebinnedSpec_2026AP.
"""

import os
import sys
import numpy as np
from pathlib import Path
from astropy.io import fits

# Paths
REPO_ROOT   = Path(__file__).parent
TREVOR      = Path("c:/Users/Alexandros Pratsos/PycharmProjects/Summer 2024/Richards/Trevor Code")
DATA_DIR    = TREVOR          # SulenticAllData/ and HSLA_coadds_wCIV/ live here
REF_DIR     = TREVOR / "RebinnedSpec_2026AP"
OUTPUT_DIR  = REPO_ROOT / "RebinnedSpec_Phase3_Test"
SDSS_DIR    = TREVOR / "SDSS_spec" / "lite"
CATALOG     = REPO_ROOT / "Data" / "HST_CIV_Sulentic2007_HSLA2018_finalprops.csv"

os.makedirs(OUTPUT_DIR, exist_ok=True)

from rebinning import coadd
import pandas as pd

sul = pd.read_csv(CATALOG)

# Test objects: one FOS, one STIS, one HSLA
TEST_OBJECTS = [
    ("FOS",  "J04232-0120"),   # no SDSS
    ("STIS", "J07086-4933"),   # no SDSS
    ("HSLA", "1H1613-097"),    # no SDSS
]

results = []

for instrument, spec_name in TEST_OBJECTS:
    row = sul[sul["Spec_Name"] == spec_name].iloc[0]
    z = row["z"]
    sdss_name = row["SDSS_NAME"]
    fn_sdss = None
    if sdss_name == sdss_name:  # not NaN
        fn_sdss = "%04d/spec-%04d-%05d-%04d.fits" % (
            int(row["PLATE"]), int(row["PLATE"]), int(row["MJD"]), int(row["FIBERID"]))

    if instrument in ("FOS", "STIS"):
        data_path = str(DATA_DIR / "SulenticAllData" / instrument) + os.sep
        ident = spec_name
        out_name = spec_name + "_" + instrument + ".fits"
    else:  # HSLA
        import glob
        hsla_path = str(DATA_DIR / "HSLA_coadds_wCIV") + os.sep
        matches = glob.glob(os.path.join(hsla_path, "original", "%s*.fits" % spec_name))
        if not matches:
            print(f"SKIP HSLA {spec_name}: no file found")
            continue
        ident = os.path.basename(matches[0])
        data_path = hsla_path
        out_name = ident

    try:
        coadd.rebin(ident, z, instrument, fn_sdss,
                    data_path=data_path,
                    output_dir=str(OUTPUT_DIR) + os.sep,
                    sdss_spec_dir=str(SDSS_DIR) + os.sep)
        print(f"OK   {instrument:4s} {spec_name}")
        results.append((instrument, spec_name, out_name, "OK"))
    except Exception as e:
        print(f"FAIL {instrument:4s} {spec_name}: {type(e).__name__}: {e}")
        results.append((instrument, spec_name, out_name, f"FAIL: {e}"))

print()
print("=== Comparison vs RebinnedSpec_2026AP ===")

for instrument, spec_name, out_name, status in results:
    if "FAIL" in status:
        print(f"  {spec_name}: skipped (rebin failed)")
        continue

    new_path = OUTPUT_DIR / out_name
    ref_path = REF_DIR / out_name

    if not new_path.exists():
        print(f"  {spec_name}: NEW FILE MISSING {new_path}")
        continue
    if not ref_path.exists():
        print(f"  {spec_name}: REF FILE MISSING {ref_path}")
        continue

    try:
        new = fits.open(new_path)
        ref = fits.open(ref_path)
        new_wav = new[1].data
        ref_wav = ref[1].data
        new_flux = new[2].data
        ref_flux = ref[2].data

        # Align on common wavelength range
        w_min = max(new_wav.min(), ref_wav.min())
        w_max = min(new_wav.max(), ref_wav.max())
        new_mask = (new_wav >= w_min) & (new_wav <= w_max)
        ref_mask = (ref_wav >= w_min) & (ref_wav <= w_max)
        nf = new_flux[new_mask]
        rf = ref_flux[ref_mask]
        min_len = min(len(nf), len(rf))
        nf = nf[:min_len]
        rf = rf[:min_len]

        good = (rf != 0) & np.isfinite(nf) & np.isfinite(rf)
        if good.sum() == 0:
            print(f"  {spec_name}: no valid pixels to compare")
            continue
        ratio = nf[good] / rf[good]
        med_ratio = np.median(ratio)
        rms = np.sqrt(np.mean((ratio - 1)**2))
        print(f"  {spec_name}: n_pix={good.sum():4d}  median(new/ref)={med_ratio:.6f}  RMS_ratio-1={rms:.6f}  wl_new=[{new_wav.min():.1f},{new_wav.max():.1f}]  wl_ref=[{ref_wav.min():.1f},{ref_wav.max():.1f}]")
        new.close()
        ref.close()
    except Exception as e:
        print(f"  {spec_name}: comparison error: {e}")
