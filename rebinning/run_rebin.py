"""
rebinning/run_rebin.py
Merged from: Trevor Code/run_rebin_2026AP_FOS_STIS.py
             Trevor Code/run_rebin_2026AP_HSLA.py
Migration date: 2026-04-30

Merger notes:
  - FOS/STIS loop and HSLA loop combined into one script with --instruments flag.
  - `import HST_coadd2_2026AP` → `from rebinning import coadd`
  - All hardcoded paths replaced with DATA_DIR-relative paths (via --data-dir
    argument or HST_PAPER_DATA_DIR environment variable).
  - Output directory configurable via --output-dir (default: DATA_DIR/RebinnedSpec).
  - SDSS spec directory configurable via --sdss-spec-dir
    (default: DATA_DIR/SDSS_spec/lite).
  - Catalog path configurable via --catalog
    (default: Data/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv relative to repo root).
No algorithmic changes to the rebinning logic.

Usage:
  # Process all instruments (FOS, STIS, HSLA):
  python -m rebinning.run_rebin --data-dir /path/to/data

  # Process only FOS and STIS:
  python -m rebinning.run_rebin --data-dir /path/to/data --instruments FOS STIS

  # Use environment variable instead of --data-dir:
  export HST_PAPER_DATA_DIR=/path/to/data
  python -m rebinning.run_rebin

Expected data directory layout (matching Trevor's SulenticAllData structure):
  DATA_DIR/
    SulenticAllData/
      FOS/<ObjectName>/NecessaryParams.csv + exposure subdirs
      STIS/<ObjectName>/NecessaryParams.csv + exposure subdirs
    HSLA_coadds_wCIV/original/<filename>.fits
    SDSS_spec/lite/<plate>/spec-<plate>-<mjd>-<fiber>.fits
"""

import os
import sys
import glob
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

from rebinning import coadd

warnings.filterwarnings('ignore')

# Default catalog location (ships with the repo)
_REPO_ROOT = Path(__file__).parent.parent
_DEFAULT_CATALOG = _REPO_ROOT / "Data" / "HST_CIV_Sulentic2007_HSLA2018_finalprops.csv"


def isNaN(val):
    return val != val


def run_fos_stis(sul, data_dir, output_dir, sdss_spec_dir):
    """Process all FOS and STIS objects from the catalog."""
    sul_hst = sul[sul["Inst_final"].isin(["FOS", "STIS"])].reset_index(drop=True)

    z         = sul_hst["z"].values
    sdss_name = sul_hst["SDSS_NAME"].values
    plate     = sul_hst["PLATE"].values
    mjd       = sul_hst["MJD"].values
    fiber     = sul_hst["FIBERID"].values
    inst      = sul_hst["Inst_final"].values
    spec_name = sul_hst["Spec_Name"].values

    n_good = 0
    n_fail = 0
    failed = []

    for i in range(len(sul_hst)):
        identifier = spec_name[i]
        instrument = inst[i]
        data_path  = os.path.join(data_dir, "SulenticAllData", instrument) + os.sep

        if not isNaN(sdss_name[i]):
            fn_sdss = "%04d/spec-%04d-%05d-%04d.fits" % (
                plate[i], plate[i], mjd[i], fiber[i])
        else:
            fn_sdss = None

        try:
            coadd.rebin(
                identifier, z[i], instrument, fn_sdss,
                data_path=data_path,
                output_dir=output_dir,
                sdss_spec_dir=sdss_spec_dir)
            print("OK %s: %s" % (instrument, identifier))
            n_good += 1
        except KeyError as e:
            print("FAIL %s %s: KeyError %s" % (instrument, identifier, e))
            n_fail += 1
            failed.append((instrument, identifier, "KeyError", str(e)))
        except TypeError as e:
            print("FAIL %s %s: TypeError %s" % (instrument, identifier, e))
            n_fail += 1
            failed.append((instrument, identifier, "TypeError", str(e)))
        except Exception as e:
            print("FAIL %s %s: %s %s" % (instrument, identifier, type(e).__name__, e))
            n_fail += 1
            failed.append((instrument, identifier, type(e).__name__, str(e)))

    return n_good, n_fail, failed


def run_hsla(sul, data_dir, output_dir, sdss_spec_dir):
    """Process all HSLA objects from the catalog."""
    sul_hsla = sul[sul["Inst_final"] == "HSLA"].reset_index(drop=True)

    z         = sul_hsla["z"].values
    sdss_name = sul_hsla["SDSS_NAME"].values
    plate     = sul_hsla["PLATE"].values
    mjd       = sul_hsla["MJD"].values
    fiber     = sul_hsla["FIBERID"].values
    spec_name = sul_hsla["Spec_Name"].values

    hsla_data_path = os.path.join(data_dir, "HSLA_coadds_wCIV") + os.sep

    n_good = 0
    n_fail = 0
    failed = []

    for i in range(len(sul_hsla)):
        if not isNaN(sdss_name[i]):
            fn_sdss = "%04d/spec-%04d-%05d-%04d.fits" % (
                plate[i], plate[i], mjd[i], fiber[i])
        else:
            fn_sdss = None

        matches = glob.glob(os.path.join(hsla_data_path, "original", "%s*.fits" % spec_name[i]))
        if len(matches) == 0:
            print("SKIP HSLA %s: no file found for prefix '%s'" % (spec_name[i], spec_name[i]))
            n_fail += 1
            failed.append(("HSLA", spec_name[i], "FileNotFound", "no file matching %s*.fits" % spec_name[i]))
            continue
        if len(matches) > 1:
            print("WARN HSLA %s: multiple files found, using first: %s" % (spec_name[i], matches))

        fn_basename = os.path.basename(matches[0])

        try:
            coadd.rebin(
                fn_basename, z[i], "HSLA", fn_sdss,
                data_path=hsla_data_path,
                output_dir=output_dir,
                sdss_spec_dir=sdss_spec_dir)
            print("OK HSLA: %s" % fn_basename)
            n_good += 1
        except KeyError as e:
            print("FAIL HSLA %s: KeyError %s" % (fn_basename, e))
            n_fail += 1
            failed.append(("HSLA", fn_basename, "KeyError", str(e)))
        except TypeError as e:
            print("FAIL HSLA %s: TypeError %s" % (fn_basename, e))
            n_fail += 1
            failed.append(("HSLA", fn_basename, "TypeError", str(e)))
        except Exception as e:
            print("FAIL HSLA %s: %s %s" % (fn_basename, type(e).__name__, e))
            n_fail += 1
            failed.append(("HSLA", fn_basename, type(e).__name__, str(e)))

    return n_good, n_fail, failed


def main():
    parser = argparse.ArgumentParser(
        description="Run the 2026AP rebinning pipeline on all HST objects.")
    parser.add_argument(
        "--data-dir", default=os.environ.get("HST_PAPER_DATA_DIR"),
        help="Root directory containing raw spectral data (FOS, STIS, HSLA, SDSS). "
             "Defaults to HST_PAPER_DATA_DIR environment variable.")
    parser.add_argument(
        "--output-dir", default=None,
        help="Directory to write rebinned FITS output. "
             "Default: DATA_DIR/RebinnedSpec")
    parser.add_argument(
        "--sdss-spec-dir", default=None,
        help="Root directory containing SDSS spec files (plate subdirs live here). "
             "Default: DATA_DIR/SDSS_spec/lite")
    parser.add_argument(
        "--catalog", default=str(_DEFAULT_CATALOG),
        help="Path to HST_CIV_Sulentic2007_HSLA2018_finalprops.csv catalog. "
             "Default: Data/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv in repo.")
    parser.add_argument(
        "--instruments", nargs="+", default=["FOS", "STIS", "HSLA"],
        choices=["FOS", "STIS", "HSLA"],
        help="Which instruments to process. Default: all (FOS STIS HSLA).")
    args = parser.parse_args()

    if args.data_dir is None:
        parser.error(
            "No data directory specified. Use --data-dir or set HST_PAPER_DATA_DIR.")

    data_dir = args.data_dir
    output_dir = args.output_dir or os.path.join(data_dir, "RebinnedSpec")
    sdss_spec_dir = args.sdss_spec_dir or os.path.join(data_dir, "SDSS_spec", "lite")

    os.makedirs(output_dir, exist_ok=True)

    sul = pd.read_csv(args.catalog)

    total_good = 0
    total_fail = 0
    all_failed = []

    if any(i in args.instruments for i in ["FOS", "STIS"]):
        # Filter to whichever of FOS/STIS were requested
        sul_filtered = sul.copy()
        if "FOS" not in args.instruments:
            sul_filtered = sul_filtered[sul_filtered["Inst_final"] != "FOS"]
        if "STIS" not in args.instruments:
            sul_filtered = sul_filtered[sul_filtered["Inst_final"] != "STIS"]
        n_good, n_fail, failed = run_fos_stis(sul_filtered, data_dir, output_dir, sdss_spec_dir)
        total_good += n_good
        total_fail += n_fail
        all_failed.extend(failed)

    if "HSLA" in args.instruments:
        n_good, n_fail, failed = run_hsla(sul, data_dir, output_dir, sdss_spec_dir)
        total_good += n_good
        total_fail += n_fail
        all_failed.extend(failed)

    print("\n--- Run complete ---")
    print("Good: %d   Failed/Skipped: %d" % (total_good, total_fail))
    if all_failed:
        print("Failed/skipped objects:")
        for inst, name, etype, msg in all_failed:
            print("  [%s] %s: %s - %s" % (inst, name, etype, msg))


if __name__ == "__main__":
    main()
