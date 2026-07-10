# Originally: Trevor Code/Run_ICA_All_Objects.py
# Copied: 2026-04-29 (Phase 2 migration)
# Updated (Phase 5.4): 2026-05-02 — added --rebin-dir CLI arg
# Updated (Phase 5.8): 2026-06-18 — added --master mode (folder-driven object
#   selection for RebinnedSpec_master / master_catalog_v21.csv spectra).
# Changes from source:
#   - Removed sys.path.append hacks
#   - from ica_manual_fix_refactored import -> from ica.manual_fix import
#   - catalog path: "HST_CIV..." -> str(Path(__file__).parent.parent / "Data" / "HST_CIV...")
"""
This script runs ICA on all the objects in the dataset.
It processes each object, applies ICA, and saves the results.
It also logs the progress and any errors encountered during processing.

It uses ica_manual_fix_refactored.py for ICA processing.

There are two object-selection modes:

  * Default (Sulentic) mode: iterates the 207 objects in
    HST_CIV_Sulentic2007_HSLA2018_finalprops.csv and looks each up in the
    rebin directory by Spec_Name/Inst_final. Use this to reproduce the
    original Sulentic-sample ICA run.

  * --master mode: iterates every FITS file present in the rebin directory
    (named {common_name}_{inst}.fits, e.g. from master_catalog_v21.csv),
    independent of the Sulentic catalog. Use this to run ICA over ALL
    rebinned spectra, including COS objects absent from the Sulentic catalog.
    Results -> ICA_Results_Rebin_master.csv, plots -> ICA_Plots_Rebin_master/All/.

Usage:
  # Run ICA on the default Sulentic sample (RebinnedSpec_2026AP):
  python -m ica.run_all_objects

  # Run ICA on newly downloaded + rebinned spectra (Sulentic catalog lookup):
  python -m ica.run_all_objects --rebin-dir /path/to/RebinnedSpec_master

  # Run ICA over ALL rebinned spectra in the master folder (no catalog filter):
  python -m ica.run_all_objects --master

  # Master mode with an explicit rebin directory:
  python -m ica.run_all_objects --master --rebin-dir /path/to/RebinnedSpec_master

  # Override via environment variable:
  export HST_PAPER_REBIN_DIR=/path/to/RebinnedSpec_master
  python -m ica.run_all_objects
"""

import argparse
import os
from pathlib import Path

from ica.manual_fix import ICAManualFixProcessor
import pandas as pd

# Default master rebin directory (where rebinning/run_rebin.py --master-catalog writes).
_DEFAULT_MASTER_REBIN_DIR = Path(__file__).parent.parent / "raw_data" / "RebinnedSpec_master"


def main():
    parser = argparse.ArgumentParser(
        description="Run ICA on all objects in the dataset.")
    parser.add_argument(
        "--rebin-dir", default=None,
        help="Directory containing rebinned FITS files to run ICA on. "
             "Default: HST_PAPER_REBIN_DIR env var, then RebinnedSpec_2026AP "
             "(or raw_data/RebinnedSpec_master when --master is set). "
             "Use RebinnedSpec_master for spectra from master_catalog_v21.csv.")
    parser.add_argument(
        "--master", action="store_true",
        help="Drive object selection from the FITS files in the rebin directory "
             "instead of the Sulentic catalog. Runs ICA over every rebinned "
             "spectrum (including COS objects not in the Sulentic catalog). "
             "Results -> ICA_Results_Rebin_master.csv, "
             "plots -> ICA_Plots_Rebin_master/All/.")
    args = parser.parse_args()

    if args.master:
        rebin_dir = args.rebin_dir if args.rebin_dir is not None else str(_DEFAULT_MASTER_REBIN_DIR)
        processor = ICAManualFixProcessor(
            rebin_path=rebin_dir, master_mode=True,
            output_path="ICA_Plots_Rebin_master")
        print("Initialized ICAManualFixProcessor (master mode).")
        print("Rebinned spec path: %s" % processor.rebin_path)

        names = processor.list_master_objects()
        print("Found %d rebinned spectra to process." % len(names))

        processor.batch_process(
            names, output_file="ICA_Results_Rebin_master.csv",
            plot=False, save_plots=True)
        print("ICA processing completed for all master rebinned spectra.")
        return

    processor = ICAManualFixProcessor(rebin_path=args.rebin_dir)
    print("Initialized ICAManualFixProcessor.")
    print("Rebinned spec path: %s" % processor.rebin_path)

    dat = pd.read_csv(str(Path(__file__).parent.parent / "Data" / "HST_CIV_Sulentic2007_HSLA2018_finalprops.csv"))
    names = dat["Final_Name"].values
    for i in range(len(names)):
        names[i] = "".join(names[i].split(" "))

    processor.batch_process(names, output_file="ICA_Results_Rebin_AP.csv", plot=False, save_plots=True)
    print("ICA processing completed for all objects.")


if __name__ == "__main__":
    main()