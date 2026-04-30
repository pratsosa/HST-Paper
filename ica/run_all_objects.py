# Originally: Trevor Code/Run_ICA_All_Objects.py
# Copied: 2026-04-29 (Phase 2 migration)
# Changes from source:
#   - Removed sys.path.append hacks
#   - from ica_manual_fix_refactored import -> from ica.manual_fix import
#   - catalog path: "HST_CIV..." -> str(Path(__file__).parent.parent / "Data" / "HST_CIV...")
"""
This script runs ICA on all the objects in the dataset.
It processes each object, applies ICA, and saves the results.
It also logs the progress and any errors encountered during processing.

It uses ica_manual_fix_refactored.py for ICA processing.
"""

from pathlib import Path
import os

from ica.manual_fix import ICAManualFixProcessor
import pandas as pd

def main():
    # Initialize the processor
    processor = ICAManualFixProcessor()
    print('Initialized ICAManualFixProcessor.')

    # Load the dataset, sanitize names
    dat = pd.read_csv(str(Path(__file__).parent.parent / "Data" / "HST_CIV_Sulentic2007_HSLA2018_finalprops.csv"))
    names = dat["Final_Name"].values
    for i in range(len(names)):
        names[i] = "".join(names[i].split(" "))

    # Batch processing using ICAManualFixProcessor
    # This will save the results to ICA_Results_All_Objects.csv and the plots to the 'ICA_Plots_Refactored' directory
    processor.batch_process(names, output_file="ICA_Results_Rebin_AP.csv", plot=False, save_plots=True)
    print("ICA processing completed for all objects.")

if __name__ == "__main__":
    main()