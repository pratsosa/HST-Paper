"""
rebinning/run_rebin.py
Merged from: Trevor Code/run_rebin_2026AP_FOS_STIS.py
             Trevor Code/run_rebin_2026AP_HSLA.py
Migration date: 2026-04-30
Updated (Phase 5.4): 2026-05-02 — added --master-catalog mode via catalog_bridge

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
  - --master-catalog flag: when set, uses pipeline.catalog_bridge.build_object_list()
    instead of the Sulentic catalog layout. Default output dir becomes
    DATA_DIR/RebinnedSpec_master in this mode.
No algorithmic changes to the rebinning logic.

Usage:
  # Process all instruments (FOS, STIS, HSLA) using old Sulentic catalog:
  python -m rebinning.run_rebin --data-dir /path/to/data

  # Process only FOS and STIS:
  python -m rebinning.run_rebin --data-dir /path/to/data --instruments FOS STIS

  # Use master_catalog_v21.csv + downloaded FITS files (Phase 5):
  python -m rebinning.run_rebin --master-catalog --data-dir /path/to/downloads

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

Expected layout when --master-catalog is used:
  DATA_DIR/
    {common_name}/COS/*.fits
    {common_name}/STIS/*.fits
    {common_name}/FOS/*.fits
"""

import os
import sys
import glob
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from rebinning import coadd

warnings.filterwarnings('ignore')

# Default catalog location (ships with the repo)
_REPO_ROOT = Path(__file__).parent.parent
_DEFAULT_CATALOG = _REPO_ROOT / "Data" / "HST_CIV_Sulentic2007_HSLA2018_finalprops.csv"


def isNaN(val):
    return val != val


def output_filename(identifier, data_origin):
    """
    Basename that coadd.rebin() will write for this object.

    NOTE: This duplicates the naming convention in coadd.rebin() (see
    coadd.py, the final `t.write(...)` block). If that naming ever changes,
    update this helper to match, or the skip-existing check will look for the
    wrong file.
    """
    if data_origin == "HSLA" and identifier.endswith(".fits"):
        identifier = identifier[:-len(".fits")]
    return "%s_%s.fits" % (identifier, data_origin)


def _process_one(task):
    """
    Worker: rebin a single object.

    `task` is a plain dict of picklable primitives (so it survives the
    spawn-based ProcessPoolExecutor on Windows). Returns a result tuple
    (status, instrument, identifier, etype, msg) where status is "ok" or
    "fail". All exceptions are caught here so a single bad object never
    kills the pool.
    """
    try:
        coadd.rebin(
            task["identifier"], task["z"], task["instrument"], task["fn_sdss"],
            data_path=task["data_path"],
            output_dir=task["output_dir"],
            sdss_spec_dir=task["sdss_spec_dir"],
            **task.get("extra_kwargs", {}))
        return ("ok", task["instrument"], task["identifier"], None, None)
    except Exception as e:
        return ("fail", task["instrument"], task["identifier"],
                type(e).__name__, str(e))


def _dispatch(tasks, workers):
    """
    Run a list of rebin tasks serially (workers <= 1) or across a process
    pool (workers > 1), then print per-object OK/FAIL lines in input order
    and tally results.

    Returns (n_good, n_fail, failed) where failed is a list of
    (instrument, identifier, etype, msg) tuples.
    """
    if workers and workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(_process_one, tasks))
    else:
        results = [_process_one(t) for t in tasks]

    n_good = 0
    n_fail = 0
    failed = []
    for status, inst, name, etype, msg in results:
        if status == "ok":
            print("OK %s: %s" % (inst, name))
            n_good += 1
        else:
            print("FAIL %s %s: %s %s" % (inst, name, etype, msg))
            n_fail += 1
            failed.append((inst, name, etype, msg))
    return n_good, n_fail, failed


def run_from_catalog_bridge(objects, output_dir, sdss_spec_dir, skip_existing=True,
                            workers=1):
    """
    Process objects built by pipeline.catalog_bridge.build_object_list().

    Each entry in `objects` is a dict with keys:
        name, redshift, instrument, data_path

    All objects use fn_sdss=None (master catalog has no SDSS plate/MJD/fiber).
    """
    n_skip = 0
    tasks = []

    for obj in objects:
        name       = obj['name']
        z          = obj['redshift']
        instrument = obj['instrument']
        data_path  = str(obj['data_path']) + os.sep

        out_path = os.path.join(output_dir, output_filename(name, instrument))
        if skip_existing and os.path.exists(out_path):
            print("SKIP %s %s: output exists (%s)" % (
                instrument, name, os.path.basename(out_path)))
            n_skip += 1
            continue

        tasks.append({
            "identifier": name, "z": z, "instrument": instrument,
            "fn_sdss": None, "data_path": data_path,
            "output_dir": output_dir, "sdss_spec_dir": sdss_spec_dir,
            "extra_kwargs": {"flat": True},
        })

    n_good, n_fail, failed = _dispatch(tasks, workers)
    return n_good, n_fail, n_skip, failed


def run_fos_stis(sul, data_dir, output_dir, sdss_spec_dir, skip_existing=True,
                 workers=1):
    """Process all FOS and STIS objects from the catalog."""
    sul_hst = sul[sul["Inst_final"].isin(["FOS", "STIS"])].reset_index(drop=True)

    z         = sul_hst["z"].values
    sdss_name = sul_hst["SDSS_NAME"].values
    plate     = sul_hst["PLATE"].values
    mjd       = sul_hst["MJD"].values
    fiber     = sul_hst["FIBERID"].values
    inst      = sul_hst["Inst_final"].values
    spec_name = sul_hst["Spec_Name"].values

    n_skip = 0
    tasks = []

    for i in range(len(sul_hst)):
        identifier = spec_name[i]
        instrument = inst[i]
        data_path  = os.path.join(data_dir, "SulenticAllData", instrument) + os.sep

        out_path = os.path.join(output_dir, output_filename(identifier, instrument))
        if skip_existing and os.path.exists(out_path):
            print("SKIP %s %s: output exists (%s)" % (
                instrument, identifier, os.path.basename(out_path)))
            n_skip += 1
            continue

        if not isNaN(sdss_name[i]):
            fn_sdss = "%04d/spec-%04d-%05d-%04d.fits" % (
                plate[i], plate[i], mjd[i], fiber[i])
        else:
            fn_sdss = None

        tasks.append({
            "identifier": identifier, "z": z[i], "instrument": instrument,
            "fn_sdss": fn_sdss, "data_path": data_path,
            "output_dir": output_dir, "sdss_spec_dir": sdss_spec_dir,
        })

    n_good, n_fail, failed = _dispatch(tasks, workers)
    return n_good, n_fail, n_skip, failed


def run_hsla(sul, data_dir, output_dir, sdss_spec_dir, skip_existing=True,
             workers=1):
    """Process all HSLA objects from the catalog."""
    sul_hsla = sul[sul["Inst_final"] == "HSLA"].reset_index(drop=True)

    z         = sul_hsla["z"].values
    sdss_name = sul_hsla["SDSS_NAME"].values
    plate     = sul_hsla["PLATE"].values
    mjd       = sul_hsla["MJD"].values
    fiber     = sul_hsla["FIBERID"].values
    spec_name = sul_hsla["Spec_Name"].values

    hsla_data_path = os.path.join(data_dir, "HSLA_coadds_wCIV") + os.sep

    n_fail = 0
    n_skip = 0
    failed = []
    tasks = []

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

        out_path = os.path.join(output_dir, output_filename(fn_basename, "HSLA"))
        if skip_existing and os.path.exists(out_path):
            print("SKIP HSLA %s: output exists (%s)" % (
                fn_basename, os.path.basename(out_path)))
            n_skip += 1
            continue

        tasks.append({
            "identifier": fn_basename, "z": z[i], "instrument": "HSLA",
            "fn_sdss": fn_sdss, "data_path": hsla_data_path,
            "output_dir": output_dir, "sdss_spec_dir": sdss_spec_dir,
        })

    n_good, n_fail_dispatch, failed_dispatch = _dispatch(tasks, workers)
    n_fail += n_fail_dispatch
    failed.extend(failed_dispatch)
    return n_good, n_fail, n_skip, failed


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
             "Default (Sulentic mode): DATA_DIR/RebinnedSpec  "
             "Default (--master-catalog mode): DATA_DIR/RebinnedSpec_master")
    parser.add_argument(
        "--sdss-spec-dir", default=None,
        help="Root directory containing SDSS spec files (plate subdirs live here). "
             "Default: DATA_DIR/SDSS_spec/lite")
    parser.add_argument(
        "--catalog", default=str(_DEFAULT_CATALOG),
        help="Path to HST_CIV_Sulentic2007_HSLA2018_finalprops.csv catalog. "
             "Default: Data/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv in repo. "
             "Ignored when --master-catalog is set.")
    parser.add_argument(
        "--master-catalog", action="store_true",
        help="Use pipeline.catalog_bridge.build_object_list() instead of the "
             "Sulentic catalog. Reads master_catalog_v21.csv and download_manifest.csv "
             "from pipeline_output/. Output defaults to DATA_DIR/RebinnedSpec_master.")
    parser.add_argument(
        "--instruments", nargs="+", default=["FOS", "STIS", "HSLA", "COS"],
        choices=["FOS", "STIS", "HSLA", "COS"],
        help="Which instruments to process. Default: all (FOS STIS HSLA COS). "
             "Ignored when --master-catalog is set (catalog_bridge controls the list).")
    parser.add_argument(
        "--no-skip-existing", dest="skip_existing", action="store_false",
        help="Reprocess every object even if its rebinned output already exists. "
             "By default, objects whose output FITS file is already present in the "
             "output directory are skipped.")
    parser.set_defaults(skip_existing=True)
    parser.add_argument(
        "-j", "--workers", type=int, default=1,
        help="Number of objects to rebin in parallel using a process pool. "
             "Default 1 (serial). Each object is independent and writes its own "
             "output file, so values up to your CPU core count are safe.")
    args = parser.parse_args()

    if args.data_dir is None:
        parser.error(
            "No data directory specified. Use --data-dir or set HST_PAPER_DATA_DIR.")

    data_dir = args.data_dir
    sdss_spec_dir = args.sdss_spec_dir or os.path.join(data_dir, "SDSS_spec", "lite")

    total_good = 0
    total_fail = 0
    total_skip = 0
    all_failed = []

    if args.master_catalog:
        from pipeline.catalog_bridge import build_object_list
        output_dir = args.output_dir or os.path.join(data_dir, "RebinnedSpec_master")
        os.makedirs(output_dir, exist_ok=True)
        print("Master-catalog mode: output → %s" % output_dir)
        objects = build_object_list(data_dir, verbose=True)
        n_good, n_fail, n_skip, failed = run_from_catalog_bridge(
            objects, output_dir, sdss_spec_dir, skip_existing=args.skip_existing,
            workers=args.workers)
        total_good += n_good
        total_fail += n_fail
        total_skip += n_skip
        all_failed.extend(failed)
    else:
        output_dir = args.output_dir or os.path.join(data_dir, "RebinnedSpec")
        os.makedirs(output_dir, exist_ok=True)
        sul = pd.read_csv(args.catalog)

        if any(i in args.instruments for i in ["FOS", "STIS"]):
            sul_filtered = sul.copy()
            if "FOS" not in args.instruments:
                sul_filtered = sul_filtered[sul_filtered["Inst_final"] != "FOS"]
            if "STIS" not in args.instruments:
                sul_filtered = sul_filtered[sul_filtered["Inst_final"] != "STIS"]
            n_good, n_fail, n_skip, failed = run_fos_stis(
                sul_filtered, data_dir, output_dir, sdss_spec_dir,
                skip_existing=args.skip_existing, workers=args.workers)
            total_good += n_good
            total_fail += n_fail
            total_skip += n_skip
            all_failed.extend(failed)

        if "HSLA" in args.instruments:
            n_good, n_fail, n_skip, failed = run_hsla(
                sul, data_dir, output_dir, sdss_spec_dir,
                skip_existing=args.skip_existing, workers=args.workers)
            total_good += n_good
            total_fail += n_fail
            total_skip += n_skip
            all_failed.extend(failed)

    print("\n--- Run complete ---")
    print("Good: %d   Skipped (already existed): %d   Failed: %d" % (
        total_good, total_skip, total_fail))
    if all_failed:
        print("Failed objects:")
        for inst, name, etype, msg in all_failed:
            print("  [%s] %s: %s - %s" % (inst, name, etype, msg))


if __name__ == "__main__":
    main()
