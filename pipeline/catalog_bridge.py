"""
pipeline/catalog_bridge.py
Written: 2026-05-02

Reads master_catalog_v21.csv and download_manifest.csv, then produces the
per-object input list that rebinning/run_rebin.py expects.

One dict per (common_name, inst_family) pair present in the manifest:

    {
        'name':       common_name,          # str
        'redshift':   best_z,               # float
        'instrument': inst_family,          # 'COS', 'STIS', or 'FOS'
        'data_path':  Path(raw_data_dir) / common_name / inst_family,
    }

Multi-instrument objects (e.g. FOS+STIS) produce one dict per instrument,
matching how the existing pipeline handles them.

Note: master_catalog_v21.csv carries no SDSS plate/MJD/fiber columns, so
fn_sdss will be None for all objects when run_rebin.py calls coadd.rebin().

Usage (library):
    from pipeline.catalog_bridge import build_object_list
    objects = build_object_list('/path/to/raw_data')
    for obj in objects:
        print(obj['name'], obj['instrument'], obj['redshift'])

Usage (verification):
    python -m pipeline.catalog_bridge --data-dir /path/to/raw_data
"""

import argparse
import os
import warnings
from pathlib import Path

import pandas as pd

# ── Paths (relative to repo root, matching download_spectra.py) ──────────────
_REPO_ROOT    = Path(__file__).resolve().parent.parent
_PIPELINE_OUT = _REPO_ROOT / 'pipeline_output'

CATALOG_PATH  = _PIPELINE_OUT / 'master_catalog_v21.csv'
MANIFEST_PATH = _PIPELINE_OUT / 'download_manifest.csv'


# ── Loaders ───────────────────────────────────────────────────────────────────

def _load_catalog():
    """
    Load master_catalog_v21.csv and return only has_civ=True rows.
    Applies the same whitespace-stripping used in download_spectra.py.
    """
    if not CATALOG_PATH.exists():
        raise FileNotFoundError(
            f"Master catalog not found: {CATALOG_PATH}\n"
            "Expected pipeline_output/master_catalog_v21.csv."
        )
    catalog = pd.read_csv(CATALOG_PATH, skipinitialspace=True)
    catalog.columns = catalog.columns.str.strip()
    for col in catalog.select_dtypes(include='object').columns:
        catalog[col] = catalog[col].str.strip()
    mask = catalog['has_civ'].astype(str).str.lower() == 'true'
    return catalog[mask].reset_index(drop=True)


def _load_manifest():
    """
    Load download_manifest.csv and return only status='ok' rows.
    Raises FileNotFoundError if the manifest does not exist (i.e.
    download_spectra.py has not been run yet).
    """
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(
            f"Download manifest not found: {MANIFEST_PATH}\n"
            "Run pipeline/download_spectra.py first."
        )
    manifest = pd.read_csv(MANIFEST_PATH)
    ok = manifest[manifest['status'] == 'ok'].reset_index(drop=True)
    n_bad = len(manifest) - len(ok)
    if n_bad:
        warnings.warn(
            f"{n_bad} manifest row(s) with status != 'ok' excluded.",
            stacklevel=2,
        )
    return ok


# ── Public API ────────────────────────────────────────────────────────────────

def build_object_list(raw_data_dir, verbose=True):
    """
    Build the input list for rebinning/run_rebin.py.

    Parameters
    ----------
    raw_data_dir : str or Path
        Root directory where download_spectra.py wrote the FITS files.
        Must match the --data-dir argument used during download.
        Expected layout: raw_data_dir/{common_name}/{inst_family}/*.fits
    verbose : bool
        Print summary counts and per-object warnings.

    Returns
    -------
    list of dict
        One entry per (common_name, inst_family) pair with keys:
            name, redshift, instrument, data_path
        Pairs whose data_path does not exist or contains no .fits files
        are silently skipped (with a warning when verbose=True).
    """
    raw_data_dir = Path(raw_data_dir)

    catalog  = _load_catalog()
    manifest = _load_manifest()

    # Redshift lookup: common_name → best_z
    z_lookup = dict(zip(catalog['common_name'], catalog['best_z']))

    # One entry per (common_name, inst_family) — ignore per-file duplicates
    pairs = (
        manifest[['common_name', 'inst_family']]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    results  = []
    n_no_z   = 0
    n_no_dir = 0
    n_no_fits = 0

    for _, row in pairs.iterrows():
        name = row['common_name']
        inst = row['inst_family']

        # Redshift lookup
        if name not in z_lookup:
            if verbose:
                warnings.warn(
                    f"Skipping {name}/{inst}: not found in master catalog.",
                    stacklevel=2,
                )
            n_no_z += 1
            continue

        z = z_lookup[name]
        if pd.isna(z):
            if verbose:
                warnings.warn(
                    f"Skipping {name}/{inst}: best_z is NaN.",
                    stacklevel=2,
                )
            n_no_z += 1
            continue

        # Directory validation
        data_path = raw_data_dir / name / inst
        if not data_path.exists():
            if verbose:
                warnings.warn(
                    f"Skipping {name}/{inst}: directory not found: {data_path}",
                    stacklevel=2,
                )
            n_no_dir += 1
            continue

        fits_files = list(data_path.glob('*.fits'))
        if not fits_files:
            if verbose:
                warnings.warn(
                    f"Skipping {name}/{inst}: no .fits files in {data_path}",
                    stacklevel=2,
                )
            n_no_fits += 1
            continue

        results.append({
            'name':       name,
            'redshift':   float(z),
            'instrument': inst,
            'data_path':  data_path,
        })

    if verbose:
        n_skip = n_no_z + n_no_dir + n_no_fits
        print(
            f"catalog_bridge: {len(results)} (object, instrument) pairs ready, "
            f"{n_skip} skipped "
            f"(no_z={n_no_z}, no_dir={n_no_dir}, no_fits={n_no_fits})."
        )

    return results


# ── CLI (verification) ────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Print the object list that catalog_bridge would hand to run_rebin.py. "
            "Use this to verify the catalog/manifest join before running the full pipeline."
        )
    )
    parser.add_argument(
        '--data-dir',
        default=os.environ.get('HST_PAPER_DATA_DIR'),
        help=(
            "Root directory containing downloaded FITS files. "
            "Default: $HST_PAPER_DATA_DIR env var."
        ),
    )
    parser.add_argument(
        '--show-paths', action='store_true',
        help="Print data_path for each entry.",
    )
    args = parser.parse_args()

    if args.data_dir is None:
        parser.error(
            "No data directory specified. Use --data-dir or set HST_PAPER_DATA_DIR."
        )

    objects = build_object_list(args.data_dir, verbose=True)

    print(f"\n{'Name':<35} {'Instrument':<10} {'Redshift':>8}")
    print("-" * 58)
    for obj in objects:
        line = f"{obj['name']:<35} {obj['instrument']:<10} {obj['redshift']:8.4f}"
        if args.show_paths:
            line += f"  {obj['data_path']}"
        print(line)
    print(f"\nTotal: {len(objects)} entries")


if __name__ == '__main__':
    main()
