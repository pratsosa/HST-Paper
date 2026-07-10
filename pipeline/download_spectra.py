"""
pipeline/download_spectra.py

Original source: Richards/HST_SDSS_Master_Pipeline_v21_script.py (Cells 1–4 logic)
Copied/adapted: 2026-05-01

Downloads calibrated HST spectra from MAST for every has_civ=True object in
master_catalog_v21.csv. Uses Richards' MAST observation cache files to recover
obs_ids without re-querying MAST, guaranteeing exactly the observations Richards
identified in his v21 pipeline.

Filtering mirrors the v21 script:
  - Coordinate validity:  Cell 4's exact mask (ra in [0,360], dec in [-90,90],
                          no NaN, no MAST sentinel -1.0)
  - Object selection:     master_catalog_v21.csv has_civ=True rows only
                          (all of Richards' Cells 1–8 filtering is already applied)
  - Spatial cross-match:  DEDUP_SEP = 2.0 arcsec (same as Cell 4)
  - Product filtering:    productType='SCIENCE', calib_level in [2,3] (excludes
                          HASP/HLSP calib_level=4), per-instrument subgroup filter

Run under the hasp-env conda environment (astroquery required).

Usage:
    python -m pipeline.download_spectra [--data-dir /path/to/raw]
"""

import argparse
import os

import numpy as np
import pandas as pd
from pathlib import Path
from astropy.coordinates import SkyCoord
from astropy.table import Table as AstropyTable
import astropy.units as u
from astroquery.mast import Observations

# ── Constants (from Richards' v21 Cell 1) ─────────────────────────────────────
DEDUP_SEP = 2.0  # arcsec — sky deduplication radius (same as Cell 4 DEDUP_SEP)

# Product filters per instrument (calib_level in [2,3] applied to all)
# FOS uses non-standard productSubGroupDescription values; filter by type+level only.
PRODUCT_SUBGROUPS = {
    'COS':  ['X1D', 'X1DSUM'],
    'STIS': ['X1D', 'SX1'],
    'FOS':  None,
}
CALIB_LEVELS = [2, 3]  # exclude HASP/HLSP calib_level=4

# ── Paths (relative to repo root) ─────────────────────────────────────────────
_REPO_ROOT    = Path(__file__).resolve().parent.parent
_PIPELINE_OUT = _REPO_ROOT / 'pipeline_output'

CATALOG_PATH         = _PIPELINE_OUT / 'master_catalog_v21.csv'
CACHE_COS            = _PIPELINE_OUT / 'v21_cache_mast_cos.csv'
CACHE_STIS           = _PIPELINE_OUT / 'v21_cache_mast_stis.csv'
CACHE_FOS            = _PIPELINE_OUT / 'v21_cache_mast_fos.csv'
PROD_CACHE           = _PIPELINE_OUT / 'v21_cache_download_products.csv'
HSLA_RECOVERY_CACHE  = _PIPELINE_OUT / 'v21_cache_hsla_recovery.csv'
MANIFEST_PATH        = _PIPELINE_OUT / 'download_manifest.csv'

# Cone-search radius for HSLA-only recovery (Section 3.5).  Matches Cell 9b's
# 30-arcsec radius in Richards' v21 pipeline; the cluster-mean RA/Dec in the
# catalog is sub-arcsec, so 5" would be enough, but 30" stays safe against
# any HSLA-vs-individual-exposure pointing offsets.
HSLA_RECOVERY_RADIUS_ARCSEC = 30.0


# ── Section 2 — Load catalog ───────────────────────────────────────────────────

def load_has_civ_catalog():
    """
    Load master_catalog_v21.csv and return only has_civ=True rows.
    Handles pandas reading 'True'/'False' as strings or booleans.
    """
    # The master catalog CSV has space-padded fields.  skipinitialspace=True
    # strips the leading space after each comma so that quoted multi-value fields
    # like proposal_ids="16679,2578" are recognised by the C parser (a leading
    # space before the opening quote breaks quoting detection).  Column names and
    # string values are stripped afterward to remove remaining trailing spaces.
    catalog = pd.read_csv(CATALOG_PATH, skipinitialspace=True)
    catalog.columns = catalog.columns.str.strip()
    for col in catalog.select_dtypes(include='object').columns:
        catalog[col] = catalog[col].str.strip()
    mask = catalog['has_civ'].astype(str).str.lower() == 'true'
    return catalog[mask].reset_index(drop=True)


# ── Section 3 — Observation cache loading and cross-match ─────────────────────

def _norm_inst_family(instrument_str):
    """Normalize MAST instrument_name (e.g. 'FOS/BL') to family tag COS/STIS/FOS."""
    s = str(instrument_str).upper()
    if s.startswith('COS'):
        return 'COS'
    if s.startswith('STIS'):
        return 'STIS'
    if s.startswith('FOS'):
        return 'FOS'
    return s


def load_obs_cache():
    """
    Load and stack Richards' three MAST observation cache CSVs.

    Applies the EXACT coordinate validity filter from Cell 4 of the v21 script:

        valid = (ra.notna() & dec.notna() &
                 (ra >= 0) & (ra <= 360) &
                 (dec >= -90) & (dec <= 90))

    This drops calibration-lamp observations (MAST target_name='WAVE', no sky
    coordinates) and any rows with the MAST sentinel value -1.0.

    calib_level is NOT filtered here — Cell 4 does not filter it either.
    Calib_level filtering is applied later at the product level (Section 4).
    """
    dfs = []
    for path, family in [(CACHE_COS, 'COS'), (CACHE_STIS, 'STIS'), (CACHE_FOS, 'FOS')]:
        df = pd.read_csv(path)
        # Tag instrument family; override from 'instrument' column if present
        df['inst_family'] = family
        if 'instrument' in df.columns:
            df['inst_family'] = df['instrument'].apply(_norm_inst_family)
        dfs.append(df)

    stack = pd.concat(dfs, ignore_index=True)

    # Coordinate validation — exact Cell 4 mask
    stack['ra_deg']  = pd.to_numeric(stack['ra_deg'],  errors='coerce')
    stack['dec_deg'] = pd.to_numeric(stack['dec_deg'], errors='coerce')

    valid = (
        stack['ra_deg'].notna()   & stack['dec_deg'].notna()  &
        (stack['ra_deg']  >= 0)   & (stack['ra_deg']  <= 360) &
        (stack['dec_deg'] >= -90) & (stack['dec_deg'] <= 90)
    )
    n_dropped = (~valid).sum()
    if n_dropped:
        print(f'      Dropped {n_dropped} rows with invalid/sentinel coordinates '
              f'(mirrors Cell 4 filter)')

    return stack[valid].reset_index(drop=True)


def cross_match_catalog_to_cache(catalog, obs_df):
    """
    Recover obs_ids for each has_civ object by spatial cross-match.

    catalog RA/Dec values are the mean cluster positions computed by Richards'
    Cell 4 greedy friends-of-friends algorithm, so a radius search with the
    same DEDUP_SEP recovers exactly the observations grouped into each cluster.

    Returns obs_df rows with a 'common_name' column added (one row per
    obs_id matched to a catalog object; the same obs_id can appear at most
    once per catalog object).
    """
    cat_coords = SkyCoord(
        ra=catalog['ra_deg'].values  * u.deg,
        dec=catalog['dec_deg'].values * u.deg,
    )
    obs_coords = SkyCoord(
        ra=obs_df['ra_deg'].values  * u.deg,
        dec=obs_df['dec_deg'].values * u.deg,
    )
    sep_limit = DEDUP_SEP * u.arcsec

    matched_rows = []
    for name, cat_c in zip(catalog['common_name'], cat_coords):
        seps = cat_c.separation(obs_coords)
        hits = obs_df[seps <= sep_limit].copy()
        hits['common_name'] = name
        matched_rows.append(hits)

    if not matched_rows:
        return pd.DataFrame()
    return pd.concat(matched_rows, ignore_index=True)


# ── Section 3.5 — HSLA-only recovery ──────────────────────────────────────────

def recover_hsla_only_exposures(catalog, obs_df_matched):
    """
    For catalog objects whose only matched MAST observations are HSLA HLSP
    co-adds (obs_id starts with 'hst_hsla_'), re-query MAST by sky position
    to recover the underlying individual COS/STIS exposures.

    Why this is needed: Richards' Cell 2 / Cell 3 MAST queries filter by
    target_classification AGN keywords.  HSLA HLSP entries have rich
    classifications that match the filter, but the *underlying* individual
    exposures of the same target often have sparser classifications
    (e.g. 'Galaxy' or empty) and therefore never enter the cache.  The
    download step then fails on the HSLA-only obs_ids because
    Observations.get_product_list('hst_hsla_*') is unreliable.

    This function mirrors Cell 9b of the v21 pipeline (cone search by sky
    position, no classification filter) but is triggered for HSLA-only
    objects rather than legacy-CSV-only objects.

    Returns a DataFrame of additional obs rows in the same schema as
    obs_df_matched (with 'common_name' and 'inst_family' columns).
    HLSP/HSLA rows are filtered out so we do not re-introduce the
    same problematic obs_ids.
    """
    if obs_df_matched.empty:
        return pd.DataFrame()

    is_hsla = obs_df_matched['obs_id'].astype(str).str.startswith('hst_hsla_')
    by_name = obs_df_matched.assign(_is_hsla=is_hsla).groupby('common_name')['_is_hsla'].all()
    hsla_only_names = sorted(by_name[by_name].index.tolist())

    if not hsla_only_names:
        print('      No HSLA-only objects detected — recovery skipped.')
        return pd.DataFrame()

    print(f'      {len(hsla_only_names)} HSLA-only objects detected:')
    for nm in hsla_only_names:
        print(f'        {nm}')

    if HSLA_RECOVERY_CACHE.exists():
        print(f'      Loading HSLA recovery cache from {HSLA_RECOVERY_CACHE.name}')
        return pd.read_csv(HSLA_RECOVERY_CACHE)

    print(f'      Cone-searching MAST (radius={HSLA_RECOVERY_RADIUS_ARCSEC}") '
          f'for individual exposures ...')
    radius_deg = HSLA_RECOVERY_RADIUS_ARCSEC / 3600.0
    new_rows = []
    for nm in hsla_only_names:
        cat_row = catalog[catalog['common_name'] == nm].iloc[0]
        ra, dec = float(cat_row['ra_deg']), float(cat_row['dec_deg'])
        # Search both COS and STIS — a few HSLA targets have STIS exposures too.
        # FOS is decommissioned and not relevant to HSLA-only objects.
        found_total = 0
        for inst_prefix in ('COS', 'STIS'):
            try:
                obs = Observations.query_criteria(
                    obs_collection='HST',
                    instrument_name=f'{inst_prefix}*',
                    dataproduct_type='spectrum',
                    s_ra=[ra - radius_deg, ra + radius_deg],
                    s_dec=[dec - radius_deg, dec + radius_deg],
                )
                if len(obs) == 0:
                    continue
                df = obs.to_pandas()
            except Exception as exc:
                print(f'        WARNING: cone search failed for {nm} / {inst_prefix}: {exc}')
                continue

            # Drop HLSP / HSLA rows; we want individual exposures only.
            obs_id_str = df['obs_id'].astype(str)
            df = df[~obs_id_str.str.startswith('hst_hsla_')]
            df = df[~obs_id_str.str.startswith('hst_hlsp_')]
            if len(df) == 0:
                continue

            # Normalise to the same column schema as obs_df_matched.
            ren = {'s_ra': 'ra_deg', 's_dec': 'dec_deg',
                   'target_name': 'mast_target_name',
                   'instrument_name': 'instrument'}
            df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
            df['common_name'] = nm
            df['inst_family'] = inst_prefix
            new_rows.append(df)
            found_total += len(df)

        if found_total == 0:
            print(f'        WARNING: no individual exposures recovered for {nm}')
        else:
            print(f'        {nm}: {found_total} individual exposures recovered')

    if not new_rows:
        print('      No HSLA-only recovery rows found.')
        recovered = pd.DataFrame()
    else:
        recovered = pd.concat(new_rows, ignore_index=True)
        recovered.to_csv(HSLA_RECOVERY_CACHE, index=False)
        print(f'      Saved HSLA recovery cache ({len(recovered)} rows) → '
              f'{HSLA_RECOVERY_CACHE.name}')
    return recovered


# ── Section 4 — Product lists (with caching) ──────────────────────────────────

def _filter_products(all_prods):
    """
    Apply per-instrument subgroup + SCIENCE + calib_level filtering to a
    DataFrame of MAST products.  Factored out so both the bulk and the
    incremental code paths use identical filtering logic.
    """
    if all_prods.empty:
        return all_prods

    all_prods = all_prods[
        (all_prods['productType'] == 'SCIENCE') &
        (all_prods['calib_level'].isin(CALIB_LEVELS))
    ]
    keep = []
    for inst, subgroups in PRODUCT_SUBGROUPS.items():
        inst_mask = all_prods['inst_family'] == inst
        if subgroups is None:
            keep.append(all_prods[inst_mask])
        else:
            sg_mask = all_prods['productSubGroupDescription'].isin(subgroups)
            keep.append(all_prods[inst_mask & sg_mask])
    if not keep:
        return pd.DataFrame()
    return pd.concat(keep, ignore_index=True)


def _query_products_for_obs(unique_obs):
    """
    Loop over (obs_id, common_name, inst_family) rows and call
    Observations.get_product_list for each.  Returns a list of per-obs
    pandas DataFrames (each tagged with common_name + inst_family).
    Failed obs_ids are logged and skipped.
    """
    prod_rows = []
    n = len(unique_obs)
    for i, row in unique_obs.iterrows():
        if i % 50 == 0:
            print(f'      ... {i}/{n}')
        try:
            products = Observations.get_product_list(row['obs_id'])
        except Exception:
            # Legacy FOS obs_ids trigger "Error converting data type varchar to
            # bigint" on MAST when passed as strings.  Fetching the observation
            # table first gives back a row with the numeric obsid; passing that
            # table to get_product_list uses the integer ID and avoids the error.
            try:
                obs_table = Observations.query_criteria(obs_id=row['obs_id'])
                if len(obs_table) == 0:
                    print(f'      Warning: no MAST record for obs_id {row["obs_id"]}')
                    continue
                products = Observations.get_product_list(obs_table)
            except Exception as exc2:
                print(f'      Warning: get_product_list failed for {row["obs_id"]}: {exc2}')
                continue
        pdf = products.to_pandas()
        pdf['common_name'] = row['common_name']
        pdf['inst_family']  = row['inst_family']
        prod_rows.append(pdf)
    return prod_rows


def get_filtered_products(obs_df_matched):
    """
    For each unique obs_id, retrieve the MAST product list (cached to disk),
    then filter to SCIENCE products at calib_level in [2, 3] with the
    per-instrument productSubGroupDescription filter.

    Products at calib_level=4 (HASP/HLSP co-adds) are excluded.
    FOS is filtered by productType + calib_level only (non-standard subgroup names).

    Cache is incremental: if PROD_CACHE exists but obs_df_matched contains
    obs_ids that are not represented in it (e.g. new rows from the HSLA
    recovery step in Section 3.5), only those new obs_ids are queried and
    their filtered products appended to the existing cache.
    """
    unique_obs = (
        obs_df_matched[['obs_id', 'common_name', 'inst_family']]
        .drop_duplicates('obs_id')
        .reset_index(drop=True)
    )

    cached = pd.DataFrame()
    if PROD_CACHE.exists():
        cached = pd.read_csv(PROD_CACHE)
        cached_ids = set(cached['obs_id'].astype(str)) \
                     if 'obs_id' in cached.columns else set()
        already = unique_obs['obs_id'].astype(str).isin(cached_ids)
        new_obs = unique_obs[~already].reset_index(drop=True)
        print(f'      Product cache: {len(cached)} cached products covering '
              f'{len(cached_ids)} obs_ids')
        if len(new_obs) == 0:
            print('      No new obs_ids to query — using cache as-is.')
            return cached
        print(f'      {len(new_obs)} new obs_ids not in cache — querying ...')
    else:
        print('      No product cache found — querying MAST ...')
        new_obs = unique_obs

    prod_rows = _query_products_for_obs(new_obs)
    if not prod_rows:
        print('      No new products returned from MAST.')
        return cached

    new_prods = pd.concat(prod_rows, ignore_index=True)
    new_filtered = _filter_products(new_prods)
    print(f'      {len(new_filtered)} new products survived filtering '
          f'(from {len(new_prods)} raw products)')

    if cached.empty:
        filtered = new_filtered
    else:
        # Align columns before concat (MAST may add/drop columns over time)
        all_cols = sorted(set(cached.columns) | set(new_filtered.columns))
        filtered = pd.concat(
            [cached.reindex(columns=all_cols), new_filtered.reindex(columns=all_cols)],
            ignore_index=True,
        )

    filtered.to_csv(PROD_CACHE, index=False)
    print(f'      Saved product cache ({len(filtered)} products) → {PROD_CACHE.name}')
    return filtered


# ── Section 5 — Download ───────────────────────────────────────────────────────

def download_all(filtered_products, raw_data_dir):
    """
    Download science products to RAW_DATA_DIR/{common_name}/{inst_family}/.
    Uses flat=True so astroquery does not create nested obs_id subdirectories.
    Returns a manifest DataFrame.
    """
    raw_data_dir = Path(raw_data_dir)
    manifest_rows = []

    names = sorted(filtered_products['common_name'].unique())
    print(f'      {len(names)} objects to download')

    for name in names:
        obj_prods = filtered_products[filtered_products['common_name'] == name]
        obj_files = 0
        obj_mb    = 0.0

        for inst in ['COS', 'STIS', 'FOS']:
            inst_prods = obj_prods[obj_prods['inst_family'] == inst]
            if inst_prods.empty:
                continue

            dest = raw_data_dir / name / inst
            dest.mkdir(parents=True, exist_ok=True)

            # download_products expects an astropy Table
            prod_table = AstropyTable.from_pandas(inst_prods)

            try:
                result = Observations.download_products(
                    prod_table,
                    download_dir=str(dest),
                    flat=True,
                )
                result_df = (result.to_pandas()
                             if hasattr(result, 'to_pandas')
                             else pd.DataFrame(result))

                for _, r in result_df.iterrows():
                    local_path = Path(str(r.get('Local Path', '')))
                    size_mb = local_path.stat().st_size / 1e6 if local_path.exists() else 0.0
                    status  = 'ok' if local_path.exists() else 'error'
                    manifest_rows.append({
                        'common_name':  name,
                        'inst_family':  inst,
                        'obs_id':       r.get('obs_id', ''),
                        'filename':     local_path.name,
                        'file_size_mb': round(size_mb, 3),
                        'status':       status,
                    })
                    obj_files += 1
                    obj_mb    += size_mb

            except Exception as exc:
                print(f'      Warning: download failed for {name}/{inst}: {exc}')
                manifest_rows.append({
                    'common_name':  name,
                    'inst_family':  inst,
                    'obs_id':       '',
                    'filename':     '',
                    'file_size_mb': 0.0,
                    'status':       f'error: {exc}',
                })

        print(f'      {name}: {obj_files} files, {obj_mb:.1f} MB')

    return pd.DataFrame(manifest_rows)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            'Download calibrated HST spectra from MAST for all has_civ=True '
            'objects in master_catalog_v21.csv. Run under hasp-env.'
        )
    )
    parser.add_argument(
        '--data-dir',
        default=os.environ.get('HST_PAPER_DATA_DIR', str(_REPO_ROOT / 'raw_data')),
        help=(
            'Root directory for downloaded FITS files. '
            'Default: $HST_PAPER_DATA_DIR env var, or ./raw_data'
        ),
    )
    args = parser.parse_args()

    # ── Section 2 — Load catalog ───────────────────────────────────────────────
    print(f'[1/5] Catalog: {CATALOG_PATH.name}')
    catalog = load_has_civ_catalog()
    print(f'      {len(catalog)} has_civ=True objects')

    # ── Section 3 — Recover obs_ids from Richards' cache files ────────────────
    print('\n[2/5] Loading MAST observation cache files ...')
    obs_df = load_obs_cache()
    print(f'      {len(obs_df)} observations with valid coordinates')

    print(f'      Cross-matching to catalog (DEDUP_SEP={DEDUP_SEP}") ...')
    obs_df_matched = cross_match_catalog_to_cache(catalog, obs_df)

    # Verification A — object-level completeness (pre-download)
    # Every has_civ object's RA/Dec is the mean cluster position computed by
    # Richards' Cell 4 FoF directly from these cache rows, so matching back
    # with the same DEDUP_SEP must always recover at least one obs_id per
    # object.  Any miss indicates a mismatch between catalog and cache versions
    # and must be investigated before downloading.
    catalog_names = set(catalog['common_name'])
    matched_names = set(obs_df_matched['common_name']) if not obs_df_matched.empty else set()
    missing       = catalog_names - matched_names
    extra         = matched_names - catalog_names

    print(f'\n      Verification A (pre-download):')
    print(f'        {len(matched_names):4d} matched')
    print(f'        {len(missing):4d} missing  (should always be 0)')
    print(f'        {len(extra):4d} extra    (should always be 0)')
    if extra:
        raise SystemExit(
            f'Cross-match produced {len(extra)} extra objects not in the catalog. '
            f'This should never happen. Investigate: {sorted(extra)}'
        )
    if missing:
        print(f'      Missing objects:')
        for nm in sorted(missing):
            print(f'        {nm}')
        raise SystemExit(
            f'{len(missing)} has_civ objects have no matched obs_ids in the cache. '
            'The catalog and cache files may be from different pipeline runs. '
            'Investigate before downloading.'
        )

    # ── Section 3.5 — HSLA-only recovery ──────────────────────────────────────
    # See Migration_Log.md "HSLA-only objects missed by Phase 5 download" for
    # background.  Some has_civ objects match only HSLA HLSP co-add obs_ids
    # (`hst_hsla_*`) in the cache because Richards' classification-filtered
    # query did not return their underlying individual exposures.  We re-query
    # MAST by sky position to recover those individual exposures here.
    print('\n[2.5/5] HSLA-only recovery ...')
    recovered = recover_hsla_only_exposures(catalog, obs_df_matched)
    if not recovered.empty:
        # Drop the HSLA-only rows for these objects from obs_df_matched and
        # replace with the recovered individual-exposure rows.
        recovered_names = set(recovered['common_name'])
        is_old_hsla_only = (
            obs_df_matched['common_name'].isin(recovered_names)
            & obs_df_matched['obs_id'].astype(str).str.startswith('hst_hsla_')
        )
        n_dropped = int(is_old_hsla_only.sum())
        obs_df_matched = obs_df_matched[~is_old_hsla_only].reset_index(drop=True)
        # Align columns before concat
        all_cols = sorted(set(obs_df_matched.columns) | set(recovered.columns))
        obs_df_matched = pd.concat(
            [obs_df_matched.reindex(columns=all_cols),
             recovered.reindex(columns=all_cols)],
            ignore_index=True,
        )
        print(f'      Replaced {n_dropped} HSLA-only rows with '
              f'{len(recovered)} individual-exposure rows '
              f'(covering {len(recovered_names)} objects)')

        # Re-check: any HSLA-only objects still without recovered exposures?
        still_hsla_only = sorted(
            nm for nm in (matched_names & {nm for nm in catalog['common_name']
                                            if nm not in recovered_names})
            if obs_df_matched[obs_df_matched['common_name'] == nm]
                ['obs_id'].astype(str).str.startswith('hst_hsla_').all()
            and len(obs_df_matched[obs_df_matched['common_name'] == nm]) > 0
        )
        if still_hsla_only:
            print(f'      WARNING: {len(still_hsla_only)} object(s) remain HSLA-only '
                  f'after recovery: {still_hsla_only}')

    # ── Section 4 — Get product lists ─────────────────────────────────────────
    print('\n[3/5] Getting product lists from MAST ...')
    filtered_products = get_filtered_products(obs_df_matched)
    print(f'      {len(filtered_products)} science products to download')

    # ── Section 5 — Download ───────────────────────────────────────────────────
    print(f'\n[4/5] Downloading to: {args.data_dir}')
    manifest = download_all(filtered_products, args.data_dir)

    # ── Section 6 — Manifest and Verification B ───────────────────────────────
    print(f'\n[5/5] Writing manifest and post-download verification ...')
    manifest.to_csv(MANIFEST_PATH, index=False)
    print(f'      Manifest written → {MANIFEST_PATH.name}')

    downloaded_names = set(manifest[manifest['status'] == 'ok']['common_name'])
    missing_after    = catalog_names - downloaded_names

    print(f'\n      Verification B (post-download):')
    print(f'        {len(downloaded_names):4d} objects with at least one file downloaded OK')
    print(f'        {len(missing_after):4d} objects missing after download')
    if missing_after:
        print(f'      Missing after download:')
        for nm in sorted(missing_after):
            print(f'        {nm}')


if __name__ == '__main__':
    main()
