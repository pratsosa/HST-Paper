"""
rebinning — HST spectral rebinning and coadding pipeline.

Migrated from Trevor Code/ in Phase 1 (2026-04-30).
All science logic is identical to the validated 2026AP run.

Key public API
--------------
coadd.rebin(Identifier, z, data_origin, fn_sdss, data_path, output_dir, sdss_spec_dir)
    Main entry point: rebin + coadd a single object.

run_rebin.main()
    CLI runner for the full sample (all instruments).

Instrument-specific spectral cuts
----------------------------------
spec_cuts_fos, spec_cuts_stis, spec_cuts_cos, spec_cuts_hsla
    BlueEdges / RedEdges dicts used by cut_edge_pix.

Internal modules (not normally called directly)
------------------------------------------------
read_spec_data   — FITS readers for FOS, STIS, COS, HSLA, SDSS-RM
cut_edge_pix     — Edge-pixel quality cutting (NoDQ / 2026AP version)
lower_res_rebin  — Rebin HST pixels to SDSS log-lambda scale
small_pix_filter — Median pixel filter for noise spike rejection
spec_morph       — Continuum fitting and spectral morphing (rebinning version)
bal_regions      — BAL trough masks used by spec_morph
"""

from rebinning import coadd
from rebinning import run_rebin
from rebinning import read_spec_data
from rebinning import spec_morph
from rebinning import small_pix_filter
from rebinning import lower_res_rebin
from rebinning import cut_edge_pix
from rebinning import bal_regions
from rebinning import spec_cuts_fos
from rebinning import spec_cuts_stis
from rebinning import spec_cuts_cos
from rebinning import spec_cuts_hsla
