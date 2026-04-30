"""
ica — ICA analysis pipeline for the HST CIV paper.

Migrated from Trevor Code/ICA Scripts/ in Phase 2 (2026-04-29).
See .claude/Migration_Log.md for the full change log.

Key entry points:
    run_ica.main_ICA       — run the full ICA pipeline on one spectrum
    manual_fix.ICAManualFixProcessor — batch processor with manual overrides
    run_all_objects.main   — run all objects from the catalog
"""

from ica import (
    run_ica,
    fit_composites,
    plot_ica,
    spec_morph,
    civ_bal_regions,
    morphing_edges,
    manual_fix,
    manual_fix_config,
    run_all_objects,
)
