# Migration Log

**Project:** Rebinning + ICA Code Migration into HST Paper Repo
**Started:** 2026-04-29

All changes made during each phase are logged here with dates, file names,
and descriptions of what was done.

---

## Phase 0 — Project Scaffolding

**Date:** 2026-04-29

### Files Created

| File | Description |
|------|-------------|
| `CLAUDE.md` | Codebase guide: repo structure, principles, data flow, migration status table |
| `pyproject.toml` | Installable package `hst_civ_paper`; Python >=3.12; all dependencies listed |
| `rebinning/__init__.py` | Empty package init (to be populated in Phase 1) |
| `ica/__init__.py` | Empty package init (to be populated in Phase 2) |
| `Data/components/.gitkeep` | Placeholder so git tracks the empty components directory |

### Files Modified

| File | Change |
|------|--------|
| `.gitignore` | Added `RebinnedSpec_*/`, `ICA_Plots_*/`, `ICA_full/` to ignore large pipeline output directories |

### Notes
- `Migration_Log.md` header was already present (created alongside `Migration_Plan.md`)
- `.gitignore` already had `*.egg-info/`, `dist/`, `build/` from previous scaffold; no duplication needed

---

## Phase 1 — Rebinning Pipeline Migration

*(not yet started)*

---

## Phase 2 — ICA Analysis Pipeline Migration

*(not yet started)*

---

## Phase 3 — Integration Testing

*(not yet started)*

---

## Phase 4 — CLAUDE.md and Documentation

*(not yet started)*

---

## Phase 5 — Connect to Master Pipeline

*(not yet started)*
