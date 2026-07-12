# Manual-fix GUI (`ica/manual_fix_gui.py`)

Interactive tool for tweaking ICA CIV fits: mask noisy/absorbed pixels and force
a component set, re-fit, and save the result. Replaces the ad-hoc
`manual_object_fix` notebooks with a point-and-click workflow.

## Install

One prebuilt wheel, no compiler or system Qt needed:

```bash
pip install PySide6
```

(matplotlib, numpy, astropy, lmfit are already project deps.)

## Run

```bash
python -m ica.manual_fix_gui --rebin-dir "F:\Richards-data\HST Paper\RebinnedSpec_master"
```

The object list is every FITS stem in `--rebin-dir` (master mode, same as
`python -m ica.run_all_objects --master`). `--rebin-dir` also reads the
`HST_PAPER_REBIN_DIR` env var if omitted.

## Using it

- **Pick an object**: dropdown (type to search) or ◀ / ▶. Selecting an object
  auto-fits it (~5–10 s; the window stays responsive — fits run off-thread).
- **Zoom/pan**: the standard matplotlib toolbar (left-drag box zoom, pan, home).
- **Mask a region**: **right-drag** on the full-spectrum or CIV panel, or type in
  the text box (`1548.5` = nearest pixel; `1548.5 1549.2` = several; `1530-1545` =
  a range). Entries appear in the mask list (remove individually or clear all).
  Masks are shaded orange but only take effect on the next **Re-fit**.
- **Unmask (force-include) a pixel the pipeline masked** — e.g. a spurious NAL
  flag: tick **Unmask mode** and use the same text box / right-drag. Unmask entries
  show green (dashed line / light-green band). On Re-fit the pipeline runs normally,
  then the fit is re-run once more with those pixels cleared from the NAL/BAL mask.
  Unmasking is skipped entirely when the list is empty, so ordinary fits stay
  byte-identical to the batch pipeline.
- **Force components**: `auto` (χ²-selected, default) / `mod` / `low` / `high`.
- **Re-fit**: re-runs the pipeline with the current mask + components. Your zoom
  is preserved across the refit.
- **Save override**: writes the object's `{mask_ranges, forced_components}` to the
  override JSON (default `ica/manual_fix_overrides.json`, or `--overrides PATH`).
  An object with no mask and `auto` components is removed from the store.

## Design (for future edits)

Three separated layers so science, UI, and storage don't entangle:

1. **Engine** — `GuiFixProcessor.fit_for_gui()` subclasses the batch processor
   and reproduces its numerical path exactly (`setup_object → main_ICA →
   get_CIV_parameters`); only the override source differs, and it returns arrays
   instead of a PNG. No matplotlib, so it is safe on a worker thread. **No
   existing module is modified.**
2. **UI** — `ManualFixWindow` holds only `{current object, mask_ranges,
   comps_use}`; buttons mutate that, Re-fit runs the engine via `QThreadPool`.
3. **Storage** — plain JSON, one entry per object, diffable.

**Adding a feature (e.g. sigma clipping)** = add a field to the override dict,
teach `fit_for_gui` to apply it, add one sidebar widget. The science code and the
plotting code don't change.

## Feeding fixes back into the batch run (not yet wired)

`feed_into_batch()` documents the mapping from the GUI JSON to the shape
`ICAManualFixProcessor.apply_manual_fix` expects. The GUI stores masks as
**ranges**; the existing batch path masks by **pixel-wavelength arrays**. Making
`run_all_objects` honor the saved overrides is a small, backward-compatible change
to `apply_manual_fix` (read `mask_ranges` in addition to `custom_mask_pixels`, and
load+merge the JSON over `MANUAL_FIX_CONFIG`). Deliberately left out of the MVP so
no working file is touched until you want it.

## Headless self-test

```bash
python -m ica.manual_fix_gui --rebin-dir "<dir>" --selftest 0115+027_FOS
# -> manual_fix_gui_selftest.png
```

Builds the window offscreen, fits one object, saves a screenshot. Useful in CI /
over SSH where no display is available.
