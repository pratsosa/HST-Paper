# Manual-fix GUI for interactive ICA CIV fits (Phase 2 tooling).
# New file (2026-07-10). Modifies NO existing module: GuiFixProcessor subclasses
# ica.manual_fix.ICAManualFixProcessor and adds a non-plotting fit method that
# returns arrays for the embedded matplotlib canvas.
#
# What it does
# ------------
#   * Loads a rebinned spectrum, runs the *exact* ICA pipeline used by the batch
#     runner (setup_object -> main_ICA -> get_CIV_parameters), and shows the
#     usual diagnostic panels in a Qt window with the real matplotlib toolbar.
#   * Sidebar lets you add wavelength ranges to mask (right-drag on a plot) and
#     force a component set (auto/mod/low/high), then Re-fit (runs off the UI
#     thread so the ~5-10 s fit never freezes the window).
#   * Save writes the per-object override to a JSON file (mask ranges + forced
#     components). This is the artifact you edit; it is designed to be merged
#     back into the batch pipeline (see feed_into_batch() docstring below).
#
# Run
# ---
#   python -m ica.manual_fix_gui --rebin-dir "F:\Richards-data\HST Paper\RebinnedSpec_master"
#
# The GUI operates in master mode (object list = FITS stems in --rebin-dir),
# matching `python -m ica.run_all_objects --master`.

from __future__ import annotations

import argparse
import json
import os
import sys
import traceback
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import SpanSelector
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)

from PySide6 import QtCore, QtWidgets

from ica.manual_fix import ICAManualFixProcessor
from ica import run_ica
from ica import civ_bal_regions as CIV_BAL_regions

# Default location of the override store (kept next to the config it complements).
DEFAULT_OVERRIDES_JSON = Path(__file__).resolve().parent / "manual_fix_overrides.json"

# Right mouse button for the SpanSelector so it never collides with the toolbar's
# left-drag pan/zoom. Left-drag = zoom (toolbar); right-drag = add mask range.
_SPAN_BUTTON = 3

# Status banner styles (busy / success / error).
_STATUS_BUSY = "QLabel{background:#fff3cd; color:#7a5b00; padding:6px; border:1px solid #e0c060;}"
_STATUS_OK = "QLabel{background:#d4edda; color:#155724; padding:6px; border:1px solid #9ad0a5;}"
_STATUS_ERR = "QLabel{background:#f8d7da; color:#721c24; padding:6px; border:1px solid #d69098;}"


# ---------------------------------------------------------------------------
# Fit engine (headless; safe to call from a worker thread)
# ---------------------------------------------------------------------------
class GuiFixProcessor(ICAManualFixProcessor):
    """ICAManualFixProcessor that fits from *explicit* overrides passed by the
    GUI instead of the global MANUAL_FIX_CONFIG dict, and returns arrays for
    plotting rather than writing a PNG.

    fit_for_gui() reproduces process_object()'s numerical path exactly
    (setup_object -> main_ICA -> get_CIV_parameters); only the source of the
    mask/comps overrides differs, so a GUI fit == the batch fit for the same
    override. No matplotlib is touched here, so it is safe off the UI thread.
    """

    def fit_for_gui(self, name, mask_ranges=None, mask_pixels=None, comps_use=None):
        mask_ranges = mask_ranges or []
        mask_pixels = mask_pixels or []

        # Reload from disk every fit -> guaranteed-clean mask state, no cross-fit
        # mutation bugs. Cheap next to the ~5-10 s fit itself.
        wave, flux, z, errs, mask, spec_name = self.setup_object(name)
        mask = np.asarray(mask, dtype=float).copy()

        for lo, hi in mask_ranges:
            mask[(wave >= lo) & (wave <= hi)] = 1
        for wl in mask_pixels:
            mask[np.argmin(np.abs(wave - wl))] = 1

        cu = None if comps_use in (None, "", "auto") else comps_use

        (wave_arb, flux_arb, errs_arb, mask_arb,
         wave_ica, flux_ica, f2500) = run_ica.main_ICA(
            wave, flux, errs, mask, z, name="", ica_path=None, comps_use=cu)

        civ_blue, civ_ew = self.get_CIV_parameters(
            wave_arb, flux_arb, wave_ica, flux_ica, name)

        return dict(
            name=name, spec_name=spec_name, z=float(z),
            wave_arb=wave_arb, flux_arb=flux_arb, errs_arb=errs_arb, mask_arb=mask_arb,
            wave_ica=wave_ica, flux_ica=flux_ica,
            f2500=float(f2500), civ_blue=float(civ_blue), civ_ew=float(civ_ew),
        )


# ---------------------------------------------------------------------------
# Worker (runs one fit on a background thread)
# ---------------------------------------------------------------------------
class _FitSignals(QtCore.QObject):
    done = QtCore.Signal(dict)
    failed = QtCore.Signal(str)


class _FitTask(QtCore.QRunnable):
    def __init__(self, proc, name, mask_ranges, mask_pixels, comps_use):
        super().__init__()
        self.proc = proc
        self.name = name
        self.mask_ranges = mask_ranges
        self.mask_pixels = mask_pixels
        self.comps_use = comps_use
        self.signals = _FitSignals()

    @QtCore.Slot()
    def run(self):
        try:
            res = self.proc.fit_for_gui(
                self.name, self.mask_ranges, self.mask_pixels, self.comps_use)
            self.signals.done.emit(res)
        except Exception:
            self.signals.failed.emit(traceback.format_exc())


# ---------------------------------------------------------------------------
# Override persistence
# ---------------------------------------------------------------------------
def load_overrides(path):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save_overrides(path, overrides):
    with open(path, "w") as f:
        json.dump(overrides, f, indent=2)


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------
class ManualFixWindow(QtWidgets.QMainWindow):
    def __init__(self, proc, names, overrides_path):
        super().__init__()
        self.proc = proc
        self.names = names
        self.overrides_path = overrides_path
        self.overrides = load_overrides(overrides_path)
        self.pool = QtCore.QThreadPool.globalInstance()

        # per-object working state
        self.mask_ranges = []      # list of [lo, hi]
        self.comps_use = "auto"
        self._saved_lims = None    # (xlim, ylim) per panel, to preserve zoom on refit
        self._busy = False

        self.setWindowTitle("ICA Manual Fix")
        self.resize(1500, 950)
        self._build_ui()

        if self.names:
            self.obj_combo.setCurrentIndex(0)
            self._load_object(self.names[0])

    # ---- UI construction -------------------------------------------------
    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        root = QtWidgets.QHBoxLayout(central)

        # --- left: figure + toolbar ---
        left = QtWidgets.QVBoxLayout()
        self.fig = Figure(figsize=(11, 9), constrained_layout=True)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        left.addWidget(self.toolbar)
        left.addWidget(self.canvas, stretch=1)
        root.addLayout(left, stretch=4)

        gs = GridSpec(7, 12, figure=self.fig)
        self.ax_full = self.fig.add_subplot(gs[:3, :8])
        self.ax_civ = self.fig.add_subplot(gs[3:6, :8])
        self.ax_err = self.fig.add_subplot(gs[6, :8], sharex=self.ax_civ)
        self.ax_fit = self.fig.add_subplot(gs[3:, 8:])
        self._panels = [self.ax_full, self.ax_civ, self.ax_err, self.ax_fit]

        # right-drag on the two spectrum panels adds a mask range
        # useblit=False: blitting caches an axes-background bitmap that goes
        # stale after a refit clears+replots, leaving the canvas showing the old
        # fit until a forced draw. Negligible perf cost for this occasional drag.
        self._spans = [
            SpanSelector(ax, self._on_span, "horizontal", useblit=False,
                         button=_SPAN_BUTTON,
                         props=dict(alpha=0.25, facecolor="orange"))
            for ax in (self.ax_full, self.ax_civ)
        ]

        # --- right: sidebar ---
        side = QtWidgets.QVBoxLayout()
        root.addLayout(side, stretch=1)

        # object picker
        side.addWidget(QtWidgets.QLabel("<b>Object</b>"))
        row = QtWidgets.QHBoxLayout()
        self.prev_btn = QtWidgets.QPushButton("◀")
        self.next_btn = QtWidgets.QPushButton("▶")
        self.prev_btn.setFixedWidth(36)
        self.next_btn.setFixedWidth(36)
        self.prev_btn.clicked.connect(lambda: self._step(-1))
        self.next_btn.clicked.connect(lambda: self._step(+1))
        self.obj_combo = QtWidgets.QComboBox()
        self.obj_combo.setEditable(True)
        self.obj_combo.addItems(self.names)
        self.obj_combo.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.obj_combo.completer().setCompletionMode(
            QtWidgets.QCompleter.PopupCompletion)
        self.obj_combo.activated.connect(
            lambda _i: self._load_object(self.obj_combo.currentText()))
        row.addWidget(self.prev_btn)
        row.addWidget(self.obj_combo, stretch=1)
        row.addWidget(self.next_btn)
        side.addLayout(row)

        # mask ranges
        side.addWidget(QtWidgets.QLabel("<b>Mask ranges</b> (right-drag a plot)"))
        self.mask_list = QtWidgets.QListWidget()
        self.mask_list.setMaximumHeight(160)
        side.addWidget(self.mask_list)
        mrow = QtWidgets.QHBoxLayout()
        rm = QtWidgets.QPushButton("Remove selected")
        clr = QtWidgets.QPushButton("Clear all")
        rm.clicked.connect(self._remove_selected_mask)
        clr.clicked.connect(self._clear_masks)
        mrow.addWidget(rm)
        mrow.addWidget(clr)
        side.addLayout(mrow)

        # components
        side.addWidget(QtWidgets.QLabel("<b>Components</b>"))
        self.comp_group = QtWidgets.QButtonGroup(self)
        crow = QtWidgets.QHBoxLayout()
        for label in ("auto", "mod", "low", "high"):
            rb = QtWidgets.QRadioButton(label)
            rb.toggled.connect(self._on_comp_changed)
            self.comp_group.addButton(rb)
            crow.addWidget(rb)
            if label == "auto":
                rb.setChecked(True)
        side.addLayout(crow)

        # actions
        self.refit_btn = QtWidgets.QPushButton("Re-fit")
        self.refit_btn.setStyleSheet("font-weight:bold; padding:8px;")
        self.refit_btn.clicked.connect(self._refit)
        side.addWidget(self.refit_btn)

        self.save_btn = QtWidgets.QPushButton("Save override")
        self.save_btn.clicked.connect(self._save)
        side.addWidget(self.save_btn)

        self.status = QtWidgets.QLabel("Ready.")
        self.status.setWordWrap(True)
        side.addWidget(self.status)

        self.results = QtWidgets.QLabel("")
        self.results.setWordWrap(True)
        self.results.setStyleSheet("font-family: monospace;")
        side.addWidget(self.results)

        side.addStretch(1)

    # ---- object / state handling ----------------------------------------
    def _step(self, delta):
        if self._busy or not self.names:
            return
        i = (self.obj_combo.currentIndex() + delta) % len(self.names)
        self.obj_combo.setCurrentIndex(i)
        self._load_object(self.names[i])

    def _load_object(self, name):
        if name not in self.names:
            return
        self.current = name
        ov = self.overrides.get(name, {})
        self.mask_ranges = [list(r) for r in ov.get("mask_ranges", [])]
        self.comps_use = ov.get("forced_components") or "auto"
        self._sync_mask_list()
        for rb in self.comp_group.buttons():
            if rb.text() == self.comps_use:
                rb.setChecked(True)
        self._saved_lims = None  # new object -> autoscale
        self.results.setText("")
        self._refit()

    def _sync_mask_list(self):
        self.mask_list.clear()
        for lo, hi in self.mask_ranges:
            self.mask_list.addItem("%.2f – %.2f Å" % (lo, hi))

    def _on_span(self, xmin, xmax):
        if xmax - xmin <= 0:
            return
        self.mask_ranges.append([round(float(xmin), 3), round(float(xmax), 3)])
        self._sync_mask_list()
        self._draw_mask_overlays()
        self.canvas.draw_idle()
        self.status.setText("Mask range added. Re-fit to apply.")

    def _remove_selected_mask(self):
        for item in self.mask_list.selectedItems():
            idx = self.mask_list.row(item)
            del self.mask_ranges[idx]
        self._sync_mask_list()
        self._draw_mask_overlays()
        self.canvas.draw_idle()

    def _clear_masks(self):
        self.mask_ranges = []
        self._sync_mask_list()
        self._draw_mask_overlays()
        self.canvas.draw_idle()

    def _on_comp_changed(self, checked):
        if checked:
            btn = self.comp_group.checkedButton()
            if btn:
                self.comps_use = btn.text()

    # ---- fitting ---------------------------------------------------------
    def _refit(self):
        if self._busy:
            return
        self._busy = True
        self._set_busy(True)
        # preserve current zoom across the refit (None on first draw of object)
        self._saved_lims = [(ax.get_xlim(), ax.get_ylim()) for ax in self._panels] \
            if self._has_drawn() else None

        task = _FitTask(self.proc, self.current,
                        [list(r) for r in self.mask_ranges], [], self.comps_use)
        task.signals.done.connect(self._on_fit_done)
        task.signals.failed.connect(self._on_fit_failed)
        self.pool.start(task)

    def _has_drawn(self):
        return any(ax.lines for ax in self._panels)

    @QtCore.Slot(dict)
    def _on_fit_done(self, res):
        try:
            self._draw(res)
            self.results.setText(
                "EW      = %7.2f Å\n"
                "blueshift = %7.1f km/s\n"
                "f2500   = %.3e\n"
                "z       = %.4f" % (res["civ_ew"], res["civ_blue"],
                                    res["f2500"], res["z"]))
            self._set_status("✓  Fit complete.", _STATUS_OK)
        finally:
            self._busy = False
            self._set_busy(False)

    @QtCore.Slot(str)
    def _on_fit_failed(self, tb):
        self._busy = False
        self._set_busy(False)
        self._set_status("✗  Fit FAILED (see console).", _STATUS_ERR)
        print(tb, file=sys.stderr)

    def _set_controls_enabled(self, on):
        for w in (self.refit_btn, self.prev_btn, self.next_btn,
                  self.obj_combo, self.save_btn):
            w.setEnabled(on)

    def _set_status(self, text, style):
        self.status.setText(text)
        self.status.setStyleSheet(style)

    def _set_busy(self, busy):
        """Toggle the fitting-in-progress visuals: disable controls, relabel the
        Re-fit button, colour the status banner, and set a wait cursor."""
        self._set_controls_enabled(not busy)
        if busy:
            self.refit_btn.setText("Fitting…")
            self._set_status("⏳  Fitting %s …" % self.current, _STATUS_BUSY)
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        else:
            self.refit_btn.setText("Re-fit")
            QtWidgets.QApplication.restoreOverrideCursor()

    # ---- drawing ---------------------------------------------------------
    def _draw(self, res):
        wave, flux = res["wave_arb"], res["flux_arb"]
        errs, mask = res["errs_arb"], res["mask_arb"]
        wave_ica, flux_ica = res["wave_ica"], res["flux_ica"]
        spec_name = res["spec_name"]

        for ax in self._panels:
            ax.clear()

        c4 = wave > 1400
        ylow = max(0, np.nanpercentile(flux, 1))
        yhigh = np.nanpercentile(flux[c4], 99) + np.nanmedian(flux[c4])

        # full spectrum
        self.ax_full.plot(wave, flux, "-k", alpha=0.6)
        self.ax_full.plot(wave_ica, flux_ica, "-r")
        self.proc.plot_HST(wave, flux, mask, self.ax_full)
        self.ax_full.set_xlim(max(min(wave), 1250), min(max(wave), max(wave_ica)))
        self.ax_full.set_ylim(ylow, yhigh + 5)
        self.ax_full.set_ylabel("Flux (arb.)")
        self.ax_full.set_title("%s  —  ICA fit" % res["name"])

        # CIV region
        self.ax_civ.plot(wave, flux, "-k", alpha=0.6)
        self.ax_civ.plot(wave_ica, flux_ica, "-r")
        self.proc.plot_HST(wave, flux, mask, self.ax_civ)
        self.ax_civ.set_xlim(1500, 1600)
        self.ax_civ.set_ylim(ylow, yhigh + 5)
        self.ax_civ.set_ylabel("Flux (arb.)")

        # errors
        self.ax_err.plot(wave, errs, "-k", alpha=0.6)
        self.ax_err.set_xlim(1500, 1600)
        self.ax_err.set_ylim(max(0, np.nanpercentile(errs, 1)),
                             np.nanpercentile(errs, 99))
        self.ax_err.set_xlabel("Rest wavelength (Å)")
        self.ax_err.set_ylabel("Error")

        # CIV fit panel (reuse processor helper)
        self.proc.plot_CIV_analysis(self.ax_fit, wave, flux, wave_ica, flux_ica,
                                    res["name"], res["civ_blue"], res["civ_ew"])

        self._draw_mask_overlays()

        # restore zoom if we saved it before this refit
        if self._saved_lims is not None:
            for ax, (xl, yl) in zip(self._panels, self._saved_lims):
                ax.set_xlim(xl)
                ax.set_ylim(yl)

        # Force a full synchronous repaint. draw_idle() can be coalesced away
        # when the restored zoom leaves the axis limits unchanged, so the new
        # fit wouldn't show until a toolbar action forced a draw.
        self.canvas.draw()

    def _draw_mask_overlays(self):
        # remove previous overlay spans, then re-add current mask ranges
        for ax in (self.ax_full, self.ax_civ):
            for patch in list(getattr(ax, "_mask_patches", [])):
                patch.remove()
            ax._mask_patches = [
                ax.axvspan(lo, hi, color="orange", alpha=0.18, zorder=0)
                for lo, hi in self.mask_ranges
            ]

    # ---- save ------------------------------------------------------------
    def _save(self):
        entry = {}
        if self.mask_ranges:
            entry["mask_ranges"] = [list(r) for r in self.mask_ranges]
        if self.comps_use and self.comps_use != "auto":
            entry["forced_components"] = self.comps_use
        if entry:
            self.overrides[self.current] = entry
        else:
            self.overrides.pop(self.current, None)  # no-op fix -> don't store
        save_overrides(self.overrides_path, self.overrides)
        self.status.setText("Saved override for %s → %s"
                            % (self.current, os.path.basename(self.overrides_path)))


# ---------------------------------------------------------------------------
# feed_into_batch: how the saved JSON is meant to reach ica.run_all_objects
# ---------------------------------------------------------------------------
def feed_into_batch(overrides_path=DEFAULT_OVERRIDES_JSON):
    """Return {name: {'custom_mask_pixels': None|array, 'forced_components': ...}}
    in the shape ICAManualFixProcessor.apply_manual_fix expects, so the GUI's
    JSON can be merged over MANUAL_FIX_CONFIG.

    Note: the batch path masks by *pixel wavelength arrays* while the GUI stores
    *ranges*. Ranges are the cleaner representation; wiring them into the batch
    means teaching apply_manual_fix to also read 'mask_ranges'. That is a small,
    backward-compatible change tracked separately. This helper documents the
    contract and lets you inspect the mapping now.
    """
    ov = load_overrides(overrides_path)
    out = {}
    for name, entry in ov.items():
        out[name] = {
            "mask_ranges": entry.get("mask_ranges"),
            "custom_mask_pixels": None,
            "forced_components": entry.get("forced_components"),
        }
    return out


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------
def _make_processor(rebin_dir):
    return GuiFixProcessor(
        rebin_path=rebin_dir, master_mode=True,
        output_path="ICA_Plots_Rebin_master")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rebin-dir", default=os.environ.get("HST_PAPER_REBIN_DIR"),
        help="Directory of rebinned FITS files (master mode). "
             "Falls back to HST_PAPER_REBIN_DIR env var.")
    parser.add_argument(
        "--overrides", default=str(DEFAULT_OVERRIDES_JSON),
        help="Path to the override JSON store (read + written by Save).")
    parser.add_argument(
        "--selftest", metavar="OBJECT", default=None,
        help="Headless self-test: build the window offscreen, fit OBJECT, "
             "save a screenshot to manual_fix_gui_selftest.png, exit.")
    args = parser.parse_args(argv)

    if not args.rebin_dir:
        parser.error("--rebin-dir is required (or set HST_PAPER_REBIN_DIR).")

    if args.selftest:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    proc = _make_processor(args.rebin_dir)
    names = proc.list_master_objects()
    if not names:
        parser.error("No FITS files found in %s" % args.rebin_dir)

    win = ManualFixWindow(proc, names, args.overrides)

    if args.selftest:
        target = args.selftest
        if target not in names:
            # allow passing a bare stem present in the dir
            print("Self-test object %r not in list; using first (%s)."
                  % (target, names[0]))
            target = names[0]
        win.obj_combo.setCurrentText(target)
        win._load_object(target)
        # _load_object triggers an async fit; run the fit synchronously here so
        # the self-test is deterministic.
        res = proc.fit_for_gui(target, win.mask_ranges, [], win.comps_use)
        win._draw(res)
        out = "manual_fix_gui_selftest.png"
        win.fig.savefig(out, dpi=110)
        print("Self-test OK -> %s (EW=%.2f, blue=%.1f)"
              % (out, res["civ_ew"], res["civ_blue"]))
        return 0

    win.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
