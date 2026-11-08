# Pipeline changes, September 2026 — summary for downstream analysis

**AP, 2026-09-29.** Self-contained wrap-up of one round of pipeline fixes, written so
that comparison figures can be made without re-reading the investigation logs.

---

## 1. Use these files

| what | file |
|---|---|
| **Rebinned spectra** | `RebinnedSpec_master_v2/` (594 spectra) |
| **C IV measurements** | `ICA_Results_Rebin_master_v2_fixed.csv` (517 fits) |
| **Diagnostic plots** | `ICA_Plots_Rebin_master_v2_fixed/All/` (517 plots) |
| **Objects safe to use now** | `diagnostics/usable_objects_v2_fixed.txt` (**478 objects**) |
| **Objects to exclude** | `diagnostics/affected_unnormalised.csv` (36, see §5) |

**Do not use** `ICA_Results_Rebin_master.csv` or `ICA_Results_Rebin_master_v2.csv`
except as deliberate before/after comparisons. They are kept only for that:

| file | rebinning | C IV equivalent width | zero-error pixel guard |
|---|---|---|---|
| `ICA_Results_Rebin_master.csv` | original | biased | none |
| `ICA_Results_Rebin_master_v2.csv` | fixed | biased | none |
| `ICA_Results_Rebin_master_v2_fixed.csv` | fixed | **fixed** | **yes** |

The 478 usable objects break down as FOS 176, STIS 157, COS 145.

---

## 2. What changed, in one table

| # | stage | defect | effect of fixing it |
|---|---|---|---|
| 1 | rebinning | COS reader kept only detector row 0 | 24 spectra regained C IV; 91% of COS spectra changed |
| 2 | rebinning | continuum returned 0 inside an emission line | 15 spectra regained C IV |
| 3 | rebinning | FOS files with 2 readout groups crashed | 13 objects went from no output to output |
| 4 | rebinning | one bad exposure aborted the whole object | 5 objects likewise |
| 5 | ICA | C IV equivalent width used the continuum from the wrong wavelengths | every EW rises, median factor **1.49** |
| 6 | ICA | pixels with near-zero error pinned the fit at zero | 17 genuine collapsed fits addressed |
| 7 | rebinning | **regression from fix 2**: continuum fell back to 1.0 | 36 spectra still broken — see §5 |

Net on the spectra: 594 rebinned vs 576 before (18 new), 53 gained usable C IV,
**none lost it**.

---

## 3. The rebinning fixes (2026-09-10)

Full detail: `diagnostics/CIV_Data_Loss_Investigation.md`.

The trigger was objects with no data across C IV in the rebinned spectrum even though
MAST holds C IV coverage. A census of all 615 (object, instrument) pairs found four
causes, three of them silent.

**Fix 1 — read every detector row.** `read_cos_flat`/`read_stis_flat` in
`rebinning/read_spec_data.py` did `data['WAVELENGTH'][0]`, keeping only the first row of
the x1d table. A COS FUV exposure has two rows (FUVA, FUVB) and a NUV exposure has three
(NUVA/B/C); each is a physically separate piece of the spectrum. Row 0 alone discarded
half to two thirds of the coverage. Every row is now its own entry in the exposure
stack. Inherited from Trevor's `read_cos()`, which was never exercised on raw
two-segment files because his COS data arrived as pre-merged HSLA co-adds.

> Exception added the same day: **COS G140L segment FUVB** (cenwaves 1230, 1280) is not
> illuminated in first order. CalCOS writes the row anyway with a placeholder wavelength
> solution starting at ~0.02 Å, and flux *is* recorded against those nonsense
> wavelengths. Those rows are skipped, plus a general floor rejecting any row whose
> wavelengths start below 500 Å (the bluest genuine value in this dataset is 703 Å).

**Fix 2 — continuum of zero.** `cont_filtered()` in `rebinning/spec_morph.py` masks
emission-line windows, medians what is left, and interpolates back. The C IV window is
rest 1460–1685 Å and a COS G160M exposure covers ~190 Å, so at low redshift an entire
exposure can sit *inside* the window: nothing is left to median, and the function filled
the continuum with `0.0`. `coadd.rebin()` divides each exposure by that continuum, so
the exposure became ±inf and vanished silently. Now falls back to a flat continuum at
the exposure's own flux scale.

**Fix 3 — FOS two-group files.** `accum_flag()` tested
`mean(flux[-1]) > mean(flux[-2]) > mean(flux[-3])`; files with `NAXIS2 = 2` have no
third group and raised `IndexError`, aborting the object. Degrades to the single
comparison available.

**Fix 4 — one bad exposure no longer aborts the object.** `Cut_Edge_Pix` raises when an
exposure has no usable pixel; that exception escaped the reader and discarded every good
exposure alongside it. Such exposures are now skipped and reported.

### What this did to the spectra

- 213 of 234 COS pairs (91%) changed; FOS and STIS are unchanged apart from the
  objects fixes 3 and 4 rescued, and floating-point noise at the 1e-5 level in FOS
  (a per-pixel Python loop was vectorised).
- Exposure rows entering the co-add: 4,619 → 9,395.
- Rest-frame coverage gain, median +174 Å (+78 Å inside the ICA range 1260–3000 Å).
- S/N in already-covered wavelengths: median ×1.00. **The gain is coverage, not depth** —
  FUVA and FUVB barely overlap.

---

## 4. The ICA measurement fixes (2026-09-23)

### Fix 5 — C IV equivalent width was biased low

Full detail: `diagnostics/CIV_EW_Indexing_Investigation.md`.

In `get_CIV_parameters` the continuum was built full-length from the ICA grid but
indexed with the loop counter, which runs over the 1500–1600 Å subset:

```python
continuum = wave_r * m + b                       # full grid, starts at 1260 A
EW = (wave_r >= 1500) & (wave_r <= 1600)
for i in range(len(wave_r[EW])):
    CIV_EW += max(((flux_r[EW][i] - continuum[i]) / continuum[i]) * dlam, 0)
                    # masked subset      full array
```

The continuum was therefore read ~240 Å too blue. Quasar continua fall with wavelength,
so the divisor was always too large and **every EW was biased low**. Fixed by masking
once (`cont_EW = continuum[EW]`) in `get_CIV_parameters` and `plot_CIV_analysis`
(`ica/manual_fix.py`) and `get_CIV` (`ica/run_ica.py`).

**This is Trevor's code, unmodified.** The identical loop is in every copy of `get_CIV`
in the untouched `HST/` archive across 2022 and 2023. His published `CIV_EW_tvm` column
reproduces the buggy value to a median of 0.50 Å, while the corrected value differs by
19.16 Å — so **his published equivalent widths carry the same bias**. His plots show the
correct shaded region (`fill_between` used `continuum[EW]`) beside an incorrect number in
the legend.

Size of the correction, over 485 objects:

| | |
|---|---|
| corrected ÷ previous EW | median **1.493** |
| EW increased for | 91% of objects |
| blueshift unchanged for | 35% of objects |
| median blueshift change where it moves | 137 km/s |

Strongly dependent on line strength — weak lines move most:

```
previous EW    n    median ratio   fraction whose blueshift moves
  <10 A       49       2.91                 90%
  10-25       91       1.84                 88%
  25-50       98       1.58                 80%
  50-100     102       1.35                 55%
  >100 A     140       1.11                 28%
```

A side effect: 12 objects previously recorded as `EW = 0` now measure a real equivalent
width. Their line had been sitting just below a continuum read too high.

### Fix 6 — pixels with near-zero error pinned the fit

Full detail: `diagnostics/Collapsed_Fits_Investigation.md`.

47 rows in the previous run read `CIV_blueshift = 9574.1 km/s` with `CIV_EW` of 0 or
blank. That number is a fingerprint: when EW is zero the half-flux index lands on the
first pixel of the window, and 1500 Å converts to 9574.1 km/s. It means the ICA
reconstruction came back as a flat line at zero.

Cause: the fit weights pixels by `ivar = 1/error**2`. A pixel with an error of ~0 gets
`ivar` ≈ 1e29 and acts as an infinitely strong constraint — any model not passing exactly
through it costs ~1e30, so the only survivor is all-zero weights. Healthy spectra have
`max(ivar)/median(ivar)` of order 5; affected ones order 1e29.

Fixed by a guard inside `ICA_fit` (`ica/run_ica.py`), which every call site funnels
through: drop non-finite `ivar`, and any `ivar` above `_IVAR_MAX_FACTOR` (1e4) × median.

Negative error values, which look like the obvious culprit, are **not** the cause —
`ivar = 1/error**2` squares the sign away, so a large negative error merely down-weights
its own pixel. They are more common among healthy fits than failing ones.

---

## 5. Outstanding — the 36 spectra that still need regenerating

**This is the one thing that must be excluded from any figure made today.**

Fix 2 above introduced a regression. `cont_filtered()` ended with

```python
med_flux[bad] = np.median(med_flux[good]) if good.any() else 1.0
```

The `else 1.0` fires when an exposure's continuum is NaN at *every* pixel — which happens
for exposures almost, but not entirely, inside an emission window, so the all-emission
guard at the top of the function does not catch them. `coadd.rebin()` then divides that
exposure by **1.0**, leaving it in raw physical units while its neighbours come back
continuum-normalised. The exposure enters the co-add about 1e13 too small.

Measured on one affected spectrum:

```
MR 2251-178_COS   1260-1450 A   median |flux| 3.67        median err 0.107
                  1500-1600 A            2.67e-13               6.73e-15
```

| | spectra with an un-normalised block |
|---|---|
| `RebinnedSpec_master` (before fix 2) | 0 of 576 |
| `RebinnedSpec_master_v2` (after) | **36 of 594**, 18 overlapping C IV |

It is also where most of the collapsed fits came from: of the 47, **30 trace to this
regression** and only 17 are genuine (the pre-fix run had 11 of 418).

**Status:** the code is fixed (`_flat_continuum_level()` in `rebinning/spec_morph.py`
returns the exposure's own flux scale, never 1.0), but regenerating the 36 co-adds
requires the raw data on the F: drive, which is not attached. **Until then exclude the
36 objects in `diagnostics/affected_unnormalised.csv`.**

When F: is available:

```bash
# delete the 36 stale outputs, then re-run only those (skip_existing is on by default)
python -m rebinning.run_rebin --master-catalog --data-dir "F:/Richards-data/HST Paper/raw_data" \
    --output-dir RebinnedSpec_master_v2 -j 6
python -m ica.run_all_objects --master --rebin-dir RebinnedSpec_master_v2 \
    --results-file ICA_Results_Rebin_master_v2_fixed.csv \
    --plots-dir ICA_Plots_Rebin_master_v2_fixed
```

Note `run_rebin.py` prints a `→` character that raises `UnicodeEncodeError` when stdout
is redirected to a file on Windows; prefix with `PYTHONIOENCODING=utf-8`.

---

## 6. Caveats when making comparison figures

1. **Exclude the 36.** Several would look convincing by eye — `MR 2251-178_COS` shows a
   clean emission line — so they cannot be caught visually.
2. **Six rows have `EW = 0` because there is no data at C IV**, not because the line is
   absent: `PG0953+415_STIS`, `LEDA 138607_COS`, `Ton 488_COS`,
   `2MASS J09460395+0139235_COS`, `QSO J1422-0001_COS`, `87GB 121807.5+173455_COS`. All
   have zero good pixels in 1500–1600 Å. They are already excluded from the usable list.
3. **A quality cut independent of fit success is still needed.** Some objects return a
   number from a spectrum with no data in the C IV window — "the fit ran" is not a
   quality test. Coverage inside 1500–1600 Å is the metric to cut on.
4. **Do not compare new EW against Trevor's published `CIV_EW_tvm` directly** — his
   carry the old bias, so the comparison shows the fix, not a change in the data. To
   compare like with like, use `ICA_Results_Rebin_master.csv`.
5. **Measure S/N with care.** `median(|flux/error|)` is unreliable on affected spectra
   because flux and error are both raw there. Prefer an empirical estimate from
   pixel-to-pixel scatter (DER_SNR).
6. **Mrk 817 and NGC 5548 remain provisional.** A temporary cap
   (`_COS_EXPOSURE_CAP` in `rebinning/read_spec_data.py`) keeps only the 100
   highest-S/N exposures of each; they have ~1,500 and ~800. The threshold is a runtime
   convenience, not a science choice, and was flagged for discussion with G. Richards.
7. **Manual ICA overrides have not been re-reviewed.** 36 of the 92 entries are on COS
   spectra, whose shape changed. The overrides pass (`--overrides-only`) has not been
   re-run against the fixed code.

---

## 7. Files changed

| file | date | what |
|---|---|---|
| `rebinning/read_spec_data.py` | 09-10 | flat readers rewritten (fixes 1, 3, 4) |
| `rebinning/spec_morph.py` | 09-10, 09-29 | `cont_filtered()` (fix 2, then fix 7) |
| `ica/spec_morph.py` | 09-10 | same all-NaN continuum guard, narrower |
| `ica/manual_fix.py` | 09-23 | EW continuum indexing (fix 5) |
| `ica/run_ica.py` | 09-23 | EW indexing + `ICA_fit` ivar guard (fixes 5, 6) |
| `ica/run_all_objects.py` | 09-10 | added `--results-file` / `--plots-dir` |

Backups of the two modules edited on 09-23: `diagnostics/_backup_manual_fix.py`,
`diagnostics/_backup_run_ica.py`.

### Supporting analysis

| file | what |
|---|---|
| `diagnostics/CIV_Data_Loss_Investigation.md` | the rebinning defects, with the census |
| `diagnostics/CIV_EW_Indexing_Investigation.md` | the EW bias, and the trace to Trevor's code |
| `diagnostics/Collapsed_Fits_Investigation.md` | the collapsed fits, and §7 the regression |
| `diagnostics/civ_coverage_audit.csv` | 615 rows, C IV coverage raw vs rebinned |
| `diagnostics/v2_vs_master.csv` | per-spectrum diff, old vs new rebin |
| `diagnostics/ew_bias_master.csv`, `ew_bias_v2.csv` | per-object EW before/after |
| `diagnostics/affected_unnormalised.csv` | **the 36 to exclude** |
| `diagnostics/usable_objects_v2_fixed.txt` | **the 478 safe to use** |

---

## 8. Corrections made during this work

Recorded so they are not re-derived from the earlier logs, which contain superseded
statements:

- The C IV equivalent-width bias was first estimated at 15–35% from a 16-object
  sub-sample. The sample-wide value is a **median factor of 1.49**.
- The blueshift was first reported as unaffected by the EW fix. It is unchanged for
  only ~35% of objects; the median change where it moves is 137 km/s.
- The collapsed fits were first attributed to negative error values, then to a
  near-zero continuum in the ICA morphing step. Both are wrong; the cause is
  near-zero *errors* producing an enormous inverse variance.
- The collapsed fits were first described as a pre-existing defect affecting Trevor's
  results too. **30 of the 47 were introduced by fix 2's regression**; only 17 are
  genuine.
- A result that "28 of 34 collapsed fits were rescued by the ivar guard" is **void**:
  28 of those 29 objects sit on the 36 un-normalised spectra, so the fits ran against
  corrupted data.
- A coverage statistic quoted in those logs measured the span between the first and last
  good pixel rather than the fraction actually filled, and overstated coverage badly
  (`PG0953+415_STIS` read 100% where the true figure is 12%).
