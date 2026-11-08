# The C IV equivalent-width continuum indexing defect

**AP, 2026-09-20**

**Question.** `get_CIV_parameters()` divides the C IV line flux by a continuum read
from the wrong part of the array. Is this real, did Trevor make the same mistake, are
his saved values affected, were his plots right even if the numbers were wrong, or are
we simply using his code incorrectly?

**Answer, in one line.** The defect is real, it is Trevor's own code unchanged, his
published `CIV_EW_tvm` values carry it, we are calling the function exactly as he did,
and his plots contain the same wrong number in the legend while the shaded region
beside it is drawn correctly.

---

## 1. What the code does

`ica/manual_fix.py`, `ICAManualFixProcessor.get_CIV_parameters()` — the function that
produces every row of `ICA_Results_Rebin_master*.csv`:

```python
m, b = np.polyfit(np.concatenate((wave_r[cont1], wave_r[cont2])),
                  np.concatenate((flux_r[cont1], flux_r[cont2])), 1)
continuum = wave_r * m + b                       # full-length array

EW = ((wave_r>=EW_region[0]) & (wave_r<=EW_region[1]))   # boolean, 1500-1600 A

for i in range(len(wave_r[EW])):
    CIV_EW += max(((flux_r[EW][i] - continuum[i]) / continuum[i]) * dlam, 0)
```

`flux_r[EW][i]` and `wave_r[EW][i]` index the **selected window**. `continuum[i]`
indexes the **full array**. They are not the same pixels.

`wave_r` is `wave_ica`, the ICA reconstruction grid, which begins at 1260 Å (`mod`
components), 1265 Å (`high`) or 1275 Å (`low`) — never at 1500 Å. So the C IV window
starts around index 740–760, and the loop reads the continuum from the blue end of the
grid instead:

```
object            grid          continuum read at      flux pixels at
NGC 985_COS       1260-3000 A   1260.1-1344.0 A        1500.0-1599.9 A
0115+027_FOS      1265-3000 A   1265.0-1349.3 A        1500.0-1599.9 A
LEDA 3095694_COS  1260-3000 A   1260.1-1344.0 A        1500.0-1599.9 A
```

**The continuum is sampled roughly 240 Å too blue.**

### Why the bias has a consistent sign

The continuum is the straight line `m*wave + b` fitted through 1445–1465 Å and
1700–1705 Å. Quasar continua fall with wavelength in f_lambda, so `m < 0` and the line
is **higher** at 1260–1344 Å than at 1500–1600 Å. Dividing by a continuum that is too
high makes `(flux - continuum)/continuum` too small, so **EW is always biased low**.
That is why the error is a systematic offset rather than scatter.

### The blueshift is affected too, but only slightly

The blueshift is the wavelength at which the running EW passes half its total. Because
the continuum error is a smooth linear ramp rather than a local distortion, the half-way
point usually lands on the same pixel, and when it moves it moves by whole pixels of the
69 km/s ICA grid. Measured over 16 objects: unchanged for 7, and otherwise shifted by
exactly one to three grid pixels.

```
|change in blueshift|: median 69 km/s, mean 52 km/s, max 207 km/s; nonzero in 9 of 16
```

(An earlier note in this investigation said the blueshift was completely unaffected.
That was based on three objects and is wrong — the figures above supersede it.)

---

## 2. Is it Trevor's, or ours?

**Trevor's.** The identical loop appears in every copy of `get_CIV` in the untouched
`HST/` archive, across three different notebooks and two different years:

| file | date context |
|---|---|
| `HST/Notebooks_final_MarchApril2023/ICA_HST_2023code.ipynb` | March 2023, his most recent |
| `HST/Notebooks_final_MarchApril2023/HST_CIV_March2023.ipynb` | March 2023 |
| `HST/Notebooks/HST_CIV_Sulentic2007_HSLA2015_CIVFitregions_wfullICA_newSpec_2022Sept.ipynb` | Sept 2022 |
| `HST/Notebooks/HST_CIV_Sulentic2007_HSLA2015_CIVFitregions_wfullICA.ipynb` | earlier |

All four contain `(flux_r[EW][i] - continuum[i]) / continuum[i]` verbatim. Nothing in
`Trevor Code/` or `HST Paper/` changed it; the migration copied it faithfully.

**We are not using it incorrectly.** His call site in
`ICA_HST_2023code.ipynb` is

```python
blue, ew = get_CIV(wave_arb, flux_arb, wave_ica, flux_ica, name, ax=ax2)
```

and ours is `get_CIV_parameters(wave_arb, flux_arb, wave_ica, flux_ica, name)` — the
same four arrays in the same order, with the same default `EW_region` and
`cont_region`. There is no calling convention under which `continuum[i]` would be
correct, because the function builds `continuum` itself, at full length, from the
`wave_r` it was handed.

---

## 3. Do his saved values carry the bias?

**Yes.** Test: run the pipeline over Trevor's own `HST/RebinnedSpec_2022Aug11/` spectra,
compute the EW both ways, and compare with the `CIV_EW_tvm` column he saved in
`HST/final_tables/HST_CIV_Sulentic2007_HSLA2018_finalprops.csv` (207 objects).

```
                 median |difference from CIV_EW_tvm|
as-coded (buggy)              0.50 A
corrected                    19.16 A
```

Several reproduce to three decimal places — 3C215 154.116 vs 154.116, 3C336 88.011 vs
88.011, 3C345 138.927 vs 138.927, 3C273 13.626 vs 13.625. His blueshifts reproduce
exactly too (1231.168, −217.002, −78.781, −9.681, 749.226 …), which also shows the
blueshift is quantised to the 69 km/s ICA grid.

The residual ~0.5 Å scatter is expected: our rebinned spectra are not pixel-for-pixel
identical to his (a known, separately documented difference), so exact agreement is not
available. The 38× gap between the two candidates leaves no ambiguity about which
formula produced his numbers.

`diagnostics/tvm_ew_check.csv` has the per-object table.

---

## 4. Were his plots right?

Partly, and inconsistently. In the same function:

```python
ax.fill_between(wave_r[EW], continuum[EW], flux_r[EW], ...,
                label="EW = %.1f A" % CIV_EW)
```

The shaded region uses `continuum[EW]` — **the correct pixels**. The legend prints
`CIV_EW` — **the incorrect number**. So in every one of Trevor's C IV panels the blue
shaded area is the true line-over-continuum region, while the equivalent width quoted
next to it is 15–35% too small. The two disagree, and the plot alone would not reveal it
unless someone integrated the shading by eye.

The magenta continuum line drawn on the plot (`ax.plot(wave_r, continuum, "-m")`) is
also correct across the full grid. Only the arithmetic in the loop reads the wrong
element.

---

## 5. Size of the effect

Corrected ÷ as-coded EW, over 16 objects from Trevor's own sample:

```
median x1.338    range x1.002 - x2.567
```

So the correction **raises** every EW, typically by about a third. Objects with small
EW move proportionally most, because the fixed continuum offset is a larger fraction of
a weak line.

Sample-wide statistics over all objects in both rebin runs are produced by
`diagnostics/ew_bias_scan.py` (outputs `ew_bias_master.csv`, `ew_bias_v2.csv`).

---

## 6. What this means

Sample-wide figures (`ew_bias_master.csv`, 480 objects; `ew_bias_v2.csv`, 485) —
these supersede the 16-object estimates in sections 1 and 5 above, which came from
Trevor's clean sub-sample and understate the effect:

| | old run | v2 run |
|---|---|---|
| corrected ÷ as-coded EW, median | **1.457** | **1.493** |
| EW increased for | 88% | 91% |
| underestimated by >25% | 73% | 78% |
| underestimated by >50% | 45% | 49% |
| blueshift unchanged | 38% | 35% |
| \|Δblueshift\| where it moves | median 137, max 10359 km/s | same |

The effect is strongly EW-dependent — weak lines move most, because the fixed
continuum offset is a larger fraction of a weak line:

```
as-coded EW    n    median ratio   median |dblue|   fraction whose blueshift moves
  <10 A       49       2.91           277 km/s              90%
  10-25       91       1.84           138                   88%
  25-50       98       1.58            69                   80%
  50-100     102       1.35            69                   55%
  >100 A     140       1.11             0                   28%
```

- Every `CIV_EW` ever produced by this pipeline — Trevor's published values, the
  `ICA_Results_Rebin_master.csv` and `ICA_Results_Rebin_master_v2.csv` — is low, by a
  median factor of ~1.46.
- `CIV_blueshift` is affected more than first thought: unchanged for only ~38% of
  objects, and 6% move by more than 500 km/s.
- The fix is to index the continuum with the same mask as the flux
  (`continuum[EW][i]`, matching what `fill_between` already does). It requires no
  re-rebinning: only the ICA measurement step re-runs.
- Because the bias is systematic and one-directional, any published correlation
  involving C IV EW has a shifted normalisation. Whether it changes a *conclusion*
  depends on whether that analysis used EW in absolute terms or only as a ranking.

 The decision is AP's.

---

## 7. Fixed (2026-09-23)

Applied as `cont_EW = continuum[EW]` masked once outside the loop (indexing
`continuum[EW][i]` inside the loop would re-mask the array on every iteration), at all
three sites: `get_CIV_parameters` and `plot_CIV_analysis` in `ica/manual_fix.py`, and
`get_CIV` in `ica/run_ica.py`.

Verified against the independent measurements in `ew_bias_v2.csv`:

```
0115+027_FOS        new 295.943  expected 295.943   (was 258.687)
LEDA 3095694_COS    new 131.969  expected 131.969   (was 104.756)
Mrk 1298_COS        new  19.924  expected  19.924   (was  12.318)
1401+0951_FOS       new  91.646  expected  91.646   (was  59.912)
```

Sample-wide over the 470 objects measurable before and after: EW ratio median 1.495,
increased for 92%; blueshift unchanged for 34%, median change 137 km/s where it moves.

A side effect worth noting: 12 objects that previously recorded `EW = 0` now measure a
real equivalent width. Their line had been sitting just below a continuum that was
read too high. See `Collapsed_Fits_Investigation.md` section 6.

Results: `ICA_Results_Rebin_master_v2_fixed.csv`.
