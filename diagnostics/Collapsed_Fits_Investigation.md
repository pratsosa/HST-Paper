# Why some ICA fits collapse to a zero reconstruction

**AP, 2026-09-20**

**Question.** 133 of 594 rebinned spectra carry unmasked negative flux errors, but only
~20 ICA fits collapse. Which objects fail, which survive despite the same symptom, and
is a sweeping fix justified?

**Answer.** Negative errors are not the cause and do not predict failure. Of the 47
collapsed fits, 16 are not failures at all (a genuine `EW = 0`), leaving 34. Those are
caused by a handful of pixels with *near-zero* errors — the opposite extreme — which
claim infinite precision and pin the fit to zero. De-weighting only those pixels
recovers **28 of the 34**, several at high signal-to-noise. Only **6** are genuinely
bad data and should be dropped; they are listed in §4.

---

## 1. Negative errors do not predict collapse

Scan of 219 objects (all 47 degenerate + all 133 with negative errors + 60 healthy
controls), `diagnostics/collapse_scan.csv`:

```
                 DEGENERATE  fit OK  no fit
negative errors?
  no                     27      52       7
  yes                    23     108       2
```

Negative errors are **more** common among healthy fits (108/160) than failing ones
(23/50). 26 of the 47 degenerate fits have none at all.

The reason is arithmetic: the fit uses `ivar = 1/errs**2`. Squaring destroys the sign,
so a large negative error produces a *small* ivar and simply **down-weights** its own
pixel. A negative error is an unphysical value and should still be masked for
cleanliness, but it is not what breaks the fit.

Two other candidate explanations also failed:

| hypothesis | result |
|---|---|
| `cont_spec` going non-positive inside `ica/spec_morph.morph2` | 5/50 degenerate vs 12/160 healthy — no discrimination |
| large \|flux\| reaching the fit | degenerate median 7.1, healthy median 11.3 — failing objects are *less* extreme |
| the iterative BAL mask covering everything | degenerate 36% of pixels flagged vs healthy 50% — the opposite |

(An earlier note in this investigation attributed the collapses to a −1.2e6 outlier in
NGC 985 produced by near-zero continuum division. That is real for NGC 985 but is not
the general mechanism — the scan above supersedes it.)

---

## 2. What actually happens: pixels with near-zero error

`get_ICA` passes `ivar = 1/errs**2` to `ICA_fit`. A pixel whose error is ~0 gets an
enormous ivar:

```
                 median max(ivar)   implied smallest error
collapsed (n=12)     1.4e29              0.00000
healthy   (n=12)     8.4                 0.375

max(ivar) > 1e6 :  12 of 12 collapsed,  0 of 12 healthy
```

That is a clean separator. Such a pixel behaves as an infinitely strong constraint: at
weights = 0 the model is zero and chi-squared stays finite, but moving any weight away
from zero costs ~1e30. So the all-zero reconstruction genuinely **is** the least-squares
optimum under that weighting — the optimiser is not failing, it is obeying the data.

Confirmed by experiment: restarting L-BFGS-B from non-zero weights (0.05) does not
rescue these fits, it makes them worse (redchi 1e33–1e44). Only de-weighting the
offending pixels changes the outcome.

`CIV_blueshift = 9574.1` km/s is the fingerprint of the resulting zero reconstruction:
with `CIV_EW == 0`, the half-EW index lands on pixel 0 of the window, and
`blueshift(1500 A)` = 9574.1.

---

## 3. The three causes, counted

Of the 50 degenerate fits in the scan (`diagnostics/degenerate_causes.csv`):

| cause | n | median usable pixels | median ivar | median redchi |
|---|---|---|---|---|
| **A** starved — fewer than 100 usable pixels | 4 | 12 | 0.0005 | 0.000 |
| **B** pinned by near-zero-error pixels | 30 | 1112 | 0.738 | 6.53 |
| **C** fit ran fine, C IV genuinely never exceeds the continuum | 16 | 844 | 1.002 | 0.431 |

Class **C is not a fault**: 39 of the 40 healthy controls fall in the same class. Those
16 objects have a real, converged fit whose reconstruction simply stays at or below the
local continuum across 1500–1600 Å, so `EW = 0` is a legitimate measurement. They should
be reported as EW = 0, not as failures.

So the honest count is **34 genuine failures (A + B)**, not 47, and certainly not 133.

---

## 4. Is it just bad data? The decisive test

Each of the 34 real failures (classes A and B) was refitted with **only** the
near-zero-error pixels de-weighted, and nothing else changed
(`diagnostics/rescue_test.csv`). A "rescue" means the refit returned a finite,
non-zero equivalent width instead of the zero reconstruction.

```
            rescued   not rescued
class A         1          3
class B        27          3
            -------------------
total          28          6
```

**28 of 34 recover.** An initial visual check of six hand-picked panels
(`diagnostics/collapsed_fits_visual_check.png`) suggested most were simply bad data;
that impression did not survive testing the full set, and this section supersedes it.
Several of the recovered objects are high quality:

```
MR 2251-178_COS      S/N 34.2   EW  65.2   blueshift  197.5
Mrk 1513_COS         S/N 22.9   EW  37.4   blueshift  128.4
QSO B1011-040_COS    S/N 21.5   EW   1.0   blueshift  473.5
NGC 985_COS          S/N 14.5   EW  10.1   blueshift  -78.8
NGC 5548_STIS        S/N 12.6   EW 117.3   blueshift   -9.7   (100% coverage)
```

Of the 28 rescued: 13 have median S/N > 3, 6 have 2-3, and 9 have S/N < 2. The
low-S/N end still needs visual review before the numbers are trusted -- a fit that
runs is not a fit that is right -- but these are not data-quality failures.

### The 6 that are genuinely bad data

| object | median S/N | coverage of 1260-3000 A | reason |
|---|---|---|---|
| `2MASS J10051419+5302398_COS` | 0.73 | 0% | no usable pixels in the fitting range |
| `2MASS J09460395+0139235_COS` | 0.72 | 1% | no usable pixels in the fitting range |
| `PG 1216+069_STIS` | 4.16 | 2% | far too little coverage to constrain 10 components |
| `3C 95_COS` | 7.16 | 3% | same |
| `QSO J1422-0001_COS` | 1.01 | 27% | S/N ~ 1, noise dominated |
| `87GB 121807.5+173455_COS` | 2.90 | 54% | still returns zero; cause not isolated |

The first four have a quotable justification: the spectrum covers almost none of the
1260-3000 A range the ICA components span. The last two are lower-confidence calls.

---

## 5. Recommendation

1. **Guard the weights.** Where `ivar = 1/errs**2` is formed in `get_ICA`, ignore
   pixels whose error is implausibly small or non-positive. Two lines, in the ICA
   stage only, **no re-rebinning**, and it recovers all 28.
2. **Drop the 6** listed above, citing coverage and S/N.
3. **Report class C (16 objects) as a measured `EW = 0`**, not as a failure --
   39 of 40 healthy controls sit in the same class.
4. Log every dropped object with its reason, so the gap between 594 spectra and the
   final fit count is explained rather than silent.

---

## 6. Outcome (2026-09-23) — both fixes applied and re-run

Both defects were fixed and the ICA stage re-run over `RebinnedSpec_master_v2/`.
Results: `ICA_Results_Rebin_master_v2_fixed.csv`, plots
`ICA_Plots_Rebin_master_v2_fixed/`. The two earlier result sets are untouched:

| file | rebin | EW indexing | ivar guard |
|---|---|---|---|
| `ICA_Results_Rebin_master.csv` | original | buggy | none |
| `ICA_Results_Rebin_master_v2.csv` | fixed | buggy | none |
| `ICA_Results_Rebin_master_v2_fixed.csv` | fixed | **fixed** | **yes** |

### Code changed

| file | function | change |
|---|---|---|
| `ica/manual_fix.py` | `get_CIV_parameters` | `cont_EW = continuum[EW]` masked once, indexed in the loop |
| `ica/manual_fix.py` | `plot_CIV_analysis` | same |
| `ica/run_ica.py` | `get_CIV` | same |
| `ica/run_ica.py` | `ICA_fit` | drop non-finite ivar, and ivar above `_IVAR_MAX_FACTOR` (1e4) x median |

The guard lives inside `ICA_fit`, which every call site funnels through, so all three
component sets and the f2500 path are covered without touching them individually.

### Result

```
degenerate rows (EW = 0 / NaN):  47  ->  6      41 recovered, 0 regressions
fits produced:                  517 -> 517      exception list byte-identical
EW on the 470 objects good in both: median ratio 1.495, increased for 92%
blueshift: unchanged for 34%, median change 137 km/s where it moves
```

Which fix rescued what, over the 50 previously-degenerate fits:

| cause | n | now give a real measurement |
|---|---|---|
| B pinned by zero-error pixels | 30 | **27** (the ivar guard) |
| C fit fine, EW was 0 | 16 | **12** (the EW fix) |
| A starved | 4 | 2 (but see caveat below) |

The class-C result matters scientifically: **12 objects previously recorded as
`EW = 0` were an artefact of the continuum being read too high**, not real
non-detections. With the continuum taken from the right wavelengths the line rises
above it.

### Caveat — a quality cut is still needed

`2MASS J10051419+5302398_COS` now returns `EW = 0.013 A`, `blueshift = -562.8`, from a
spectrum with **0% coverage** of the 1260-3000 A fitting range. The guard let a
meaningless fit through. "Produced a number" is not a usability test; the explicit
drop-list in section 4 is still required.

### The 6 rows still showing `EW = 0`

| object | S/N | coverage | assessment |
|---|---|---|---|
| `2MASS J09460395+0139235_COS` | 0.72 | 1% | bad data, drop |
| `QSO J1422-0001_COS` | 1.01 | 27% | bad data, drop |
| `87GB 121807.5+173455_COS` | 2.90 | 54% | marginal |
| `Ton 488_COS` | 7.71 | 10% | **needs eyeballing** - decent S/N |
| `PG0953+415_STIS` | 9.35 | 100% | **needs eyeballing** - good S/N, full coverage; an EW of 0 here is suspicious |
| `LEDA 138607_COS` | 11.97 | 23% | **needs eyeballing** - good S/N |

The last three are good-quality spectra whose C IV reconstruction still does not rise
above the local continuum. That can be genuine (an absorbed or intrinsically weak line)
or a continuum-window problem, and only visual inspection will tell.

`PG 1216+069_STIS` and `3C 95_COS` raise exceptions rather than producing rows; they
did so before these fixes too (the exception list is unchanged).

---

## 7. CORRECTION (2026-09-29) — most of "Defect 2" was self-inflicted

AP spotted that `MR 2251-178_COS` and `Mrk 1513_COS`, which section 6 listed as
high-quality rescues, have visibly bad data across C IV. Investigating that overturned
a large part of this document.

### What was actually wrong

`rebinning/spec_morph.py` `cont_filtered()` ended with a guard added on 2026-09-10:

```python
med_flux[bad] = np.median(med_flux[good]) if good.any() else 1.0
```

The `else 1.0` fires when an exposure's continuum is NaN at every pixel — every pixel
emission-masked, so the interpolation has no finite neighbour. `coadd.rebin()` then
divides that exposure by **1.0**, leaving it in raw physical units (~1e-14) while its
neighbours come back continuum-normalised (~1). The exposure enters the co-add about
1e13 too small.

Measured directly:

```
MR 2251-178_COS   1260-1450 A   median |flux| 3.67        median err 0.107
                  1500-1600 A            2.67e-13               6.73e-15
NGC 985_COS       1700-2000 A                 2.89                   0.314
```

### Blast radius

```
RebinnedSpec_master     (before the 2026-09-10 fix):   0 of 576 spectra
RebinnedSpec_master_v2  (after):                      36 of 594 spectra, 18 over C IV
```

### This explains Defect 2 itself

Those blocks read as ~0 flux with ~0 error once the ICA normalises by the whole
spectrum's median. `ivar = 1/err**2` then reaches ~1e29 — which is exactly the
"pixels with near-zero error" mechanism described in section 2. So the causal chain is:

1. the `else 1.0` fallback leaves a block un-normalised;
2. that block reads as zero flux with zero error;
3. ivar ~1e29 pins the ICA fit at zero weights.

Splitting the 47 collapsed fits by whether the spectrum carries an un-normalised block:

```
caused by this bug : 30
genuine            : 17     (the old run, with no affected spectra, had 11 of 418)
```

**Sections 1-5 are still correct about the mechanism** (near-zero-error pixels pin the
fit; negative errors are irrelevant). What they got wrong is the origin: this was
overwhelmingly introduced on 2026-09-10, not inherited.

### Section 6's rescue list is not trustworthy

**28 of the 29 objects listed as rescued by the ivar guard sit on affected spectra.**
The guard stopped the pinning, but the fit then ran against a spectrum with a corrupted
block, producing plausible-looking numbers from bad data. The single unaffected one,
`2MASS J10051419+5302398_COS`, has 0% coverage of the fitting range and is useless for a
different reason. **The keep/drop triage in section 6 should be disregarded.**

The S/N figures quoted there are also unreliable: `median(|flux/error|)` computed inside
an un-normalised block divides raw flux by raw error, so it looks healthy. Empirical S/N
from pixel-to-pixel scatter gives 6.2 rather than 34.1 for `MR 2251-178_COS`, and 3.6
rather than 16.2 for `Mrk 1513_COS`.

### Fixed, but the affected spectra cannot be regenerated yet

`cont_filtered()` now falls back to the exposure's own flux scale via a new
`_flat_continuum_level()` helper, never to 1.0. **The 36 affected spectra need
re-rebinning, which requires the raw data on F: — that drive is no longer attached**,
so they cannot be regenerated in this session.

Until then: exclude the 36 objects in `diagnostics/affected_unnormalised.csv` from
`ICA_Results_Rebin_master_v2_fixed.csv`. The ivar guard stays justified for the 17
genuine cases. The C IV equivalent-width fix is independent of all this and unaffected.
