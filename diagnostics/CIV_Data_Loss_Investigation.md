# Where the C IV data went — investigation of `RebinnedSpec_master/`

**AP, 2026-09-10**

**Question.** Some objects in `RebinnedSpec_master/` have no data across C IV
(rest 1549.48 Å) even though MAST clearly holds C IV coverage for them. Where in the
pipeline is that data lost?

**Answer.** Three independent bugs in `rebinning/`, all of them silent, plus one that
crashes loudly. Together they cost **57 spectra** their C IV: 24 to a reader that only
ever reads the first row of each COS x1d table, 15 to a continuum fit that returns zero,
and 18 to two crashes that produce no output file at all. All 57 have been
addressed with a patched pipeline (`diagnostics/rebin_patched.py`): 47 now have C IV, 8
recover everything the telescope actually recorded, and 2 turn out to have no C IV at all.
Nothing in `rebinning/` has been changed.

Related existing documentation: `.claude/Migration_Log.md` (Phase 5 — the MAST download
and master-catalog rebin, including the earlier `_x1dsum` double-counting and empty-x1d
fixes), `Trevor Code/SpectralCoverage_CrossCheck_Log.md` (rest- vs observed-frame
conventions per instrument), `CLAUDE.md` → Known Issues (the two diverged `spec_morph.py`
copies).

---

## 1. The census

`diagnostics/audit_civ_coverage.py` walks every (object, instrument) pair the rebin run
was asked to process — the same 615 pairs `pipeline.catalog_bridge.build_object_list()`
produces — and compares three rest-frame coverages:

| | what it is |
|---|---|
| `raw_all` | every row of every exposure file the reader globs |
| `raw_read` | only the rows the reader actually keeps |
| `output` | pixels in the rebinned FITS with finite, non-zero flux |

"C IV is present" means the median per-pixel flux/error over rest 1520–1580 Å exceeds
0.3 in the raw data, and that more than half the rest 1500–1600 Å pixels are good in the
output. The measured signal-to-noise is written to the CSV so that threshold can be
re-judged without re-running.

Result — `diagnostics/civ_coverage_audit.csv`, 615 rows:

| verdict | COS | FOS | STIS | total | meaning |
|---|---|---|---|---|---|
| `OK` | 110 | 161 | 130 | **401** | C IV in the output |
| `NO_CIV_ON_SKY` | 64 | 6 | 44 | 114 | C IV outside the wavelengths observed |
| **`LOST_ROWS`** | 24 | 0 | 0 | **24** | C IV in a detector segment the reader never reads |
| `NO_CIV_SIGNAL` | 21 | 0 | 1 | 22 | C IV nominally covered but the raw flux is noise |
| `NO_RAW_FILES` | 10 | 2 | 9 | 21 | catalogued but nothing downloaded |
| **`NO_OUTPUT`** | 5 | 13 | 0 | **18** | rebin crashed; no file written |
| **`LOST_COADD`** | 10 | 0 | 5 | **15** | C IV read in, then destroyed inside `coadd.py` |

The three bold rows are the recoverable loss: **57 spectra**.

`NO_CIV_SIGNAL` is worth a note because it is easy to mistake for a bug. Most of these are
COS G140L/800 snapshots of z ≈ 0.28 quasars: the wavelength array nominally reaches past
C IV at ~1980 Å observed, but the detector's throughput there is gone — median flux is
exactly 0.0 with only 42% of pixels flagged good. `Cut_Edge_Pix` trims that tail, correctly.
There is nothing to recover.

`NO_RAW_FILES` is a download-stage gap, outside the scope of this investigation. All 21 are
objects whose folder is missing or empty under `F:\Richards-data\HST Paper\raw_data`.

---

## 2. Bug A — only the first row of each COS exposure is read (`LOST_ROWS`, 24 spectra)

`rebinning/read_spec_data.py`, `read_cos_flat()` and `read_stis_flat()`:

```python
wavelength = data['WAVELENGTH'][0]
flux       = data['FLUX'][0]
```

A COS x1d table has **one row per detector segment**, and each row is a physically separate
piece of the spectrum, not a repeat of the same one:

- FUV (G130M, G160M, G140L): 2 rows — `FUVA`, `FUVB`
- NUV (G185M, G225M, G285M, G230L): 3 rows — `NUVA`, `NUVB`, `NUVC`

Reading `[0]` keeps `FUVA`/`NUVA` and discards the rest. For a G130M+G160M object that
leaves a hole in the middle of the coverage: G130M FUVA reaches ~1450 Å and G160M FUVA
starts at ~1563 Å, while the two FUVB rows that would have filled 1375–1587 Å are thrown
away.

The audit checked, for each of the 24 `LOST_ROWS` objects, which row actually carries the
C IV signal. In **all 24** it is a discarded one:

```
2MASS J10423110+2213003    NUVB row1 snr=0.91      Mrk 1044      FUVB row1 snr=3.53
2MASS J11493175+2230285    NUVB row1 snr=2.65      Mrk 352       FUVB row1 snr=1.48
7C 081438.10+224634.00     NUVB row1 snr=1.95      NGC 5548      FUVB row1 snr=6.66
CLASS-B1152+199-ABS01      NUVB row1 snr=1.67      NGC-4051      FUVB row1 snr=5.29
ESO 462-9                  FUVB row1 snr=1.67      NGC-7469      FUVB row1 snr=6.52
ESO 499-41                 FUVB row1 snr=0.73      NGC4395       FUVB row1 snr=0.78
ICRF J154743.5+205216      NUVB row1 snr=3.27      PG 1216+069   NUVC row2 snr=2.20
ICRF J213745.1-143255      NUVB row1 snr=3.38      QSO B2128-123 NUVC row2 snr=4.47
J04561-2159                NUVC row2 snr=3.07      QSO J1016+4706 NUVB row1 snr=4.32
LBQS 1340-0038             NUVC row2 snr=2.26      RX J1131.6+1557 NUVB row1 snr=0.59
Mrk 1179                   FUVB row1 snr=0.78      SDSS J132921.75+144431.2 NUVB row1 snr=1.29
UGC 12163                  FUVB row1 snr=1.04      US 1867       NUVB row1 snr=1.22
```

Note that this also degrades objects the audit calls `OK`: any multi-segment COS spectrum
in the sample is currently built from roughly half the photons available, C IV or not.

**Where it came from.** The `[0]` indexing was inherited verbatim from Trevor's
`read_cos()`, which is still in the file above the flat readers. Trevor's COS objects came
in as pre-built HSLA co-adds — a single already-merged array — so `read_cos()` was
essentially never exercised on raw two-segment x1d files and the bug never showed.
`read_stis_flat()` carries the same line, but STIS first-order x1d tables have one row, and
the one multi-row STIS mode present (E140M) has its own branch, so no STIS object is
affected by this bug.

---

## 3. Bug B — the continuum comes back as zero (`LOST_COADD`, 15 spectra)

`coadd.rebin()` normalises **each exposure separately** before co-adding:

```python
continuum = spec_morph.cont_filtered(new_binned_waves, new_binned_fluxes, z, Identifier)
new_binned_fluxes /= continuum
new_binned_errs   /= continuum
```

`cont_filtered()` builds that continuum by masking out the emission-line windows defined in
`emission_lines2()`, taking a 100 Å running median of everything left over, and interpolating
back across the masked gaps.

The C IV window is **rest 1460–1685 Å** — 225 Å wide. A single COS G160M exposure covers
about 190 Å. So for a low-redshift object a G160M exposure can sit *entirely inside the C IV
window*: for LEDA 3095694 (z = 0.069) the G160M exposures span rest 1475–1660 Å, wholly
enclosed. Every pixel is masked as emission, the running median has nothing to work with, and
`med_flux` is NaN from end to end.

`cont_filtered()` then hits this branch (added during the migration as a defensive fix):

```python
if starts_at_zero and ends_at_last:
    valid_flux = med_flux[~np.isnan(med_flux)]
    fill_value = np.nanmedian(valid_flux) if len(valid_flux) > 0 else 0.0
    med_flux[jp_list[i]:jf_list[i]+1] = fill_value
```

`valid_flux` is empty, so the continuum is **filled with 0.0**. Dividing by it turns the whole
exposure into ±inf. In `coadd.py` those get weight `1/inf² = 0`, the weighted median returns
NaN, the NaN-recovery branch finds `nansum(new_errs) == 0`, and the pixel is written out as
NaN with mask = 1. The exposure vanishes without a single warning.

`diagnostics/trace_exposures.py` shows it per exposure — `cont<=0` counts pixels where the
division destroys the flux, `surv` counts pixels that survive it:

```
LEDA 3095694 COS  z=0.0690  5 exposures  (C IV rest 1549.48 = 1656.4 A observed)
  i read (rest)         rebin (rest)            npix  cont<=0     surv  CIV?
  0   1206.4-  1348.3    1206.2-  1348.1      484        0      452
  1   1215.3-  1357.4    1215.1-  1357.8      483        0      449
  2   1475.6-  1649.4    1475.4-  1649.8      486      486        0  yes
  3   1475.7-  1649.4    1475.4-  1649.8      486      486        0  yes
  4   1486.4-  1660.3    1486.4-  1660.4      482      482        0  yes
```

Every exposure containing C IV is annihilated; the two that do not are fine.

**It is not only C IV.** Any exposure that falls entirely inside any window in
`emission_lines2()` is destroyed the same way. For 4C 11.72 (STIS) two exposures at rest
1199–1240 Å are lost inside the Lyα window (1180–1285 Å); for Ton S 210 two exposures at rest
2488–2527 Å are lost inside the small blue bump window (2220–2640 Å).

**Where it came from.** Before the defensive branch was added, this case raised `IndexError`
(the original code does `med_flux[jf_list[i]+1]` with `jf_list[i] == len(wave)-1`). The fix
stopped the crash but replaced a loud failure with a silent one. The equivalent function in
`ica/spec_morph.py` still has the original unguarded code; it operates on full co-added
spectra, which are always far wider than any single emission window, so it is not affected —
but see `CLAUDE.md` → Known Issues about keeping the two copies in step.

---

## 4. Bugs C and D — objects that crash and write nothing (`NO_OUTPUT`, 18 spectra)

Re-running `coadd.rebin()` on the 18 objects with no output file gives two error signatures.

### C. FOS exposures with only two readout groups — 13 objects

```
IndexError: index -3 is out of bounds for axis 0 with size 2
  at read_spec_data.py:575
```

`read_fos_flat()`'s `accum_flag()` decides whether a FOS file's readout groups are cumulative
by testing `mean(flux[-1]) > mean(flux[-2]) > mean(flux[-3])`. Files with `NAXIS2 = 2` have no
third group. One such file anywhere in an object's folder aborts the whole object, including
all its other exposures. Affected: **3C 273, AKN120, LEDA 29208, MR 2251-178, Mrk 1018,
PG1322+659, PG2308+098, PKS1055+20, PKS2216-03, QSO B1103-006, QSO B1626+5529, Ton 1565,
Ton 202 (FOS)**.

### D. One dead exposure aborts the object — 5 objects

```
IndexError: No good pixels found for NGC 3516 - lf3f01axq - all pixels failed basic quality cut
  at cut_edge_pix.py:77
```

`Cut_Edge_Pix` raises when an exposure has no pixel with non-zero flux, positive error and
positive wavelength. That is a reasonable thing to notice, but the exception propagates out of
the reader and discards every good exposure alongside the bad one. Affected: **2MASS
J14314604+0301444, NGC 3516, NGC 3783, Ton 202 (COS), Ton 469**.

This is the same family as the empty-`x1d` problem already handled in `read_cos_flat`
(Migration_Log, 2026-06-01) — MAST archives failed CalCOS products, and bulk downloads take
them. That fix caught files with a zero-row table; these files have rows, but the rows are
unusable.

---

## 5. Verification — the data is recoverable

`diagnostics/rebin_patched.py` monkey-patches `read_spec_data.read_data_flat` and
`spec_morph.cont_filtered`, then calls the unmodified `coadd.rebin`. **Nothing in
`rebinning/` was edited.** The four candidate fixes:

| fix | change |
|---|---|
| A | every detector row becomes its own entry in the exposure stack (segments are not contiguous or sorted, so they are kept separate rather than concatenated; `coadd.py` already handles independent overlapping exposures) |
| B | an exposure that fails `Cut_Edge_Pix` is skipped and reported, not fatal |
| C | `accum_flag` falls back to the single comparison it can make when there are only two groups |
| D | when an exposure lies entirely inside an emission window, the continuum becomes a constant — the exposure's own median flux — instead of 0.0 |

Fix D deserves scrutiny: it is the only one that changes what a number *means* rather than
which data is used. Normalising by a constant is safe here because the absolute scale of each
exposure is not preserved anyway — `coadd.py` morphs the co-added spectrum afterwards — but a
flat continuum is a placeholder, not a fit. If an exposure sits entirely inside C IV, there is
genuinely no line-free wavelength in it from which to measure a continuum.

Recovery on a sample, before vs after (`diagnostics/compare_rebins.py`):

```
object                       patched good (rest)   civ frac    original
4C 11.72_STIS                1199-1553              0.99       all-NaN
ESO 462-9_COS                1120-1746              0.92       1262-1746, civ 0.34
LEDA 3095694_COS             1120-1648              1.00       1206-1347, civ 0.00
Mrk 1298_COS                 1120-1656              1.00       1120-1388, civ 0.00
PG2308+098_FOS               1120-2284              1.00       no file
PKS2216-03_FOS               1120-2514              1.00       no file
```

### Full run over every flagged pair

All 79 flagged pairs (`LOST_ROWS`, `LOST_COADD`, `NO_OUTPUT`, plus `NO_CIV_SIGNAL` for
completeness) were re-rebinned into `RebinnedSpec_master_patched/`. **79 written, 0
failed.** Per-object results are in `diagnostics/civ_recovery_summary.csv`.

Of the 57 spectra with recoverable C IV:

| verdict | flagged | C IV core now >50% good | at the instrument's coverage limit | genuinely no C IV |
|---|---|---|---|---|
| `LOST_ROWS` | 24 | 17 | 7 | 0 |
| `LOST_COADD` | 15 | 14 | 1 | 0 |
| `NO_OUTPUT` | 18 | 16 | 0 | 2 |

The 8 "coverage limit" cases are not still broken — the pipeline is now extracting
everything the telescope recorded. They are COS NUV medium-resolution observations, where
each stripe is only ~35 Å wide, so 30–40% of the 1500–1600 Å window is all there is. The
patched good fraction matches the raw coverage to within a few percent in every case:

```
object                  patched good   raw coverage of 1500-1600   gratings
ICRF J154743.5+205216       0.33            0.36                   G130M G160M G185M
ICRF J213745.1-143255       0.33            0.39                   G130M G160M G185M
J04561-2159                 0.36            0.31                   G160M G185M G225M
QSO B2128-123               0.39            0.28                   G130M G160M G185M G225M
RX J1131.6+1557             0.34            0.40                   G130M G185M
US 1867                     0.30            0.33                   G130M G160M G185M G225M
Ton S 210 (STIS)            0.50            0.48                   E140M E230M
PG 1216+069                 0.31            1.00 (see below)       G130M G140L G160M G185M
```

`PG 1216+069` looks like an exception but is not: its nominal 100% coverage comes from a
G140L FUVA segment whose median flux/error across the window is **−0.01** — dead. The only
real data is a G185M NUVC stripe at rest 1528–1562 Å (S/N ≈ 2.2), which is 34% of the
window, and the patched run recovers 31% of it.

The two genuinely-no-C IV cases are `Mrk 1018` FOS (its FOS spectrum is rest 6013–8839 Å;
its C IV comes from its COS spectrum, which was already fine) and
`2MASS J14314604+0301444` COS (raw coverage of the window is 3%).

Fixes C and D also rescued 5 of the 22 `NO_CIV_SIGNAL` objects, which turned out to have
usable C IV after all. The other 17 remain correctly empty — mostly COS G140L snapshots
whose red throughput is gone.

### Regression check

Six objects the audit already calls `OK` (two per instrument):

```
object                        A good (rest)   B good (rest)   civ_A civ_B  B/A med scatter
2MASS J08265888+0611428_STIS  1120-2049       1120-2049        1.00  1.00    1.000   0.000
2MASS J11243917+4201450_FOS   1299-1886       1299-1886        1.00  1.00    1.000   0.000
2MASS J12253479-0247571_COS   1120-1666       1120-1665        1.00  1.00    0.979   0.370
Mrk 1018_COS                  1120-1760       1120-1760        1.00  1.00    0.997   0.151
Mrk 1502_FOS                  1120-3106       1120-3106        1.00  1.00    1.000   0.000
NGC 4593_STIS                 1128-10135      1128-10135       1.00  1.00    1.000   0.000
```

FOS and STIS come out bit-identical, as expected: fix C only fires on two-group FOS files and
fix A does nothing to single-row STIS tables. **COS does not**, and this is the important
caveat: the patched COS spectra now include the FUVB/NUVB/NUVC segments, so the co-add changes
wherever those overlap the segment that was already being used. The median ratio stays at
0.98–1.00 and the 51-pixel-smoothed ratio at 0.99, so there is no systematic shift — but
per-pixel scatter reaches 0.15–0.37 because there is genuinely more data going in. Adopting
fix A means **every COS spectrum in the sample changes**, including the ones that look fine
now, and all COS ICA fits would need re-running.

---

## 6. Files added

All new; nothing existing was modified.

| file | purpose |
|---|---|
| `diagnostics/audit_civ_coverage.py` | the census; writes `civ_coverage_audit.csv` |
| `diagnostics/civ_coverage_audit.csv` | 615 rows, one per (object, instrument) |
| `diagnostics/civ_recovery_summary.csv` | before/after C IV fraction and raw coverage for the 79 re-rebinned pairs |
| `diagnostics/patched_run.log` | log of the full patched run |
| `diagnostics/trace_exposures.py` | per-exposure trace for one object: where each exposure lands and whether it survives the continuum division |
| `diagnostics/rebin_patched.py` | the four candidate fixes, applied by monkey-patching; can re-rebin any subset |
| `diagnostics/compare_rebins.py` | side-by-side comparison of two rebin directories |
| `diagnostics/CIV_Data_Loss_Investigation.md` | this document |

```bash
# census
python -m diagnostics.audit_civ_coverage --raw-data-dir "F:/Richards-data/HST Paper/raw_data"

# why did one object lose its C IV?
python -m diagnostics.trace_exposures "LEDA 3095694" COS --raw-data-dir "F:/Richards-data/HST Paper/raw_data"

# re-rebin everything the audit flags, with the fixes applied
python -m diagnostics.rebin_patched --raw-data-dir "F:/Richards-data/HST Paper/raw_data" \
    --verdicts LOST_ROWS,LOST_COADD,NO_OUTPUT

# compare
python -m diagnostics.compare_rebins RebinnedSpec_master RebinnedSpec_master_patched
```

---

## 7. Open decisions

1. **Adopt fix A sample-wide, or only for the 24 objects that need it?** Sample-wide is the
   scientifically correct choice — using one of two detector segments is not defensible — but
   it changes every COS spectrum and invalidates the existing COS ICA fits. A partial
   application would leave the sample internally inconsistent.
2. **Is a flat continuum acceptable for exposures wholly inside an emission line (fix D)?**
   The alternative is to not normalise those exposures at all and let `coadd.py` handle the
   scale mismatch, or to fit the continuum on the co-added spectrum rather than per exposure.
   Worth raising with G. Richards.
3. **`NO_RAW_FILES` (21 pairs)** — a download-stage gap not investigated here.
4. **The `_COS_EXPOSURE_CAP` stopgap** for Mrk 817 and NGC 5548 (Migration_Log, 2026-06-17)
   is still in force and still arbitrary. It interacts with fix A: with all rows read, the cap
   of 100 files now yields ~200 exposure rows. NGC 5548 is in the `LOST_ROWS` list, so its
   recovered C IV is subject to that cap.
