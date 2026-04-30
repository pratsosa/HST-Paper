"""
rebinning/cut_edge_pix.py
Migrated from: Trevor Code/Cut_Edge_Pix_TVM_NoDQ.py
Migration date: 2026-04-30

Import changes only:
  - Instrument-specific SpecCuts imports updated to use package-relative imports
    (e.g., `from rebinning import spec_cuts_fos as SpecCuts`)
No algorithmic changes.

Original file header:
  Cut_Edge_Pix_TVM_NoDQ.py
  AP: 2026-03-08
  This is Cut_Edge_Pix_TVM_DQOnly.py with the DQ==0 condition removed entirely.
  Used by the 2026AP run to test the impact of DQ cuts on the final rebinned spectra.

  Version lineage:
    Python 2 era (Cut_Edge_Pix.py, Cut_Edge_Pix0.py): DQ-walking loops, no SNR condition.
    _TVM Aug 2022 (hypothesized): good_mask with & (DQ==0)
    _TVM Jan 18 2023: replaced (DQ==0) with (all_snrs>1.5)
    PreSNR (AP 2026-03-07): removed SNR condition but did NOT restore DQ==0 — no edge quality cut at all
    DQOnly (AP 2026-03-07): restores & (DQ==0), the hypothesized Aug 2022 condition
    NoDQ   (AP 2026-03-08): removes DQ==0 entirely — no DQ quality cut applied

  Purpose: quantify the effect of DQ cuts by comparing RebinnedSpec_2026AP/ against
  RebinnedSpec_2026E/ (DQ==0 applied). All other settings (FOS error formula, edge
  wavelength cuts) are identical to 2026E.

  See Rebin_Implementation_Log.md for full rationale.
"""

import numpy as np


def maskEdge(wave, flux, errs, mask):
    '''
    mask the edges of (mainly SDSS) spectrum so we ignore it in the ICA fitting
    '''
    s2n = flux / errs
    good_mask = ( (flux!=0) & (~np.isnan(flux)) & (wave>0) & (errs!=0) )
    istart = (np.abs(wave - wave[good_mask][0])).argmin()
    iend   = (np.abs(wave - wave[good_mask][-1])).argmin()
    #while s2n[iend-30:iend+30]


def Cut_Edge_Pix(DQ, wavelength, flux, fluxerr, input_data, array_len, data_type, Cut_Spikes, z=0, exp_id=None, Instrument=None):
    if Instrument == "FOS":
        from rebinning import spec_cuts_fos as SpecCuts
    elif Instrument == "STIS":
        from rebinning import spec_cuts_stis as SpecCuts
    elif Instrument == "COS":
        from rebinning import spec_cuts_cos as SpecCuts
    else:
        from rebinning import spec_cuts_hsla as SpecCuts

    all_snrs = np.zeros(len(flux))
    for i in range(len(fluxerr)):
        all_snrs[i] = np.nan if fluxerr[i] == 0 else flux[i]/fluxerr[i]

    if exp_id in SpecCuts.BlueEdges and exp_id in SpecCuts.RedEdges:
        wavecut_blue = SpecCuts.BlueEdges[exp_id] * (1+z)
        wavecut_red  = SpecCuts.RedEdges[exp_id] * (1+z)
        good_mask = ( (flux!=0) & (~np.isnan(flux)) & (wavelength>0) & (fluxerr!=0) &
                        (wavelength>wavecut_blue) & (wavelength<wavecut_red) )  # AP: no DQ cut (2026AP)
    elif exp_id in SpecCuts.BlueEdges:
        wavecut_blue = SpecCuts.BlueEdges[exp_id] * (1+z)
        good_mask = ( (flux!=0) & (~np.isnan(flux)) & (wavelength>0) & (fluxerr!=0) & (wavelength>wavecut_blue) )  # AP: no DQ cut
    elif exp_id in SpecCuts.RedEdges:
        wavecut_red = SpecCuts.RedEdges[exp_id] * (1+z)
        good_mask = ( (flux!=0) & (~np.isnan(flux)) & (wavelength>0) & (fluxerr!=0) & (wavelength<wavecut_red) )  # AP: no DQ cut
    else:
        good_mask = ( (flux!=0) & (~np.isnan(flux)) & (wavelength>0) & (fluxerr!=0) )  # AP: no DQ cut

    #good_mask = ( (flux!=0) & (~np.isnan(flux)) & (wavelength>0) & (fluxerr!=0) & (all_snrs>2.5) ) #& (DQ==0)

    if good_mask.sum() == 0:
        raise IndexError(f"No good pixels found for {exp_id} - all pixels failed basic quality cut")

    begin_stop_index = (np.abs(wavelength - wavelength[good_mask][0])).argmin()
    end_stop_index = (np.abs(wavelength - wavelength[good_mask][-1])).argmin()
    #FIXME: maybe change so that most pixels have S2N > 2 or something?
    while (all_snrs[begin_stop_index:begin_stop_index+65]<=0).sum()>20:
        begin_stop_index += 1
    while (all_snrs[end_stop_index-65:end_stop_index]<=0).sum()>20:
        end_stop_index -= 1

    if data_type == "wavelength":
        trunc_data = np.array(list(np.zeros(begin_stop_index)) + list(input_data[begin_stop_index:]))
        return trunc_data

    if ~Cut_Spikes:
        #Save SNRs for each pixel
        all_snrs = np.zeros(len(flux))
        for i in range(len(fluxerr)):
            all_snrs[i] = np.nan if fluxerr[i] == 0 else flux[i]/fluxerr[i]

        if data_type == "flux":
            trunc_data = np.concatenate((np.nan*np.ones(begin_stop_index), input_data[begin_stop_index:]))
            trunc_data[end_stop_index:] = np.nan if end_stop_index is not None else trunc_data[end_stop_index:]

        elif data_type == "flux error":
            trunc_data = np.concatenate((np.zeros(begin_stop_index), input_data[begin_stop_index:]))
            trunc_data[end_stop_index:] = 0. if end_stop_index is not None else trunc_data[end_stop_index:]

        if data_type == "masks":
            trunc_data = np.concatenate((np.zeros(begin_stop_index), input_data[begin_stop_index:]))
            trunc_data[end_stop_index:] = 0. if end_stop_index is not None else trunc_data[end_stop_index:]

        return trunc_data
