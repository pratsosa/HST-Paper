# Originally: Trevor Code/ICA Scripts/plot_ICA.py
# Copied: 2026-04-29 (Phase 2 migration)
# Changes from source: header added only; no algorithm or import changes needed
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def spec_wreconstruction(spec_wave, spec_flux, spec_errs, spec_mask, rec_wave, rec_flux):
    good_flux = (~np.isnan(spec_flux) & ~np.isinf(spec_flux))
    good_errs = (~np.isnan(spec_errs) & ~np.isinf(spec_errs))
    ylims1 = np.percentile(spec_flux[good_flux],1), np.percentile(spec_flux[good_flux],99)
    #ylims1 = np.nanmin(flux[good_flux]), np.nanmax(flux[good_flux])
    ylims2 = np.percentile(spec_errs[good_errs],1), np.percentile(spec_errs[good_errs],99)

    #Plot spectrum with error specrtrum under
    fig, [ax1,ax2] = plt.subplots(2, 1, figsize=(16,8), sharex=True, gridspec_kw={'height_ratios': [4, 1]})
    fig.subplots_adjust(hspace=0)

    ax1.plot(spec_wave, spec_flux, color="k", label="Spectrum", zorder=1)
    ax1.scatter(spec_wave[spec_mask>0], spec_flux[spec_mask>0], color="r", label="Masked Pixels - found", zorder=2)
    #ax1.scatter(spec_wave[spec_mask==2], spec_flux[spec_mask==2], color="g", label="Masked Pixels - bordering found", zorder=2)
    ax1.plot(rec_wave, rec_flux, color="b", linewidth=2.5, label="Reconstruction", zorder=2)
    ax2.plot(spec_wave, spec_errs, color="k", label="Spectrum Errors", zorder=1)

    #red_expxlim = 1875
    ax1.set_xlim(min(spec_wave)-20, max(spec_wave)+20)
    #ax1.set_xlim(min(wave)-2, red_expxlim)
    ax1.set_ylim(ylims1[0], ylims1[1])
    #ax1.set_ylim(10.**np.nanpercentile(np.log10(spec_flux[good_flux], 1)), 10.**np.nanpercentile(np.log10(spec_flux[good_flux], 99)))
    ax2.set_ylim(ylims2)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ticks = np.array([1200,1500,1750,2000,2200,2500,3000,3500,4000,4500,5000,6000,7000,8000])
    ticks_use = ticks[(ticks>min(spec_wave))&(ticks<max(spec_wave))]
    #ticks_use = ticks[(ticks>min(wave))&(ticks<red_expxlim)]
    ax2.set_xticks(ticks_use)
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.get_xaxis().set_tick_params(which='minor', size=0)
    ax2.get_xaxis().set_tick_params(which='minor', width=0)
    plt.minorticks_off()
    ax2.set_xlabel("Wavelength (Å)", fontsize=17.5)
    ax1.set_ylabel("Flux (Arbitrary units)", fontsize=17.5)
    ax2.set_ylabel("Error", fontsize=17.5)
    ax1.legend(loc=1)
    #plt.suptitle(Identifier)
    #plt.show()
