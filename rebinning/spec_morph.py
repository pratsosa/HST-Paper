"""
rebinning/spec_morph.py
Migrated from: Trevor Code/spec_morph.py
Migration date: 2026-04-30

This is the REBINNING version of spec_morph. It differs from the ICA version
(Trevor Code/ICA Scripts/spec_morph.py), which will live in ica/spec_morph.py.
Use this module for the coadd/rebinning pipeline; use ica.spec_morph for ICA.

!! TWO-VERSION WARNING !!
There are two diverged copies of spec_morph in this repo:
  rebinning/spec_morph.py  ← this file (rebinning pipeline)
  ica/spec_morph.py        ← Phase 2 (ICA pipeline)
If you fix a bug or make a change here, check whether ica/spec_morph.py needs
the same fix. See Migration_Log.md → "Known Issues" for the full diff summary
and future consolidation options.

Import/path changes:
  - `import Small_Pix_Filter_TVM` → `from rebinning import small_pix_filter as Small_Pix_Filter_TVM`
  - `import BAL_regions` → `from rebinning import bal_regions as BAL_regions`
  - `np.load("qsomod_z250_ng_continuumSDSSrebin.npy")` → path relative to repo root
    via Data/components/qsomod_z250_ng_continuumSDSSrebin.npy
  - `np.load("/Users/Trevor1/...")` (vandenberk, used only by `morph()`) → path
    relative to repo root via Data/components/vandenberk01_continuum.npy.
    Note: `morph()` is NOT called by the rebinning pipeline (only `morph2()` is used).
No algorithmic changes.
"""

import matplotlib.pyplot as plt  #for testing
import numpy as np
from pathlib import Path

from rebinning import small_pix_filter as Small_Pix_Filter_TVM
from rebinning import bal_regions as BAL_regions

# Paths to reference data files (ship with the repo in Data/components/)
_DATA_COMPONENTS = Path(__file__).parent.parent / "Data" / "components"
_VANDENBERK_FILE = _DATA_COMPONENTS / "vandenberk01_continuum.npy"
_QSOMOD_FILE     = _DATA_COMPONENTS / "qsomod_z250_ng_continuumSDSSrebin.npy"


def morph(wave, flux, z, median, unique_markers=None, return_markers=False):
    #vb01 = np.load("/Volumes/MyPassportforMac/HST/AllData/vandenberk01_continuum.npy")
    #vb01 = np.load("/Volumes/Macintosh HD 1/Users/Trevor1/Desktop/HST/AllData/vandenberk01_continuum.npy")
    vb01 = np.load(str(_VANDENBERK_FILE))
    vb01_wave = vb01[:,0] * (1+z) #haven't de-redshifted wave yet
    vb01_flux = vb01[:,1] #* 10**(-17)

    mask_exp = (vb01_wave>=wave[0])&(vb01_wave<=wave[-1])
    wave_morph = vb01_wave[mask_exp]
    cont_morph = vb01_flux[mask_exp]

    #Now fit continuum of this spectrum
    if unique_markers is None:
        markers = np.array([1145,1290,1450,1690,1810,1985,2200,2660,2915,3550,4200,4720,5150,6125,6850,7000]) * (1+z) #these general markers better
    else:
        markers = unique_markers * (1+z)
    if min(wave)<10:
        print(wave)
    good_mask = ( (~np.isnan(flux)) & (flux!=0) )
    argstart = max(5, np.abs(wave-wave[good_mask][0]).argmin())
    argfin   = min(len(wave)-5, np.abs(wave-wave[good_mask][-1]).argmin())
    use_markers = markers[(markers>=wave[argstart])&(markers<=wave[argfin])]
    cont_region = {}

    for i in range(len(use_markers)-1):
        ind1 = abs(wave-use_markers[i]).argmin()
        ind2 = abs(wave-use_markers[i+1]).argmin()
        fit_region = np.polyfit([np.nanmedian(wave[(ind1-5):(ind1+5)]), np.nanmedian(wave[(ind2-5):(ind2+5)])],
                                [np.nanmedian(flux[(ind1-5):(ind1+5)]), np.nanmedian(flux[(ind2-5):(ind2+5)])],
                                 1) #do a linear fit on this section

        if len(use_markers)-2==0:
            x = wave
        elif i==0:
            x = wave[:ind2]
        elif i==(len(use_markers)-2):
            x = wave[ind1:]
        else:
            x = wave[ind1:ind2]
        y = fit_region[0]*x + fit_region[1]
        cont_region[i] = [x, y]

    cont_spec = cont_region[0][1]
    for i in range(1, len(use_markers)-1):
        cont_spec = np.append(cont_spec, cont_region[i][1])

    if len(cont_spec)-len(cont_morph)==1:
        cont_morph = np.append(cont_morph, cont_morph[-1]) #fix an annoying index error below, coming from definition of mask_exp

    cont_morph = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave, cont_morph, np.nanmin([len(wave)-1,300]))
    cont_spec = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave, cont_spec, np.nanmin([len(wave)-1,300]))
    morph_coeff = cont_morph / cont_spec

    if return_markers:
        return flux * morph_coeff, morph_coeff, use_markers/(1+z)
    return flux * morph_coeff, morph_coeff


"""
Changed morphing a bit on 4 August 2022;
now pads edges before fitting continuum;
means don't need to manually set markers for
individual exposures.
"""
def morph2(wave, flux, errs, z, Identifier):
    #Load reference spectrum
    qso = np.load(str(_QSOMOD_FILE))
    qso_wave = qso[:,0] * (1+z) #haven't de-redshifted wave yet
    qso_flux = qso[:,1] #* 10**(-17)

    expStart = np.argmin(np.abs(qso_wave-min(wave)))
    wave_morph = qso_wave[expStart:expStart+len(wave)]
    cont_morph = qso_flux[expStart:expStart+len(wave)]

    cont_spec = cont_filtered(wave, flux, z, Identifier)
    cont_spec = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave, cont_spec, np.nanmin([len(wave)-1, 300]))
    cont_morph = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave, cont_morph, np.nanmin([len(wave)-1, 300]))
    morph_coeff = cont_morph / cont_spec

    return flux * morph_coeff, morph_coeff


def pad_spec(wave, flux, errs, z):
    #return original spec with padded edges

    filterval = len(wave)//3  #number of pixels to pad onto each side - large fraction for smaller exposures
    sdss_c1 = 0.0001 #logarithmic wavelength scale of SDSS

    #pad wave array
    loglam = np.log10(wave)
    loglam_leftpad  = np.flip(np.arange(min(loglam), min(loglam)-filterval*sdss_c1-sdss_c1, -sdss_c1))
    loglam_rightpad = np.arange(max(loglam), max(loglam)+filterval*sdss_c1+sdss_c1, sdss_c1)
    loglam_wpad = np.concatenate((loglam_leftpad, loglam, loglam_rightpad))
    wave_wpad = 10.**loglam_wpad

    #"initialize" padded flux+error array
    flux_wpad = np.concatenate((np.zeros(len(loglam_leftpad)), flux, np.zeros(len(loglam_rightpad))))
    errs_wpad = np.concatenate((np.zeros(len(loglam_leftpad)), errs, np.zeros(len(loglam_rightpad))))

    em = emission_lines(wave/(1+z))
    em_wpad = emission_lines(wave_wpad/(1+z))

    #pad onto edges - should ignore emission line regions
    for k in range(filterval):
        #reflect fluxes
        df_left   = flux_wpad[~em_wpad][filterval+k] - np.nanmedian(flux[~em][:filterval//6+1])
        df_right  = flux_wpad[~em_wpad][-filterval-k] - np.nanmedian(flux[~em][-filterval//6:])
        flux_wpad[filterval-k-1] = -2*df_left + flux_wpad[~em_wpad][filterval+k]
        flux_wpad[-filterval+k]  = -2*df_right + flux_wpad[~em_wpad][-filterval-k]
        #same for errors
        de_left   = errs_wpad[~em_wpad][filterval+k] - np.nanmedian(errs[~em][:filterval//6+1])
        de_right  = errs_wpad[~em_wpad][-filterval-k] - np.nanmedian(errs[~em][-filterval//6:])
        errs_wpad[filterval-k-1] = -2*de_left + errs_wpad[~em_wpad][filterval+k]
        errs_wpad[-filterval+k]  = -2*de_right + errs_wpad[~em_wpad][-filterval-k]

    #padded well; now apply ~601-pix median filter
    median_flux = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave_wpad, flux_wpad, len(wave)//6)
    median_errs = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave_wpad, errs_wpad, len(wave)//6)
    for i in range(len(flux_wpad)):
        if abs(flux_wpad[i]-median_flux[i]) > 5.*median_errs[i]:
            flux_wpad[i] = median_flux[i]

    return wave_wpad, flux_wpad


def continuum_fit(wave, flux, z):
    #fit continuum to a spectrum; assumes SDSS wavelength resolution

    #will interpolate across emission lines using these markers
    markers = np.array([1145,1290,1450,1690,1810,1985,2200,2660,2915,3550,4200,4720,5150,6125,6850,7000]) * (1+z) #these general markers better
    good_mask = ( (~np.isnan(flux)) & (flux!=0) )

    #only use the markers within spec wavelength range
    argstart = max(5, np.abs(wave-wave[good_mask][0]).argmin())
    argfin   = min(len(wave)-5, np.abs(wave-wave[good_mask][-1]).argmin())
    use_markers = markers[(markers>=wave[argstart])&(markers<=wave[argfin])]

    #save [wave,flux] for continuum in dictionary
    cont_region = {}
    for i in range(len(use_markers)-1):
        #interpolate between markers on spectrum
        ind1 = abs(wave-use_markers[i]).argmin()
        ind2 = abs(wave-use_markers[i+1]).argmin()
        fit_region = np.polyfit([np.nanmedian(wave[(ind1-5):(ind1+5)]), np.nanmedian(wave[(ind2-5):(ind2+5)])],
                                [np.nanmedian(flux[(ind1-5):(ind1+5)]), np.nanmedian(flux[(ind2-5):(ind2+5)])],
                                 1) #do a linear fit on this section

        if len(use_markers)<=1:
            print("Wave coverage not big enough; <= 1 marker")

        elif len(use_markers)==2:
            x = wave

        elif i==0:
            x = wave[:ind2]

        elif i==(len(use_markers)-2):
            x = wave[ind1:]

        else:
            x = wave[ind1:ind2]
        y = fit_region[0]*x + fit_region[1]
        cont_region[i] = [x, y]

    cont_spec = cont_region[0][1]
    for i in range(1, len(use_markers)-1):
        cont_spec = np.append(cont_spec, cont_region[i][1])

    cont_spec = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave, cont_spec, np.nanmin([len(wave)-1, 300]))
    return cont_spec


def continuum_fit2(wave, flux, z):
    #Taken from default SDSS continuum fitting windows - should be better than individual markers in v1
    window_all = np.array(
            [[980,1000], [1150., 1170.], [1275., 1290.], [1350., 1360.], [1435., 1465.], [1695, 1705], [1730., 1810.],
             [1970., 2400.], [2480., 2700.], [2925., 3400.], [3575., 3832.], [4000., 4050.], [4200., 4230.],
             [4435., 4640.], [5100., 5535.], [6005., 6035.], [6110., 6250.], [6800., 7000.], [7160., 7180.],
             [7500., 7800.], [8050., 8150.]])
    rfwave = wave / (1+z)

    #locate which windows to use
    ind_use = []
    for i,window in enumerate(window_all):
        if ( (rfwave>window[0]).sum() > 0 ) and \
           ( (rfwave<window[1]).sum() > 0 ):
            ind_use.append(i)
    ind_use = np.array(ind_use)

    #fit the spectrum piecewise
    use_windows = [] #these can extend beyond actual window to fit whole spec
    cont_windows = []
    for n,i in enumerate(ind_use):
        ifit = (rfwave>window_all[i][0]) & (rfwave<window_all[i][1]) & (~np.isnan(flux))
        m, b = np.polyfit(rfwave[ifit], flux[ifit], 1)

        if len(ind_use)==1:
            ixfit = rfwave
            iyfit = m*ixfit + b
        elif n==0:
            #if first fit, use out to left edge
            ixfit = rfwave[0:np.where(ifit)[0][-1]] #exclude right edge pixel
        elif n==len(ind_use)-1:
            #if final fit, extend out to right edge
            ixfit = rfwave[np.where(ifit)[0][0]:]
        else:
            #if central fit, start at window's left edge,
            #out to window's right edge exclusive
            ixfit = rfwave[np.where(ifit)[0][0]:np.where(ifit)[0][-1]]
        iyfit = m*ixfit + b

        #now need to connect edges of individual window fits
        use_windows.append(ixfit)
        cont_windows.append(iyfit)

    #fill in the gaps between window fits
    use_windows = np.array(use_windows, dtype=object)
    cont_windows = np.array(cont_windows, dtype=object)
    cont_final_windows = []
    if len(cont_windows) > 1:
        #then interpolate between gaps
        for k in range(len(cont_windows)-1):
            cont_final_windows.append(cont_windows[k])
            #now fit from right edge of this window fit
            #to just before left edge of next window
            indstart = np.where(rfwave>use_windows[k][-1])[0][0]
            indend   = np.where(rfwave<use_windows[k+1][0])[0][-1]

            x1, y1 = use_windows[k][-1], cont_windows[k][-1]
            x2, y2 = use_windows[k+1][0], cont_windows[k+1][0]
            mgap, bgap = np.polyfit([x1,x2], [y1,y2], 1)
            xgap = rfwave[indstart:indend+1]
            ygap = mgap*xgap + bgap
            cont_final_windows.append(ygap)

        cont_final_windows.append(cont_windows[-1])
        continuum = cont_final_windows[0]
        for i in range(1, len(cont_final_windows)):
            continuum = np.append(continuum, cont_final_windows[i])
    else:
        continuum = np.array(cont_windows[0], dtype=float)

    continuum = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave, continuum, np.nanmin([len(wave)-1, 300]))
    return continuum


def emission_lines(wavelength):
    LyA     = ( (wavelength>1180) & (wavelength<1270) )
    SiIV    = ( (wavelength>1370) & (wavelength<1430) )
    CIV     = ( (wavelength>1460) & (wavelength<1600) )
    CIII    = ( (wavelength>1820) & (wavelength<1935) )
    MgII    = ( (wavelength>2720) & (wavelength<2835) )
    HG      = ( (wavelength>4320) & (wavelength<4370) )
    HB_OIII = ( (wavelength>4835) & (wavelength<5030) )
    HA      = ( (wavelength>6530) & (wavelength<6600) )
    return ( LyA | CIV | CIII | MgII | HB_OIII | HA )


def emission_lines2(wavelength):
    OVI     = ( (wavelength>1010) & (wavelength<1050) )
    LyA     = ( (wavelength>1180) & (wavelength<1285) )
    OI      = ( (wavelength>1295) & (wavelength<1320) )
    CII     = ( (wavelength>1325) & (wavelength<1350) )
    SiIV    = ( (wavelength>1360) & (wavelength<1440) )
    CIV     = ( (wavelength>1460) & (wavelength<1685) )
    CIII    = ( (wavelength>1820) & (wavelength<1960) ) #contains CII
    SBB     = ( (wavelength>2220) & (wavelength<2640) )
    MgII    = ( (wavelength>2695) & (wavelength<2900) )

    NeV     = ( (wavelength>3325) & (wavelength<3365) )
    NeVI    = ( (wavelength>3390) & (wavelength<3455) )
    OII     = ( (wavelength>3700) & (wavelength<3755) )
    HeI     = ( (wavelength>3835) & (wavelength<3925) )
    HD      = ( (wavelength>4050) & (wavelength<4150) )
    HG      = ( (wavelength>4270) & (wavelength<4415) )
    HB_OIII = ( (wavelength>4770) & (wavelength<5050) )
    HA      = ( (wavelength>6410) & (wavelength<6790) )

    return ( OVI | LyA | OI | CII | SiIV | CIV | CIII | SBB | MgII |
           NeV | NeVI | OII | HeI | HD | HG | HB_OIII | HA )


def cont_filtered(wave, flux, z, Identifier):
    # AP: Fixed to handle edge cases where emission lines extend to array boundaries
    #Ignore BAL troughs
    bal = np.zeros(len(wave))
    if Identifier in BAL_regions.Troughs:
        bal_region = BAL_regions.Troughs[Identifier]
        bal[((wave/(1+z))>bal_region[0])&((wave/(1+z))<bal_region[1])] = 1
    #compute median of spectrum outside emission line regions
    em = emission_lines2(wave/(1+z))
    med_flux = np.zeros(len(wave))
    for k in range(len(wave)):
        if em[k] or bal[k]:
            med_flux[k] = np.nan
        else:
            waverange = np.abs(wave-wave[k])<100 #use pixels within 100Å
            med_flux[k] = np.nanmedian(flux[(~em)&(waverange)])

    #locate indices where nan (i.e. emission line) appears
    jp_list = []
    jf_list = []
    myiter = iter(range(len(wave)))
    for j in myiter:
        if np.isnan(med_flux[j]):
            jp = j #index where branch of nans starts
            while np.isnan(med_flux[j]):
                j += 1
                next(myiter, None)
                if j >= len(wave): break
            jf = j-1 #index where branch of nans ends
            jp_list.append(jp)
            jf_list.append(jf)

    #interpolate across emission-line regions

    # AP: Original code below - replaced with defensive version to handle edge cases
    # for i in range(len(jp_list)):
    #     if jp_list[i]==0:
    #         med_flux[jp_list[i]:jf_list[i]+1] = med_flux[jf_list[i]+1]
    #     elif jf_list[i]==len(wave)-1:
    #         med_flux[jp_list[i]:jf_list[i]+1] = med_flux[jp_list[i]-1]
    #     else:
    #         x = np.array([wave[jp_list[i]-1], wave[jf_list[i]+1]])
    #         y = np.array([med_flux[jp_list[i]-1], med_flux[jf_list[i]+1]])
    #         m, b = np.polyfit(x, y, 1)
    #         w    = wave[jp_list[i]:jf_list[i]+1]
    #         med_flux[jp_list[i]:jf_list[i]+1] = m*w + b
    # End original code

    # AP: Defensive version - adds bounds checking to prevent IndexError
    for i in range(len(jp_list)):
        # Check if emission line region extends to edges (prevents index out of bounds)
        starts_at_zero = (jp_list[i] == 0)
        ends_at_last = (jf_list[i] == len(wave) - 1)

        if starts_at_zero and ends_at_last:
            # Edge case: emission spans entire spectrum - use median of any valid flux
            valid_flux = med_flux[~np.isnan(med_flux)]
            fill_value = np.nanmedian(valid_flux) if len(valid_flux) > 0 else 0.0
            med_flux[jp_list[i]:jf_list[i]+1] = fill_value
        elif starts_at_zero:
            # Emission starts at first pixel - use value after emission region
            med_flux[jp_list[i]:jf_list[i]+1] = med_flux[jf_list[i]+1]
        elif ends_at_last:
            # Emission extends to last pixel - use value before emission region
            med_flux[jp_list[i]:jf_list[i]+1] = med_flux[jp_list[i]-1]
        else:
            # Normal case: emission in middle - interpolate between edges
            x = np.array([wave[jp_list[i]-1], wave[jf_list[i]+1]])
            y = np.array([med_flux[jp_list[i]-1], med_flux[jf_list[i]+1]])
            m, b = np.polyfit(x, y, 1)
            w = wave[jp_list[i]:jf_list[i]+1]
            med_flux[jp_list[i]:jf_list[i]+1] = m*w + b
    # End AP

    return med_flux


def cont_filtered_AP(wave, flux, z, Identifier):
    """
    Improved version of cont_filtered with defensive fixes for edge cases.
    Fits continuum using local 100Å median windows and interpolates across emission lines.
    """
    #Ignore BAL troughs
    bal = np.zeros(len(wave))
    if Identifier in BAL_regions.Troughs:
        bal_region = BAL_regions.Troughs[Identifier]
        bal[((wave/(1+z))>bal_region[0])&((wave/(1+z))<bal_region[1])] = 1
    #compute median of spectrum outside emission line regions
    em = emission_lines2(wave/(1+z))
    med_flux = np.zeros(len(wave))
    for k in range(len(wave)):
        if em[k] or bal[k]:
            med_flux[k] = np.nan
        else:
            waverange = np.abs(wave-wave[k])<100 #use pixels within 100Å
            med_flux[k] = np.nanmedian(flux[(~em)&(waverange)])

    #locate indices where nan (i.e. emission line) appears
    jp_list = []
    jf_list = []
    myiter = iter(range(len(wave)))
    for j in myiter:
        if np.isnan(med_flux[j]):
            jp = j #index where branch of nans starts
            # Defensive fix: ensure we don't go beyond array bounds
            while j < len(wave) and np.isnan(med_flux[j]):
                j += 1
                next(myiter, None)
            jf = j-1 #index where branch of nans ends
            jp_list.append(jp)
            jf_list.append(jf)

    #interpolate across emission-line regions
    for i in range(len(jp_list)):
        # Defensive fix: add bounds checking to prevent index errors
        if jp_list[i]==0 or jf_list[i] >= len(wave)-1:
            #if first or extends to last pixel, use nearest valid value
            if jf_list[i] < len(wave)-1:
                med_flux[jp_list[i]:jf_list[i]+1] = med_flux[jf_list[i]+1]
            elif jp_list[i] > 0:
                med_flux[jp_list[i]:jf_list[i]+1] = med_flux[jp_list[i]-1]
            else:
                # Edge case: entire spectrum is NaN - use median of non-NaN values if any exist
                valid_flux = med_flux[~np.isnan(med_flux)]
                if len(valid_flux) > 0:
                    med_flux[jp_list[i]:jf_list[i]+1] = np.nanmedian(valid_flux)
        else:
            #if middle, interpolate
            x = np.array([wave[jp_list[i]-1], wave[jf_list[i]+1]])
            y = np.array([med_flux[jp_list[i]-1], med_flux[jf_list[i]+1]])
            m, b = np.polyfit(x, y, 1)
            w    = wave[jp_list[i]:jf_list[i]+1]
            med_flux[jp_list[i]:jf_list[i]+1] = m*w + b

    return med_flux
