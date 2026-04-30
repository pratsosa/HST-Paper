# Originally: Trevor Code/ICA Scripts/spec_morph.py
# Copied: 2026-04-29 (Phase 2 migration)
# Changes from source:
#   - Removed sys.path.append("/Users/Trevor1/Dropbox/HST/HSTCode/")
#   - import Small_Pix_Filter_TVM -> from rebinning import small_pix_filter as Small_Pix_Filter_TVM
#   - import BAL_regions -> from rebinning import bal_regions as BAL_regions
#   - np.load('components/qsomod_...') -> np.load(str(_QSOMOD_FILE)) [both morph() and morph2()]
# NOTE: This is the ICA version of spec_morph.py.  It differs from rebinning/spec_morph.py:
#   - continuum_fit() uses np.polynomial.polynomial.Polynomial.fit (not np.polyfit)
#   - continuum_fit() takes extra 'value' and 'markers_arr' parameters
#   - cont_filtered() here does NOT have AP's bounds-checking fixes present in rebinning/spec_morph.py
#   See Migration_Log.md -> Known Issues for details.
"""
Keep morphing code in it's own file.
"""

from pathlib import Path
import matplotlib.pyplot as plt #for testing
import numpy as np
from rebinning import small_pix_filter as Small_Pix_Filter_TVM
from rebinning import bal_regions as BAL_regions

_QSOMOD_FILE = Path(__file__).resolve().parent.parent / "Data" / "components" / "qsomod_z250_ng_continuumSDSSrebin.npy"

def morph(wave, flux, z, median, unique_markers=None, return_markers=False):
    #vb01 = np.load("/Volumes/MyPassportforMac/HST/AllData/vandenberk01_continuum.npy")
    #vb01 = np.load("/Volumes/Macintosh HD 1/Users/Trevor1/Desktop/HST/AllData/vandenberk01_continuum.npy")
    #vb01 = np.load("/Users/Trevor1/Dropbox/HST/ICA/vandenberk01_continuum.npy")
    vb01 = np.load(str(_QSOMOD_FILE))

    vb01_wave = vb01[:,0] * (1+z) #haven't de-redshifted wave yet
    vb01_flux = vb01[:,1] #* 10**(-17)

    mask_exp = (vb01_wave>=wave[0])&(vb01_wave<=wave[-1])
    wave_morph = vb01_wave[mask_exp]
    cont_morph = vb01_flux[mask_exp]

    #Now fit continuum of this spectrum
    #markers = np.array([vb01_wave[5],950,1000,1145,1285,1322,1355,1450,1490,1590,1695,1740,1755,1800,1830,1975,2030,2160,2200,2235,2300,2500,2660,2900,3075,3550,3802,3934,4200,4440,4760,4793,5050,5085,5405,6000,6200,6350,6800,7000,7600,8000,8300,8700,vb01_wave[-5]]) * (1+z)
    if unique_markers is None:
        markers = np.array([1145,1290,1450,1690,1810,1985,2200,2660,2915,3550,4200,4720,5150,6125,6850,7000]) * (1+z) #these general markers better
    else:
        markers = unique_markers * (1+z)
    if min(wave)<10:
        print(wave)
    #print(wave, markers)
    good_mask = ( (~np.isnan(flux)) & (flux!=0) )
    argstart = max(5, np.abs(wave-wave[good_mask][0]).argmin())
    argfin   = min(len(wave)-5, np.abs(wave-wave[good_mask][-1]).argmin())
    use_markers = markers[(markers>=wave[argstart])&(markers<=wave[argfin])]
    #use_markers = np.append(use_markers, wave[-5])
    cont_region = {}

    for i in range(len(use_markers)-1):
        ind1 = abs(wave-use_markers[i]).argmin()
        ind2 = abs(wave-use_markers[i+1]).argmin()
        #print([wave[(ind1-5):(ind1+5)], wave[(ind2-5):(ind2+5)]],
        #                        [flux[(ind1-5):(ind1+5)], flux[(ind2-5):(ind2+5)]])
        #print(f"ind1: {ind1}, ind2: {ind2}.") - AP
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
    morph_coeff = cont_morph / cont_spec  #/ np.nanmedian(cont_morph) #* (np.nanmedian(cont_spec) / np.nanmedian(cont_morph))

    if return_markers:
        return flux * morph_coeff, morph_coeff, use_markers/(1+z)
    return flux * morph_coeff, morph_coeff #np.nanmedian(morph_coeff)


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

    #Now fit continuum of this spectrum
    """
    #markers = np.array([qso_wave[5],950,1000,1145,1285,1322,1355,1450,1490,1590,1695,1740,1755,1800,1830,1975,2030,2160,2200,2235,2300,2500,2660,2900,3075,3550,3802,3934,4200,4440,4760,4793,5050,5085,5405,6000,6200,6350,6800,7000,7600,8000,8300,8700,qso_wave[-5]]) * (1+z)
    if unique_markers is None:
        markers = np.array([1145,1290,1450,1690,1810,1985,2200,2660,2915,3550,4200,4720,5150,6125,6850,7000]) * (1+z) #these general markers better
    else:
        markers = unique_markers * (1+z)
    """

    #pad spectrum
    #wave_wpad, flux_wpad = pad_spec(wave, flux, errs, z)
    #fit continuum on padded spectrum
    #cont_spec_padded = cont_filtered(wave_wpad, flux_wpad, z)
    #cont_spec_padded = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave_wpad, cont_spec_padded, np.nanmin([len(wave_wpad)-1, 300]))
    #cut to only original spectrum; only care about shape
    #nopadStart = np.argmin(np.abs(min(wave)-wave_wpad))
    #cont_spec = cont_spec_padded[nopadStart:nopadStart+len(wave)]
    cont_spec = cont_filtered(wave, flux, z, Identifier)
    #cont_spec = continuum_fit2(wave, flux, z)
    cont_spec = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave, cont_spec, np.nanmin([len(wave)-1, 300]))
    cont_morph = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave, cont_morph, np.nanmin([len(wave)-1, 300]))
    #cont_spec = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave, cont_spec, 300)
    morph_coeff = cont_morph / cont_spec  #/ np.nanmedian(cont_morph) #* (np.nanmedian(cont_spec) / np.nanmedian(cont_morph))

    """
    fig = plt.figure(figsize=(8,6))
    #plt.plot(wave/(1+z), flux/np.nanmedian(flux), "-b", alpha=0.5)
    plt.plot(wave/(1+z), flux, "-b")
    plt.plot(wave/(1+z), cont_spec, "-r")
    #plt.plot(wave/(1+z), cont_spec/np.nanmedian(flux), "-b", alpha=0.5)
    #plt.plot(wave_wpad/(1+z), flux_wpad/np.nanmedian(flux), c="purple")
    #plt.plot(wave/(1+z), flux*morph_coeff / np.nanmedian(flux*morph_coeff))
    #plt.xlim(2100,3100)
    plt.ylim(0,4)
    plt.show()
    print(min(wave)/(1+z))
    """
    return flux * morph_coeff, morph_coeff #np.nanmedian(morph_coeff)

def pad_spec(wave, flux, errs, z):
    #return original spec with padded edges

    filterval = len(wave)//3  #number of pixels to pad onto each side - large fraction for smaller exposures
    sdss_c1 = 0.0001 #logarithmic wavelength scale of SDSS

    #pad wave array
    loglam = np.log10(wave)
    #loglam_wpad = np.arange(np.nanmin(loglam)-filterval*sdss_c1,
    #                        np.nanmax(loglam)+filterval*sdss_c1+sdss_c1,
    #                        sdss_c1)
    loglam_leftpad  = np.flip(np.arange(min(loglam), min(loglam)-filterval*sdss_c1-sdss_c1, -sdss_c1))
    loglam_rightpad = np.arange(max(loglam), max(loglam)+filterval*sdss_c1+sdss_c1, sdss_c1)
    loglam_wpad = np.concatenate((loglam_leftpad, loglam, loglam_rightpad))
    wave_wpad = 10.**loglam_wpad

    #"initialize" padded flux+error array
    flux_wpad = np.concatenate((np.zeros(len(loglam_leftpad)), flux, np.zeros(len(loglam_rightpad))))
    errs_wpad = np.concatenate((np.zeros(len(loglam_leftpad)), errs, np.zeros(len(loglam_rightpad))))

    #subtract emission lines from orig spec for padding? (maybe no)

    #print(len(loglam_wpad), len(wave_wpad), len(flux_wpad), filterval)
    em = emission_lines(wave/(1+z))
    em_wpad = emission_lines(wave_wpad/(1+z))

    #pad onto edges - should ignore emission line regions
    #without messing up reflection somehow
    for k in range(filterval):
        #reflect fluxes
        df_left   = flux_wpad[~em_wpad][filterval+k] - np.nanmedian(flux[~em][:filterval//6+1])
        df_right  = flux_wpad[~em_wpad][-filterval-k] - np.nanmedian(flux[~em][-filterval//6:])
        flux_wpad[filterval-k-1] = -2*df_left + flux_wpad[~em_wpad][filterval+k] #change to median here?
        flux_wpad[-filterval+k]  = -2*df_right + flux_wpad[~em_wpad][-filterval-k]
        #same for errors - maybe turn into function
        de_left   = errs_wpad[~em_wpad][filterval+k] - np.nanmedian(errs[~em][:filterval//6+1])
        de_right  = errs_wpad[~em_wpad][-filterval-k] - np.nanmedian(errs[~em][-filterval//6:])
        errs_wpad[filterval-k-1] = -2*de_left + errs_wpad[~em_wpad][filterval+k] #change to median here?
        errs_wpad[-filterval+k]  = -2*de_right + errs_wpad[~em_wpad][-filterval-k]

    #padded well; now apply ~601-pix median filter
    median_flux = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave_wpad, flux_wpad, len(wave)//6)
    median_errs = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave_wpad, errs_wpad, len(wave)//6)
    for i in range(len(flux_wpad)):
        if abs(flux_wpad[i]-median_flux[i]) > 5.*median_errs[i]:
            flux_wpad[i] = median_flux[i]


    #plt.plot(wave_wpad, flux_wpad)
    #plt.plot(wave_wpad, median_flux)
    #plt.plot(wave_wpad, errs_wpad)
    #plt.ylim(0,18)
    #plt.xlim(None,5000)
    return wave_wpad, flux_wpad

def continuum_fit(wave, flux, z, value, markers_arr):
    #fit continuum to a spectrum; assumes SDSS wavelength resolution

    #will interpolate across emission lines using these markers
    markers = np.array([1145,1290,1450,1690,1810,1985,2200,2660,2915,3550,4200,4720,5150,6125,6850,7000]) * (1+z) #these general markers better

    #Above is Trevor's original markers array - AP


    markers = np.array([1145,1200, 1290,1350, 1450,1550, 1600, 1690,1810,1880, 1985,2100, 2200,2300, 2660,2915,3250, 3550, 3700, 4200,4720,5150,6125,6850,7000]) * (1+z)

    #Here is a markers array I was testing out (not the one I ended up using), just added more markers more or less between the previous ones - AP

    # I ended up using the following array as the parameter: [1145,1290,1450,1690,1810,1985,2200,2400,2660,2915,3550,4200,4720,5150,6125,6850,7000]
    # Note that I only added the value 2400 as this was where the green spectrum looked weird.
    markers = np.array(markers_arr) * (1+z)

    #This just allowed me to test specific markers arrays more easily by passing them into the function through the jupyter notebook - AP

    good_mask = ( (~np.isnan(flux)) & (flux!=0) )

    #only use the markers within spec wavelength range
    val = value

    # I created this variable to more easily change the "window" of fitting for the continuum_fit, Trevor had val set to 5 - AP

    argstart = max(val, np.abs(wave-wave[good_mask][0]).argmin())
    argfin   = min(len(wave)-val, np.abs(wave-wave[good_mask][-1]).argmin())
    use_markers = markers[(markers>=wave[argstart])&(markers<=wave[argfin])]

    print(f"min good wave: {wave[good_mask][0]}, max good wave: {wave[good_mask][-1]}")
    print(f"use_markers: {use_markers}")

    #save [wave,flux] for continuum in dictionary
    cont_region = {}
    for i in range(len(use_markers)-1):
        #interpolate between markers on spectrum
        ind1 = abs(wave-use_markers[i]).argmin()
        ind2 = abs(wave-use_markers[i+1]).argmin()
        # fit_region = np.polyfit([np.nanmedian(wave[(ind1-5):(ind1+5)]), np.nanmedian(wave[(ind2-5):(ind2+5)])],
        #                         [np.nanmedian(flux[(ind1-5):(ind1+5)]), np.nanmedian(flux[(ind2-5):(ind2+5)])],
        #                          1) #do a linear fit on this section
        # print(fit_region)
        # Above is Trevor version


        fit_region = np.polynomial.polynomial.Polynomial.fit([np.nanmedian(wave[(ind1-val):(ind1+val)]), np.nanmedian(wave[(ind2-val):(ind2+val)])],
                                [np.nanmedian(flux[(ind1-val):(ind1+val)]), np.nanmedian(flux[(ind2-val):(ind2+val)])],
                                 1) #do a linear fit on this section

        # print(fit_region.convert().coef)
        # print(type(fit_region))
        fit_region = np.flip(fit_region.convert().coef)

        if len(use_markers)<=1:
            print("Wave coverage not big enough; <= 1 marker")

        #if only two markers in wavelength coverage, use
        #slope for whole continuum
        elif len(use_markers)==2:
            x = wave

        #if >2 markers and first, include indices <marker 1
        #up to second marker (inclusive)
        elif i==0:
            x = wave[:ind2] #double check this!

        #if >2 markers, this last interpolation; skip over
        #final marker index to end of wave
        elif i==(len(use_markers)-2):
            x = wave[ind1:]

        #otherwise interpolate from marker 1 to index
        #preceding marker 2; start inclusive and ends on ind1 so OK
        else:
            x = wave[ind1:ind2]
        y = fit_region[0]*x + fit_region[1]
        cont_region[i] = [x, y]

    #concatenate all the separate fits
    cont_spec = cont_region[0][1]
    for i in range(1, len(use_markers)-1):
        cont_spec = np.append(cont_spec, cont_region[i][1])

    #if len(cont_spec)-len(cont_morph)==1:
        #cont_morph = np.append(cont_morph, cont_morph[-1]) #fix an annoying index error below, coming from definition of mask_exp
    #print(len(cont_spec))
    #cont_spec = Small_Pix_Filter_TVM.SDSS_pixel_filter(wave, cont_spec, np.nanmin([len(wave)-1, 300]))
    #print(len(cont_spec))
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
    #Ignore BAL troughs
    bal = np.zeros(len(wave))
    if Identifier in BAL_regions.Troughs:
        bal_region = BAL_regions.Troughs[Identifier]
        bal[((wave/(1+z))>bal_region[0])&((wave/(1+z))<bal_region[1])] = 1
    #print((bal==1).sum())
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

    #print(jp_list)
    #print(jf_list)
    #print(len(wave))
    #interpolate across emission-line regions
    for i in range(len(jp_list)):
        if jp_list[i]==0:
            #if first or last
            med_flux[jp_list[i]:jf_list[i]+1] = med_flux[jf_list[i]+1]
        elif jf_list[i]==len(wave)-1:
            med_flux[jp_list[i]:jf_list[i]+1] = med_flux[jp_list[i]-1]
        else:
            #if middle, interpolate
            x = np.array([wave[jp_list[i]-1], wave[jf_list[i]+1]])
            y = np.array([med_flux[jp_list[i]-1], med_flux[jf_list[i]+1]])
            #idx = ~np.isnan(y) #for np.polyfit only
            m, b = np.polyfit(x, y, 1)
            w    = wave[jp_list[i]:jf_list[i]+1]
            med_flux[jp_list[i]:jf_list[i]+1] = m*w + b

    return med_flux
