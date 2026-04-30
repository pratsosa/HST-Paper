# Originally: Trevor Code/manual_fix_config.py
# Copied: 2026-04-29 (Phase 2 migration)
# Changes from source: header added only
"""
Configuration file for manual ICA fixes.
Each object can have:
1. custom_mask_pixels: array of wavelengths to mask
2. forced_components: force ICA to use specific component set ('mod', 'low', 'high')
3. Both can be specified, or neither (in which case no manual intervention)
"""

import numpy as np

# Manual fix configurations for individual objects
MANUAL_FIX_CONFIG = {
    'RBS1763': {
        'custom_mask_pixels': np.concatenate((
            np.linspace(1530, 1540, 50),  # filter_low
            np.array([1544.20174858, 1545.98060034, 1546.33661651, 1546.69271467, 
                     1547.40515702, 1548.47443591, 1548.83102638, 1549.18769896, 
                     1549.54445368, 1552.75894544, 1553.47418203])  # bad_values
        )),
        'forced_components': None
    },
    
    '1H1613-097': {
        'custom_mask_pixels': np.concatenate((
            np.linspace(1544, 1545.84769975, 100),  # filter_low
            np.array([1546.91590242, 1547.62844757, 1549.41124673, 1549.76805293, 
                     1550.8389646, 1551.55331683]),  # filter_waves
            np.linspace(1552, 1563, 100)  # filter_high
        )),
        'forced_components': 'high'
    },
    
    '3C057': {
        'custom_mask_pixels': np.concatenate((
            np.linspace(1880, 1900, 1000),  # filter_low
            np.linspace(1390, 1440, 1000),  # filter_2
            np.linspace(1300, 1320, 1000),  # filter_3
            np.linspace(1575, 1610, 1000),  # filter_4
            np.array([1546.12850417, 1546.4845544, 1546.84068663, 1547.19690087, 
                     1554.33843728, 1554.69637814, 1555.05440143, 1555.41250717, 
                     1555.77069537])
        )),
        'forced_components': None
    },
    
    '3C095': {
        'custom_mask_pixels': None,  # No masking
        'forced_components': None
    },
    
    '3C110': {
        'custom_mask_pixels': np.array([1575.24631086, 1575.60906649, 1575.9719056, 
                                       1576.33482838, 1579.24122051, 1579.60489611, 
                                       1579.96865545, 1580.33249857, 1580.69642547]),
        'forced_components': None
    },
    
    '3C390': {
        #'custom_mask_pixels': np.linspace(1525, 1540, 1000),
        'custom_mask_pixels': None,
        'forced_components': 'high'
    },
    
    '3C207': {
        'custom_mask_pixels': np.array([1555.20386586, 1556.6369214, 1557.35394428, 
                                       1557.71257957, 1558.78898103, 1559.8661263]),
        'forced_components': None
    },
    
    '3C215': {
        'custom_mask_pixels': None,
        'forced_components': None
    },
    
    '3C232': {
        'custom_mask_pixels': None,
        'forced_components': None
    },
    
    '3C263': {
        'custom_mask_pixels': None,
        'forced_components': None
    },
    
    '3C273': {
        'custom_mask_pixels': None,
        'forced_components': None
    },
    # ---------------------------------------------------------------------------------------------
    # Here begin the Probably Good Objects
    '3C254': {
        'custom_mask_pixels': None,
        'forced_components': 'high'
    },
    
    'FBQSJ1010+3003': {
        # Can't improve this one without better data
        'custom_mask_pixels': None,
        'forced_components': None
    },
    
    'HE0132-4313': {
        'custom_mask_pixels': None,
        'forced_components': None
    },
    
    'Mrk231': {
        'custom_mask_pixels': np.linspace(1540, 1545, 1000),
        'forced_components': 'high'
    },
    
    'FBQSJ1251+2404': {
        'custom_mask_pixels': np.linspace(1580, 1650, 1000),
        'forced_components': None
    },
    
    'PDS456': {
        # Famous high outflow object - no changes needed
        'custom_mask_pixels': None,
        'forced_components': None
    },
    
    'SDSSJ165958.94+620218.1': {
        'custom_mask_pixels': None,
        'forced_components': 'high'
    },
    
    'RXJ1230.8+0115': {
        'custom_mask_pixels': np.concatenate((
            np.array([1548.5, 1549]),  # filter_pix
            np.linspace(1580, 1605, 1000)  # filter_high
        )),
        'forced_components': 'mod'
    },
    
    'Mrk478': {
        # GTR: Still needs work?  Definitely don't mask 1800-1900 or >2200.
        # TODO: Fix this
        'custom_mask_pixels': None,  # Complex masking logic needed - handle separately
        'forced_components': None
    },
    
    'SDSSJ163013.56+375821.7': {
        # GTR: Is redshift of SDSS part wrong?  Why is MgII line in the wrong place?
        # TODO: Check this
        'custom_mask_pixels': np.concatenate((
            np.array([1541.9897,1542.344, 1542.700 , 1544.477, 1544.832 , 1545.188 , 1546.968, 1547.324, 1547.681, 1548.750,
                      1555.89928562, 1556.25758592, 1559.84513017, 1560.20433914, 1562.00162517, 1562.36133074, 1562.72111916 ]), 
            np.linspace(1551.24888161, 1556.61596873  , 1000)
        
        )), 
        'forced_components': None
    },
    'ESO499-41': {
        # GTR: But I don't think that it is a type-1 quasar.  Maybe type-2, should check.
        # TODO: Check this
        'custom_mask_pixels': np.linspace(1547, 1551.5, 1000),
        'forced_components': None
    },
    'PKS2300-68': {
        'custom_mask_pixels': np.concatenate((
            np.linspace(1543, 1547, 1000),  # filter_low,
            np.linspace(1551, 1554 , 1000)  # filter_med
        )),
        'forced_components': None
    },
    'NGC1566':{
        'custom_mask_pixels': np.linspace(1250, 1500, 1000),
        'forced_components': 'mod'
    },
    'NGC4253':{
        # GTR: I don't think that it is type-1.
        # TODO: Check this
        'custom_mask_pixels': None,
        'forced_components': 'high'
    },
    'NGC4395':{
        'custom_mask_pixels': np.linspace(1550, 1560, 1000),
        'forced_components': None
    },
    'PG1415+451':{
        'custom_mask_pixels': None,
        'forced_components': 'high'
    },     
    'WPVS7':{
        # GTR: Could mark as good, except that the lines could be narrow, so leave with some uncertainty.
        'custom_mask_pixels': None,
        'forced_components': None
    },
    'RBS9':{
        'custom_mask_pixels': None,
        'forced_components': None
    },
    '3C351.0':{
        'custom_mask_pixels': None,
        'forced_components': 'mod'
    },
    '2MASSJ10053274-2417161':{
        # GTR: Insufficient spectral coverage to include in analysis.
        # TODO: Remove / Flag
        'custom_mask_pixels': None,
        'forced_components': None
    }, 
    'HS0033+4300':{
        # GTR: Not "bad", just "marginal".  The values you get are probably about right.
        # TODO: Move to Marginal
        'custom_mask_pixels': np.linspace(1545, 1554, 1000),
        'forced_components': None
    },
    '2E1644':{
        # GTR: Not "bad", just "marginal".  The values you get are probably about right.
        # TODO: Move to Marginal
        'custom_mask_pixels': None,
        'forced_components': 'high'
    },
    'SDSSJ132059.41+295728.1':{
        # GTR: Not "bad", just "marginal".  The values you get are probably about right.
        # TODO: Move to Marginal
        'custom_mask_pixels': None,
        'forced_components': 'high'
    },
    '6C134441+623604':{
        # GTR: Not "bad", just "marginal".  The values you get are probably about right.
        # TODO: Move to Marginal
        'custom_mask_pixels': np.linspace(1546, 1551, 1000),
        'forced_components': None
    },
    'SDSSÂ\xa0J093653.84+533126.9':{
        # GTR: Not "bad", just "marginal".  The values you get are probably about right.
        # TODO: Move to Marginal
        'custom_mask_pixels': np.linspace(1542, 1552.5, 1000),
        'forced_components': None
    },
    'NGC1275':{
        # GTR: Both of those are bad.  Not sure why.  Redshift might be wrong for starters.  But I'd think that it would do better.  
        # Much higher EW and smaller blueshift is what I'd expect.
        # TODO Check this
        'custom_mask_pixels': None,
        'forced_components': 'high'
    },
    'PG1351+640':{
        # GTR: Agreed, though all bad.  I'd guess that it would help to temporarily fake some data points (using the ICA fit) shortward of 1510A, then fit again?  
        # Though the one you picked looks like it is reasonable in terms of blueshift and EW.
        # TODO: ICA extrapolation option?
        'custom_mask_pixels':None,
        'forced_components': None
    },
    'PG1404+226':{
        # Was in Bad, got moved to Probably Good.
        'custom_mask_pixels': None,
        'forced_components': None
    },
    'SDSSJ015536.03+311518.0':{
        # GTR: No, fit pulled up arbitrarily because of spike at 1585.  
        # As with PG1351, it might help to extrapolate past 1580A with the ICA, then fit again.
        # TODO: ICA extrapolation option?
        'custom_mask_pixels': None,
        'forced_components': None
    },
    'SDSSJ005346.15+223222.3':{
        # GTR: "marginal" because those two fits give very different predictions longward of 1550A.
        # TODO: Move to Marginal
        'custom_mask_pixels': None,
        'forced_components': None
    },
    'SDSS J115758.72-002220.8':{
        # GTR: Bad.  Too noisy and not enough spectrum to work with.  I wouldn't include in any analysis.
        # TODO: Flag / Remove
        'custom_mask_pixels': None,
        'forced_components': None 
    },
    'SDSS J152139.66+033729.2':{
        # GTR: Marginal?  Gap at 1400 makes it hard to tell if CIV fit is good or not.
        # TODO: Move to Marginal
        'custom_mask_pixels': None,
        'forced_components': None  
    }
}