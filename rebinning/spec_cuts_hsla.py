"""
rebinning/spec_cuts_hsla.py
Migrated from: Trevor Code/SpecCuts_HSLA.py
Migration date: 2026-04-30
No algorithmic changes — only header added.

Original description:
In some cases, the edges of the co-added spectra from the HSLA are just
a complete mess.  We can automate the process for most, but it will in the
end be easier and less time consuming to just look at the final re-bin for
each, check the edges, then write down where to cut them here.
"""

BlueEdges = {
    "NEWQZ003_coadd_G140L_final": 1180,
    "QSO-003916-511701_coadd_G140L_final": 1190,
    "QSO-B0923+201_coadd_G140L_final": 1105
    #"SDSSJ031027.82-004950.7_coadd_FUVM_final": 1615
}

RedEdges = {
    "J084349.75+261910.7_coadd_G140L_final": 1610,
    "MARK509_coadd_G140L_final": 2000,
    "MRK-231_coadd_G140L_final": 1800,
    "MRK-335_coadd_G140L_final": 2000,
    "NEWQZ003_coadd_G140L_final": 1650,
    "NEWQZ004_coadd_G140L_final": 1700,
    "NEWQZ013_coadd_G140L_final": 1660,
    "NEWQZ015_coadd_G140L_final": 1690,
    "NEWQZ027_coadd_G140L_final": 1670,
    "NGC-3516_coadd_FUVM_final": 1700,
    "NGC-985_coadd_G140L_final": 1990,
    "PDS456_coadd_G140L_final": 1680,
    "PG1011-040_coadd_FUVM_final": 1660,
    "PG1211+143_coadd_G140L_final": 1750,
    "Q1659+6202_coadd_G140L_final": 1600,
    "QSO-003916-511701_coadd_G140L_final": 2000,
    "QSO-B0923+201_coadd_G140L_final": 1623,
    "QSO-B1351+6400_coadd_G140L_final": 1800,
    "SDSSJ001224.01-102226.5_coadd_G140L_final": 1620,
    "SDSSJ015530.02-085704.0_coadd_G140L_final": 1695,
    "SDSSJ031027.82-004950.7_coadd_FUVM_final": 1615,
    "SDSSJ075620.08+304535.3_coadd_G140L_final": 1585,
    "SDSSJ115758.72-002220.8_coadd_G140L_final": 1580,
    "SDSSJ122534.79-024757.1_coadd_G140L_final": 1620,
    "SDSSJ132144.97+033055.7_coadd_G140L_final": 1580,  # AP 2026-04-22: noisy red tail from ~1582 A onwards; G140L limit for this object at z=0.269
    "SDSSJ230845.59-091123.9_coadd_G140L_final": 1690,
    #These aren't as obvious as the above, but will still benefit from an additional cut
    "PG1309+355_coadd_G140L_final": 1650,
    "NGC-5548_coadd_FUVM_final": 1690
}
