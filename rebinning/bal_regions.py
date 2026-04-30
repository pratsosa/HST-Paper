"""
rebinning/bal_regions.py
Migrated from: Trevor Code/BAL_regions.py
Migration date: 2026-04-30
No algorithmic changes — only header added.

Note: BAL_regions.py is identical between Trevor Code/ and Trevor Code/ICA Scripts/.
Placed here (rebinning/) because spec_morph.py (rebinning version) imports it.
Phase 2 ICA code can import from this module: `from rebinning import bal_regions`.

Original description:
Name objects here that are affected by BAL troughs, and mask
appropriate region.  Following the Rankine+20 recipe exactly
won't quite work here, since we are limited in wavelength
coverage already, and don't want to mask anything unnecessary.
"""

#10 August
BlueEdges = {
    "B21425+26 - o65637010": 1180,
    "J13122+3515 - o65630010": 1370,
    "J13253-3824 - o56j01010": 1535, "J13253-3824 - o56j01020": 1535, "J13253-3824 - o56j01030": 1535, "J13253-3824 - o56j01040": 1535,
    "J21377-1432 - o65631010": 1375,
    "Mrk509 - odjh01040": 1570, "Mrk509 - odjh01050": 1570,
    "NGC3227 - o5kp01010": 1680, "NGC3227 - o5kp01020": 2950,
    "NGC3516 - o4st02020": 2900, "NGC3516 - o56c01050": 2900,
    "PG1351+640 - o65616010": 1505,
    "TONS180 - o58p01020": 1498
}

Troughs = {
    "J21148+0607": [1439, 1510]
}
