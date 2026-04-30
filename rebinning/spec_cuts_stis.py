"""
rebinning/spec_cuts_stis.py
Migrated from: Trevor Code/SpecCuts_STIS.py
Migration date: 2026-04-30
No algorithmic changes — only header added.

Note: SpecCuts_STIS.py was not listed in the original Phase 1 source→destination
map but is required because cut_edge_pix.py and coadd.py have STIS code paths,
and the FOS+STIS runner processes STIS objects.
"""

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

RedEdges = {

}
