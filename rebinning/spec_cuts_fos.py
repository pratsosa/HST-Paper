"""
rebinning/spec_cuts_fos.py
Migrated from: Trevor Code/SpecCuts_FOS.py
Migration date: 2026-04-30
No algorithmic changes — only header added.

Original description:
In some cases, the edges of the co-added spectra from the HSLA are just
a complete mess.  We can automate the process for most, but it will in the
end be easier and less time consuming to just look at the final re-bin for
each, check the edges, then write down where to cut them here.
"""

BlueEdges = {
    "3C215 - y0pe0602t": 1155, "3C215 - y0pe0603t": 1155, "3C215 - y0pe0604t": 1155, "3C215 - y0pe0607t": 2300,
    "IRAS04505-2958 - y3gs0302t": 1272,
    "J00392-5117 - y3790105t": 1692,
    "J10309+3102 - y34w5102t": 1375,
    "J12047+2754 - y38o0603t": 1405, "J12047+2754 - y17j0203t": 1405, "J12047+2754 - y17j0204t": 1405,
    "J14052+2555 - y1hk0e02t": 1390, "J14052+2555 - y38o1005t": 1390, "J14052+2555 - y29c0c02t": 1390,
    "J14297+4747 - y38o1308t": 1330,
    "J14467+4035 - y10q0304t": 1295, "J14467+4035 - y10q0305t": 1295,
    "J15455+4846 - y27o0102t": 1166, "J15455+4846 - y27o0103t": 1166, "J15455+4846 - y27o0104t": 1166,
    "J16142+2604 - y1hk0p02t": 1444,
    "J16279+5522 - y38o1605t": 1465,
    "J21377-1432 - y10f0104t": 1360,
    "J23519-0109 - y1hk1002t": 1395,
    "LBQS0003+0146 - y29c0102t": 1360,
    "LBQS0017+0209 - y29c0202t": 1170,
    "LBQS1132-0302 - y29c0502t": 1330,
    "LBQS1138+0204 - y29c0602t": 1180,
    "LBQS1317-0142 - y29c0a02t": 1350, "LBQS1317-0142 - y29c0a03t": 1350,
    "Mrk0478 - y1hk0l02t": 1495, "Mrk0478 - y38o1403t": 1495,
    "Mrk205 - y1hi0x02t": 1144, "Mrk205 - y0na0105t": 1510, "Mrk205 - y0na0106t": 1510,
    "Mrk813 - y1hk0f02t": 1468,
    "NGC1566 - y0h7510at": 1390,
    "NGC5548 - y0ya0204t": 1565,
    "NGC7469 - y3b60106t": 1155, "NGC7469 - y3b60107t": 1155,
    "PG0026+129 - y2jk0109t": 1435,
    "PG0052+251 - y1hh0202t": 1445,
    "PG0947+396 - y38o0105t": 1345,
    "PG0953+414 - y0ml0103t": 1315,
    "PG1001+054 - y38o0208t": 1435,
    "PG1114+445 - y38o0309t": 1435, "PG1114+445 - y3ai100bt": 1435, "PG1114+445 - y3ai100ct": 1435,
    "PG1115+407 - y38o0405t": 1430, "PG1115+407 - y38o0407t": 3150,
    "PG1116+215 - y38o0503t": 1385, "PG1116+215 - y17j0103t": 1385, "PG1116+215 - y38o0505t": 2755,
    "PG1121+422 - y29c0302t": 1345,
    "PG1211+143 - y0iz0406t": 1490,
    "PG1302-102 - y1020103t": 1265,
    "PG1307+085 - y1hk0102t": 1440,
    "PG1322+659 - y38o0804t": 1390,
    "PG1352+183 - y38o0903t": 1435,
    "PG1404+226 - y33s0204t": 1485,
    "PG1415+451 - y38o1103t": 1445,
    "PG1545+210 - y0pe0d02t": 1315,
    "SBS1704+608 - y0rv0g03t": 1190,
    "UGC12163 - y3790405t": 1572, "UGC12163 - y3790402t": 1190, "UGC12163 - y3790402t": 1190
}

RedEdges = {
    "J10309+3102 - y34w5102t": 1900,
    "PG1116+215 - y38o0502t": 2600, "PG1116+215 - y17j0104t": 2600
}
