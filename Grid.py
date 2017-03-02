SSOCatalog = ['Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus']
SSOscore   = [
[ 5.0,     4.0,       1.0,        1.0,        0.1,     0.2,  0.3,      0.2,     0.1,      0.1,     0.1,    10.0,   10.0,     10.0],   # Band 1
[ 6.0,     5.0,       1.0,        1.0,        0.1,     0.6,  0.5,      0.3,     0.2,      0.2,     0.2,    10.0,   10.0,     10.0],   # Band 2
[ 7.0,     6.0,       1.0,        1.0,        0.2,     0.7,  0.7,      0.5,     0.3,      0.3,     0.3,    10.0,   10.0,     10.0],   # Band 3
[ 8.0,     7.0,       2.0,        2.0,        0.3,     0.9,  0.9,      0.6,     0.4,      0.4,     0.4,     8.0,   10.0,      8.0],   # Band 4
[ 9.0,     8.0,       3.0,        3.0,        0.5,     1.0,  1.0,      0.8,     0.5,      0.5,     0.5,     4.0,    4.0,      4.0],   # Band 5
[10.0,     9.0,       4.0,        4.0,        5.0,     3.0,  3.0,      1.0,     0.6,      0.6,     0.6,     2.0,    2.0,      2.0],   # Band 6
[10.0,     9.0,       5.0,        5.0,        7.0,     4.0,  4.0,      3.0,     0.7,      0.7,     0.7,     1.0,    1.0,      1.0]]   # Band 7
ELshadow = np.pi* 40.0 / 180.0
sourceDic = {
'J0237+288':'J0237+2848',
'J0238+166':'J0238+1636',
'3c84':'J0319+4130',
'0334-401':'J0334-4008',
'J0423-013':'J0423-0120',
'J0510+180':'J0510+1800',
'J0519-454':'J0519-4546',
'J0522-364':'J0522-3627',
'J1037-295':'J1037-2934',
'J1058+015':'J1058+0133',
'J1107-448':'J1107-4449',
'J1146+399':'J1146+3958',
'3c273':'J1229+0203',
'3c279':'J1256-0547',
'J1337-129':'J1337-1257',
'J1427-421':'J1427-4206',
'J1517-243':'J1517-2422',
'J1550+054':'J1550+0527',
'J1613-586':'J1617-5848',
'3c345':'J1642+3948',
'J1733-130':'J1733-1304',
'J1751+096':'J1751+0939',
'J1924-292':'J1924-2914',
'J2025+337':'J2025+3343',
'J2056-472':'J2056-4714',
'J2148+069':'J2148+0657',
'J2232+117':'J2232+1143',
'3c454.3':'J2253+1608',
'J2258-279':'J2258-2758'}

def sourceRename(sourceList):
    renameList = []
    for srcname in sourceList: renameList = renameList + [sourceDic.get(srcname, srcname)]
    return renameList
#
catalogStokesI = {
'J1256-0547':   7.0,
'J2253+1608':  12.6,
'J0006-0623':   3.1,
'J0522-3627':   4.5,
'J1229+0203':   8.2,
'J0854+2006':   3.2,
'J1924-2914':   3.4,
'J1058+0133':   3.0,
'J1751+0939':   2.0,
'J1337-1257':   2.0,
'J2232+1143':   2.2,
'J0238+1636':   1.3,
'J1037-2934':   0.8,
'J0334-4008':   0.5,
'J2056-4714':   1.1,
'J1517-2422':   2.0,
'J1642+3948':   1.9,
'J0750+1231':   0.8,
'J2258-2758':   0.9,
'J1427-4206':   2.0,
'J1107-4449':   0.8,
'J0237+2848':   1.2,
'J1550+0527':   0.6,
'J0538-4405':   1.5,
'J0510+1800':   3.2,
'J2148+0657':   0.9,
'J1733-1304':   2.0,
'J1617-5848':   0.5,
'J0423-0120':   0.7,
'J2357-5311':   0.4,
'J1146+3958':   0.3,
'J2025+3343':   0.6,
'J0519-4546':   0.9,
'J0635-7516':   0.4}
#
catalogStokesQ = {
'J1256-0547':  0.010,
'J2253+1608': -0.033,
'J0006-0623':  0.089,
'J0522-3627': -0.035,
'J1229+0203':  0.005,
'J0854+2006':  0.014,
'J1924-2914': -0.011,
'J1058+0133':  0.018,
'J1751+0939':  0.018,
'J1337-1257':  0.026,
'J2232+1143':  0.018,
'J0238+1636':  0.007,
'J1037-2934': -0.057,
'J0334-4008': -0.029,
'J2056-4714':  0.021,
'J1517-2422':  0.011,
'J1642+3948':  0.002,
'J0750+1231':  0.010,
'J2258-2758': -0.006,
'J1427-4206': -0.000,
'J1107-4449':  0.007,
'J0237+2848':  0.012,
'J1550+0527':  0.021,
'J0538-4405': -0.001,
'J0510+1800': -0.005,
'J2148+0657':  0.016,
'J1733-1304':  0.006,
'J1617-5848': -0.020,
'J0423-0120':  0.005,
'J2357-5311': -0.006,
'J1146+3958': -0.010,
'J2025+3343':  0.003,
'J0519-4546': -0.001,
'J0635-7516': -0.004}
#
catalogStokesU = {
'J1256-0547':  0.105,
'J2253+1608':  0.006,
'J0006-0623':  0.060,
'J0522-3627': -0.029,
'J1229+0203': -0.020,
'J0854+2006': -0.043,
'J1924-2914': -0.016,
'J1058+0133':  0.001,
'J1751+0939': -0.020,
'J1337-1257': -0.005,
'J2232+1143':  0.014,
'J0238+1636':  0.033,
'J1037-2934': -0.006,
'J0334-4008': -0.046,
'J2056-4714':  0.013,
'J1517-2422': -0.005,
'J1642+3948':  0.012,
'J0750+1231':  0.025,
'J2258-2758': -0.024,
'J1427-4206': -0.010,
'J1107-4449': -0.024,
'J0237+2848':  0.008,
'J1550+0527':  0.019,
'J0538-4405':  0.011,
'J0510+1800': -0.001,
'J2148+0657': -0.003,
'J1733-1304':  0.002,
'J1617-5848':  0.019,
'J0423-0120': -0.016,
'J2357-5311': -0.021,
'J1146+3958': -0.008,
'J2025+3343':  0.004,
'J0519-4546': -0.002,
'J0635-7516': -0.001}
#