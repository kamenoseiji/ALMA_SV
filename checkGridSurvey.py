import sys
import re
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
#---- Definitions
BANDPA = [0.0, 45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, 90.0, 0.0]   # X-pol orientation for Band-1, 2, 3, 4, 5, 6, 7, 8, 9, and 10
ELshadow = np.pi* 40.0 / 180.0
SSOCatalog = ['Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus']
SSOscore   = [[ 5.0,     4.0,       1.0,        1.0,        0.1,     0.2,  0.3,      0.2,     0.1,      0.1,     0.1,    10.0,   10.0,     10.0],   # Band 1
              [ 6.0,     5.0,       1.0,        1.0,        0.1,     0.6,  0.5,      0.3,     0.2,      0.2,     0.2,    10.0,   10.0,     10.0],   # Band 2
              [ 7.0,     6.0,       1.0,        1.0,        0.2,     0.7,  0.7,      0.5,     0.3,      0.3,     0.3,    10.0,   10.0,     10.0],   # Band 3
              [ 8.0,     7.0,       2.0,        2.0,        0.3,     0.9,  0.9,      0.6,     0.4,      0.4,     0.4,     8.0,   10.0,      8.0],   # Band 4
              [ 9.0,     8.0,       3.0,        3.0,        0.5,     1.0,  1.0,      0.8,     0.5,      0.5,     0.5,     4.0,    4.0,      4.0],   # Band 5
              [10.0,     9.0,       4.0,        4.0,        5.0,     3.0,  3.0,      1.0,     0.6,      0.6,     0.6,     2.0,    2.0,      2.0],   # Band 6
              [10.0,     9.0,       5.0,        5.0,        7.0,     4.0,  4.0,      3.0,     0.7,      0.7,     0.7,     1.0,    1.0,      1.0]]   # Band 7
#-------- Procedures
msfile = wd + prefix + '.ms'
#-------- Check Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
blNum = antNum* (antNum - 1) / 2
#-------- Check SPWs of atmCal
msmd.open(msfile)
print prefix
print '---Checking spectral windows'
spw = list(set(msmd.tdmspws()) & set(msmd.spwsforintent("CALIBRATE_ATMOSPHERE*"))); spw.sort()
spwNum = len(spw)
spwNames = msmd.namesforspws(spw)
BandNames, pattern = [], r'RB_..'
for spwName in spwNames: BandNames = BandNames + re.findall(pattern, spwName)
UniqBands = unique(BandNames).tolist(); NumBands = len(UniqBands)
spwLists, BandScans, BandPA = [], [], []
for band_index in range(NumBands):
    BandPA = BandPA + [(BANDPA[int(UniqBands[band_index][3:5])] + 90.0)*pi/180.0]
    spwLists = spwLists + [np.array(spw)[indexList( np.array([UniqBands[band_index]]), np.array(BandNames))].tolist()]
    BandScans = BandScans + [msmd.scansforspw(spwLists[band_index][0])]
    print ' ',
    print UniqBands[band_index] + ': SPW=' + `spwLists[band_index]`
#
#-------- Check source list
print '---Checking source list'
sourceList, posList = GetSourceList(msfile) 
numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList) ) # Find Solar System Objects in the source list
if len(SSOList) == 0: print '  No Solar System Object was observed.'; sys.exit()
QSOList   = list(set(range(numSource)) - set(SSOList))
for source in sourceList: print source,
print ''; print '  Solar System Objects:',
for index in SSOList: print sourceList[index],
print ''
#-------- Check MJD for Ambient Load
print '---Checking time for ambient and hot load'
timeOFF, timeAMB, timeHOT = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT")
if len(timeAMB) == 0:
    timeON  = msmd.timesforintent("CALIBRATE_ATMOSPHERE#ON_SOURCE")
    timeAMB = timeON[(np.where( np.diff(timeON) > 20.0)[0] - 8).tolist()]
    timeHOT = timeON[(np.where( np.diff(timeON) > 20.0)[0] - 2).tolist()]
#
#-------- Configure Array
print '---Checking array configulation'
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
blMap, blInv= range(UseBlNum), [False]* UseBlNum
try:
    refantID = np.where(antList[UseAnt] == refant )[0][0]
except:
    ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
    for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
    timeStamp, UVW = GetUVW(msfile, spw[0], msmd.scansforspw(spw[0])[0])
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = bestRefant(uvDist)
#
print '  Use ' + antList[UseAnt[refantID]] + ' as the refant.'
#-------- Baseline Mapping
print '---Baseline Mapping'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
UseAntNum = len(antMap)
UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
antDia = np.ones([UseAntNum])
for ant_index in range(UseAntNum): antDia[ant_index] = msmd.antennadiameter(antList[antMap[ant_index]])['value']
#-------- Check Scans of BandPass, EQualization, and FluxScaling
try:
    FCScans = msmd.scansforintent("CALIBRATE_FLUX#ON_SOURCE")
except:
    FCScans = msmd.scansforintent("CALIBRATE_AMPLI#ON_SOURCE")
#
if not BPcal: BPScans = msmd.scansforintent("CALIBRATE_BANDPASS#ON_SOURCE")
#
ONScans = msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE")
print '---SPWs and Scans for each receiver band'
msmd.done()
#for band_index in range(1):
for band_index in range(NumBands):
    bandID = int(UniqBands[band_index][3:5])-1
    msmd.open(msfile)
    #-------- Select Flux Calibrator
    #if not FLcal in sourceList: FLcal = sourceList[SSOList[0]]
    #-------- Check Calibrators
    FCScan = BandScans[band_index][indexList( FCScans, BandScans[band_index] )]
    BPScan = BandScans[band_index][indexList( BPScans, BandScans[band_index] )][0]
    ONScan = BandScans[band_index][indexList( ONScans, BandScans[band_index] )]
    EQScan = BPScan
    onsourceScans = unique([BPScan] + FCScan.tolist() + ONScan.tolist()).tolist()
    scanNum = len(onsourceScans)
    #-------- Check AZEL
    azelTime, AntID, AZ, EL = GetAzEl(msfile)
    azelTime_index = np.where( AntID == UseAnt[refantID] )[0].tolist() 
    azel = np.r_[AZ[azelTime_index], EL[azelTime_index]].reshape(2, len(azelTime_index))
    OnEL, sourceIDscan, FLscore = [], [], np.zeros(scanNum)
    for scan_index in range(scanNum):
        sourceIDscan.append( msmd.sourceidforfield(msmd.fieldsforscan(onsourceScans[scan_index])[0]))
        refTime = np.median(msmd.timesforscan(onsourceScans[scan_index]))
        OnEL.append(EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - refTime))]])
        print 'Scan%d : %s EL=%4.1f' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0*OnEL[scan_index]/np.pi)
        if sourceIDscan[scan_index] in SSOList: FLscore[scan_index] = np.exp(np.log(math.sin(OnEL[scan_index])-0.34))* SSOscore[bandID][SSOCatalog.index(sourceList[sourceIDscan[scan_index]])]
    #
    #if FLcal in sourceList: FCScan = list(set(msmd.scansforfield(FLcal)) & set(onsourceScans))[0]
    FCScan = onsourceScans[np.argmax(FLscore)]; FLcal = sourceList[sourceIDscan[np.argmax(FLscore)]]
    print 'Use %s [EL = %4.1f] as Flux Equalizer' % (FLcal, 180.0* OnEL[np.argmax(FLscore)]/np.pi)
    if BPcal in sourceList: BPScan = list(set(msmd.scansforfield(BPcal)) & set(onsourceScans))[0]
    if EQcal in sourceList: EQScan = list(set(msmd.scansforfield(EQcal)) & set(onsourceScans))[0]
    #-------- SSO in observed source list
    BandSSOList = list( set(SSOList) & set(sourceIDscan) )
    #-------- Avoid EQ == FL
    if EQScan == FCScan or OnEL[onsourceScans.index(EQScan)] < ELshadow :
        QSOscanIndex = indexList(QSOList, np.array(sourceIDscan))
        QSOEL = np.array(OnEL)[QSOscanIndex]
        EQScan = onsourceScans[QSOscanIndex[np.argmax(QSOEL)]]
        BPScan = EQScan
    #
    EQcal  = sourceList[sourceIDscan[onsourceScans.index(EQScan)]]
    print 'Use %s [EL = %4.1f] as Gain Equalizer' % (EQcal, 180.0* OnEL[onsourceScans.index(EQScan)]/np.pi)
    #-------- Polarization setup
    spw = spwLists[band_index]; spwNum = len(spw); polNum = msmd.ncorrforpol(msmd.polidfordatadesc(spw[0]))
    pPol, cPol = [0,1], []  # parallel and cross pol
    PolList = ['X', 'Y']
    if polNum == 4: pPol, cPol = [0,3], [1,2]  # parallel and cross pol
    ppolNum, cpolNum = len(pPol), len(cPol)
    msmd.done()
    if polNum == 2: pPol, cPol = [0,1], []   ; ppolNum, cpolNum = len(pPol), len(cPol); execfile(SCR_DIR + 'checkSEFD.py')
    if polNum == 4:
        pPol, cPol = [0,3], [1,2]; ppolNum, cpolNum = len(pPol), len(cPol)
        if Apriori: execfile(SCR_DIR + 'aprioriFlux.py')
        else:  execfile(SCR_DIR + 'checkSEFDStokes.py')
    #
#
