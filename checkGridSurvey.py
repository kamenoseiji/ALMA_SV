import sys
import re
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
#
ELshadow = np.pi* 40.0 / 180.0
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
#spwNames = ['ALMA_RB_03#BB_1#SW-01#FULL_RES', 'ALMA_RB_03#BB_2#SW-01#FULL_RES', 'ALMA_RB_03#BB_3#SW-01#FULL_RES', 'ALMA_RB_03#BB_4#SW-01#FULL_RES']
BandNames, pattern = [], r'RB_..'
for spwName in spwNames: BandNames = BandNames + re.findall(pattern, spwName)
UniqBands = unique(BandNames).tolist(); NumBands = len(UniqBands)
spwLists, BandScans = [], []
for band_index in range(NumBands):
    spwLists = spwLists + [np.array(spw)[indexList( np.array([UniqBands[band_index]]), np.array(BandNames))].tolist()]
    BandScans = BandScans + [msmd.scansforspw(spwLists[band_index][0])]
    print ' ',
    print UniqBands[band_index] + ': SPW=' + `spwLists[band_index]`
#
#-------- Check source list
print '---Checking source list'
sourceList, posList = GetSourceList(msfile) 
numSource = len(sourceList)
SSOList   = np.where( (np.array(posList)[:,0] == 0.0) & (np.array(posList)[:,1] == 0.0) )[0].tolist()   # Solar System Objects
QSOList   = list(set(range(numSource)) - set(SSOList))
if FLcal in sourceList:
    print ' ',
    for source in sourceList: print source,
    print ''; print '  Solar System Objects:',
    for index in SSOList: print sourceList[index],
    print ''
    if len(SSOList) == 0: print '  No Solar System Object was observed.'; sys.exit()
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
    for bl_index in range(UseBlNum):
        blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
    #
    print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
    antDia = np.ones([UseAntNum])
    for ant_index in range(UseAntNum): antDia[ant_index] = msmd.antennadiameter(antList[antMap[ant_index]])['value']
    #-------- Check Scans of BandPass, EQualization, and FluxScaling
    if FCScans == []: FCScans = msmd.scansforintent("CALIBRATE_FLUX#ON_SOURCE")
    if BPScans == []: BPScans = msmd.scansforintent("CALIBRATE_BANDPASS#ON_SOURCE")
    ONScans = msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE")
    print '---SPWs and Scans for each receiver band'
    for band_index in range(NumBands):
        #-------- Check Calibrators
        FCScan = BandScans[band_index][indexList( FCScans, BandScans[band_index] )][0]
        BPScan = BandScans[band_index][indexList( BPScans, BandScans[band_index] )][0]
        ONScan = BandScans[band_index][indexList( ONScans, BandScans[band_index] )]
        EQScan = BPScan
        onsourceScans = unique([BPScan] + [FCScan] + ONScan.tolist()).tolist()
        scanNum = len(onsourceScans)
        #-------- Check AZEL
        azelTime, AntID, AZ, EL = GetAzEl(msfile)
        azelTime_index = np.where( AntID == UseAnt[refantID] )[0].tolist() 
        OnEL, sourceIDscan = [], []
        for scan_index in range(scanNum):
            sourceIDscan.append( msmd.fieldsforscan(onsourceScans[scan_index])[0])
            refTime = np.median(msmd.timesforscan(onsourceScans[scan_index]))
            OnEL.append(EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - refTime))]])
            print 'Scan%d : %s EL=%4.1f' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0*OnEL[scan_index]/np.pi)
        #
        if BPcal in sourceList: BPScan = list(set(msmd.scansforfield(sourceList.index(BPcal))) & set(onsourceScans))[0]
        if FLcal in sourceList: FCScan = list(set(msmd.scansforfield(sourceList.index(FLcal))) & set(onsourceScans))[0]
        if EQcal in sourceList: EQScan = list(set(msmd.scansforfield(sourceList.index(EQcal))) & set(onsourceScans))[0]
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
        execfile(SCR_DIR + 'checkSEFD.py')
    #
    msmd.done()
#
