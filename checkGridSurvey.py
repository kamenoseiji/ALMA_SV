import sys
import re
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
#
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
spwLists, BandScans = [], []
for band_index in range(NumBands):
    spwLists = spwLists + [np.array(spw)[indexList( np.array([UniqBands[band_index]]), np.array(BandNames))].tolist()]
    BandScans = BandScans + [msmd.scansforspw(spwLists[band_index][0])]
#
#-------- Check source list
print '---Checking source list'
sourceList, posList = GetSourceList(msfile) 
numSource = len(sourceList)
SSOList   = np.where( (np.array(posList)[:,0] == 0.0) & (np.array(posList)[:,1] == 0.0) )[0].tolist()   # Solar System Objects
for source in sourceList: print source,
print ''; print 'Solar System Objects:',
for index in SSOList: print sourceList[index],
print ''
if len(SSOList) == 0: print 'No Solar System Object was observed.'; sys.exit()
#-------- Check MJD for Ambient Load
print '---Checking time for ambient and hot load'
timeOFF = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE")
timeAMB = msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT")
timeHOT = msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT")
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
print 'Use ' + antList[UseAnt[refantID]] + ' as the refant.'
#-------- Baseline Mapping
print '---Baseline Mapping'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
UseAntNum = len(antMap)
UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
for bl_index in range(UseBlNum):
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#
print `len(np.where( blInv )[0])` + ' baselines are inverted.'
antDia = np.ones([UseAntNum])
for ant_index in range(UseAntNum): antDia[ant_index] = msmd.antennadiameter(antList[antMap[ant_index]])['value']
#-------- Check Scans of BandPass, EQualization, and FluxScaling
FCScans, BPScans, ONScans = msmd.scansforintent("CALIBRATE_FLUX#ON_SOURCE"), msmd.scansforintent("CALIBRATE_BANDPASS#ON_SOURCE"), msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE")
print '---SPWs and Scans for each receiver band'
for band_index in range(NumBands):
    print '%s: SPW = ' % (UniqBands[band_index]),
    print spwLists[band_index]
    FCScan = BandScans[band_index][indexList( FCScans, BandScans[band_index] )][0]
    BPScan = BandScans[band_index][indexList( BPScans, BandScans[band_index] )][0]
    ONScan = BandScans[band_index][indexList( ONScans, BandScans[band_index] )]
    EQScan = BPScan
    onsourceScans = [BPScan] + [FCScan] + ONScan.tolist()
    scanNum = len(onsourceScans)
    if BPcal in sourceList: BPScan = list(set(msmd.scansforfield(sourceList.index(BPcal))) & set(onsourceScans))[0]
    if FLcal in sourceList: FCScan = list(set(msmd.scansforfield(sourceList.index(FLcal))) & set(onsourceScans))[0]
    if EQcal in sourceList: EQScan = list(set(msmd.scansforfield(sourceList.index(EQcal))) & set(onsourceScans))[0]
    #-------- Polarization setup
    spw = spwLists[band_index]; spwNum = len(spw); polNum = msmd.ncorrforpol(msmd.polidfordatadesc(spw[0]))
    pPol, cPol = [0,1], []  # parallel and cross pol
    PolList = ['X', 'Y']
    if polNum == 4: pPol, cPol = [0,3], [1,2]  # parallel and cross pol
    ppolNum, cpolNum = len(pPol), len(cPol)
    execfile(SCR_DIR + 'checkSEFD.py')
#
msmd.done()
