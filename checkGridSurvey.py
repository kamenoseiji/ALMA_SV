import sys
import re
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Grid.py')
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
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList); numSource = len(sourceList)
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
#-------- Check D-term files
polNum = msmd.ncorrforpol(msmd.polidfordatadesc(spw[0]))
if polNum == 4:
    print '---Checking D-term files'
    DantList, noDlist = [], []
    Dpath = SCR_DIR + 'DtermB' + `int(UniqBands[band_index][3:5])` + '/'
    for ant_index in UseAnt:
        Dfile = Dpath + 'B' + `int(UniqBands[band_index][3:5])` + '-SPW0-' + antList[ant_index] + '.DSpec.npy'
        if os.path.exists(Dfile): DantList += [ant_index]
        else: noDlist += [ant_index]
    #
    DantNum, noDantNum = len(DantList), len(noDlist)
    print 'Antennas with D-term file (%d):' % (DantNum),
    for ant_index in DantList: print '%s ' % antList[ant_index],
    print ''
    if noDantNum > 0:
        print 'Antennas without D-term file (%d) : ' % (noDantNum),
        for ant_index in noDlist: print '%s ' % antList[ant_index],
        sys.exit(' Run DtermTransfer first!!')
    #
#
#-------- Antenna and baseline mapping
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
antDia = np.ones(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
#-------- Check Scans of BandPass, EQualization, and FluxScaling
try:
    FCScans = np.append(msmd.scansforintent("CALIBRATE_FLUX#ON_SOURCE"), msmd.scansforintent("OBSERVE_CHECK_SOURCE*"))
except:
    FCScans = np.append(msmd.scansforintent("CALIBRATE_AMPLI#ON_SOURCE"), msmd.scansforintent("OBSERVE_CHECK_SOURCE*"))
#
ONScans = msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE")
try:
    BPScans = msmd.scansforintent("CALIBRATE_BANDPASS#ON_SOURCE")
except:
    BPScans = ONScans
print '---SPWs and Scans for each receiver band'
msmd.done()
#for band_index in range(1):
for band_index in range(NumBands):
    bandID = int(UniqBands[band_index][3:5])-1
    msmd.open(msfile)
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
    refTime, OnAZ, OnEL, OnPA, BPquality, EQquality, sourceIDscan, FLscore = [], [], [], [], [], [], [], np.zeros(scanNum)
    for scan_index in range(scanNum):
        sourceIDscan.append( msmd.sourceidforfield(msmd.fieldsforscan(onsourceScans[scan_index])[0]))
        interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spwLists[band_index][0], onsourceScans[scan_index])
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; dPA = abs(np.sin(max(PA) - min(PA)))
        OnAZ.append(np.median(AzScan)); OnEL.append(np.median(ElScan)); OnPA.append(np.median(PA))
        refTime = refTime + [np.median(timeStamp)]
        catalogIQUV = np.array([catalogStokesI.get(sourceList[sourceIDscan[scan_index]], 0.0), catalogStokesQ.get(sourceList[sourceIDscan[scan_index]], 0.0), catalogStokesU.get(sourceList[sourceIDscan[scan_index]], 0.0), 0.0])
        CS, SN = np.cos(2.0* (OnPA[scan_index] + BandPA[band_index])), np.sin(2.0* (OnPA[scan_index] + BandPA[band_index]))
        QCpUS = catalogIQUV[1]*CS + catalogIQUV[2]*SN   # Qcos + Usin
        UCmQS = catalogIQUV[2]*CS - catalogIQUV[1]*SN   # Ucos - Qsin
        BPquality = BPquality + [1000.0* abs(UCmQS)* catalogIQUV[0]* dPA * np.sin(OnEL[scan_index])]
        EQquality = EQquality + [catalogIQUV[0]* np.sin(OnEL[scan_index] - ELshadow) / (0.001 + QCpUS**2)]
        print 'Scan%02d : %10s AZ=%6.1f EL=%4.1f PA=%6.1f dPA=%5.2f pRes=%5.2f BPQuality=%7.4f EQquality=%6.0f' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0*OnAZ[scan_index]/np.pi, 180.0*OnEL[scan_index]/np.pi, 180.0*OnPA[scan_index]/np.pi, 180.0*dPA/np.pi, UCmQS, BPquality[scan_index], EQquality[scan_index]) 
        if sourceIDscan[scan_index] in SSOList: FLscore[scan_index] = np.exp(np.log(math.sin(OnEL[scan_index])-0.34))* SSOscore[bandID][SSOCatalog.index(sourceList[sourceIDscan[scan_index]])]
    #
    if not 'FLcal' in locals(): FLcal = sourceList[sourceIDscan[np.argmax(FLscore)]]
    if FLcal in sourceList: FCScan = list(set(msmd.scansforfield(FLcal)) & set(onsourceScans))[0]; FLsel = FLcal
    else: FCScan = onsourceScans[np.argmax(FLscore)]; FLsel = sourceList[sourceIDscan[np.argmax(FLscore)]]
    FLScaleText = 'Use %s [EL = %4.1f deg] as Flux Scaler' % (FLsel, 180.0* OnEL[onsourceScans.index(FCScan)]/np.pi); print FLScaleText
    BPcal = sourceList[sourceIDscan[np.argmax(BPquality)]]; BPScan = onsourceScans[np.argmax(BPquality)]
    EQcal = sourceList[sourceIDscan[np.argmax(EQquality)]]; EQScan = onsourceScans[np.argmax(EQquality)]
    #-------- SSO in observed source list
    BandSSOList = list( set(SSOList) & set(sourceIDscan) )
    BPcalText = 'Use %s [EL = %4.1f deg] as Bandpass Calibrator' % (BPcal, 180.0* OnEL[onsourceScans.index(BPScan)]/np.pi); print BPcalText
    EQcalText = 'Use %s [EL = %4.1f deg] as Gain Equalizer' % (EQcal, 180.0* OnEL[onsourceScans.index(EQScan)]/np.pi); print EQcalText
    #-------- Polarization setup
    spw = spwLists[band_index]; spwNum = len(spw)
    pPol, cPol = [0,1], []  # parallel and cross pol
    PolList = ['X', 'Y']
    if polNum == 4: pPol, cPol = [0,3], [1,2]  # parallel and cross pol
    ppolNum, cpolNum = len(pPol), len(cPol)
    msmd.done()
    if polNum == 2: pPol, cPol = [0,1], []   ; ppolNum, cpolNum = len(pPol), len(cPol); execfile(SCR_DIR + 'checkSEFD.py')
    if polNum == 4:
        pPol, cPol = [0,3], [1,2]; ppolNum, cpolNum = len(pPol), len(cPol)
        execfile(SCR_DIR + 'checkSEFDStokes.py')
    #
#
