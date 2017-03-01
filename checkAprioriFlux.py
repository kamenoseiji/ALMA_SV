import sys
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
SSOCatalog = ['Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus']
msfile = prefix + '.ms'
#-------- Check Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
blNum = antNum* (antNum - 1) / 2
#-------- Check SPWs
print '---Checking spectral windows for ' + prefix
msmd.open(msfile)
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
sourceList, posList = GetSourceList(msfile); numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList))
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
    sys.exit(' Run DtermTransfer, or flag these antennas out.')
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
antDia = np.ones([UseAntNum])
for ant_index in range(UseAntNum): antDia[ant_index] = msmd.antennadiameter(antList[antMap[ant_index]])['value']
#-------- Scan Intents
ONScans = msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE")
BPScans = msmd.scansforintent("CALIBRATE_BANDPASS#ON_SOURCE")
try:
    FCScans = msmd.scansforintent("CALIBRATE_FLUX#ON_SOURCE")
except:
    FCScans = msmd.scansforintent("CALIBRATE_AMPLI#ON_SOURCE")
#
PolList = ['X', 'Y']
#-------- Loop for Bands
for band_index in range(NumBands):
    bandID = int(UniqBands[band_index][3:5])-1
    ONScan = BandScans[band_index][indexList( ONScans, BandScans[band_index] )]
    BPScan = BandScans[band_index][indexList( BPScans, BandScans[band_index] )][0]
    FCScan = BandScans[band_index][indexList( FCScans, BandScans[band_index] )]
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
    #
    if BPcal in sourceList: BPScan = list(set(msmd.scansforfield(BPcal)) & set(onsourceScans))[0]
    #-------- Polarization setup 
    spw = spwLists[band_index]; spwNum = len(spw); polNum = msmd.ncorrforpol(msmd.polidfordatadesc(spw[0]))
    if polNum == 4: pPol, cPol = [0,3], [1,2]   # Full polarizations
    else:   pPol, cPol = [0,1], []              # Only parallel polarizations
    ppolNum, cpolNum = len(pPol), len(cPol)
    execfile(SCR_DIR + 'aprioriFlux.py')
#
msmd.done()
#del msfile, UniqBands
