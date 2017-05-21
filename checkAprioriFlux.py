import sys
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Grid.py')
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
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList); numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList))
#-------- Check MJD for Ambient Load
print '---Checking time for ambient and hot load'
timeOFF, timeAMB, timeHOT = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT")
if len(timeAMB) == 0:
    timeON  = msmd.timesforintent("CALIBRATE_ATMOSPHERE#ON_SOURCE")
    timeAMB = timeON[(np.where( np.diff(timeON) > 20.0)[0] - 8).tolist()]
    timeHOT = timeON[(np.where( np.diff(timeON) > 20.0)[0] - 2).tolist()]
#
#-------- Scan Intents
try:
    FCScans = np.append(msmd.scansforintent("CALIBRATE_FLUX#ON_SOURCE"), msmd.scansforintent("OBSERVE_CHECK_SOURCE*"))
except:
    FCScans = np.append(msmd.scansforintent("CALIBRATE_AMPLI#ON_SOURCE"), msmd.scansforintent("OBSERVE_CHECK_SOURCE*"))
try:
    ONScans = msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE")
except:
    ONScans = FCScans
try:
    BPScans = msmd.scansforintent("CALIBRATE_BANDPASS#ON_SOURCE")
except:
    BPScans = ONScans
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
    azelTime_index = np.where( AntID == 0 )[0].tolist()
    azel = np.r_[AZ[azelTime_index], EL[azelTime_index]].reshape(2, len(azelTime_index))
    OnAZ, OnEL, OnPA, BPquality, EQquality, PQquality, sourceIDscan, FLscore, refTime = [], [], [], [], [], [], [], np.zeros(scanNum), []
    for scan_index in range(scanNum):
        sourceIDscan.append( msmd.sourceidforfield(msmd.fieldsforscan(onsourceScans[scan_index])[0]))
        interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spwLists[band_index][0], onsourceScans[scan_index])
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, 0, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; np.std(np.sin(PA))
        OnAZ.append(np.median(AzScan)); OnEL.append(np.median(ElScan)); OnPA.append(np.median(PA))
        refTime = refTime + [np.median(timeStamp)]
        catalogIQUV = np.array([catalogStokesI.get(sourceList[sourceIDscan[scan_index]], 0.0), catalogStokesQ.get(sourceList[sourceIDscan[scan_index]], 0.0), catalogStokesU.get(sourceList[sourceIDscan[scan_index]], 0.0), 0.0])
        CS, SN = np.cos(2.0* (OnPA[scan_index] + BandPA[band_index])), np.sin(2.0* (OnPA[scan_index] + BandPA[band_index]))
        QCpUS = catalogIQUV[1]*CS + catalogIQUV[2]*SN   # Qcos + Usin
        UCmQS = catalogIQUV[2]*CS - catalogIQUV[1]*SN   # Ucos - Qsin
        BPquality = BPquality + [1000.0* abs(UCmQS)* catalogIQUV[0]* dPA * np.sin(OnEL[scan_index])]
        EQquality = EQquality + [catalogIQUV[0]* np.sin(OnEL[scan_index] - ELshadow) / (0.001 + QCpUS**2)]
        print 'Scan%02d : %10s AZ=%6.1f EL=%4.1f PA=%6.1f dPA=%5.2f pRes=%5.2f BPquality=%7.4f EQquality=%6.0f' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0*OnAZ[scan_index]/np.pi, 180.0*OnEL[scan_index]/np.pi, 180.0*OnPA[scan_index]/np.pi, 180.0*dPA/np.pi, UCmQS, BPquality[scan_index], EQquality[scan_index])
        if sourceIDscan[scan_index] in SSOList:
            FLscore[scan_index] = np.exp(np.log(math.sin(OnEL[scan_index])-0.34))* SSOscore[bandID][SSOCatalog.index(sourceList[sourceIDscan[scan_index]])]
        #
    #
    BPcal = sourceList[sourceIDscan[np.argmax(BPquality)]]; BPScan = onsourceScans[np.argmax(BPquality)]; timeLabelBP = qa.time('%fs' % (refTime[np.argmax(BPquality)]), form='ymd')[0]
    if max(FLscore) > 0.1:
        EQcal = sourceList[sourceIDscan[np.argmax(FLscore)]]; EQScan = onsourceScans[np.argmax(FLscore)]; timeLabelEQ = qa.time('%fs' % (refTime[np.argmax(FLscore)]), form='ymd')[0]
    else: 
        EQcal = sourceList[sourceIDscan[np.argmax(EQquality)]]; EQScan = onsourceScans[np.argmax(EQquality)]; timeLabelEQ = qa.time('%fs' % (refTime[np.argmax(EQquality)]), form='ymd')[0]
    #
    BPcalText = 'Use %s [EL = %4.1f deg] %s as Bandpass Calibrator' % (BPcal, 180.0* OnEL[onsourceScans.index(BPScan)]/np.pi, timeLabelBP); print BPcalText
    EQcalText = 'Use %s [EL = %4.1f deg] %s as Gain Equalizer' % (EQcal, 180.0* OnEL[onsourceScans.index(EQScan)]/np.pi, timeLabelEQ); print EQcalText
    #-------- Polarization setup 
    spw = spwLists[band_index]; spwNum = len(spw); polNum = msmd.ncorrforpol(msmd.polidfordatadesc(spw[0]))
    if polNum == 4: pPol, cPol = [0,3], [1,2]   # Full polarizations
    else:   pPol, cPol = [0,1], []              # Only parallel polarizations
    ppolNum, cpolNum = len(pPol), len(cPol)
    execfile(SCR_DIR + 'aprioriFlux.py')
#
del msfile, UniqBands
msmd.done()
