import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Grid.py')
class END(Exception):
    pass
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
if not 'spwList' in locals():
    spwList = list(set(msmd.tdmspws()) & set(msmd.spwsforintent("CALIBRATE_ATMOSPHERE*"))); spwList.sort()
#
if not 'spwNames' in locals(): spwNames = msmd.namesforspws(spwList)
if len(spwList) < 1: sys.exit('No available SPWs found.\n')
#
BandNames, pattern = [], r'RB_..'
for spwName in spwNames: BandNames = BandNames + re.findall(pattern, spwName)
UniqBands = unique(BandNames).tolist(); NumBands = len(UniqBands)
spwLists, BandScans, BandPA = [], [], []
for band_index in range(NumBands):
    BandPA = BandPA + [(BANDPA[int(UniqBands[band_index][3:5])] + 90.0)*pi/180.0]
    spwLists = spwLists + [np.array(spwList)[indexList( np.array([UniqBands[band_index]]), np.array(BandNames))].tolist()]
    BandScans = BandScans + [msmd.scansforspw(spwLists[band_index][0])]
    print ' ',
    print UniqBands[band_index] + ': SPW=' + `spwLists[band_index]`
#
#-------- Check source list
print '---Checking source list'
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList); numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList) ) # Find Solar System Objects in the source list
try:
    if len(SSOList) == 0: raise END
except END:
    print '  No Solar System Object was observed.'
#
QSOList   = list(set(range(numSource)) - set(SSOList))
for source in sourceList: print source,
print ''; print '  Solar System Objects:',
for index in SSOList: print sourceList[index],
print ''
#-------- Check Scans of BandPass, EQualization, and FluxScaling
polNum = msmd.ncorrforpol(msmd.polidfordatadesc(spwList[0]))
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
print '---SPWs and Scans for each receiver band'
msmd.done()
#for band_index in range(1):
for band_index in range(NumBands):
    bandName = UniqBands[band_index]; bandID = int(UniqBands[band_index][3:5])
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
    azelTime_index = np.where( AntID == 0 )[0].tolist() 
    azel = np.r_[AZ[azelTime_index], EL[azelTime_index]].reshape(2, len(azelTime_index))
    refTime, OnAZ, OnEL, OnPA, BPquality, EQquality, sourceIDscan, FLscore = [], [], [], [], [], [], [], np.zeros(scanNum)
    #-------- Check QU catalog
    if QUMODEL:
        os.system('rm -rf CalQU.data')
        text_sd = R_DIR + 'Rscript %spolQuery.R -D%s -F%f' % (SCR_DIR, qa.time('%fs' % (azelTime[0]), form='ymd')[0], BANDFQ[bandID])
        for source in sourceList: text_sd = text_sd + ' ' + source
    #
    os.system(text_sd)
    fp = open('CalQU.data')
    lines = fp.readlines()
    fp.close()
    for eachLine in lines:
        catalogStokesI[eachLine.split()[0]] = float(eachLine.split()[1])
        catalogStokesQ[eachLine.split()[0]] = float(eachLine.split()[2])
        catalogStokesU[eachLine.split()[0]] = float(eachLine.split()[3])
    #
    for scan_index in range(scanNum):
        sourceIDscan.append( msmd.sourceidforfield(msmd.fieldsforscan(onsourceScans[scan_index])[0]))
        interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spwLists[band_index][0], onsourceScans[scan_index])
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, 0, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; dPA = np.std(np.sin(PA)) #dPA = abs(np.sin(max(PA) - min(PA)))
        OnAZ.append(np.median(AzScan)); OnEL.append(np.median(ElScan)); OnPA.append(np.median(PA))
        refTime = refTime + [np.median(timeStamp)]
        catalogIQUV = np.array([catalogStokesI.get(sourceList[sourceIDscan[scan_index]], 0.0), catalogStokesQ.get(sourceList[sourceIDscan[scan_index]], 0.0), catalogStokesU.get(sourceList[sourceIDscan[scan_index]], 0.0), 0.0])
        #if sourceIDscan[scan_index] in SSOList:
        #    BPquality = BPquality + [-1.0]
        #    EQquality = EQquality + [-1.0]
        #    QCpUS, UCmQS = 0.0, 0.0
        #else:
        if catalogIQUV[0] > 0.1:
            CS, SN = np.cos(2.0* OnPA[scan_index]), np.sin(2.0* OnPA[scan_index])
            QCpUS = catalogIQUV[1]*CS + catalogIQUV[2]*SN   # Qcos + Usin
            UCmQS = catalogIQUV[2]*CS - catalogIQUV[1]*SN   # Ucos - Qsin
            if QUMODEL:
                BPquality = BPquality + [1000.0* abs(UCmQS)* np.sin(OnEL[scan_index] - 0.5*ELshadow) / np.sqrt(catalogIQUV[0])]
            else:
                BPquality = BPquality + [1000.0* abs(UCmQS)* dPA* np.sin(OnEL[scan_index] - 0.5*ELshadow) / np.sqrt(catalogIQUV[0])]
            #
            EQquality = EQquality + [catalogIQUV[0]**2 * np.sin(OnEL[scan_index] - ELshadow) / (1.0e-4 + abs(QCpUS))]
        else:
            QCpUS, UCmQS = 0.0, 0.0
            BPquality = BPquality + [-1.0]
            EQquality = EQquality + [-1.0]
        #
        print 'Scan%02d : %10s AZ=%6.1f EL=%4.1f PA=%6.1f dPA=%5.2f pRes=%5.2f BPquality=%7.4f EQquality=%6.0f' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0*OnAZ[scan_index]/np.pi, 180.0*OnEL[scan_index]/np.pi, 180.0*OnPA[scan_index]/np.pi, 180.0*dPA/np.pi, UCmQS, BPquality[scan_index], EQquality[scan_index]) 
        if sourceIDscan[scan_index] in SSOList: FLscore[scan_index] = np.exp(np.log(math.sin(OnEL[scan_index])-0.34))* SSOscore[bandID-1][SSOCatalog.index(sourceList[sourceIDscan[scan_index]])]
    #
    if not 'FLcal' in locals(): FLcal = sourceList[sourceIDscan[np.argmax(FLscore)]]
    if FLcal in sourceList: FCScan = list(set(msmd.scansforfield(FLcal)) & set(onsourceScans))[0]; FLsel = FLcal
    else: FCScan = onsourceScans[np.argmax(FLscore)]; FLsel = sourceList[sourceIDscan[np.argmax(FLscore)]]
    FLEL = 180.0* OnEL[onsourceScans.index(FCScan)]/np.pi
    #FLScaleText = 'Use %s [EL = %4.1f deg] as Flux Scaler' % (FLsel, FLEL); print FLScaleText
    if 'BPcal' in locals():
        BPScan = onsourceScans[sourceIDscan.index(sourceList.index(BPcal))]
    else:
        BPcal = sourceList[sourceIDscan[np.argmax(BPquality)]]; BPScan = onsourceScans[np.argmax(BPquality)]
    #
    EQcal = sourceList[sourceIDscan[np.argmax(EQquality)]]; EQScan = onsourceScans[np.argmax(EQquality)]
    #-------- SSO in observed source list
    BandSSOList = list( set(SSOList) & set(sourceIDscan) )
    BPcalText = 'Use %s [EL = %4.1f deg] as Bandpass Calibrator' % (BPcal, 180.0* OnEL[onsourceScans.index(BPScan)]/np.pi); print BPcalText
    EQcalText = 'Use %s [EL = %4.1f deg] as Gain Equalizer' % (EQcal, 180.0* OnEL[onsourceScans.index(EQScan)]/np.pi); print EQcalText
    #-------- Polarization setup
    spwList = spwLists[band_index]; spwNum = len(spwList)
    pPol, cPol = [0,1], []  # parallel and cross pol
    PolList = ['X', 'Y']
    if polNum == 4: pPol, cPol = [0,3], [1,2]  # parallel and cross pol
    ppolNum, cpolNum = len(pPol), len(cPol)
    msmd.done()
    if Apriori:
        if polNum == 2: pPol, cPol = [0,1], []   ; ppolNum, cpolNum = len(pPol), len(cPol); execfile(SCR_DIR + 'aprioriFlux.py')
        if polNum == 4: pPol, cPol = [0,3], [1,2]; ppolNum, cpolNum = len(pPol), len(cPol); execfile(SCR_DIR + 'aprioriStokes.py')
    else:
        try:
            if polNum == 2: pPol, cPol = [0,1], []   ; ppolNum, cpolNum = len(pPol), len(cPol); execfile(SCR_DIR + 'checkSEFD.py')
            if polNum == 4: pPol, cPol = [0,3], [1,2]; ppolNum, cpolNum = len(pPol), len(cPol); execfile(SCR_DIR + 'checkSEFDStokes.py')
        except:
            if polNum == 2: pPol, cPol = [0,1], []   ; ppolNum, cpolNum = len(pPol), len(cPol); execfile(SCR_DIR + 'aprioriFlux.py')
            if polNum == 4: pPol, cPol = [0,3], [1,2]; ppolNum, cpolNum = len(pPol), len(cPol); execfile(SCR_DIR + 'aprioriStokes.py')
        #
    #
    del(FLcal)
    del(BPcal)
#
del(spwList)
