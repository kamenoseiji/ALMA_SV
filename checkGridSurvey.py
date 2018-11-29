import sys
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Grid.py')
msfile = wd + prefix + '.ms'
execfile(SCR_DIR + 'TsysCal.py')
class END(Exception):
    pass
#
#-------- Get Bandpass SPWs
def GetBPcalSPWs(msfile):
    msmd.open(msfile)
    bpSPWs  = msmd.spwsforintent("CALIBRATE_BANDPASS*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_FLUX*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_DELAY*").tolist(); bpSPWs.sort()
    msmd.close()
    return bpSPWs
#
#-------- Check Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
blNum = antNum* (antNum - 1) / 2
#-------- Check SPWs of atmCal
print '---Checking spectral windows with atmCal for ' + prefix
atmSPWs = GetAtmSPWs(msfile)
bpSPWs  = GetBPcalSPWs(msfile)
msmd.open(msfile)
atmspwNames, bpspwNames = msmd.namesforspws(atmSPWs), msmd.namesforspws(bpSPWs)
bpSPWs = np.array(bpSPWs)[indexList(np.array(atmspwNames), np.array(bpspwNames))].tolist(); bpspwNames = msmd.namesforspws(bpSPWs)
atmBandNames, atmPattern = [], r'RB_..'
for spwName in atmspwNames : atmBandNames = atmBandNames + re.findall(atmPattern, spwName)
UniqBands = unique(atmBandNames).tolist(); NumBands = len(UniqBands)
atmspwLists, bpspwLists, atmscanLists, bpscanLists, BandPA = [], [], [], [], []
for band_index in range(NumBands):
    BandPA = BandPA + [(BANDPA[int(UniqBands[band_index][3:5])] + 90.0)*pi/180.0]
    atmspwLists = atmspwLists + [np.array(atmSPWs)[indexList(np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()]
    bpspwLists  = bpspwLists  + [np.array(bpSPWs)[indexList( np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()]
    atmscanLists= atmscanLists+ [msmd.scansforspw(atmspwLists[band_index][0]).tolist()]
    bpscanLists = bpscanLists + [msmd.scansforspw(bpspwLists[band_index][0]).tolist()]
    print ' ',
    print UniqBands[band_index] + ': atmSPW=' + `atmspwLists[band_index]` + ' bpSPW=' + `bpspwLists[band_index]`
#
#-------- Check source list
print '---Checking source list'
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList); numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList))
try: ONScans = msmd.scansforintent("CALIBRATE_BANDPASS")
try: ONScans = np.array(set(ONScans) + set(msmd.scansforintent("CALIBRATE_FLUX")))
try: ONScans = np.array(set(ONScans) + set(msmd.scansforintent("CALIBRATE_PHASE")))
msmd.close()
msmd.done()
#-------- Loop for Bands
for band_index in range(NumBands):
    msmd.open(msfile)
    bandName = UniqBands[band_index]; bandID = int(UniqBands[band_index][3:5])
    ONScan = np.array(bpscanLists[band_index])[indexList( ONScans, np.array(bpscanLists[band_index]))]
    ATMScan= np.array(atmscanLists[band_index])[indexList( ONScans, np.array(atmscanLists[band_index]))]
    onsourceScans, atmScans = ONScan.tolist(), ATMScan.tolist()
    scanNum, atmscanNum = len(onsourceScans), len(atmScans)
    SSOScanIndex = []
    #-------- Check AZEL
    azelTime, AntID, AZ, EL = GetAzEl(msfile)
    azelTime_index = np.where( AntID == 0 )[0].tolist() 
    azel = np.r_[AZ[azelTime_index], EL[azelTime_index]].reshape(2, len(azelTime_index))
    OnAZ, OnEL, OnPA, BPquality, EQquality, PQquality, sourceIDscan, FLscore, refTime = [], [], [], [], [], [], [], np.zeros(scanNum), []
    #-------- Check QU catalog
    if QUMODEL: # A priori QU model from Rdata
        os.system('rm -rf CalQU.data')
        text_sd = R_DIR + 'Rscript %spolQuery.R -D%s -F%f' % (SCR_DIR, qa.time('%fs' % (azelTime[0]), form='ymd')[0], BANDFQ[bandID])
        for source in sourceList: text_sd = text_sd + ' ' + source
        os.system(text_sd)
        fp = open('CalQU.data')
        lines = fp.readlines()
        fp.close()
        for eachLine in lines:
            catalogStokesI[eachLine.split()[0]] = float(eachLine.split()[1])
            catalogStokesQ[eachLine.split()[0]] = float(eachLine.split()[2])
            catalogStokesU[eachLine.split()[0]] = float(eachLine.split()[3])
        #
    #
    for scan_index in range(scanNum):
        sourceIDscan.append( msmd.sourceidforfield(msmd.fieldsforscan(onsourceScans[scan_index])[0]))
        interval, timeStamp = GetTimerecord(msfile, 0, 0, bpspwLists[band_index][0], onsourceScans[scan_index])
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, 0, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; dPA = np.std(np.sin(PA)) #dPA = abs(np.sin(max(PA) - min(PA)))
        OnAZ.append(np.median(AzScan)); OnEL.append(np.median(ElScan)); OnPA.append(np.median(PA))
        refTime = refTime + [np.median(timeStamp)]
        catalogIQUV = np.array([catalogStokesI.get(sourceList[sourceIDscan[scan_index]], 0.0), catalogStokesQ.get(sourceList[sourceIDscan[scan_index]], 0.0), catalogStokesU.get(sourceList[sourceIDscan[scan_index]], 0.0), 0.0])
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
            BPquality = BPquality + [-9999.9]
            EQquality = EQquality + [-9999.9]
        #
        print 'Scan%02d : %10s AZ=%6.1f EL=%4.1f PA=%6.1f dPA=%5.2f pRes=%5.2f BPquality=%7.4f EQquality=%6.0f' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0*OnAZ[scan_index]/np.pi, 180.0*OnEL[scan_index]/np.pi, 180.0*OnPA[scan_index]/np.pi, 180.0*dPA/np.pi, UCmQS, BPquality[scan_index], EQquality[scan_index]) 
        if sourceIDscan[scan_index] in SSOList: FLscore[scan_index] = np.exp(np.log(math.sin(OnEL[scan_index])-0.34))* SSOscore[bandID-1][SSOCatalog.index(sourceList[sourceIDscan[scan_index]])]
    #
    #-------- Select Bandpass Calibrator
    if 'BPScans' not in locals():
        BPscanIndex = np.argmax(BPquality)
    else:
        if len(BPScans) < NumBands:
            BPscanIndex = np.argmax(BPquality)
        else:
            BPscanIndex = onsourceScans.index(BPScans[band_index])
        #
    #
    BPScan = onsourceScans[BPscanIndex]; BPcal = sourceList[sourceIDscan[BPscanIndex]]; timeLabelBP = qa.time('%fs' % (refTime[BPscanIndex]), form='ymd')[0]
    BPEL = OnEL[onsourceScans.index(BPScan)]
    #-------- Select Equalization Calibrator
    if 'EQScans' not in locals():
        EQscanIndex = np.argmax(EQquality)
    else:
        if len(EQScans) < NumBands:
            EQscanIndex = np.argmax(EQquality)
        else:
            EQscanIndex = onsourceScans.index(EQScans[band_index])
        #
    #
    EQScan = onsourceScans[EQscanIndex]; EQcal = sourceList[sourceIDscan[EQscanIndex]]; timeLabelEQ = qa.time('%fs' % (refTime[EQscanIndex]), form='ymd')[0]
    BPcalText = 'Use %s [Scan%d EL=%4.1f deg] %s as Bandpass Calibrator' % (BPcal, BPScan, 180.0* OnEL[onsourceScans.index(BPScan)]/np.pi, timeLabelBP); print BPcalText
    EQcalText = 'Use %s [Scan%d EL=%4.1f deg] %s as Gain Equalizer' % (EQcal, EQScan, 180.0* OnEL[onsourceScans.index(EQScan)]/np.pi, timeLabelEQ); print EQcalText
    #-------- SSO in observed source list
    BandSSOList = list( set(SSOList) & set(sourceIDscan) )
    if len(BandSSOList) == 0: Apriori = True
    #-------- Polarization setup
    atmspw = atmspwLists[band_index]; spwNum = len(atmspw)
    scnspw = bpspwLists[band_index]; scnspwNum = len(scnspw)
    polNum = msmd.ncorrforpol(msmd.polidfordatadesc(scnspw[0]))
    PolList = ['X', 'Y']
    msmd.done()
    if polNum == 4:
        pPol, cPol = [0,3], [1,2];  ppolNum, cpolNum = len(pPol), len(cPol)
        if not Apriori:
            try:
                execfile(SCR_DIR + 'checkSEFDStokes.py')
            except:
                print '  --SSO-based flux calibration falied. Switch to a priori (SEFD) calibration.'
                execfile(SCR_DIR + 'aprioriStokes.py')
            #
        else:
            execfile(SCR_DIR + 'aprioriStokes.py')
        #
    #
    if polNum == 2:
        cPol = [0,1], []; ppolNum, cpolNum = len(pPol), len(cPol)
        if not Apriori:
            try:
                execfile(SCR_DIR + 'checkSEFD.py')
            except:
                print  ' --SSO-based flux calibration falied. Switch to a priori (SEFD) calibration.'
                execfile(SCR_DIR + 'aprioriFlux.py')
        else:
            execfile(SCR_DIR + 'aprioriFlux.py')
    #
#
del msfile, UniqBands, UseAnt
if 'flagAnt' in locals(): del flagAnt
if 'BPScans' in locals(): del BPScans
if 'EQScans' in locals(): del EQScans
