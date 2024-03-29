import sys
import pickle
import analysisUtils as au
import xml.etree.ElementTree as ET
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Grid.py')
msfile = wd + prefix + '.ms'
#
#-------- Check Correlator Type
BLCORR = True
tree = ET.parse(prefix + '/CorrelatorMode.xml')
root = tree.getroot()
correlatorName = root.find("row").find("correlatorName").text
if 'ACA' in correlatorName: BLCORR = False
#-------- Check Receivers
tree = ET.parse(prefix + '/Receiver.xml')
root = tree.getroot()
rowList = root.findall(".//row[frequencyBand]")
BandLists, BandList = [], []
for row in rowList:
    BandName = row.find("frequencyBand").text
    if 'ALMA' in BandName: BandLists = BandLists + [BandName.replace('ALMA_', '')]
#
BandList = unique(BandLists).tolist()
#-------- Tsys measurement
execfile(SCR_DIR + 'TsysCal.py')
class END(Exception):
    pass
#
#-------- Check Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
blNum = antNum* (antNum - 1) / 2
#-------- Check SPWs of atmCal
bpspwLists, bpscanLists, BandPA = [], [], []
msmd.open(msfile)
for band_index in range(NumBands):
    if atmBandNames == []:
        bpspwLists = bpspwLists + [np.array(atmSPWs)]
    else:
        bpspwLists  = bpspwLists  + [np.array(atmSPWs)[indexList( np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()]
    if 'spwFlag' in locals():
        flagIndex = indexList(np.array(spwFlag), np.array(bpspwLists[band_index]))
        for index in flagIndex: del bpspwLists[band_index][index]
    #
    bpscanLists = bpscanLists + [msmd.scansforspw(bpspwLists[band_index][0]).tolist()]
    BandPA = BandPA + [(BANDPA[int(UniqBands[band_index][3:5])] + 90.0)*pi/180.0]
#
msmd.close()
#-------- Check source list
print('---Checking source list')
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList); numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList))
msmd.open(msfile)
ONScans = sort(np.array(list(set(msmd.scansforintent("*CALIBRATE_POLARIZATION*")) | set(msmd.scansforintent("*CALIBRATE_AMPLI*")) | set(msmd.scansforintent("*CALIBRATE_BANDPASS*")) | set(msmd.scansforintent("*CALIBRATE_FLUX*")) | set(msmd.scansforintent("*CALIBRATE_PHASE*")) | set(msmd.scansforintent("*OBSERVE_CHECK*")) | set(msmd.scansforintent("*CALIBRATE_DELAY*")) )))
msmd.close()
msmd.done()
#-------- Loop for Bands
for band_index in range(NumBands):
    msmd.open(msfile)
    timeSum = 0
    bandName = UniqBands[band_index]; bandID = int(UniqBands[band_index][3:5])
    if Tau0Max[band_index] < 0.01: print('Failed in Trx/Tsys measurements...'); os.system('touch ' + prefix + '-' + UniqBands[band_index] + '-Flux.log'); continue
    if Tau0Max[band_index] > TauLimit[bandID]: print('Too high optical depth...'); os.system('touch ' + prefix + '-' + UniqBands[band_index] + '-Flux.log'); continue
    ONScan = np.array(bpscanLists[band_index])[indexList( ONScans, np.array(bpscanLists[band_index]))]
    ATMScan= np.array(atmscanLists[band_index])[indexList( ONScans, np.array(atmscanLists[band_index]))]
    onsourceScans, atmScans = ONScan.tolist(), ATMScan.tolist()
    scanNum, atmscanNum = len(onsourceScans), len(atmScans)
    SSOScanIndex = []
    StokesDic = dict(zip(sourceList, [[]]*numSource))   # Stokes parameters for each source
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
        print(text_sd)
        os.system(text_sd)
        fp = open('CalQU.data')
        lines = fp.readlines()
        fp.close()
        for eachLine in lines:
            sourceName = eachLine.split()[0]
            StokesDic[sourceName] = [float(eachLine.split()[1]), float(eachLine.split()[2]), float(eachLine.split()[3]), 0.0]
        #
    #
    for scan_index in range(scanNum):
        sourceID = msmd.sourceidforfield(msmd.fieldsforscan(onsourceScans[scan_index])[0])
        #sourceIDscan.append( msmd.sourceidforfield(msmd.fieldsforscan(onsourceScans[scan_index])[0]))
        sourceIDscan = sourceIDscan + [sourceID]
        SunAngle = SunAngleSourceList[sourceID]
        interval, timeStamp = GetTimerecord(msfile, 0, 0, bpspwLists[band_index][0], onsourceScans[scan_index])
        timeSum += len(timeStamp)
        try:
            AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, 0, AZ, EL)
        except:
            AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, 1, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; dPA = np.std(np.sin(PA)) #dPA = abs(np.sin(max(PA) - min(PA)))
        OnAZ.append(np.median(AzScan)); OnEL.append(np.median(ElScan)); OnPA.append(np.median(PA))
        refTime = refTime + [np.median(timeStamp)]
        sourceName = sourceList[sourceID]
        if len(StokesDic[sourceName]) == 4:
            CS, SN = np.cos(2.0* OnPA[scan_index]), np.sin(2.0* OnPA[scan_index])
            QCpUS = StokesDic[sourceName][1]*CS + StokesDic[sourceName][2]*SN   # Qcos + Usin
            UCmQS = StokesDic[sourceName][2]*CS - StokesDic[sourceName][1]*SN   # Ucos - Qsin
            if QUMODEL:
                BPquality = BPquality + [1000.0* StokesDic[sourceName][0]* abs(UCmQS)* np.sin(OnEL[scan_index] - 0.5*ELshadow) / np.sqrt(StokesDic[sourceName][0] + 0.5)]
            else:
                BPquality = BPquality + [1000.0* abs(UCmQS)* dPA* np.sin(OnEL[scan_index] - 0.5*ELshadow) / np.sqrt(StokesDic[sourceName][0])]
            #
            EQquality = EQquality + [(StokesDic[sourceName][0] - abs(QCpUS))**2 * np.sin(OnEL[scan_index] - ELshadow) ]
        else:
            QCpUS, UCmQS = 0.0, 0.0
            BPquality = BPquality + [-1]
            EQquality = EQquality + [-1]
        #
        if SunAngle < SunAngleThresh: BPquality[scan_index], EQquality[scan_index] = -1, -1
        print('Scan%02d : %10s AZ=%6.1f EL=%4.1f SA=%5.1f PA=%6.1f dPA=%5.2f pRes=%5.2f BPquality=%7.4f EQquality=%6.0f') % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0*OnAZ[scan_index]/np.pi, 180.0*OnEL[scan_index]/np.pi, SunAngle, 180.0*OnPA[scan_index]/np.pi, 180.0*dPA/np.pi, UCmQS, BPquality[scan_index], EQquality[scan_index]) 
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
    BPcalText = 'Use %s [Scan%d EL=%4.1f deg] %s as Bandpass Calibrator' % (BPcal, BPScan, 180.0* OnEL[onsourceScans.index(BPScan)]/np.pi, timeLabelBP); print(BPcalText)
    EQcalText = 'Use %s [Scan%d EL=%4.1f deg] %s as Gain Equalizer' % (EQcal, EQScan, 180.0* OnEL[onsourceScans.index(EQScan)]/np.pi, timeLabelEQ); print(EQcalText)
    #-------- Check Baseline Limit
    if 'uvLimit' in locals(): antFlag = angFlagBL(msfile, uvLimit, bpspwLists[band_index][0], BPScan, antFlag)
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
        #execfile(SCR_DIR + 'checkSEFDStokes.py')
        #execfile(SCR_DIR + 'aprioriStokes.py')
        if Apriori:
            try:
                execfile(SCR_DIR + 'aprioriStokes.py')
            except:
                print('  --A priori flux calibration falied.')
        else:
            try:
                execfile(SCR_DIR + 'checkSEFDStokes.py')
            except:
                print('  --SSO-based flux calibration falied. Switch to a priori (SEFD) calibration.')
                try:
                    execfile(SCR_DIR + 'aprioriStokes.py')
                except:
                    print('  --A priori flux calibration falied.')
            #
        #
    #
    if polNum == 2:
        cPol = [0,1], []; ppolNum, cpolNum = len(pPol), len(cPol)
        if not Apriori:
            try:
                execfile(SCR_DIR + 'checkSEFD.py')
            except:
                print(' --SSO-based flux calibration falied. Switch to a priori (SEFD) calibration.')
                execfile(SCR_DIR + 'aprioriFlux.py')
        else:
            execfile(SCR_DIR + 'aprioriFlux.py')
    #
#---- end of band loop
del msfile, UniqBands, useAnt, atmspwLists, atmSPWs
if 'spwFlag' in locals(): del spwFlag
if 'BPScans' in locals(): del BPScans
if 'EQScans' in locals(): del EQScans
if 'antFlag' in locals(): del antFlag
if 'flagAnt' in locals(): del flagAnt
