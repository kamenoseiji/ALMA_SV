#---- Script for Band-3 Astroholograpy Data
import sys
from scipy import stats
import matplotlib.pyplot as plt
execfile(SCR_DIR + 'interferometry.py')
#-------- Initial Settings
if 'SNR_THRESH' not in locals(): SNR_THRESH = 1.0
msfile = wd + prefix + '.ms'; msmd.open(msfile)
antList = GetAntName(msfile)
antNum = len(antList)
blNum = antNum* (antNum - 1) / 2
spwName = msmd.namesforspws(spw)[0]
BandName = re.findall(r'RB_..', spwName)[0]; BandID = int(BandName[3:5])
#-------- Array Configuration
print '---Checking array configuration'
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
text_sd = '  Usable antennas (%d) : ' % (len(UseAnt))
for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
print text_sd
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
timeStamp, UVW = GetUVW(msfile, spw, msmd.scansforspw(spw)[0])
uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
if 'refant' in locals():    refantID = indexList(np.array([refant]), antList[UseAnt])[0]
else: refantID = bestRefant(uvDist)
print '  Use ' + antList[UseAnt[refantID]] + ' as the refant.'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
msmd.done(); msmd.close()
#-------- Az El
azelTime, AntID, AZ, EL = GetAzEl(msfile)
azelTime_index = np.where( AntID == refantID )[0].tolist()
#-------- Bandpass Table
if 'BPprefix' in locals():
    BPfileName = BPprefix + '-REF' + antList[UseAnt[refantID]] +'-SPW' + `spw` + '-BPant.npy'
    print '---Loading bandpass table : ' + BPfileName
    BP_ant = np.load(BPfileName)
#
chNum, chWid, Freq = GetChNum(msfile, spw); Freq *= 1.0e-9
spwIndex = spwList.index(spw)
#-------- Tau Table
if 'TAUprefix' in locals():
    scanEL = []
    Tau0 = np.load(TAUprefix + '.Tau0.npy')
    interval, BPtime = GetTimerecord(msfile, refantID, refantID, 0, spw, BPscan)
    BPEL = EL[azelTime_index[argmin( abs(azelTime[azelTime_index] - median(BPtime)))]]
    TauBP = Tau0[spwIndex] / np.sin(BPEL)
    BP_ant *= np.exp(0.5*TauBP)
    print ' Bandpass Scan %d EL=%.1f Tau=%.3f' %(BPscan, 180.0*BPEL/pi, np.median(Tau0[spwList.index(spw)]/np.sin(BPEL)))
    for scan in scanList:
        interval, scanTime = GetTimerecord(msfile, refantID, refantID, 0, spw, scan)
        scanEL = scanEL + [EL[azelTime_index[argmin( abs(azelTime[azelTime_index] - median(scanTime)))]]]
    #
#
#-------- Loop for Scan
GainAP, timeList, SNRList, flagList, scanVis, scanWT = [], [], [], [], [], []
for scan in scanList:
    scanIndex = scanList.index(scan)
    print ' Processing Scan %d EL=%.1f Tau=%.3f' %(scan, 180.0*scanEL[scanIndex]/pi, np.median(Tau0[spwList.index(spw)]/np.sin(scanEL[scanIndex]))),
    tauSpec = np.exp(0.5* Tau0[spwList.index(spw)] / np.sin(scanEL[scanIndex]))
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
    timeNum, polNum, chNum = Xspec.shape[3], Xspec.shape[0], Xspec.shape[1]
    if polNum == 4: polIndex = [0, 3]
    if polNum == 2: polIndex = [0, 1]
    if polNum == 1: polIndex = [0]
    polNum = len(polIndex)
    tempSpec = ParaPolBL(Xspec[polIndex][:,:,blMap], blInv).transpose(3,2,0,1)  # Parallel Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    if 'BP_ant' in locals():
        BPCaledXspec = (tempSpec* tauSpec / (BP_ant[ant0]* BP_ant[ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    else:
        BPCaledXspec = tempSpec.transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    scanSNR, scanGain, scanFlag = [], [], []
    for time_index in range(timeNum):
        #------ Progress bar
        progress = (time_index + 1.0) / timeNum
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
        #------ X-pol
        gainFlag = np.ones(UseAntNum)
        for pol_index in range(2):
            tempGain, tempErr = gainComplexErr(chAvgVis[polIndex[pol_index], :, time_index])
            GainAP = GainAP + [tempGain]; tempSNR = abs(tempGain) / tempErr; SNRList = SNRList + [tempSNR]
            scanSNR = scanSNR + [tempSNR]; scanGain = scanGain + [tempGain]
            gainFlag[np.where(tempSNR  < SNR_THRESH)[0]] = 0.0
        #
        scanFlag = scanFlag + [gainFlag]; flagList = flagList + [gainFlag]
    #
    print ''
    #-------- Summary of scan
    timeList.extend(timeStamp.tolist())
    scanFlag = np.array(scanFlag).reshape(timeNum, UseAntNum)
    scanSNR = np.mean(np.array(scanSNR).reshape(timeNum, 2, UseAntNum), axis=0)
    for ant_index in range(UseAntNum):
        print '%s : SNR(med) = %.1f, %.1f  [%d/%d flagged]' % (antList[antMap[ant_index]], scanSNR[0,ant_index], scanSNR[1,ant_index], len(np.where(scanFlag[:,ant_index] == 0.0)[0]), timeNum)
    #
    #-------- Average cross power spectra for scan
    wt = (np.array([scanFlag*scanSNR[0] , scanFlag*scanSNR[1]]).transpose(0,2,1))     # [pol, ant, time]
    wtSum = np.sum(wt[:,ant0]* wt[:,ant1], axis=(1,2))                                # Weight sum
    scanGain = np.array(scanGain).reshape(timeNum, 2, UseAntNum).transpose(1,2,0)     # [pol, ant, time]
    scanPhas = scanGain / abs(scanGain)
    gCalVis = (BPCaledXspec.transpose(1,0,2,3) * wt[:,ant0]* wt[:,ant1]) / (scanPhas[:,ant0]* scanPhas[:,ant1].conjugate()) # [ch, pol, bl, time]
    timeAvgVis = np.sum(gCalVis,axis=3).transpose(1,2,0)  # [pol, bl, ch]
    if DELAYCAL:
        scanVis = scanVis + [np.sum(np.array([delayCalSpec(timeAvgVis[0],chRange), delayCalSpec(timeAvgVis[1],chRange)]), axis=1) ]
    else:
        scanVis = scanVis + [np.sum(timeAvgVis, axis=1)]
    #
    scanWT = scanWT + [wtSum]
#
scanVis = np.array(scanVis)
scanWT  = np.array(scanWT)
totalTimeNum = len(timeList)
antFG  = np.array(flagList).T                          # [ant, time]
antSNR = np.array(SNRList).reshape(totalTimeNum, 2, UseAntNum).transpose(2,1,0)  # [ant, pol, time]
Gain = np.array(GainAP).reshape(totalTimeNum, 2, UseAntNum).transpose(2,1,0)
np.save(prefix + '.Ant.npy', antList[antMap]) 
np.save(prefix + '-SPW' + `spw` + '.TS.npy', np.array(timeList)) 
np.save(prefix + '-SPW' + `spw` + '.GA.npy', Gain) 
np.save(prefix + '-SPW' + `spw` + '.SN.npy', antSNR) 
np.save(prefix + '-SPW' + `spw` + '.FG.npy', antFG) 
np.save(prefix + '-SPW' + `spw` + '.SP.npy', scanVis) 
np.save(prefix + '-SPW' + `spw` + '.WT.npy', scanWT) 
np.save(prefix + '-SPW' + `spw` + '.FQ.npy', Freq) 
