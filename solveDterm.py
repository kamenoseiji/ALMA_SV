execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Grid.py')
execfile(SCR_DIR + 'Plotters.py')
from matplotlib.backends.backend_pdf import PdfPages
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#----------------------------------------- Procedures
fileNum, spwNum = len(prefixList), len(spwList)
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
#
trkAntSet = set(range(64))
scansFile = []
pattern = r'RB_..'
timeNum = 0
sourceList = []
for file_index in range(fileNum):
    prefix = prefixList[file_index]
    msfile = wd + prefix + '.ms'
    sources, posList = GetSourceList(msfile); sourceList = sourceList + sourceRename(sources)
#
sourceList = unique(sourceList).tolist()
sourceScan = []
scanDic   = dict(zip(sourceList, [[]]*len(sourceList)))
StokesDic = dict(zip(sourceList, [[]]*len(sourceList)))
scanIndex = 0
for file_index in range(fileNum):
    prefix = prefixList[file_index]
    msfile = wd + prefix + '.ms'
    print '-- Checking %s ' % (msfile)
    sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList)
    antList = GetAntName(msfile)
    refAntID  = indexList([refantName], antList)
    flagAntID = indexList(antFlag, antList)
    if len(refAntID) < 1:
        print 'Antenna %s didn not participate in this file.' % refantName
        continue
    else:
        refAntID = refAntID[0]
    #
    msmd.open(msfile)
    if 'scanList' in locals():
        scanLS = scanList[file_index]
    else:
        scanLS = msmd.scannumbers().tolist()
    #
    spwName = msmd.namesforspws(spwList[0])[0]; BandName = re.findall(pattern, spwName)[0]; bandID = int(BandName[3:5])
    BandPA = (BANDPA[bandID] + 90.0)*pi/180.0
    for scan in scanLS:
        timeStamp = msmd.timesforscans(scan).tolist()
        timeNum += len(timeStamp)
        trkAnt, scanAnt, Time, Offset = antRefScan( msfile, [timeStamp[0], timeStamp[-1]], antFlag )
        trkAnt = list(set(trkAnt) - set(flagAntID))
        if refAntID in trkAnt:
            trkAntSet = set(trkAnt) & trkAntSet
        else:
            scanLS = list( set(scanLS) - set([scan]) )
        #
        sourceName = sourceList[msmd.sourceidforfield(msmd.fieldsforscan(scan))]
        sourceScan = sourceScan + [sourceName]
        scanDic[sourceName] = scanDic[sourceName] + [scanIndex]
        IQU = GetPolQuery(sourceName, timeStamp[0], BANDFQ[bandID], SCR_DIR)
        StokesDic[sourceName] = [IQU[0][sourceName], IQU[1][sourceName], IQU[2][sourceName], 0.0]
        print '---- Scan%3d : %d tracking antennas : %s' % (scan, len(trkAnt), sourceName)
        scanIndex += 1
    #
    scansFile.append(scanLS)
#
if not 'scanList' in locals(): scanList = scansFile
#-------- Check source list and Stokes Parameters
msmd.done()
antMap = [refAntID] + list(trkAntSet - set([refAntID]))
antNum = len(antMap); blNum = antNum * (antNum - 1)/2
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
blMap, blInv= range(blNum), [False]* blNum
for bl_index in range(blNum):
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#
if not 'bunchNum' in locals(): bunchNum = 1
def bunchVecCH(spec): return bunchVec(spec, bunchNum)
print 'Total %d integration periods.' % (timeNum)
#-------- Loop for SPW
DxList, DyList, FreqList = [], [], []
for spw_index in range(spwNum):
    spw = spwList[spw_index]
    chNum, chWid, Freq = GetChNum(msfile, spw); chRange = range(int(0.05*chNum/bunchNum), int(0.95*chNum/bunchNum)); FreqList = FreqList + [1.0e-9* Freq]
    DxSpec, DySpec = np.zeros([antNum, chNum], dtype=complex), np.zeros([antNum, chNum], dtype=complex)
    caledVis = np.ones([4,blNum, 0], dtype=complex)
    if 'FGprefix' in locals():  # Flag table
        FG = np.load(FGprefix + '-SPW' + `spw` + '.FG.npy'); FG = np.min(FG, axis=0)
        TS = np.load(FGprefix + '-SPW' + `spw` + '.TS.npy')
    #
    if 'BPprefix' in locals():  # Bandpass file
        BPantList, BP_ant = np.load(BPprefix + '-REF' + refantName + '.Ant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-BPant.npy')
        BP_ant = BP_ant[indexList(antList[antMap], BPantList)]      # BP antenna mapping
    if 'XYprefix' in locals():
        XYspec = np.load(XYprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYspec.npy')
        print 'Apply XY phase into Y-pol Bandpass.'; BP_ant[:,1] *= XYspec  # XY phase cal
    #
    BP_ant = np.apply_along_axis(bunchVecCH, 2, BP_ant)
    BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()    # Baseline-based bandpass table
    #
    chAvgVis = np.zeros([4, blNum, timeNum], dtype=complex)
    VisSpec  = np.zeros([4, chNum/bunchNum, blNum, timeNum], dtype=complex)
    mjdSec, scanST, scanET, Az, El, PA, QCpUS, UCmQS = [], [], [], [], [], [], [], []
    timeIndex, scanIndex  = 0, 0
    #-------- Store visibilities into memory
    for file_index in range(fileNum):
        prefix = prefixList[file_index]
        msfile = wd + prefix + '.ms'
        #-------- AZ, EL, PA
        azelTime, AntID, AZ, EL = GetAzEl(msfile)
        azelTime_index = np.where( AntID == refAntID )[0].tolist()
        timeThresh = np.median( np.diff( azelTime[azelTime_index]))
        for scan in scanList[file_index]:
            #-------- Load Visibilities
            print '-- Loading visibility data %s SPW=%d SCAN=%d...' % (prefix, spw, scan)
            timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan, -1, True)  # Xspec[POL, CH, BL, TIME]
            if 'FG' in locals():
                flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
            else:
                flagIndex = range(len(timeStamp))
            if chNum == 1:
                print '  -- Channel-averaged data: no BP and delay cal'
                # chAvgVis= np.c_[chAvgVis, CrossPolBL(Xspec[:,:,blMap], blInv)[:,0]]
                chAvgVis[:, :, timeIndex:timeIndex + len(timeStamp)] =  CrossPolBL(Xspec[:,:,blMap], blInv)[:,0]
            else:
                print '  -- Apply bandpass cal'
                tempSpec = CrossPolBL(np.apply_along_axis(bunchVecCH, 1, Xspec[:,:,blMap]), blInv).transpose(3, 2, 0, 1)[flagIndex]
                VisSpec[:,:,:,timeIndex:timeIndex + len(timeStamp)] = (tempSpec / BP_bl).transpose(2,3,1,0) 
            #
            scanST = scanST + [timeIndex]
            timeIndex += len(timeStamp)
            scanET = scanET + [timeIndex]
            #-------- Expected polarization responses
            scanAz, scanEl = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refAntID, AZ, EL)
            mjdSec = mjdSec + timeStamp[flagIndex].tolist()
            scanPA = AzEl2PA(scanAz, scanEl, ALMA_lat) + BandPA
            CS, SN = np.cos(2.0* scanPA), np.sin(2.0* scanPA)
            QCpUS = QCpUS + ((StokesDic[sourceScan[scanIndex]][1]*CS + StokesDic[sourceScan[scanIndex]][2]*SN)/StokesDic[sourceScan[scanIndex]][0]).tolist()
            UCmQS = UCmQS + ((StokesDic[sourceScan[scanIndex]][2]*CS - StokesDic[sourceScan[scanIndex]][1]*SN)/StokesDic[sourceScan[scanIndex]][0]).tolist()
            Az = Az + scanAz.tolist()
            El = El + scanEl.tolist()
            PA = PA + scanPA.tolist()
            scanIndex += 1
        #
    #
    scanNum = scanIndex
    if chNum > 1: chAvgVis = np.mean(VisSpec[:,chRange], axis=1)
    mjdSec, Az, El, QCpUS, UCmQS = np.array(mjdSec), np.array(Az), np.array(El), np.array(QCpUS), np.array(UCmQS)
    print '---- Antenna-based gain solution using tracking antennas'
    Gain = np.array([ gainComplexVec(chAvgVis[0]), gainComplexVec(chAvgVis[3]) ])   # Parallel-pol gain
    Gamp = np.sqrt(np.mean(abs(Gain)**2, axis=0))
    Gain = Gamp* Gain/abs(Gain)
    #-------- Gain-calibrated visibilities
    print '  -- Apply parallel-hand gain calibration'
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    PAnum = len(PA); PA = np.array(PA)
    """
    print '  -- Solution for Q and U'
    #-------- Coarse estimation of Q and U using XX and YY
    if 'QUsol' not in locals():
        QUsol   = XXYY2QU(PA, Vis[[0,3]])             # XX*, YY* to estimate Q, U
        text_sd = '[XX,YY]  Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (QUsol[0], QUsol[1], 100.0* np.sqrt(QUsol[0]**2 + QUsol[1]**2), np.arctan2(QUsol[1],QUsol[0])*90.0/pi); print text_sd
    """
    ##-------- XY phase determination
    print '  -- XY phase determination'
    XYphase = XY2Phase(UCmQS, Vis[[1,2]])    # XY*, YX* to estimate X-Y phase
    XYsign = np.sign(np.cos(XYphase))
    text_sd = '  XY Phase = %6.2f [deg]  sign = %3.0f' % (XYphase* 180.0 / pi, XYsign); print text_sd
    #-------- Gain adjustment
    print '  -- Polarized gain calibration'
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY* XYsign])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    CS, SN = np.cos(2.0* PA), np.sin(2.0* PA)
    #-------- Fine estimation of Q and U using XY and YX
    for sourceName in sourceList:
        scanList = scanDic[sourceName]
        if len(scanList) < 1 : continue
        timeIndex = []
        for scanIndex in scanList: timeIndex = timeIndex + range(scanST[scanIndex], scanET[scanIndex])
        Qsol, Usol = XY2Stokes(PA[timeIndex], Vis[[1,2]][:,timeIndex])
        text_sd = '[XY,YX] %s:  Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, Qsol, Usol, 100.0* np.sqrt(Qsol**2 + Usol**2), np.arctan2(Usol,Qsol)*90.0/pi); print text_sd
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
    #
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    #-------- XY phase correction
    XYphase, XYproduct = XY2PhaseVec(mjdSec - np.median(mjdSec), UCmQS, Vis[[1,2]])
    twiddle = np.exp((1.0j)* XYphase)
    caledVis[1] /= twiddle
    caledVis[2] *= twiddle
    #-------- Fine estimation of Q and U using XY and YX
    Vis    = np.mean(caledVis, axis=1)
    for sourceName in sourceList:
        scanList = scanDic[sourceName]
        if len(scanList) < 1 : continue
        timeIndex = []
        for scanIndex in scanList: timeIndex = timeIndex + range(scanST[scanIndex], scanET[scanIndex])
        Qsol, Usol = XY2Stokes(PA[timeIndex], Vis[[1,2]][:,timeIndex])
        text_sd = '[XY,YX] %s:  Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, Qsol, Usol, 100.0* np.sqrt(Qsol**2 + Usol**2), np.arctan2(Usol,Qsol)*90.0/pi); print text_sd
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
    #
    #QUsol[0], QUsol[1] = XY2Stokes(PA, Vis[[1,2]])
    #text_sd = '[XY,YX]  Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (QUsol[0], QUsol[1], 100.0* np.sqrt(QUsol[0]**2 + QUsol[1]**2), np.arctan2(QUsol[1],QUsol[0])*90.0/pi); print text_sd
    #GainX, GainY = polariGain(caledVis[0], caledVis[3], PA, QUsol[0], QUsol[1])
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    GainCaledVisSpec = VisSpec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    #-------- Antenna-based on-axis D-term (chAvg)
    #Dx, Dy = VisPA_solveD(caledVis, PA, np.array([1.0, QUsol[0], QUsol[1], 0.0]))
    #CS, SN = np.cos(2.0* PA), np.sin(2.0*PA)
    #QCpUS = QUsol[0]*CS + QUsol[1]*SN
    #UCmQS = QUsol[1]*CS - QUsol[0]*SN
    Dx, Dy = VisMuiti_solveD(caledVis, QCpUS, UCmQS)
    #-------- D-term-corrected Stokes parameters
    Minv = InvMullerVector(Dx[ant1], Dy[ant1], Dx[ant0], Dy[ant0], np.ones(blNum, dtype=complex))
    for sourceName in sourceList:
        scanList = scanDic[sourceName]
        if len(scanList) < 1 : continue
        timeIndex = []
        for scanIndex in scanList: timeIndex = timeIndex + range(scanST[scanIndex], scanET[scanIndex])
        timeNum = len(timeIndex)
        PS = InvPAVector(PA[timeIndex], np.ones(timeNum))
        StokesVis = PS.reshape(4, 4*timeNum).dot(Minv.reshape(4, 4*blNum).dot(caledVis[:,:,timeIndex].reshape(4*blNum, timeNum)).reshape(4*timeNum)) / (timeNum* blNum)
        Qsol, Usol = StokesVis[1].real, StokesVis[2].real; Qerr, Uerr = abs(StokesVis[1].imag), abs(StokesVis[2].imag)
        text_sd = '%s:  Q/I= %6.3f+-%6.4f  U/I= %6.3f+-%6.4f EVPA = %6.2f deg' % (sourceName, Qsol, Qerr, Usol, Uerr, np.arctan2(Usol,Qsol)*90.0/pi); print text_sd
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
    #
    #-------- get D-term spectra
    for ch_index in range(int(chNum/bunchNum)):
        DxSpec[:,ch_index], DySpec[:,ch_index] = VisMuiti_solveD(GainCaledVisSpec[ch_index], QCpUS, UCmQS, Dx, Dy)
        progress = (ch_index + 1.0) / (chNum / bunchNum)
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
    #
    sys.stderr.write('\n'); sys.stderr.flush()
    for ant_index in range(antNum):
        if Dsmooth :
            DX_real, DX_imag = UnivariateSpline(bunchVecCH(Freq), DxSpec[ant_index,range(int(chNum/bunchNum))].real), UnivariateSpline(bunchVecCH(Freq), DxSpec[ant_index,range(int(chNum/bunchNum))].imag)
            DY_real, DY_imag = UnivariateSpline(bunchVecCH(Freq), DySpec[ant_index,range(int(chNum/bunchNum))].real), UnivariateSpline(bunchVecCH(Freq), DySpec[ant_index,range(int(chNum/bunchNum))].imag)
            DxSpec[ant_index] = DX_real(Freq) + (0.0 + 1.0j)* DX_imag(Freq)
            DySpec[ant_index] = DY_real(Freq) + (0.0 + 1.0j)* DY_imag(Freq)
        #
    #
    #-------- D-term-corrected visibilities (invD dot Vis = PS)
    """
    M  = InvMullerVector(DxSpec[ant0], DySpec[ant0], DxSpec[ant1], DySpec[ant1], np.ones([blNum,chNum]))
    DcorrectedVis = np.zeros([4, chNum, blNum, PAnum], dtype=complex)
    for bl_index in range(blNum):
        for ch_index in range(chNum):
            for pa_index in range(PAnum):
                DcorrectedVis[:,ch_index,bl_index,pa_index] = M[:,:,bl_index, ch_index].dot(VisSpec[:,ch_index,bl_index,pa_index])
            #
        #
    #
    chAvgVis = np.mean(DcorrectedVis[:,chRange], axis=1)
    Gain = np.array([gainComplexVec(chAvgVis[0]), gainComplexVec(chAvgVis[3])])   # Gain[pol, ant, time]
    caledVis = DcorrectedVis.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    for sourceName in sourceList:
        scanList = scanDic[sourceName]
        if len(scanList) < 1 : continue
        timeIndex = []
        for scanIndex in scanList: timeIndex = timeIndex + range(scanST[scanIndex], scanET[scanIndex])
        timeNum = len(timeIndex)
        PS = InvPAVector(PA[timeIndex], np.ones(timeNum))
        StokesVis = PS*  np.mean(caledVis[:,:,timeIndex], axis=1)

.reshape(4*blNum, timeNum)).reshape(4*timeNum)) / (timeNum* blNum)

        Qsol, Usol = StokesVis[1].real, StokesVis[2].real; Qerr, Uerr = abs(StokesVis[1].imag), abs(StokesVis[2].imag)
        text_sd = '%s:  Q/I= %6.3f+-%6.4f  U/I= %6.3f+-%6.4f EVPA = %6.2f deg' % (sourceName, Qsol, Qerr, Usol, Uerr, np.arctan2(Usol,Qsol)*90.0/pi); print text_sd
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
    #
    #XX, YY = np.mean(GainCaledVisSpec[:,0], axis=2), np.mean(GainCaledVisSpec[:,3], axis=2)
    #BP_ant[:,0] *= gainComplexVec(XX.T); BP_ant[:,1] *= gainComplexVec(YY.T)
    #BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()    # Baseline-based bandpass table

    #-------- D-term-corrected visibilities (invD dot Vis = PS)
    print 'Update Bandpass with D-term-corrected visibilities'
    PS = (PAVector(PA, np.ones(PAnum)).transpose(2,0,1).dot(StokesVis.real))
    M  = MullerVector(DxSpec[ant0], DySpec[ant0], DxSpec[ant1], DySpec[ant1], np.ones([blNum,chNum])).transpose(0,3,2,1)
    DcorrectedVis = VisSpec / M.dot(PS.T)
    chAvgVis = np.mean(DcorrectedVis[:,chRange], axis=1)
    Gain = np.array([gainComplexVec(chAvgVis[0]), gainComplexVec(chAvgVis[3])])   # Gain[pol, ant, time]
    caledVis = DcorrectedVis.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    #-------- Updating bandpass
    XX, YY = np.mean(caledVis[:,0], axis=2), np.mean(caledVis[:,3], axis=2)
    BP_ant[:,0] *= gainComplexVec(XX.T); BP_ant[:,1] *= gainComplexVec(YY.T)
    #-------- Updating bandpass XY phase
    XY  = np.mean(caledVis[:,1:3], axis=(1,2))  # XY[ch, time]
    XYresponse = XY.dot(abs(PS[:,1]))
    XYtwiddle  = XYresponse / abs(XYresponse)
    BP_ant[:,1] *= XYtwiddle
    XYspec = np.ones(chNum, dtype=complex)
    np.save(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-BPant.npy', BP_ant); print 'Updating BP table : ' + BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-BPant.npy'
    np.save(XYprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYspec.npy', XYspec); print 'Resetting XY phase : ' + XYprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYspec.npy'
    #-------- Time-series XY phase
    XYphase = np.angle((XY.T).dot(XYresponse.conjugate()))
    #-------- Plot
    if np.mean(np.cos(PA)) < 0.0: PA = np.arctan2(-np.sin(PA), -np.cos(PA)) +  np.pi
    PArange = np.arange(min(PA), max(PA), 0.01); CSrange, SNrange = np.cos(2.0*PArange), np.sin(2.0*PArange)
    UCmQS, QCpUS = QUsol[1]*CSrange - QUsol[0]* SNrange, QUsol[0]*CSrange + QUsol[1]* SNrange
    #plt.plot(PArange,  QCpUS + np.mean(Dx).real * UCmQS, '-', color='green')              # XX* - 1.0
    plt.plot(PArange,  QCpUS, '-', color='green')              # XX* - 1.0
    plt.plot(PArange,  -QCpUS, '-', color='orange')            # YY* - 1.0
    plt.plot(PArange,  UCmQS + np.mean(Dx).real * (1.0 - QCpUS) + np.mean(Dy).real* (1.0 + QCpUS), '-', color='cyan')   # ReXY
    plt.plot(PArange,  UCmQS + np.mean(Dx).real* (1.0 - QCpUS) + np.mean(Dy).real* (1.0 + QCpUS), '-', color='magenta')# ReYX
    plt.hlines(0.0, min(PArange), max(PArange), color='grey')
    #plt.plot(PArange,  np.mean(Dx).imag* (1.0 - QCpUS) - np.mean(Dy).imag* (1.0 + QCpUS), '-', color='darkblue')       # ImXY
    #plt.plot(PArange, -np.mean(Dx).imag* (1.0 - QCpUS) + np.mean(Dy).imag* (1.0 + QCpUS), '-', color='darkred')        # ImYX
    plt.plot(PA, Vis[0].real - 1.0, '.', label = 'XX* - 1.0',   color='green')
    plt.plot(PA, Vis[1].real, '.', label = 'ReXY*', color='cyan')
    plt.plot(PA, Vis[1].imag, '.', label = 'ImXY*', color='darkblue')
    plt.plot(PA, Vis[2].real, '.', label = 'ReYX*', color='magenta')
    plt.plot(PA, Vis[2].imag, '.', label = 'ImYX*', color='darkred')
    plt.plot(PA, Vis[3].real - 1.0, '.', label = 'YY* - 1.0',   color='orange')
    plt.xlabel('X-Feed Position Angle [rad]'); plt.ylabel('Normalized cross correlations')
    polMax = np.sqrt(QUsol[0]**2 + QUsol[1]**2); plt.ylim([-1.5* polMax, 1.5* polMax])
    plt.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    text_sd = '(Q, U)/I = (%7.4f+-%6.4f, %7.4f+-%6.4f) XY-phase=%6.2f deg (Ref:%s)' % (QUsol[0], QUerr[0], QUsol[1], QUerr[1], np.angle(np.mean(XYtwiddle[chRange])*np.mean(np.exp((0.0 + 1.0j)* XYphase)))* 180.0/pi,  antList[refAntID]); plt.text(min(PA), polMax*1.2, text_sd, size='x-small')
    plt.savefig(prefixList[0] + '-SPW' + `spw` + '-' + refantName + 'QUXY.pdf', form='pdf')
    #-------- Save Results
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Ant.npy', antList[antMap])
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Azel.npy', np.array([mjdSec, Az, El, PA]))
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.QUXY.npy', QUsol )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.TS.npy', mjdSec )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.GA.npy', Gain )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.XYPH.npy', XYphase )
    """
    for ant_index in range(antNum):
        DtermFile = np.array([FreqList[spw_index], DxSpec[ant_index].real, DxSpec[ant_index].imag, DySpec[ant_index].real, DySpec[ant_index].imag])
        np.save(prefixList[0] + '-SPW' + `spw` + '-' + antList[antMap[ant_index]] + '.DSpec.npy', DtermFile)
    #
    plt.close('all')
    DxList, DyList = DxList + [DxSpec], DyList + [DySpec]
#
#-------- Plot D-term spectra
DxList, DyList = np.array(DxList).transpose(1,0,2), np.array(DyList).transpose(1,0,2)
pp = PdfPages('D_' + prefixList[0] + '-REF' + refantName + '-Dspec.pdf')
plotDSpec(pp, prefixList[0], antList[antMap], spwList, DxList, DyList)
