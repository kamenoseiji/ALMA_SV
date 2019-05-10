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
scanDic   = dict(zip(sourceList, [[]]*len(sourceList))) # Scan list index for each source
timeDic   = dict(zip(sourceList, [[]]*len(sourceList))) # Time index list for each source
StokesDic = dict(zip(sourceList, [[]]*len(sourceList))) # Stokes parameters for each source
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
        interval, timeStamp = GetTimerecord(msfile, 0, 1, spwList[0], scan)
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
        print '---- Scan%3d : %d tracking antennas : %s, %d records, expected I = %.1f p=%.1f%%' % (scan, len(trkAnt), sourceName, len(timeStamp), StokesDic[sourceName][0], 100.0*sqrt(StokesDic[sourceName][1]**2 + StokesDic[sourceName][2]**2)/StokesDic[sourceName][0])
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
    mjdSec, scanST, scanET, Az, El, PA = [], [], [], [], [], []
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
            mjdSec, scanPA = mjdSec + timeStamp[flagIndex].tolist(), AzEl2PA(scanAz, scanEl, ALMA_lat) + BandPA
            Az, El, PA = Az + scanAz.tolist(), El + scanEl.tolist(), PA + scanPA.tolist()
            scanIndex += 1
        #
    #
    scanNum = scanIndex
    if chNum > 1: chAvgVis = np.mean(VisSpec[:,chRange], axis=1)
    mjdSec, Az, El = np.array(mjdSec), np.array(Az), np.array(El)
    print '---- Antenna-based gain solution using tracking antennas'
    Gain = np.array([ gainComplexVec(chAvgVis[0]), gainComplexVec(chAvgVis[3]) ])   # Parallel-pol gain
    Gamp = np.sqrt(np.mean(abs(Gain)**2, axis=0))
    Gain = Gamp* Gain/abs(Gain)
    #-------- Gain-calibrated visibilities
    print '  -- Apply parallel-hand gain calibration'
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    PAnum = len(PA); PA = np.array(PA)
    QCpUS, UCmQS =  np.zeros(PAnum), np.zeros(PAnum)
    print '  -- Solution for Q and U'
    #-------- Coarse estimation of Q and U using XX and YY
    if 'QUmodel' not in locals(): QUmodel = False
    CS, SN = np.cos(2.0* PA), np.sin(2.0* PA)
    for sourceName in sourceList:
        scanList = scanDic[sourceName]
        if len(scanList) < 1 : continue
        timeIndex = []
        for scanIndex in scanList: timeIndex = timeIndex + range(scanST[scanIndex], scanET[scanIndex])
        timeDic[sourceName] = timeIndex
        if QUmodel:
            QUsol = np.array(StokesDic[sourceName])[[1,2]]/StokesDic[sourceName][0]
        else:
            QUsol   = XXYY2QU(PA[timeIndex], Vis[[0,3]][:,timeIndex])             # XX*, YY* to estimate Q, U
            text_sd = '[XX,YY] %s: Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, QUsol[0], QUsol[1], 100.0* np.sqrt(QUsol[0]**2 + QUsol[1]**2), np.arctan2(QUsol[1],QUsol[0])*90.0/pi); print text_sd
        #
        QCpUS[timeIndex] = QUsol[0]* CS[timeIndex] + QUsol[1]* SN[timeIndex]
        UCmQS[timeIndex] = QUsol[1]* CS[timeIndex] - QUsol[0]* SN[timeIndex]
    #
    ##-------- XY phase determination
    print '  -- Degenerating pi-ambiguity in XY phase'
    XYphase = XY2Phase(UCmQS, Vis[[1,2]])    # XY*, YX* to estimate X-Y phase
    XYsign = np.sign(np.cos(XYphase))
    text_sd = '  XY Phase = %6.2f [deg]  sign = %3.0f' % (XYphase* 180.0 / pi, XYsign); print text_sd
    #-------- Gain adjustment
    print '  -- Polarized gain calibration'
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY* XYsign])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    #-------- Fine estimation of Q and U using XY and YX
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
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
    print '  -- XY phase correction'
    Vis    = np.mean(caledVis, axis=1)
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        Qsol, Usol = XY2Stokes(PA[timeIndex], Vis[[1,2]][:,timeIndex])
        text_sd = '[XY,YX] %s:  Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, Qsol, Usol, 100.0* np.sqrt(Qsol**2 + Usol**2), np.arctan2(Usol,Qsol)*90.0/pi); print text_sd
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
    #
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    GainCaledVisSpec = VisSpec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    #-------- Antenna-based on-axis D-term (chAvg)
    StokesI = np.ones(PAnum)
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        caledVis[:,:,timeIndex] *= StokesDic[sourceName][0]
        QCpUS[timeIndex] *= StokesDic[sourceName][0]
        UCmQS[timeIndex] *= StokesDic[sourceName][0]
        StokesI[timeIndex] *= StokesDic[sourceName][0]
    #
    Dx, Dy = VisMuiti_solveD(caledVis, QCpUS, UCmQS, [], [], StokesI)
    #-------- D-term-corrected Stokes parameters
    Minv = InvMullerVector(Dx[ant1], Dy[ant1], Dx[ant0], Dy[ant0], np.ones(blNum, dtype=complex))
    print '  -- D-term-corrected visibilities'
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        timeNum = len(timeIndex)
        if timeNum < 1 : continue
        PS = InvPAVector(PA[timeIndex], np.ones(timeNum))
        StokesVis = PS.reshape(4, 4*timeNum).dot(Minv.reshape(4, 4*blNum).dot(caledVis[:,:,timeIndex].reshape(4*blNum, timeNum)).reshape(4*timeNum)) / (timeNum* blNum)
        Isol, Qsol, Usol = StokesVis[0].real, StokesVis[1].real, StokesVis[2].real
        Ierr, Qerr, Uerr = abs(StokesVis[0].imag), abs(StokesVis[1].imag), abs(StokesVis[2].imag)
        text_sd = '%s: I= %6.3f  Q= %6.3f+-%6.4f  U= %6.3f+-%6.4f EVPA = %6.2f deg' % (sourceName, Isol, Qsol, Qerr, Usol, Uerr, np.arctan2(Usol,Qsol)*90.0/pi); print text_sd
        StokesDic[sourceName] = (np.array([Isol, Qsol, Usol, 0.0])).tolist()
        StokesI[timeIndex] = Isol
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
        GainCaledVisSpec[:,:,:,timeIndex] *= Isol
    #
    #-------- get D-term spectra
    print '  -- Determining D-term spectra'
    for ch_index in range(int(chNum/bunchNum)):
        DxSpec[:,ch_index], DySpec[:,ch_index] = VisMuiti_solveD(GainCaledVisSpec[ch_index], QCpUS, UCmQS, Dx, Dy, StokesI)
        progress = (ch_index + 1.0) / (chNum / bunchNum)
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
    #
    sys.stderr.write('\n'); sys.stderr.flush()
    for ant_index in range(antNum):
        if 'Dsmooth' in locals():
            if Dsmooth :
                DX_real, DX_imag = UnivariateSpline(bunchVecCH(Freq), DxSpec[ant_index,range(int(chNum/bunchNum))].real), UnivariateSpline(bunchVecCH(Freq), DxSpec[ant_index,range(int(chNum/bunchNum))].imag)
                DY_real, DY_imag = UnivariateSpline(bunchVecCH(Freq), DySpec[ant_index,range(int(chNum/bunchNum))].real), UnivariateSpline(bunchVecCH(Freq), DySpec[ant_index,range(int(chNum/bunchNum))].imag)
                DxSpec[ant_index] = DX_real(Freq) + (0.0 + 1.0j)* DX_imag(Freq)
                DySpec[ant_index] = DY_real(Freq) + (0.0 + 1.0j)* DY_imag(Freq)
            #
        #
    #
    #-------- D-term-corrected visibilities (invD dot Vis = PS)
    print '  -- Applying D-term spectral correction'
    M  = InvMullerVector(DxSpec[ant0], DySpec[ant0], DxSpec[ant1], DySpec[ant1], np.ones([blNum,chNum])).transpose(0, 2, 3, 1)
    DcorrectedVis = np.zeros([4, blNum, chNum, PAnum], dtype=complex)
    for pa_index in range(PAnum):
        progress = (pa_index + 1.0) / PAnum
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
        DcorrectedVis[:,:,:,pa_index] = np.sum(M* GainCaledVisSpec[:,:,:,pa_index].transpose(2, 0, 1), axis=3)
    #
    sys.stderr.write('\n'); sys.stderr.flush()
    chAvgVis = np.mean(DcorrectedVis[:,:,chRange], axis=2)
    GainX, GainY = polariGain(chAvgVis[0], chAvgVis[3], QCpUS)
    Gain = np.array([GainX, GainY])
    chAvgVis = np.mean( (chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())), axis=1)
    maxP = 0.0
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        #chAvgVis[:,timeIndex] *= StokesDic[sourceName][0]
        QUsol = np.array([StokesDic[sourceName][1], StokesDic[sourceName][2]])
        maxP = max(maxP, sqrt(QUsol.dot(QUsol)))
        EVPA = 0.5* np.arctan2(QUsol[1], QUsol[0])
        ThetaPlot = PA[timeIndex] - EVPA; ThetaPlot = np.arctan(np.tan(ThetaPlot))
        ThetaMin, ThetaMax = min(ThetaPlot), max(ThetaPlot)
        PArange = np.arange(ThetaMin + EVPA, ThetaMax + EVPA, 0.01)
        ThetaRange = np.arange(ThetaMin, ThetaMax, 0.01)
        CSrange, SNrange = np.cos(2.0*PArange), np.sin(2.0*PArange)
        UCmQS, QCpUS = QUsol[1]*CSrange - QUsol[0]* SNrange, QUsol[0]*CSrange + QUsol[1]* SNrange
        ThetaRange[ThetaRange >  1.56] = np.inf
        ThetaRange[ThetaRange < -1.56] = -np.inf
        plt.plot(RADDEG* ThetaRange,  QCpUS, '-', label=sourceName + ' XX* - I')     # XX* - 1.0
        plt.plot(RADDEG* ThetaRange, -QCpUS, '-', label=sourceName + ' YY* - I')     # YY* - 1.0
        plt.plot(RADDEG* ThetaRange,  np.mean(Dy).real* (StokesDic[sourceName][0] + QCpUS) + UCmQS + np.mean(Dx).real* (StokesDic[sourceName][1] - QCpUS), '-', label=sourceName + ' ReXY*')
        plt.plot(RADDEG* ThetaRange, -np.mean(Dy).imag* (StokesDic[sourceName][1] + QCpUS) + np.mean(Dx).imag* (StokesDic[sourceName][1] - QCpUS), '-', label=sourceName + ' ImXY*')
        plt.plot(RADDEG* ThetaRange,  np.mean(Dy).imag* (StokesDic[sourceName][1] + QCpUS) - np.mean(Dx).imag* (StokesDic[sourceName][1] - QCpUS), '-', label=sourceName + ' ImYX*')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[0][timeIndex].real - StokesDic[sourceName][0], ',')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[1][timeIndex].real, ',')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[1][timeIndex].imag, ',')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[2][timeIndex].real, ',')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[2][timeIndex].imag, ',')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[3][timeIndex].real - StokesDic[sourceName][0], ',')
        #text_sd = '%s: (Q, U)/I = (%7.4f+-%6.4f, %7.4f+-%6.4f) XY-phase=%6.2f deg (Ref:%s)' % (sourceName, QUsol[0], QUerr[0], QUsol[1], QUerr[1], np.angle(np.mean(XYtwiddle[chRange])*np.mean(np.exp((0.0 + 1.0j)* XYphase)))* 180.0/pi,  antList[refAntID])
    #
    plt.xlabel('Linear polarization angle w.r.t. X-Feed [deg]'); plt.ylabel('Normalized cross correlations')
    plt.xlim([-90.0, 90.0])
    plt.ylim([-1.5* maxP, 1.5*maxP])
    plt.legend(loc = 'best', prop={'size' :6}, numpoints = 1)
    plt.savefig(prefixList[0] + '-SPW' + `spw` + '-' + refantName + 'QUXY.pdf', form='pdf')
    #-------- Save Results
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Ant.npy', antList[antMap])
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Azel.npy', np.array([mjdSec, Az, El, PA]))
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.TS.npy', mjdSec )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.GA.npy', Gain )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.XYPH.npy', XYphase )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.XYC.npy', chAvgVis[[1,2]])
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
