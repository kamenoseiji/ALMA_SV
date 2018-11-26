execfile(SCR_DIR + 'interferometry.py')
from matplotlib.backends.backend_pdf import PdfPages
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#----------------------------------------- Procedures
fileNum, spwNum = len(prefixList), len(spwList)
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
#
trkAntSet = set(range(64))
scansFile = []
pattern = r'RB_..'
for file_index in range(fileNum):
    prefix = prefixList[file_index]
    msfile = wd + prefix + '.ms'
    print '-- Checking %s ' % (msfile)
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
    spwName = msmd.namesforspws(spwList[0])[0]; BandName = re.findall(pattern, spwName)[0]; BandPA = (BANDPA[int(BandName[3:5])] + 90.0)*pi/180.0
    for scan in scanLS:
        timeStamp = msmd.timesforscans(scan).tolist()
        trkAnt, scanAnt, Time, Offset = antRefScan( msfile, [timeStamp[0], timeStamp[-1]], antFlag )
        trkAnt = list(set(trkAnt) - set(flagAntID))
        if refAntID in trkAnt:
            trkAntSet = set(trkAnt) & trkAntSet
        else:
            scanLS = list( set(scanLS) - set([scan]) )
        #
        print '---- Scan %d : %d tracking antennas' % (scan, len(trkAnt))
    #
    scansFile.append(scanLS)
#
if not 'scanList' in locals(): scanList = scansFile
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
        BPantList, BP_ant, XYspec = np.load(BPprefix + '-REF' + refantName + '.Ant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-BPant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYspec.npy')
        BP_ant = BP_ant[indexList(antList[antMap], BPantList)]      # BP antenna mapping
        BP_ant[:,1] *= XYspec                                       # XY phase cal
        BP_ant = np.apply_along_axis(bunchVecCH, 2, BP_ant)
        BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()    # Baseline-based bandpass table
    #
    #
    chAvgVis = np.zeros([4, blNum, 0], dtype=complex)
    VisSpec  = np.zeros([4, chNum/bunchNum, blNum, 0], dtype=complex)
    mjdSec, Az, El = [], [], []
    for file_index in range(fileNum):
        prefix = prefixList[file_index]
        msfile = wd + prefix + '.ms'
        #-------- AZ, EL, PA
        azelTime, AntID, AZ, EL = GetAzEl(msfile)
        azelTime_index = np.where( AntID == refAntID )[0].tolist()
        timeThresh = np.median( np.diff( azelTime[azelTime_index]))
        #for scan in scansFile[file_index]:
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
                chAvgVis= np.c_[chAvgVis, CrossPolBL(Xspec[:,:,blMap], blInv)[:,0]]
            else:
                print '  -- Apply bandpass cal'
                tempSpec = CrossPolBL(np.apply_along_axis(bunchVecCH, 1, Xspec[:,:,blMap]), blInv).transpose(3, 2, 0, 1)[flagIndex]
                BPCaledXspec = (tempSpec / BP_bl).transpose(2,3,1,0) 
                #-------- Antenna-based Gain
                print '  -- Channel-averaging'
                chAvgVis = np.c_[chAvgVis, np.mean(BPCaledXspec[:,chRange], axis=1)]
                VisSpec  = np.c_[VisSpec, BPCaledXspec]
            #
            #-------- Time index at on axis
            scanAz, scanEl = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refAntID, AZ, EL)
            mjdSec = mjdSec + timeStamp[flagIndex].tolist()
            Az = Az + scanAz.tolist()
            El = El + scanEl.tolist()
            #for time_index in range(timeNum):
            #    scanAz, scanEl = AzElMatch(timeStamp[time_index], azelTime, AntID, refAntID, AZ, EL)
            #    Az, El = np.append(Az, scanAz), np.append(El, scanEl)
            #
        #
    #
    mjdSec = np.array(mjdSec)
    Az = np.array(Az)
    El = np.array(El)
    print '---- Antenna-based gain solution using tracking antennas'
    Gain  = np.ones([2, antNum, len(mjdSec)], dtype=complex)
    Gain[0, 0:antNum] = gainComplexVec(chAvgVis[0])
    Gain[1, 0:antNum] = gainComplexVec(chAvgVis[3])
    Gamp = np.sqrt(np.mean(abs(Gain)**2, axis=0))
    Gain[0, 0:antNum] = Gain[0, 0:antNum] * Gamp / abs(Gain[0, 0:antNum])
    Gain[1, 0:antNum] = Gain[1, 0:antNum] * Gamp / abs(Gain[1, 0:antNum])
    #-------- Gain-calibrated visibilities
    print '  -- Apply gain calibration'
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    PA = AzEl2PA(Az, El, ALMA_lat) + BandPA     # Apply feed orientation
    PA = np.arctan2( np.sin(PA), np.cos(PA))    # to set in [-pi, pi]
    Vis    = np.mean(caledVis, axis=1)
    PAnum = len(PA)
    print '  -- Solution for Q and U'
    #-------- Coarse estimation of Q and U using XX and YY
    if 'QUsol' not in locals():
        QUsol   = XXYY2QU(PA, Vis[[0,3]])             # XX*, YY* to estimate Q, U
        text_sd = '[XX,YY]  Q/I= %6.3f  U/I= %6.3f  EVPA = %6.2f deg' % (QUsol[0], QUsol[1], np.arctan2(QUsol[1],QUsol[0])*90.0/pi); print text_sd
    #-------- XY phase determination
    XYphase = XY2Phase(PA, QUsol[0], QUsol[1], Vis[[1,2]])    # XY*, YX* to estimate X-Y phase
    XYsign = np.sign(np.cos(XYphase))
    text_sd = '  XY Phase = %6.2f [deg]  sign = %3.0f' % (XYphase* 180.0 / pi, XYsign); print text_sd
    #-------- Gain adjustment
    GainX, GainY = polariGain(caledVis[0], caledVis[3], PA, QUsol[0], QUsol[1])
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY* XYsign])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    #-------- Fine estimation of Q and U using XY and YX
    QUsol[0], QUsol[1] = XY2Stokes(PA, Vis[[1,2]])
    text_sd = '[XY,YX]  Q/I= %6.3f  U/I= %6.3f  EVPA = %6.2f deg' % (QUsol[0], QUsol[1], np.arctan2(QUsol[1],QUsol[0])*90.0/pi); print text_sd
    GainX, GainY = polariGain(caledVis[0], caledVis[3], PA, QUsol[0], QUsol[1])
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    #-------- XY phase correction
    XYphaseVec = XY2PhaseVec(mjdSec - np.median(mjdSec), PA, QUsol[0], QUsol[1], Vis[[1,2]])
    twiddle = np.exp((-1.0j)* XYphaseVec)
    caledVis[1] *= twiddle
    caledVis[2] /= twiddle
    #-------- Fine estimation of Q and U using XY and YX
    Vis    = np.mean(caledVis, axis=1)
    QUsol[0], QUsol[1] = XY2Stokes(PA, Vis[[1,2]])
    text_sd = '[XY,YX]  Q/I= %6.3f  U/I= %6.3f  EVPA = %6.2f deg' % (QUsol[0], QUsol[1], np.arctan2(QUsol[1],QUsol[0])*90.0/pi); print text_sd
    GainX, GainY = polariGain(caledVis[0], caledVis[3], PA, QUsol[0], QUsol[1])
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    caledVis[1] *= twiddle
    caledVis[2] /= twiddle
    Vis    = np.mean(caledVis, axis=1)
    VisSpec[:,1] *= twiddle
    VisSpec[:,2] /= twiddle
    VisSpec = VisSpec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    #
    #-------- Antenna-based on-axis D-term (chAvg)
    Dx, Dy = VisPA_solveD(caledVis, PA, np.array([1.0, QUsol[0], QUsol[1], 0.0]))
    #-------- D-term-corrected Stokes parameters
    Minv = InvMullerVector(Dx[ant1], Dy[ant1], Dx[ant0], Dy[ant0], np.ones(blNum, dtype=complex))
    PS = InvPAVector(PA, np.ones(PAnum))
    StokesVis = PS.reshape(4, 4*PAnum).dot(Minv.reshape(4, 4*blNum).dot(caledVis.reshape(4*blNum, PAnum)).reshape(4*PAnum)) / (PAnum* blNum)
    QUsol[0], QUsol[1] = StokesVis[1].real, StokesVis[2].real; QUerr = np.array([abs(StokesVis[1].imag), abs(StokesVis[2].imag)])
    text_sd = '  Q/I= %6.3f+-%6.4f  U/I= %6.3f+-%6.4f EVPA = %6.2f deg' % (QUsol[0], QUerr[0], QUsol[1], QUerr[1], np.arctan2(QUsol[1],QUsol[0])*90.0/pi); print text_sd
    #-------- get D-term spectra
    for ch_index in range(int(chNum/bunchNum)):
        DxSpec[:,ch_index], DySpec[:,ch_index] = VisPA_solveD(VisSpec[ch_index], PA, np.array([1.0, QUsol[0], QUsol[1], 0.0]), Dx, Dy)
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
    text_sd = '(Q, U)/I = (%6.3f+-%6.3f, %6.3f+-%6.3f) XY-phase=%6.2f deg (Ref:%s)' % (QUsol[0], QUerr[0], QUsol[1], QUerr[1], XYphase*180.0/pi, antList[refAntID]); plt.text(min(PA), polMax*1.2, text_sd, size='x-small')
    plt.savefig(prefixList[0] + '-SPW' + `spw` + '-' + refantName + 'QUXY.pdf', form='pdf')
    #-------- Save Results
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Ant.npy', antList[antMap])
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Azel.npy', np.array([mjdSec, Az, El, PA]))
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.QUXY.npy', QUsol )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.TS.npy', mjdSec )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.XYPH.npy', XYphaseVec )
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
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index, figsize = (11, 8))
    figAnt.suptitle(prefixList[0] + ' ' + antList[antMap[ant_index]])
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'D-term Spectra (Real and Imaginary)', rotation=90)
#
#-------- Plot D-spec
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index)
    for spw_index in range(spwNum):
        DxPL = figAnt.add_subplot( 2, spwNum, spw_index + 1 )
        DyPL = figAnt.add_subplot( 2, spwNum, spw_index + spwNum + 1 )
        #
        plotDx, plotDy = DxList[ant_index, spw_index], DyList[ant_index, spw_index]
        DxPL.plot( FreqList[spw_index], plotDx.real, ls='steps-mid', label = 'reDx')
        DxPL.plot( FreqList[spw_index], plotDx.imag, ls='steps-mid', label = 'imDx')
        DxPL.axis([np.min(FreqList[spw_index]), np.max(FreqList[spw_index]), -0.12, 0.12])
        #
        DyPL.plot( FreqList[spw_index], plotDy.real, ls='steps-mid', label = 'reDy')
        DyPL.plot( FreqList[spw_index], plotDy.imag, ls='steps-mid', label = 'imDy')
        DyPL.axis([np.min(FreqList[spw_index]), np.max(FreqList[spw_index]), -0.12, 0.12])
        #
        DxPL.text( np.min(FreqList[spw_index]), 0.09, 'SPW=' + `spwList[spw_index]`)
    #
    DxPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
    DyPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
    figAnt.savefig(pp, format='pdf')
#
#
plt.close('all')
pp.close()
