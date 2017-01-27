execfile(SCR_DIR + 'interferometry.py')
from matplotlib.backends.backend_pdf import PdfPages
#from scipy.constants import constants
#from scipy.interpolate import UnivariateSpline
#from scipy.interpolate import griddata
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#----------------------------------------- Procedures
fileNum, spwNum = len(prefixList), len(spwList)
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
#
trkAntSet = set(range(64))
scansFile = []
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
    scanList = msmd.scannumbers().tolist()
    for scan in scanList:
        timeStamp = msmd.timesforscans(scan).tolist()
        trkAnt, scanAnt, Time, Offset = antRefScan( msfile, [timeStamp[0], timeStamp[-1]] )
        trkAnt = list(set(trkAnt) - set(flagAntID))
        if refAntID in trkAnt:
            trkAntSet = set(trkAnt) & trkAntSet
        else:
            scanList = list( set(scanList) - set([scan]) )
        #
        print '---- Scan %d : %d tracking antennas' % (scan, len(trkAnt))
    #
    scansFile.append(scanList)
#
antMap = [refAntID] + list(trkAntSet - set([refAntID]))
antNum = len(antMap); blNum = antNum * (antNum - 1)/2
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
blMap, blInv= range(blNum), [False]* blNum
for bl_index in range(blNum):
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#
#-------- Loop for SPW
DxList, DyList, FreqList = [], [], []
for spw_index in range(spwNum):
    spw = spwList[spw_index]
    chNum, chWid, Freq = GetChNum(msfile, spw); chRange = range(int(0.05*chNum), int(0.95*chNum)); FreqList = FreqList + [1.0e-9* Freq]
    DxSpec, DySpec = np.zeros([antNum, chNum], dtype=complex), np.zeros([antNum, chNum], dtype=complex)
    caledVis = np.ones([4,blNum, 0], dtype=complex)
    mjdSec, Az, El = np.ones([0]), np.ones([0]), np.ones([0])
    if BPprefix != '':  # Bandpass file
        BPantList, BP_ant, XYspec, XYdelay = np.load(BPprefix + '-REF' + refantName + '.Ant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-BPant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYspec.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYdelay.npy')
        BP_ant = BP_ant[indexList(antList[antMap], BPantList)]
        BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()
    #
    #
    chAvgVis = np.zeros([4, blNum, 0], dtype=complex)
    VisSpec  = np.zeros([4, chNum, blNum, 0], dtype=complex)
    for file_index in range(fileNum):
        prefix = prefixList[file_index]
        msfile = wd + prefix + '.ms'
        #-------- AZ, EL, PA
        azelTime, AntID, AZ, EL = GetAzEl(msfile)
        azelTime_index = np.where( AntID == refAntID )[0].tolist()
        timeThresh = np.median( np.diff( azelTime[azelTime_index]))
        for scan in scansFile[file_index]:
            #-------- Load Visibilities
            print '-- Loading visibility data %s SPW=%d SCAN=%d...' % (prefix, spw, scan)
            timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)  # Xspec[POL, CH, BL, TIME]
            #chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
            if chNum == 1:
                print '  -- Channel-averaged data: no BP and delay cal'
                chAvgVis= np.c_[chAvgVis, CrossPolBL(Xspec[:,:,blMap], blInv)[:,0]]
            else:
                print '  -- Apply bandpass cal'
                tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3, 2, 0, 1)
                BPCaledXspec = (tempSpec / BP_bl).transpose(2,3,1,0) 
                #-------- XY delay cal
                print '  -- XY delay cal'
                BPCaledXspec[1] = (BPCaledXspec[1].transpose(1,2,0)* XYspec.conjugate()).transpose(2,0,1)
                BPCaledXspec[2] = (BPCaledXspec[2].transpose(1,2,0)* XYspec).transpose(2,0,1)
                #-------- Antenna-based Gain
                print '  -- Channel-averaging'
                chAvgVis = np.c_[chAvgVis, np.mean(BPCaledXspec[:,chRange], axis=1)]
                VisSpec  = np.c_[VisSpec, BPCaledXspec]
            #
            #-------- Time index at on axis
            timeNum = len(timeStamp)
            mjdSec = np.append(mjdSec, timeStamp)
            for time_index in range(timeNum):
                scanAz, scanEl = AzElMatch(timeStamp[time_index], azelTime, AntID, refAntID, timeThresh, AZ, EL)
                Az, El = np.append(Az, scanAz), np.append(El, scanEl)
            #
        #
    #
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
    VisSpec  = VisSpec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    PA = AzEl2PA(Az, El, ALMA_lat) + BANDPA     # Apply feed orientation
    PA = np.arctan2( np.sin(PA), np.cos(PA))    # to set in [-pi, pi]
    Vis    = np.mean(caledVis, axis=1)
    PAnum = len(PA)
    #-------- Solve for Stokes Parameters and XY phase
    print '  -- Solution for Q and U'
    solution = np.r_[XXYY2QU(PA, Vis[[0,3]]), np.zeros(5)]              # XX*, YY* to estimate Q, U
    solution[2] = XY2Phase(PA, solution[0], solution[1], Vis[[1,2]])    # XY*, YX* to estimate X-Y phase
    solution, solerr = XY2Stokes(PA, Vis[[1,2]], solution)
    GainX, GainY = polariGain(caledVis[0], caledVis[3], PA, solution[0], solution[1]); Gain = np.array([GainX, GainY])
    caledVis /= (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    VisSpec  /= (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    solution, solerr = XY2Stokes(PA, Vis[[1,2]], solution)
    text_sd = '  Q/I= %6.3f+-%6.4f  U/I= %6.3f+-%6.4f  X-Y phase= %6.3f+-%6.4f rad EVPA = %6.2f deg' % (solution[0], solerr[0], solution[1], solerr[1], solution[2], solerr[2], np.arctan2(solution[1],solution[0])*90.0/pi); print text_sd
    #
    #-------- Antenna-based on-axis D-term (chAvg)
    Dx, Dy = VisPA_solveD(caledVis, PA, np.array([1.0, solution[0], solution[1], 0.0]))
    #-------- D-term-corrected Stokes parameters
    Minv = InvMullerVector(Dx[ant1], Dy[ant1], Dx[ant0], Dy[ant0], np.ones(blNum))
    PS = InvPAVector(PA, np.ones(PAnum))
    StokesVis = PS.reshape(4, 4*PAnum).dot(Minv.reshape(4, 4*blNum).dot(caledVis.reshape(4*blNum, PAnum)).reshape(4*PAnum)) / (PAnum* blNum)
    solution[0], solution[1], solerr[0], solerr[1] = StokesVis[1].real, StokesVis[2].real, abs(StokesVis[1].imag), abs(StokesVis[2].imag)
    text_sd = '  Q/I= %6.3f+-%6.4f  U/I= %6.3f+-%6.4f EVPA = %6.2f deg' % (solution[0], solerr[0], solution[1], solerr[1], np.arctan2(solution[1],solution[0])*90.0/pi); print text_sd
    #-------- apply Gain cal for every channel and get D-term spectra
    for ch_index in range(chNum):
        #GainX, GainY = polariGain(VisSpec[0][ch_index], VisSpec[3][ch_index], PA, solution[0], solution[1])
        #Gain = np.array([GainX,  GainY])
        #VisSpec[:,ch_index] /= (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
        DxSpec[:,ch_index], DySpec[:,ch_index] = VisPA_solveD(VisSpec[ch_index], PA, np.array([1.0, solution[0], solution[1], 0.0]))
        progress = (ch_index + 1.0) / chNum
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
    #
    sys.stderr.write('\n'); sys.stderr.flush()
    #-------- Plot
    if np.mean(np.cos(PA)) < 0.0: PA = np.arctan2(-np.sin(PA), -np.cos(PA)) +  np.pi
    PArange = np.arange(min(PA), max(PA), 0.01);  CSrange, SNrange = np.cos(2.0*PArange), np.sin(2.0*PArange)
    plt.plot(PArange,  CSrange* solution[0] + SNrange* solution[1], '-', color='green')
    plt.plot(PArange,  np.cos(solution[2])* (-SNrange* solution[0] + CSrange* solution[1]) + solution[3], '-', color='cyan')
    plt.plot(PArange,  np.sin(solution[2])* (-SNrange* solution[0] + CSrange* solution[1]) + solution[4], '-', color='darkblue')
    plt.plot(PArange,  np.cos(solution[2])* (-SNrange* solution[0] + CSrange* solution[1]) + solution[5], '-', color='magenta')
    plt.plot(PArange, -np.sin(solution[2])* (-SNrange* solution[0] + CSrange* solution[1]) + solution[6], '-', color='darkred')
    plt.plot(PArange, -CSrange* solution[0] - SNrange* solution[1], '-', color='orange')
    plt.plot(PA, Vis[0].real - 1.0, '.', label = 'XX* - 1.0',   color='green')
    plt.plot(PA, Vis[1].real, '.', label = 'ReXY*', color='cyan')
    plt.plot(PA, Vis[1].imag, '.', label = 'ImXY*', color='darkblue')
    plt.plot(PA, Vis[2].real, '.', label = 'ReYX*', color='magenta')
    plt.plot(PA, Vis[2].imag, '.', label = 'ImYX*', color='darkred')
    plt.plot(PA, Vis[3].real - 1.0, '.', label = 'YY* - 1.0',   color='orange')
    plt.xlabel('X-Feed Position Angle [rad]'); plt.ylabel('Normalized cross correlations')
    plt.ylim([-0.15,0.15])
    plt.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    text_sd = '(Q, U)/I = (%6.3f+-%6.3f, %6.3f+-%6.3f) X-Y phase=%6.3f+-%6.3f rad (Ref:%s)' % (solution[0], solerr[0], solution[1], solerr[1], solution[2], solerr[2], antList[refAntID]); plt.text(min(PA), 0.16, text_sd, size='x-small')
    plt.savefig(prefixList[0] + '-SPW' + `spw` + '-' + refantName + 'QUXY.pdf', form='pdf')
    #-------- Save Results
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Ant.npy', antList[antMap])
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Azel.npy', np.array([mjdSec, Az, El, PA]))
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.QUXY.npy', solution )
    np.save(prefixList[0] + '-SPW' + `spw` + '-DxSpec.npy', DxSpec)
    np.save(prefixList[0] + '-SPW' + `spw` + '-DySpec.npy', DySpec)
    plt.close('all')
    DxList, DyList = DxList + [DxSpec], DyList + [DySpec]
#
#-------- Plot D-term spectra
DxList, DyList = np.array(DxList).transpose(1,0,2), np.array(DyList).transpose(1,0,2)
pp = PdfPages('D_' + prefixList[0] + '-Dspec.pdf')
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index, figsize = (11, 8))
    figAnt.suptitle(prefixList[0] + ' ' + antList[antMap[ant_index]])
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'D-term Spectra (Amplitude and Phase)', rotation=90)
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
        DxPL.axis([np.min(FreqList[spw_index]), np.max(FreqList[spw_index]), -0.1, 0.1])
        #
        DyPL.plot( FreqList[spw_index], plotDy.real, ls='steps-mid', label = 'reDy')
        DyPL.plot( FreqList[spw_index], plotDy.imag, ls='steps-mid', label = 'imDy')
        DyPL.axis([np.min(FreqList[spw_index]), np.max(FreqList[spw_index]), -0.1, 0.1])
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
