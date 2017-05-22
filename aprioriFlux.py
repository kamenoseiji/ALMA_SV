import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
from matplotlib.backends.backend_pdf import PdfPages
#-------- Configure Array
print '---Checking array configulation'
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
Tatm_OFS  = 15.0     # Ambient-load temperature - Atmosphere temperature
kb        = 1.38064852e3
#-------- Check Scans for atmCal
logfile = open(prefix + '-' + UniqBands[band_index] + '-Flux.log', 'w') 
logfile.write(BPcalText + '\n')
logfile.write(EQcalText + '\n')
print '---Checking time series in MS and atmCal scans'
tb.open(msfile); timeXY = tb.query('ANTENNA1 == 0 && ANTENNA2 == 0 && DATA_DESC_ID == '+`spw[0]`).getcol('TIME'); tb.close()
OnTimeIndex = []
for scan_index in range(scanNum): OnTimeIndex.append( indexList(msmd.timesforscan(onsourceScans[scan_index]), timeXY) )
offTime, ambTime, hotTime = sort( list(set(timeXY) & set(timeOFF)) ), sort( list(set(timeXY) & set(timeAMB)) ), sort( list(set(timeXY) & set(timeHOT)) )
offTimeIndex, ambTimeIndex, hotTimeIndex = indexList(offTime, timeXY),  indexList(ambTime, timeXY),  indexList(hotTime, timeXY)
#-------- Tsys measurements
try:
    if TSYSCAL :
        execfile(SCR_DIR + 'TsysCal.py')
    else : 
        execfile(SCR_DIR + 'TsysTransfer.py')
except:
    execfile(SCR_DIR + 'TsysCal.py')
#
######## Outputs from TsysCal.py :
#  TantN[ant, spw, pol] : Antenna noise pickup. ant order is the same with MS
#  chAvgTrx[ant, spw, pol, scan]  : Channel-averaged receiver noise temperature
#  chAvgTsky[ant, spw, pol, scan] : Channel-averaged sky noise temperature
#  chAvgTsys[ant, spw, pol, scan] : Channel-averaged system noise temperature
#  TsysFlag[ant, spw, pol, scan] : 0 = invalid, 1=valid
#  Tau0med[spw] : median-value of the zenith optical depth
#  onTau[spw, scan] : on-source optical depth
#
#  They include all of antennas (even if flagged) in MS order
########
#-------- Array Configuration
print '---Checking array configuration'
flagList = np.where(np.median(chAvgTrx.reshape(antNum, 2* spwNum), axis=1) > 2.0* np.median(chAvgTrx))[0].tolist()
flagList = unique(flagList + np.where(np.min(chAvgTrx.reshape(antNum, 2* spwNum), axis=1) < 1.0 )[0].tolist()).tolist()
flagAnt[flagList] = 0.0 # Flagging by abnormal Trx
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
print '  Usable antennas: ',
for ants in antList[UseAnt].tolist(): print ants,
print ''
print '  Flagged by Trx:  ',
for ants in antList[flagList].tolist(): print ants,
print ''
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
timeStamp, UVW = GetUVW(msfile, spw[0], msmd.scansforspw(spw[0])[0])
uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
refantID = bestRefant(uvDist)
print '  Use ' + antList[UseAnt[refantID]] + ' as the refant.'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
antDia = np.ones(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
#-------- Load Aeff file
Afile = open(SCR_DIR + 'AeB' + `int(UniqBands[band_index][3:5])` + '.data')
Alines = Afile.readlines()
Afile.close()
AeX, AeY = 0.25* np.pi* antDia**2, 0.25* np.pi* antDia**2       # antenna collecting area (100% efficiency)
AeX, AeY = [], []
for ant_index in antMap:
    for Aline in Alines:
        if antList[ant_index] in Aline:
            AeX = AeX + [(0.0025* np.pi* float(Aline.split()[1]))* antDia[ant_index]**2]
            AeY = AeY + [(0.0025* np.pi* float(Aline.split()[2]))* antDia[ant_index]**2]
        #
    #
#
AeX, AeY = np.array(AeX), np.array(AeY) # in antMap order 
#-------- Load D-term file
DxList, DyList = [], []
print '---Loading D-term table'
for ant_index in antMap:
    Dpath = SCR_DIR + 'DtermB' + `int(UniqBands[band_index][3:5])` + '/'
    for spw_index in range(4):
        Dfile = Dpath + 'B' + `int(UniqBands[band_index][3:5])` + '-SPW' + `spw_index` + '-' + antList[ant_index] + '.DSpec.npy'
        Dterm = np.load(Dfile)
        DxList = DxList + [Dterm[1] + (0.0 + 1.0j)* Dterm[2]]
        DyList = DyList + [Dterm[3] + (0.0 + 1.0j)* Dterm[4]]
    #
#
chNum = np.array(DxList).shape[1]
DxSpec, DySpec = np.array(DxList).reshape([UseAntNum, spwNum, chNum]), np.array(DyList).reshape([UseAntNum, spwNum, chNum])
#-------- Bandpass Table
print '---Generating antenna-based bandpass table'
BPList = []
for spw_index in spw:
    BP_ant, XY_BP, XYdelay, Gain = BPtable(msfile, spw_index, BPScan, blMap, blInv)
    BP_ant[:,1] *= XY_BP
    BPList = BPList + [BP_ant]
#
if PLOTBP:
    figAnt = PlotBP(msfile, antList[antMap], spw, BPList)
    fileExt = '.pdf'
    if PLOTFMT == 'png': fileExt = '.png'
    for ant_index in range(UseAntNum):
        figAnt = plt.figure(ant_index)
        plotFigFileName = 'BP_' + prefix + '_' + antList[antMap[ant_index]] + '_REF' + antList[UseAnt[refantID]] + '_Scan' + `BPScan` + fileExt
        figAnt.savefig(plotFigFileName)
    #
    plt.close('all')
#
#-------- XY phase using BP scan
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spw[0], BPScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PA = np.arctan2( np.sin(PA), np.cos(PA))
GainP, SEFD, XYphase, caledVis = [], [], [], []
relGain = np.ones([spwNum, UseAntNum, 2])
scan_index = onsourceScans.index(BPScan)
for spw_index in range(spwNum):
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], BPScan)
    timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainP = GainP + [np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])]
    SEFD = SEFD + [2.0* kb* chAvgTsys[antMap, spw_index, :,scan_index].T / np.array([AeX, AeY])]
    aprioriVis = np.mean((chAvgVis / (GainP[spw_index][polYindex][:, ant0]* GainP[spw_index][polXindex][:,ant1].conjugate())).transpose(2, 0, 1)* np.sqrt(SEFD[spw_index][polYindex][:,ant0]* SEFD[spw_index][polXindex][:,ant1]), axis=0)[[0,1]].T
    relGain[spw_index] = abs(gainComplexVec(aprioriVis))
    relGain[spw_index] = relGain[spw_index] / np.median(relGain[spw_index], axis=0)
    SEFD[spw_index] /= (relGain[spw_index]**2).T
    caledVis.append(np.mean((chAvgVis / (GainP[spw_index][polYindex][:,ant0]* GainP[spw_index][polXindex][:,ant1].conjugate())).transpose(2, 0, 1)* np.sqrt(SEFD[spw_index][polYindex][:,ant0]* SEFD[spw_index][polXindex][:,ant1]), axis=2).T)
#
SEFD  = np.array(SEFD)  # SEFD[spw, pol, ant]
GainP = np.array(GainP) # GainP[spw, pol, ant, time]
spwPhase = [0.0]* 2* spwNum
for ant_index in range(1,UseAntNum):
    for pol_index in range(2):
        spwPhase = spwPhase + [0.0]
        for spw_index in range(1,spwNum):
            spwPhase = spwPhase + [np.angle(GainP[spw_index, pol_index, ant_index].dot(GainP[0, pol_index, ant_index].conjugate()))]
        #
    #
#
spwPhase = np.array(spwPhase).reshape([UseAntNum, 2, spwNum]); spwTwiddle = exp(1.0j *spwPhase)
caledVis = np.array(caledVis)   # [spw, pol, time]
#QUsolution = XXYY2QU(PA, np.mean(caledVis[:,[0,3]], axis=0))
QUsolution = np.array([catalogStokesQ.get(BPcal), catalogStokesU.get(BPcal)])
#-------- XY phase cal in Bandpass table
XYsign = np.ones(spwNum)
for spw_index in range(spwNum):
    XYphase = XY2Phase(PA, QUsolution[0], QUsolution[1], caledVis[spw_index][[1,2]])
    XYsign[spw_index] = np.sign(np.cos(XYphase))
    BPList[spw_index] = (BPList[spw_index].transpose(2,0,1)* spwTwiddle[:,:,spw_index]).transpose(1,2,0)
    BPList[spw_index][:,1] *= XYsign[spw_index]
    print 'SPW[%d] : XYphase = %6.1f [deg] sign = %3.0f' % (spw[spw_index], 180.0*XYphase/pi, XYsign[spw_index])
#
#-------- Gain Equalization between X and Y
if 'PolEQ' in locals():
    if PolEQ: execfile(SCR_DIR + 'PolEqualize.py')
#
#-------- Flux Density
ScanFlux = np.zeros([scanNum, spwNum, 4])
ScanSlope= np.zeros([scanNum, spwNum, 4])
ErrFlux  = np.zeros([scanNum, spwNum, 4])
ScanEL     = np.zeros([scanNum])
centerFreqList = []
for spw_index in range(spwNum):
    chNum, chWid, Freq = GetChNum(msfile, spw[spw_index])
    centerFreqList.append( np.median(Freq)*1.0e-9 )
#
print '---Flux densities of sources ---'
pp = PdfPages('FL_' + prefix + '_' + UniqBands[band_index] + '.pdf')
polLabel = ['I', 'Q', 'U', 'V']
Pcolor   = ['black', 'blue', 'red', 'green']
#for scan_index in range(1):
for scan_index in range(scanNum):
    figScan = plt.figure(scan_index, figsize = (11, 8))
    figScan.suptitle(prefix + ' ' + UniqBands[band_index])
    figScan.text(0.75, 0.95, qa.time('%fs' % timeStamp[0], form='ymd')[0]) 
    figScan.text(0.45, 0.05, 'Projected baseline [m]') 
    figScan.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90) 
    ScanEL[scan_index] = np.median(OnEL[:,scan_index])
    if sourceIDscan[scan_index] in SSOList: SSO_flag = T
    else:
        SSO_flag = F 
        SAantennas, SAbl, SAblFlag, SAant0, SAant1 = range(UseAntNum), range(UseBlNum), np.ones([blNum]), ant0, ant1
        SAantMap, SAblMap, SAblInv = antMap, blMap, blInv
    if SSO_flag: continue
    text_sd = ' %02d %010s EL=%4.1f deg' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0* ScanEL[scan_index]/pi ); logfile.write(text_sd + '\n'); print text_sd
    figScan.text(0.05, 0.95, text_sd) 
    text_sd = ' SPW  Frequency    I               Q               U               V               %Pol     EVPA '; logfile.write(text_sd + '\n'); print text_sd
    text_sd = ' ------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
    #
    BPCaledXspec = []
    for spw_index in range(spwNum):
        chNum, chWid, Freq = GetChNum(msfile, spw[spw_index])
        centerFreqList.append( np.median(Freq)*1.0e-9 )
        atmCorrect = np.exp(-onTau[spw_index, scan_index])
        TA = 0.0
        #
        SEFD[spw_index][:,SAantennas] = 2.0* kb* (chAvgTsys[SAantMap, spw_index, :,scan_index] + TA).T / (np.array([AeX[SAantennas], AeY[SAantennas]])* atmCorrect)/ (relGain[spw_index][SAantennas]**2).T
        #-------- UV distance
        timeStamp, UVW = GetUVW(msfile, spw[spw_index], onsourceScans[scan_index])
        uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], onsourceScans[scan_index])
        timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); UseChNum = len(chRange)
        if np.max(abs(Xspec)) < 1.0e-9: continue
        #timeThresh = np.median(diff(timeStamp))
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
        tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        #-------- Bandpass Calibration
        BPCaledXspec = BPCaledXspec + [(tempSpec / (BPList[spw_index][SAant0][:,polYindex]* BPList[spw_index][SAant1][:,polXindex].conjugate())).transpose(2,3,1,0)] # Bandpass Cal
        #
    #
    #-------- Antenna-based Phase Solution
    SAantNum = len(SAantennas); SAblNum = len(SAblMap)
    BPCaledXspec = np.array(BPCaledXspec)   # BPCaledXspec[spw, pol, ch, bl, time]
    chAvgVis = np.mean(np.array(BPCaledXspec)[:,:,chRange], axis=(0,2))
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    pCalVis = (BPCaledXspec.transpose(0,2,1,3,4) / (GainP[polYindex][:,ant0[0:SAblNum]]* GainP[polXindex][:,ant1[0:SAblNum]].conjugate()))[:,chRange]
    for spw_index in range(spwNum):
        chAvgVis = np.mean(pCalVis[spw_index], axis=(0,3))[[0,3]]* np.sqrt(SEFD[spw_index][:,ant0[0:SAblNum]]* SEFD[spw_index][:,ant1[0:SAblNum]]) 
        indivRelGain = abs(gainComplexVec(chAvgVis.T))
        SEFD[spw_index][:,SAantennas] *= ((np.percentile(indivRelGain, 75, axis=0) / indivRelGain)**2).T
    #
    AmpCalVis = (pCalVis.transpose(1,4,0,2,3)* np.sqrt(SEFD[:,polYindex][:,:,ant0[0:SAblNum]]* SEFD[:,polXindex][:,:,ant1[0:SAblNum]])).transpose(2,4,3,0,1) # AmpCalVis[spw,bl,pol,ch,time]
    for spw_index in range(spwNum):
        StokesI_PL = figScan.add_subplot( 2, spwNum, spw_index + 1 )
        StokesP_PL = figScan.add_subplot( 2, spwNum, spwNum + spw_index + 1 )
        text_sd = ' SPW%02d %5.1f GHz' % (spw[spw_index], centerFreqList[spw_index]); logfile.write(text_sd); print text_sd,
        Stokes = np.zeros([4,SAblNum], dtype=complex)
        for bl_index in range(SAblNum):
            Minv = InvMullerVector(DxSpec[SAant1[bl_index], spw_index][chRange], DySpec[SAant1[bl_index], spw_index][chRange], DxSpec[SAant0[bl_index], spw_index][chRange], DySpec[SAant0[bl_index], spw_index][chRange], np.ones(UseChNum))
            Stokes[:,bl_index] = PS.reshape(4, 4*PAnum).dot(Minv.reshape(4, 4*UseChNum).dot(AmpCalVis[spw_index][bl_index].reshape(4*UseChNum, PAnum)).reshape(4*PAnum)) / (PAnum* UseChNum)
        #
        StokesVis, StokesErr = Stokes.real, Stokes.imag
        for pol_index in range(4):
            visFlag = np.where(abs(StokesVis[pol_index] - np.percentile(StokesVis[pol_index], 75))/np.percentile(StokesVis[0], 75) < 0.2 )[0]
            #visFlag = np.where(abs(StokesVis[pol_index] - np.median(StokesVis[pol_index]))/np.median(StokesVis[0]) < 0.2 )[0]
            #if len(visFlag) < 4:
            #    text_sd = 'Only %d vis.    ' % (len(visFlag)) ; logfile.write(text_sd + '\n'); print text_sd,
            #    continue
            #
            weight = np.zeros(SAblNum); weight[visFlag] = 1.0/np.var(StokesVis[pol_index][visFlag])
            P, W = np.c_[np.ones(SAblNum), uvDist], np.diag(weight)
            PtWP_inv = scipy.linalg.inv(P.T.dot(W.dot(P)))
            solution, solerr = PtWP_inv.dot(P.T.dot(weight* StokesVis[pol_index])),  np.sqrt(np.diag(PtWP_inv))
            if solution[1] > -2.0* solerr[1]: solution[0] = np.median(StokesVis[pol_index]); solution[1] = 0.0
            ScanFlux[scan_index, spw_index, pol_index], ScanSlope[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index] = solution[0], solution[1], solerr[0]
            text_sd = '%6.3f (%.3f) ' % (ScanFlux[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index]); logfile.write(text_sd); print text_sd,
        #
        StokesI_PL.plot( uvDist, StokesVis[0], '.', label=polLabel[0], color=Pcolor[0])
        StokesP_PL.plot( uvDist, StokesVis[1], '.', label=polLabel[1], color=Pcolor[1])
        StokesP_PL.plot( uvDist, StokesVis[2], '.', label=polLabel[2], color=Pcolor[2])
        StokesP_PL.plot( uvDist, StokesVis[3], '.', label=polLabel[3], color=Pcolor[3])
        text_sd = '%6.3f   %6.1f \n' % (100.0* np.sqrt(ScanFlux[scan_index, spw_index, 1]**2 + ScanFlux[scan_index, spw_index, 2]**2)/ScanFlux[scan_index, spw_index, 0], np.arctan2(ScanFlux[scan_index, spw_index, 2],ScanFlux[scan_index, spw_index, 1])*90.0/pi); logfile.write(text_sd); print text_sd,
        #
    #
    uvMax, IMax = max(uvDist), max(ScanFlux[scan_index,:,0])
    for spw_index in range(spwNum):
        StokesI_PL = figScan.add_subplot( 2, spwNum, spw_index + 1 )
        StokesP_PL = figScan.add_subplot( 2, spwNum, spwNum + spw_index + 1 )
        StokesI_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 0], ScanFlux[scan_index, spw_index, 0]+ uvMax* ScanSlope[scan_index, spw_index, 0]]), '-', color=Pcolor[0])
        StokesP_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 1], ScanFlux[scan_index, spw_index, 1]+ uvMax* ScanSlope[scan_index, spw_index, 1]]), '-', color=Pcolor[1])
        StokesP_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 2], ScanFlux[scan_index, spw_index, 2]+ uvMax* ScanSlope[scan_index, spw_index, 2]]), '-', color=Pcolor[2])
        StokesP_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 3], ScanFlux[scan_index, spw_index, 3]+ uvMax* ScanSlope[scan_index, spw_index, 3]]), '-', color=Pcolor[3])
        StokesI_PL.axis([0.0, uvMax, 0.0, 1.25*IMax]); StokesP_PL.axis([0.0, uvMax, -0.25*IMax, 0.25*IMax])
        StokesI_PL.text(0.0, 1.26*IMax, 'SPW%2d %5.1f GHz' % (spw[spw_index], centerFreqList[spw_index]))
    #
    StokesI_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    StokesP_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    figScan.savefig(pp, format='pdf')
    freqArray = np.array(centerFreqList)[range(spwNum)]; meanFreq = np.mean(freqArray); relFreq = freqArray - meanFreq
    text_sd = ' ------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
    text_sd = ' mean  %5.1f GHz' % (meanFreq); logfile.write(text_sd); print text_sd,
    pflux, pfluxerr = np.zeros(4), np.zeros(4)
    for pol_index in range(4):
        sol, solerr = linearRegression(relFreq, ScanFlux[scan_index, :, pol_index], ErrFlux[scan_index, :, pol_index] ); pflux[pol_index], pfluxerr[pol_index] = sol[0], solerr[0]
        text_sd = '%6.3f (%.3f) ' % (pflux[pol_index], pfluxerr[pol_index]) ; logfile.write(text_sd); print text_sd,
    #
    text_sd = '%6.3f   %6.1f \n' % (100.0* np.sqrt(pflux[1]**2 + pflux[2]**2)/pflux[0], np.arctan2(pflux[2],pflux[1])*90.0/pi); logfile.write(text_sd); print text_sd,
    print '\n'; logfile.write('\n')
    if COMPDB & (not SSO_flag) : 
        print ' -------- Comparison with ALMA Calibrator Catalog --------'
        au.searchFlux(sourcename='%s' % (sourceList[sourceIDscan[scan_index]]), band=int(UniqBands[band_index][3:5]), date=timeLabelBP[0:10], maxrows=3)
    #
    print '\n'; logfile.write('')
#
plt.close('all')
pp.close()
logfile.close()
np.save(prefix + '-' + UniqBands[band_index] + '.Flux.npy', ScanFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Ferr.npy', ErrFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Source.npy', np.array(sourceList)[sourceIDscan])
np.save(prefix + '-' + UniqBands[band_index] + '.EL.npy', ScanEL)
