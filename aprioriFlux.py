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
tb.open(msfile); timeXY = tb.query('ANTENNA1 == 0 && ANTENNA2 == 0 && DATA_DESC_ID == '+`atmspw[0]`).getcol('TIME'); tb.close()
atmTimeIndex = []
for scan_index in range(atmscanNum): atmTimeIndex.append( indexList(msmd.timesforscan(atmScans[scan_index]), timeXY) )
#for scan_index in range(scanNum): OnTimeIndex.append( indexList(msmd.timesforscan(onsourceScans[scan_index]), timeXY) )
#-------- Tsys measurements
#spwList  = atmspwLists[band_index]
#scanList = atmscanLists[band_index]
#try:
#    if TSYSCAL :
#        execfile(SCR_DIR + 'TsysCal.py')
#    else : 
#        execfile(SCR_DIR + 'TsysTransfer.py')
#except:
#    execfile(SCR_DIR + 'TsysCal.py')
#
#-------- Load Tsys table
Tau0spec = np.load(prefix +  '-' + UniqBands[band_index] + '.Tau0.npy') # Tau0spec[spw][ch]
Trxspec  = np.load(prefix +  '-' + UniqBands[band_index] + '.Trx.npy')  # Trxspec[ant*spw][pol, ch]
OnEL = np.load(prefix +  '-' + UniqBands[band_index] + '.OnEL.npy')
AtmEL = np.load(prefix +  '-' + UniqBands[band_index] + '.AtmEL.npy')
OnEL = np.median(OnEL, axis=0)
msmd.open(msfile)
#-------- Array Configuration
print '---Checking array configuration'
flagList = np.where(np.median(chAvgTrx.reshape(antNum, 2* spwNum), axis=1) > 2.0* np.median(chAvgTrx))[0].tolist()
flagList = unique(flagList + np.where(np.min(chAvgTrx.reshape(antNum, 2* spwNum), axis=1) < 1.0 )[0].tolist()).tolist()
flagAnt[flagList] = 0.0 # Flagging by abnormal Trx
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
text_sd = '  Usable antennas: '
for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
logfile.write(text_sd + '\n'); print text_sd
text_sd = '  Flagged by Trx:  '
for ants in antList[flagList].tolist(): text_sd = text_sd + ants + ' '
logfile.write(text_sd + '\n'); print text_sd
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
timeStamp, UVW = GetUVW(msfile, scnspw[0], msmd.scansforspw(scnspw[0])[0])
uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
refantID = bestRefant(uvDist)
print '  Use ' + antList[UseAnt[refantID]] + ' as the refant.'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
antDia = np.ones(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
msmd.done(); msmd.close()
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
#-------- Bandpass Table
print '---Generating antenna-based bandpass table'
BPList = []
for spw_index in scnspw:
    BP_ant, XY_BP, XYdelay, Gain = BPtable(msfile, spw_index, BPScan, blMap, blInv)
    BP_ant[:,1] *= XY_BP
    BPList = BPList + [BP_ant]
#
if PLOTBP:
    figAnt = PlotBP(msfile, antList[antMap], scnspw, BPList)
    fileExt = '.pdf'
    if PLOTFMT == 'png': fileExt = '.png'
    for ant_index in range(UseAntNum):
        figAnt = plt.figure(ant_index)
        plotFigFileName = 'BP_' + prefix + '_' + antList[antMap[ant_index]] + '_REF' + antList[UseAnt[refantID]] + '_Scan' + `BPScan` + fileExt
        figAnt.savefig(plotFigFileName)
    #
    plt.close('all')
#
##-------- Equalization using EQ scan
GainP, relGain = [], np.ones([scnspwNum, 2, UseAntNum])
polXindex, polYindex, scan_index = (arange(4)//2).tolist(), (arange(4)%2).tolist(), onsourceScans.index(EQScan)
for spw_index in range(scnspwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, scnspw[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = ParaPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)       # Parallel Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0]* BPList[spw_index][ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)[pPol]
    GainP = GainP + [np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])]
    aprioriSEFD = 2.0* kb* np.array(chAvgTsys).reshape(NumBands, scanNum, antNum, spwNum, polNum)[band_index,scan_index][antMap][:,spw_index].T / np.array([AeX, AeY])
    aprioriVisX = np.mean(chAvgVis[0] / (GainP[spw_index][0,ant0]* GainP[spw_index][0,ant1].conjugate()), axis=1) * np.sqrt(aprioriSEFD[0, ant0]* aprioriSEFD[0, ant1])
    aprioriVisY = np.mean(chAvgVis[1] / (GainP[spw_index][1,ant0]* GainP[spw_index][1,ant1].conjugate()), axis=1) * np.sqrt(aprioriSEFD[1, ant0]* aprioriSEFD[1, ant1])
    #-------- Determine Antenna-based Gain
    relGain[spw_index, 0] = abs(gainComplex(aprioriVisX)); relGain[spw_index, 0] /= np.median( abs(relGain[spw_index, 0]) ) # X-pol delta gain
    relGain[spw_index, 1] = abs(gainComplex(aprioriVisY)); relGain[spw_index, 1] /= np.median( abs(relGain[spw_index, 1]) ) # Y-pol delta gain
#
#-------- SPW-specific phase using BP scan
GainP = np.array(GainP) # GainP[spw, pol, ant, time]
spwPhase = [0.0]* 2* scnspwNum
for ant_index in range(1,UseAntNum):
    for pol_index in range(2):
        spwPhase = spwPhase + [0.0]
        for spw_index in range(1,scnspwNum): spwPhase = spwPhase + [np.angle(GainP[spw_index, pol_index, ant_index].dot(GainP[0, pol_index, ant_index].conjugate()))]
    #
#
spwPhase = np.array(spwPhase).reshape([UseAntNum, 2, scnspwNum]); spwTwiddle = exp(1.0j *spwPhase)
for spw_index in range(scnspwNum): BPList[spw_index] = (BPList[spw_index].transpose(2,0,1)* spwTwiddle[:,:,spw_index]).transpose(1,2,0)
#-------- Flux Density
ScanFlux, ScanSlope, ErrFlux, centerFreqList = np.zeros([scanNum, scnspwNum]), np.zeros([scanNum, scnspwNum]), np.zeros([scanNum, scnspwNum]), []
for spw_index in range(scnspwNum):
    chNum, chWid, Freq = GetChNum(msfile, scnspw[spw_index])
    centerFreqList.append( np.median(Freq)*1.0e-9 )
#
centerFreqList = np.array(centerFreqList)
print '---Flux densities of sources ---'
pp = PdfPages('FL_' + prefix + '_' + UniqBands[band_index] + '.pdf')
text_sd = ' Scan     Source     EL(deg) '
for spw_index in range(scnspwNum): text_sd = text_sd + ' SPW%02d %5.1f GHz   ' % (scnspw[spw_index], centerFreqList[spw_index])
text_sd = text_sd + '|  mean  %5.1f GHz' % (np.mean(centerFreqList)); logfile.write(text_sd + '\n'); print text_sd
text_sd = ' ------------------------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
for scan_index in range(scanNum):
    figScan = plt.figure(scan_index, figsize = (11, 8))
    figScan.suptitle(prefix + ' ' + UniqBands[band_index])
    figScan.text(0.75, 0.95, qa.time('%fs' % timeStamp[0], form='ymd')[0]) 
    figScan.text(0.45, 0.05, 'Projected baseline [m]') 
    figScan.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90) 
    if sourceIDscan[scan_index] in SSOList: SSO_flag = True
    else:
        SSO_flag = False 
        SAantennas, SAbl, SAblFlag, SAant0, SAant1 = range(UseAntNum), range(UseBlNum), np.ones([blNum]), ant0, ant1
        SAantMap, SAblMap, SAblInv = antMap, blMap, blInv
    if SSO_flag: continue
    SAantNum = len(SAantennas); SAblNum = len(SAblMap)
    text_sd = ' %02d %016s %4.1f' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0* OnEL[scan_index]/pi )
    figScan.text(0.05, 0.95, text_sd) 
    BPCaledXspec = []
    #-------- UV distance
    timeStamp, UVW = GetUVW(msfile, scnspw[0], onsourceScans[scan_index]); uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    for spw_index in range(scnspwNum):
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, scnspw[spw_index], onsourceScans[scan_index])
        timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); UseChNum = len(chRange)
        if np.max(abs(Xspec)) < 1.0e-9: continue
        tempSpec = ParaPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Parallel Polarization Baseline Mapping
        #-------- Bandpass Calibration
        BPCaledXspec = BPCaledXspec + [(tempSpec / (BPList[spw_index][SAant0]* BPList[spw_index][SAant1].conjugate())).transpose(2,3,1,0)] # Bandpass Cal
    #
    #-------- Antenna-based Phase Solution
    BPCaledXspec = np.array(BPCaledXspec)   # BPCaledXspec[spw, pol, ch, bl, time]
    chAvgVis = np.mean(np.array(BPCaledXspec)[:,:,chRange], axis=(0,2)) # chAvgVis[pol, bl, time]
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
    pCalVis = (BPCaledXspec.transpose(0,2,1,3,4) / (GainP[:,SAant0]* GainP[:,SAant1].conjugate()))[:,chRange]
    StokesI = []
    for spw_index in range(scnspwNum):
        chNum, chWid, Freq = GetChNum(msfile, scnspw[spw_index])
        atmCorrect, TA = np.exp(-np.mean(Tau0spec[band_index* spwNum + spw_index])/np.sin(OnEL[scan_index])), 0.0
        SEFD = 2.0* kb* (np.array(chAvgTsys).reshape(NumBands, scanNum, antNum, spwNum, polNum)[band_index,scan_index][SAantMap][:,spw_index] + TA).T / (np.array([AeX[SAantennas], AeY[SAantennas]])* atmCorrect)/ (relGain[spw_index][:,SAantennas]**2)
        AmpCalVis = np.mean(pCalVis[spw_index], axis=(0,3))* np.sqrt(SEFD[:,SAant0]* SEFD[:,SAant1])
        indivRelGain = abs(gainComplexVec(AmpCalVis.T)); indivRelGain /= np.percentile(indivRelGain, 75, axis=0)
        SEFD /= (indivRelGain**2).T
        AmpCalVis = np.mean(pCalVis[spw_index], axis=(0,3))* np.sqrt(SEFD[:,SAant0]* SEFD[:,SAant1])
        StokesI = StokesI + [np.mean(AmpCalVis, axis=0).real]
    #
    StokesI = np.array(StokesI)
    for spw_index in range(scnspwNum):
        visFlag = np.where(abs(StokesI[spw_index] - np.percentile(StokesI[spw_index], 75))/np.percentile(StokesI[spw_index], 75) < 0.2 )[0]
        weight = np.zeros(SAblNum); weight[visFlag] = 1.0/np.var(StokesI[spw_index][visFlag])
        P, W = np.c_[np.ones(SAblNum), uvDist], np.diag(weight)
        PtWP_inv = scipy.linalg.inv(np.dot(P.T, np.dot(W, P)))
        solution, solerr = PtWP_inv.dot(P.T.dot(weight* StokesI[spw_index])),  np.sqrt(np.diag(PtWP_inv))
        if abs(solution[1] < 3.0* solerr[1]): solution[0] = np.median(StokesI[spw_index]); solution[1] = 0.0
        ScanFlux[scan_index, spw_index], ScanSlope[scan_index, spw_index], ErrFlux[scan_index, spw_index] = solution[0], solution[1], solerr[0]
        text_sd = text_sd + '  %6.3f (%.3f) ' % (ScanFlux[scan_index, spw_index], ErrFlux[scan_index, spw_index])
        StokesI_PL = figScan.add_subplot( 2, scnspwNum/2, spw_index + 1 )
        StokesI_PL.plot( uvDist, StokesI[spw_index], 'k.')
        uvMax, IMax = max(uvDist), max(ScanFlux[scan_index])
        StokesI_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index], ScanFlux[scan_index, spw_index]+ uvMax* ScanSlope[scan_index, spw_index]]), '-')
        StokesI_PL.axis([0.0, uvMax, 0.0, 1.25*IMax])
        StokesI_PL.text(0.0, 1.26*IMax, 'SPW%2d %5.1f GHz' % (scnspw[spw_index], centerFreqList[spw_index]))
    #
    #-------- Statistics for all SPWs
    relFreq = centerFreqList - np.mean(centerFreqList)
    sol, solerr = linearRegression(relFreq, ScanFlux[scan_index], ErrFlux[scan_index] )
    text_sd = text_sd + '  | %7.4f (%.4f) ' % (sol[0], solerr[0])
    #
    logfile.write(text_sd + '\n'); print text_sd
    figScan.savefig(pp, format='pdf')
    if COMPDB & (not SSO_flag) : 
        print ' -------- Comparison with ALMA Calibrator Catalog --------'
        au.searchFlux(sourcename='%s' % (sourceList[sourceIDscan[scan_index]]), band=int(UniqBands[band_index][3:5]), date=qa.time('%fs' % (refTime[scan_index]), form='ymd')[0][0:10], maxrows=3)
    #
    logfile.write('')
#
text_sd = ' ------------------------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
plt.close('all')
pp.close()
logfile.close()
np.save(prefix + '-' + UniqBands[band_index] + '.Flux.npy', ScanFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Ferr.npy', ErrFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Source.npy', np.array(sourceList)[sourceIDscan])
np.save(prefix + '-' + UniqBands[band_index] + '.EL.npy', np.array(OnEL))
