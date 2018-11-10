import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
from matplotlib.backends.backend_pdf import PdfPages
RADperHzMeterArcsec = 2.0* pi / 299792458 / (180*3600/pi)
#-------- Configure Array
msmd.open(msfile)
print '---Checking array configulation'
antDia = np.ones(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
print '  -- usable antenna checking for BP scan : '
spwList = scnspw
gainFlag = np.ones([antNum])
for spw_index in range(spwNum):
    #-------- Checking usable baselines and antennas
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPScan)
    timeNum, chNum, blNum = Xspec.shape[3], Xspec.shape[1], Xspec.shape[2]; chRange, timeRange = range(int(0.05*chNum), int(0.95*chNum)), range(timeNum-4, timeNum-1)
    for polID in pPol:
        blD, blA = np.apply_along_axis(delay_search, 0, np.mean(Xspec[polID][chRange][:,:,timeRange], axis=2))
        blA = blA / np.sqrt(antDia[ANT0[0:blNum]]* antDia[ANT1[0:blNum]])
        errD, errA = np.where(abs(blD - np.median(blD)) > 4.0)[0].tolist(), np.where(abs(blA - np.median(blA)) > 0.5* np.median(blA))[0].tolist()
        errCount = np.zeros(antNum)
        for bl in set(errD) or set(errA): errCount[ list(Bl2Ant(bl)) ] += 1
        gainFlag[np.where(errCount > 2.5 )[0].tolist()] *= 0.0
    #
#-------- Load Tsys table
TrxFreq  = np.load(prefix +  '-' + UniqBands[band_index] + '.TrxFreq.npy') # TrxFreq[spw][ch]
TrxAnts  = np.load(prefix +  '-' + UniqBands[band_index] + '.TrxAnt.npy') # TrxAnts[ant]
Tau0spec = np.load(prefix +  '-' + UniqBands[band_index] + '.Tau0.npy') # Tau0spec[spw][ch]
Trxspec  = np.load(prefix +  '-' + UniqBands[band_index] + '.Trx.npy')  # Trxspec[spw, ant, pol, ch]
Tau0E    = np.load(prefix +  '-' + UniqBands[band_index] + '.TauE.npy') # Tau0E[spw, atmScan]
atmTimeRef = np.load(prefix +  '-' + UniqBands[band_index] + '.atmTime.npy') # atmTimeRef[atmScan]
TrxMap = indexList(TrxAnts, antList); TrxFlag = np.zeros([antNum]); TrxFlag[TrxMap] = 1.0
Tau0E = np.nanmedian(Tau0E, axis=0); Tau0E[np.isnan(Tau0E)] = np.nanmedian(Tau0E); Tau0E[np.isnan(Tau0E)] = 0.0
TrxMed = np.median(Trxspec, axis=3)
#-------- Tsys channel interpolation
chNum, chWid, Freq = GetChNum(msfile, spwList[0])
if TrxFreq.shape[1] != chNum:
    tmpTAU0, tmpTRX = np.zeros([spwNum, chNum]), np.zeros([spwNum, antNum, 2, chNum])
    for spw_index in range(spwNum):
        chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq *= 1.0e-9
        TAU0 = interpolate.interp1d(TrxFreq[spw_index], Tau0spec[spw_index])
        tmpTAU0[spw_index] = TAU0(Freq)
        for ant_index in range(len(TrxAnts)):
            for pol_index in range(2):
                TRX = interpolate.interp1d(TrxFreq[spw_index], Trxspec[spw_index, ant_index, pol_index])
                tmpTRX[spw_index, ant_index, pol_index] = TRX(Freq)
        #
    #
    Tau0spec = tmpTAU0
    Trxspec  = tmpTRX
#
for spw_index in range(spwNum):
    for pol_index in range(2): TrxFlag[np.where(abs(TrxMed[spw_index][:,pol_index] - np.median(TrxMed[spw_index][:,pol_index])) > 0.8* np.median(TrxMed[spw_index][:,pol_index]))[0].tolist()] *= 0.0
if np.min(np.median(Tau0spec[:,chRange], axis=1)) < 0.0: TrxFlag *= 0.0    # Negative Tau(zenith) 
#
print 'Ant:',
for ant_index in range(antNum): print antList[ant_index],
print; print 'givn',
for ant_index in range(antNum): print '   %.0f' % (flagAnt[ant_index]),
print; print 'Trx ',
for ant_index in range(antNum): print '   %.0f' % (TrxFlag[ant_index]),
print; print 'gain',
for ant_index in range(antNum): print '   %.0f' % (gainFlag[ant_index]),
print
flagAnt = flagAnt* TrxFlag* gainFlag
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
if len(UseAnt) < 5: sys.exit('Too few usable antennas. Reduction failed.')
#-------- Check Scans for atmCal
ingestFile = open(prefix + '-' + UniqBands[band_index] + '-Ingest.log', 'w')
text_sd = '#source,   RA,eRA,dec,edec,frequency,  flux,eflux,degree,edeg,EVPA,eEVPA,uvmin,uvmax,date\n'; ingestFile.write(text_sd)
logfile = open(prefix + '-' + UniqBands[band_index] + '-Flux.log', 'w')
logfile.write(BPcalText + '\n'); logfile.write(EQcalText + '\n')
scanList = onsourceScans
msmd.close()
msmd.done()
#-------- Array Configuration
print '---Determining refant'
msmd.open(msfile)
timeStamp, UVW = GetUVW(msfile, spwList[0], msmd.scansforspw(spwList[0])[0])
uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
refantID = bestRefant(uvDist, UseAnt)
print '  Use ' + antList[refantID] + ' as the refant.'
#
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
if 'gainRef' in locals(): flagRef = np.zeros([UseAntNum]); refIndex = indexList(gainRef, antList[antMap]); flagRef[refIndex] = 1.0; del(gainRef)
if 'refIndex' not in locals(): refIndex = range(UseAntNum)
if len(refIndex) == 0: refIndex = range(UseAntNum)
#-------- Load Aeff file
Afile = open(SCR_DIR + 'AeB' + `int(UniqBands[band_index][3:5])` + '.data')
Alines = Afile.readlines()
Afile.close()
AeX, AeY = 0.25* np.pi* antDia**2, 0.25* np.pi* antDia**2       # antenna collecting area (100% efficiency)
AeX, AeY, etaX, etaY = [], [], [], []
for ant_index in antMap:
    for Aline in Alines:
        if antList[ant_index] in Aline:
            etaX = etaX + [float(Aline.split()[1])]
            etaY = etaY + [float(Aline.split()[2])]
            AeX = AeX + [(0.0025* np.pi* float(Aline.split()[1]))* antDia[ant_index]**2]
            AeY = AeY + [(0.0025* np.pi* float(Aline.split()[2]))* antDia[ant_index]**2]
            #print '%s  : etaX = %.2f  etaY = %.2f' % (antList[ant_index], float(Aline.split()[1]), float(Aline.split()[2]))
        #
    #
#
Ae = np.array([AeX, AeY])
Trx2antMap = indexList( antList[antMap], antList[TrxMap] )
#AeX, AeY = np.array(AeX), np.array(AeY) # in antMap order 
#-------- Flag table
if 'FGprefix' in locals():
    print '---Checking Flag File'
    FGList = []
    for spw_index in range(spwNum): FG = np.load(FGprefix + '-SPW' + `spwList[spw_index]` + '.FG.npy'); FGList = FGList + [np.min(FG, axis=0)]
    FG = np.min( np.array(FGList), axis=0)
    TS = np.load(FGprefix + '-SPW' + `spwList[spw_index]` + '.TS.npy')
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
    flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
else :
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
    flagIndex = range(timeNum)
#
BPList = []
print '---Generating antenna-based bandpass table'
for spw_index in range(spwNum):
    BP_ant, XY_BP, XYdelay, Gain = BPtable(msfile, spwList[spw_index], BPScan, blMap, blInv)
    exp_Tau = np.exp(-Tau0spec[spw_index] / np.sin(BPEL))
    atmCorrect = 1.0 / exp_Tau
    TsysBPScan = atmCorrect* (Trxspec[spw_index][Trx2antMap] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau)) # [antMap, pol, ch]
    TsysBPShape = (TsysBPScan.transpose(2,0,1) / np.median(TsysBPScan, axis=2)).transpose(1,2,0)
    BPList = BPList + [BP_ant* np.sqrt(TsysBPShape)]
#
if PLOTBP:
    pp = PdfPages('BP_' + prefix + '_REF' + antList[UseAnt[refantID]] + '_Scan' + `BPScan` + '.pdf')
    plotBP(pp, prefix, antList[antMap], spwList, BPScan, BPList)
#
BPDone = True
##-------- Equalization using EQ scan
scanList = onsourceScans
relGain = np.ones([spwNum, 2, UseAntNum])
polXindex, polYindex, scan_index = (arange(4)//2).tolist(), (arange(4)%2).tolist(), scanList.index(EQScan)
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refantID, AZ, EL)
if len(atmTimeRef) > 5:
    exTauSP  = UnivariateSpline(atmTimeRef, Tau0E, np.ones(len(atmTimeRef)), s=0.1*np.std(Tau0E))
else:
    tempTime = np.arange(np.min(atmTimeRef) - 3600.0,  np.max(atmTimeRef) + 3600.0, 300.0)
    tempTauE = np.repeat(np.median(Tau0E), len(tempTime))
    exTauSP = UnivariateSpline(tempTime, tempTauE, np.ones(len(tempTime)), s=0.1)
#
for spw_index in range(spwNum):
    exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
    TsysEQScan = np.mean(Trxspec[spw_index,:,:,chRange].transpose(1,2,0) + Tcmb*exp_Tau[chRange] + tempAtm* (1.0 - exp_Tau[chRange]), axis=2)[Trx2antMap] # [antMap, pol]
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = ParaPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)[flagIndex]       # Cross Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0]* BPList[spw_index][ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
    pCaledVis = np.array([chAvgVis[0] / (GainP[0,ant0]* GainP[0,ant1].conjugate()), chAvgVis[1]/(GainP[1,ant0]* GainP[1,ant1].conjugate())])
    aprioriSEFD = 2.0* kb* TsysEQScan.T / Ae
    #
    aprioriVisX = np.mean(pCaledVis[0], axis=1) * np.sqrt(aprioriSEFD[0, ant0]* aprioriSEFD[0, ant1])
    aprioriVisY = np.mean(pCaledVis[1], axis=1) * np.sqrt(aprioriSEFD[1, ant0]* aprioriSEFD[1, ant1])
    #-------- Determine Antenna-based Gain
    relGain[spw_index, 0] = abs(gainComplex(aprioriVisX)); relGain[spw_index, 0] /= np.median( abs(relGain[spw_index, 0, refIndex]) ) # X-pol delta gain
    relGain[spw_index, 1] = abs(gainComplex(aprioriVisY)); relGain[spw_index, 1] /= np.median( abs(relGain[spw_index, 1, refIndex]) ) # Y-pol delta gain
#
print '---Equalized aperture efficiencies (Pol-X, Pol-Y) in %'
antRelGain = np.median(relGain, axis=0)
for ant_index in range(UseAntNum):
    text_sd = '%s  %.2f  %.2f' % (antList[antMap[ant_index]], etaX[ant_index]* antRelGain[0,ant_index]**2, etaY[ant_index]* antRelGain[1,ant_index]**2)
    logfile.write(text_sd + '\n'); print text_sd
#
#-------- SPW-specific phase using BP scan
GainP = []
scan_index = scanList.index(BPScan)
for spw_index in range(spwNum):
    exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
    atmCorrect = 1.0 / exp_Tau
    TsysSPW = (Trxspec[spw_index] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))[Trx2antMap] # [antMap, pol, ch]
    TsysBL  = np.sqrt( TsysSPW[ant0]* TsysSPW[ant1])
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPScan); timeNum = len(timeStamp)
    if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
    else : flagIndex = range(timeNum)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = ParaPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec * TsysBL/ (BPList[spw_index][ant0]* BPList[spw_index][ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainP = GainP + [np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])]
    SEFD = 2.0* kb / (np.array([AeX, AeY]) * (relGain[spw_index]**2))
    caledVis.append(np.mean((chAvgVis / (GainP[spw_index][:,ant0]* GainP[spw_index][:,ant1].conjugate())).transpose(2, 0, 1)* np.sqrt(SEFD[:,ant0]* SEFD[:,ant1]), axis=2).T)
#
caledVis = np.array(caledVis)
GainP = np.array(GainP) # GainP[spw, pol, ant, time]
spwPhase = [0.0]* 2* spwNum
for ant_index in range(1,UseAntNum):
    for pol_index in range(2):
        spwPhase = spwPhase + [0.0]
        for spw_index in range(1,spwNum): spwPhase = spwPhase + [np.angle(GainP[spw_index, pol_index, ant_index].dot(GainP[0, pol_index, ant_index].conjugate()))]
    #
#
spwPhase = np.array(spwPhase).reshape([UseAntNum, 2, spwNum]); spwTwiddle = exp(1.0j *spwPhase)
for spw_index in range(spwNum): BPList[spw_index] = (BPList[spw_index].transpose(2,0,1)* spwTwiddle[:,:,spw_index]).transpose(1,2,0)
ScanFlux, ScanSlope, ErrFlux, centerFreqList, FreqList = np.zeros([scanNum, spwNum, 4]), np.zeros([scanNum, spwNum, 4]), np.zeros([scanNum, spwNum, 4]), [], []
#-------- Flux Density
for spw_index in range(spwNum):
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index])
    centerFreqList.append( np.median(Freq)*1.0e-9 )
    FreqList = FreqList + [Freq* 1.0e-9]
#
print '---Flux densities of sources ---'
pp = PdfPages('SP_' + prefix + '_' + UniqBands[band_index] + '.pdf')
figSP = plt.figure(figsize = (11, 8))
figSP.suptitle(prefix + ' ' + UniqBands[band_index])
figSP.text(0.45, 0.05, 'Frequency [GHz]')
figSP.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90)
for scan_index in range(scanNum):
    scan = scanList[scan_index]
    if scan_index > 0:
        for SP in IList: figSP.delaxes(SP)
        for SP in PList: figSP.delaxes(SP)
    #
    timeStamp, UVW = GetUVW(msfile, spwList[0], scan);  timeNum = len(timeStamp)
    uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    #-------- Flagging
    if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
    else: flagIndex = range(timeNum)
    #-------- Parallactic Angle
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    OnAZ[scan_index], OnEL[scan_index] = np.median(AzScan), np.median(ElScan)
    OnPA[scan_index] = AzEl2PA(OnAZ[scan_index], OnEL[scan_index])
    scanTime = scanTime + [np.median(timeStamp)]
    BPCaledXspec, IList, PList = [], [], []
    text_time = qa.time('%fs' % scanTime[scan_index], form='ymd')[0]
    text_src  = ' %02d %010s EL=%4.1f deg' % (scanList[scan_index], sourceList[sourceIDscan[scan_index]], 180.0* OnEL[scan_index]/pi) # ; logfile.write(text_src + ' ' + text_time + '\n'); print text_src + ' ' + text_time
    #-------- Subarray formation
    SAantMap, SAblMap, SAblInv, SAant0, SAant1 = antMap, blMap, blInv, ant0, ant1
    bpAntMap = indexList(antList[SAantMap],antList[antMap])
    Trx2antMap = indexList( antList[SAantMap], antList[TrxMap] )
    #-------- Baseline-based cross power spectra
    for spw_index in range(spwNum):
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], scan)
        timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); UseChNum = len(chRange)
        if np.max(abs(Xspec)) < 1.0e-9: continue
        #-------- Position offset phase correction
        if 'offAxis' in locals():
            lm = np.array(offAxis[scan])
            Twiddle =  np.exp((0.0 + 1.0j)* np.outer(FreqList[spw_index]*1.0e9, uvw[0:2].transpose(1,2,0).dot(lm)).reshape([chNum, UseBlNum, timeNum])* RADperHzMeterArcsec)
            tempSpec = ParaPolBL(Xspec[:,:,SAblMap]*Twiddle, SAblInv).transpose(3,2,0,1)[flagIndex]
   # Cross Polarization Baseline Mapping
        else:
            tempSpec = ParaPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)[flagIndex]      # Cross Polarization Baseline Mapping
        #-------- Bandpass Calibration
        BPCaledXspec = BPCaledXspec + [(tempSpec / (BPList[spw_index][SAant0]* BPList[spw_index][SAant1].conjugate())).transpose(2,3,1,0)]
    #
    #-------- Antenna-based Phase Solution
    BPCaledXspec = np.array(BPCaledXspec)   # BPCaledXspec[spw, pol, ch, bl, time]
    chAvgVis = np.mean(np.array(BPCaledXspec)[:,:,chRange], axis=(0,2)) # chAvgVis[pol, bl, time]
    if 'timeBunch' in locals():
        useTimeNum = timeNum / timeBunch * timeBunch
        leapNum = timeNum - useTimeNum
        timeAvgVis = np.mean(chAvgVis[:,:,range(useTimeNum)].reshape(2, UseBlNum, timeNum / timeBunch, timeBunch), axis=3)
        GainP = np.array([np.apply_along_axis(clphase_solve, 0, timeAvgVis[0]), np.apply_along_axis(clphase_solve, 0, timeAvgVis[1])]).repeat(timeBunch, axis=2)
        if leapNum > 0: GainP = np.append(GainP,  GainP[:,:,(useTimeNum-1):(useTimeNum)].repeat(leapNum, axis=2), axis=2)
    else:
        GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
    #
    pCalVis = (BPCaledXspec.transpose(0,2,1,3,4) / (GainP[:,SAant0]* GainP[:,SAant1].conjugate()))[:,chRange]
    #-------- Stokes I
    IMax = 0.0
    for spw_index in range(spwNum):
        SpecI, SpecE = np.zeros([UseChNum]), np.zeros([UseChNum])
        StokesI_SP = figSP.add_subplot( 1, spwNum, spw_index + 1 )
        IList = IList + [StokesI_SP]
        exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
        atmCorrect = 1.0 / exp_Tau
        TsysSPW = (Trxspec[spw_index] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))[Trx2antMap]    # [ant, pol, ch]
        SEFD = 2.0* kb* (TsysSPW * atmCorrect).transpose(2,1,0) / Ae[:,bpAntMap]                # SEFD[ch,pol,antMap]
        SAantNum = len(SAantMap); SAblNum = len(SAblMap)
        AmpCalVis = np.mean(np.mean(pCalVis[spw_index], axis=3) * np.sqrt(SEFD[chRange][:,:,ant0[0:SAblNum]]* SEFD[chRange][:,:,ant1[0:SAblNum]]), axis=0)
        indivRelGain = abs(gainComplexVec(AmpCalVis.T)); indivRelGain /= np.percentile(indivRelGain, 75, axis=0)
        SEFD /= (indivRelGain**2).T
        AmpCalVis = np.mean((pCalVis[spw_index].transpose(3,0,1,2)* np.sqrt(SEFD[chRange][:,:,SAant0]* SEFD[chRange][:,:,SAant1])), axis=1)        # AmpCalVis[ch, bl]
        for ch_index in range(UseChNum):
            medVis, sdVis = np.median(AmpCalVis[ch_index].real), np.std(AmpCalVis[ch_index].imag)
            visFlag = np.where(abs(AmpCalVis[ch_index].real - medVis) < 4.0* sdVis )[0].tolist()
            weight = np.zeros(SAblNum); weight[visFlag] = 1.0/(sdVis**2)
            P, W = np.c_[np.ones(SAblNum), uvDist[SAblMap]], np.diag(weight)
            PtWP_inv = scipy.linalg.inv(P.T.dot(W.dot(P)))
            solution, solerr = PtWP_inv.dot(P.T.dot(weight* AmpCalVis[ch_index].real)),  np.sqrt(np.diag(PtWP_inv)) # solution[0]:intercept, solution[1]:slope
            SpecI[ch_index], SpecE[ch_index] = solution[0], solerr[0]
        #
        IMax = max(IMax, np.max(SpecI))
        np.save(prefix + '-' + UniqBands[band_index] + '-SPW' + `spwList[spw_index]` + '-Scan' + `scan` + '.SpecI.npy', SpecI)
        np.save(prefix + '-' + UniqBands[band_index] + '-SPW' + `spwList[spw_index]` + '-Scan' + `scan` + '.SpecE.npy', SpecE)
        np.save(prefix + '-' + UniqBands[band_index] + '-SPW' + `spwList[spw_index]` + '-Scan' + `scan` + '.Freq.npy', np.array(FreqList)[:,chRange])
        StokesI_SP.plot(FreqList[spw_index][chRange], StokesSpec[0], ls='steps-mid', label=polLabel[0], color=Pcolor[0])
    #
    for spw_index in range(spwNum):
        StokesI_SP = IList[spw_index]
        StokesI_SP.tick_params(axis='both', labelsize=6)
        StokesI_SP.axis([min(FreqList[spw_index]), max(FreqList[spw_index]), 0.0, 1.25*IMax])
    #
    IList[0].text(min(FreqList[0]), IMax*1.35, text_src)
    IList[-1].text(min(FreqList[-1]), IMax*1.35, text_time)
    StokesI_SP.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    plt.show()
    figSP.savefig(pp, format='pdf')
#
logfile.close()
plt.close('all')
pp.close()
