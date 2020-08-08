import sys
from scipy import stats
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
from matplotlib.backends.backend_pdf import PdfPages
RADperHzMeterArcsec = 2.0* pi / 299792458 / (180*3600/pi)
fluxCalText = 'SEFD'
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
        errD = np.where(abs(blD - np.median(blD)) > 4.0)[0].tolist()
        errA = np.where(blA / np.median(blA) > 2.5)[0].tolist() + np.where(blA / np.median(blA) < 0.4)[0].tolist()
        errCount = np.zeros(antNum)
        for bl in set(errD) or set(errA): errCount[ list(Bl2Ant(bl)) ] += 1
        gainFlag[np.where(errCount > 2.5 )[0].tolist()] *= 0.0
    #
#
#-------- Load D-term file
Dcat = GetDterm(TBL_DIR, antList,int(UniqBands[band_index][3:5]), np.mean(timeStamp))
#-------- Load Tsys table
Tau0spec, TrxList, TrxFreqList = [], [], []
for spw in spwList:
    TrxList = TrxList + [np.median(np.load(prefix +  '-' + UniqBands[band_index] + '-SPW' + `spw` + '.Trx.npy'), axis=3)]  # TrxList[spw][pol, ch, ant]
    Tau0spec = Tau0spec + [np.load(prefix +  '-' + UniqBands[band_index] + '-SPW' + `spw` + '.Tau0.npy')]  # Tau0spec[spw][ch]
    TrxFreqList  = TrxFreqList + [np.load(prefix +  '-' + UniqBands[band_index] + '-SPW' + `spw` + '.TrxFreq.npy')] # TrxFreq[spw][ch]
#
TrxAnts  = np.load(prefix +  '-' + UniqBands[band_index] + '.TrxAnt.npy') # TrxAnts[ant]
Tau0E    = np.load(prefix +  '-' + UniqBands[band_index] + '.TauE.npy') # Tau0E[spw, atmScan]
atmTimeRef = np.load(prefix +  '-' + UniqBands[band_index] + '.atmTime.npy') # atmTimeRef[atmScan]
TrxMap = indexList(TrxAnts, antList); TrxFlag = np.zeros([antNum]); TrxFlag[TrxMap] = 1.0
Tau0E = np.nanmedian(Tau0E, axis=0); Tau0E[np.isnan(Tau0E)] = np.nanmedian(Tau0E); Tau0E[np.isnan(Tau0E)] = 0.0
#-------- Tsys channel interpolation
'''
chNum, chWid, Freq = GetChNum(msfile, spwList[0])
for spw_index in range(spwNum):
    if len(TrxFreqList[spw_index]) != chNum: 
        tmpTAU0 = np.zeros(chNum)
        tmpTRX = np.zeros([antNum, 2, chNum])
        chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq *= 1.0e-9
        TAU0 = interpolate.interp1d(TrxFreq[spw_index], Tau0spec[spw_index])
        tmpTAU0 = TAU0(Freq)
        for ant_index in range(len(TrxAnts)):
            for pol_index in range(2):
                TRX = interpolate.interp1d(TrxFreq[spw_index], TrxList[spw_index][pol_index, :, ant_index])
                tmpTRX[ant_index, pol_index] = TRX(Freq)
        #
        TrxList[spw_index]  = tmpTRX.transpose(1,2,0)
        Tau0spec[spw_index] = tmpTAU0
    #
#
'''
for spw_index in range(spwNum):
    TrxMed = np.median(TrxList[spw_index], axis=1)  # TrxMed[pol, ant]
    for pol_index in range(2):
        Trx2MedianRatio = TrxMed[pol_index] / np.median(TrxMed[pol_index])
        TrxFlag[ np.where(Trx2MedianRatio < 0.1)[0].tolist() ] *= 0.0   # Flagged by negative Trx
        TrxFlag[ np.where(Trx2MedianRatio > 4.0)[0].tolist() ] *= 0.0   # Flagged by too-high Trx 
    if np.median(Tau0spec[spw_index][chRange]) < 0.0: TrxFlag *= 0.0    # Negative Tau(zenith)
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
UseAnt = np.where(flagAnt > 0.0)[0]; UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
print `UseAntNum` + ' / ' + `antNum` + ' usable antennas'
if len(UseAnt) < 4:
    sys.exit('Too few usable antennas. Reduction failed.')
#-------- Check Scans for atmCal
ingestFile = open(prefix + '-' + UniqBands[band_index] + '-Ingest.log', 'w')
text_sd = '#source,    RA,eRA,dec,edec, frequency,  flux,  eflux,     %P,    d%P,  EVPA, eEVPA, uvmin, uvmax,         date, fluxCal, ExecBlock\n'; ingestFile.write(text_sd)
logfile = open(prefix + '-' + UniqBands[band_index] + '-Flux.log', 'w')
if 'BPcalText' in locals(): logfile.write(BPcalText + '\n')
if 'EQcalText' in locals(): logfile.write(EQcalText + '\n')
scanList = onsourceScans
msmd.close()
msmd.done()
#-------- Array Configuration
print '---Determining refant'
if 'refant' in locals(): refantID = np.where(antList == refant)[0][0]
if 'refantID' not in locals():
    msmd.open(msfile)
    timeStamp, UVW = GetUVW(msfile, spwList[0], msmd.scansforspw(spwList[0])[0])
    uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = UseAnt[bestRefant(uvDist, UseAnt)]
print '  Use ' + antList[refantID] + ' as the refant.'
#
antMap = [refantID] + UseAnt[np.where(UseAnt != refantID)].tolist()
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
if 'gainRef' in locals(): flagRef = np.zeros([UseAntNum]); refIndex = indexList(gainRef, antList[antMap]); flagRef[refIndex] = 1.0; del(gainRef)
if 'refIndex' not in locals(): refIndex = range(UseAntNum)
if len(refIndex) == 0: refIndex = range(UseAntNum)
#-------- Load Aeff file
etaA = GetAeff(TBL_DIR, antList[antMap], int(UniqBands[band_index][3:5]), np.mean(timeStamp)).T
Ae = 0.0025* np.pi* etaA* antDia[antMap]**2
if BLCORR: Ae = 1.15* Ae         # BL Correlator correction factor
#
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
#-------- Bandpass Table
Trx2antMap = indexList( antList[antMap], antList[TrxMap] )
BPList = []
print '---Generating antenna-based bandpass table'
for spw_index in range(spwNum):
    BP_ant, XY_BP, XYdelay, Gain = BPtable(msfile, spwList[spw_index], BPScan, blMap, blInv)
    BP_ant[:,1] *= XY_BP
    exp_Tau = np.exp(-Tau0spec[spw_index] / np.sin(BPEL))
    atmCorrect = 1.0 / exp_Tau
    TsysBPScan = atmCorrect* (TrxList[spw_index].transpose(2,0,1)[Trx2antMap] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau)) # [antMap, pol, ch]
    TsysBPShape = (TsysBPScan.transpose(2,0,1) / np.median(TsysBPScan, axis=2)).transpose(1,2,0)
    BPList = BPList + [BP_ant* np.sqrt(TsysBPShape)]
#
if PLOTBP:
    pp = PdfPages('BP_' + prefix + '_REF' + antList[refantID] + '_Scan' + `BPScan` + '.pdf')
    plotBP(pp, prefix, antList[antMap], spwList, BPScan, BPList)
#
BPDone = True
##-------- Equalization using EQ scan
spwStokesDic = dict(zip(sourceList, [[]]*len(sourceList))) # StokesDic[sourceName][pol* spw]
scanList = onsourceScans
relGain = np.ones([spwNum, 2, UseAntNum])
polXindex, polYindex, scan_index = (arange(4)//2).tolist(), (arange(4)%2).tolist(), scanList.index(EQScan)
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PA = np.arctan2( np.sin(PA), np.cos(PA))
if len(StokesDic[EQcal]) > 0:
    StokesEQ = np.array(StokesDic[EQcal])
else:
    StokesEQ = np.array([1.0, 0.0, 0.0, 0.0])
#
QCpUS = (StokesEQ[1]* np.cos(2.0* PA) + StokesEQ[2]* np.sin(2.0* PA)) / StokesEQ[0]
##-------- Smooth excess Tau 
exTauSP = tauSMTH(atmTimeRef, Tau0E)
##-------- Equalize aperture efficiency
for spw_index in range(spwNum):
    exp_Tau = np.exp(-(Tau0spec[spw_index] + scipy.interpolate.splev(np.median(timeStamp), exTauSP) ) / np.mean(np.sin(ElScan)))
    TsysEQScan = np.mean(TrxList[spw_index].transpose(2,0,1)[:,:,chRange] + Tcmb*exp_Tau[chRange] + tempAtm* (1.0 - exp_Tau[chRange]), axis=2)[Trx2antMap] # [antMap, pol]
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)[flagIndex]       # Cross Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)[pPol]
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
    pCaledVis = np.array([chAvgVis[0] / (GainP[0,ant0]* GainP[0,ant1].conjugate()), chAvgVis[1]/(GainP[1,ant0]* GainP[1,ant1].conjugate())])
    aprioriSEFD = 2.0* kb* TsysEQScan.T / Ae
    #
    aprioriVisX = np.mean(pCaledVis[0] / (1.0 + QCpUS), axis=1) * np.sqrt(aprioriSEFD[0, ant0]* aprioriSEFD[0, ant1])
    aprioriVisY = np.mean(pCaledVis[1] / (1.0 - QCpUS), axis=1) * np.sqrt(aprioriSEFD[1, ant0]* aprioriSEFD[1, ant1])
    #-------- Determine Antenna-based Gain
    relGain[spw_index, 0] = abs(gainComplex(aprioriVisX)); relGain[spw_index, 0] /= np.median( abs(relGain[spw_index, 0, refIndex]) ) # X-pol delta gain
    relGain[spw_index, 1] = abs(gainComplex(aprioriVisY)); relGain[spw_index, 1] /= np.median( abs(relGain[spw_index, 1, refIndex]) ) # Y-pol delta gain
#
print '---Equalized aperture efficiencies (Pol-X, Pol-Y) in %'
antRelGain = np.median(relGain, axis=0)
for ant_index in range(UseAntNum):
    text_sd = '%s  %.2f  %.2f' % (antList[antMap[ant_index]], etaA[0,ant_index]* antRelGain[0,ant_index]**2, etaA[1,ant_index]* antRelGain[1,ant_index]**2)
    logfile.write(text_sd + '\n'); print text_sd
#
##-------- Iteration for Equalization using EQ scan
print '---XY phase determination in bandpass scan'
#-------- XY phase using BP scan
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], BPScan); timeNum = len(timeStamp)
if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
else: flagIndex = range(timeNum)
#
AzScan, ElScan = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PA = np.arctan2( np.sin(PA), np.cos(PA))
GainP, XYphase, caledVis = [], [], []
DtermDic = {'mjdSec': np.median(timeStamp)}
StokesDic['mjdSec'] = np.median(timeStamp)
Trx2antMap = indexList( antList[antMap], antList[TrxMap] )
scan_index = scanList.index(BPScan)
for spw_index in range(spwNum):
    exp_Tau = np.exp(-(Tau0spec[spw_index] + scipy.interpolate.splev(np.median(timeStamp), exTauSP) ) / np.mean(np.sin(ElScan)))
    atmCorrect = 1.0 / exp_Tau
    TsysSPW = (TrxList[spw_index].transpose(2,0,1) + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))[Trx2antMap] # [antMap, pol, ch]
    TsysBL  = np.sqrt( TsysSPW[ant0][:,polYindex]* TsysSPW[ant1][:,polXindex])
    #
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPScan); timeNum = len(timeStamp)
    if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
    else : flagIndex = range(timeNum)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec * TsysBL/ (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainP = GainP + [np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])]
    SEFD = 2.0* kb / (Ae * (relGain[spw_index]**2))
    caledVis.append(np.mean((chAvgVis / (GainP[spw_index][polYindex][:,ant0]* GainP[spw_index][polXindex][:,ant1].conjugate())).transpose(2, 0, 1)* np.sqrt(SEFD[polYindex][:,ant0]* SEFD[polXindex][:,ant1]), axis=2).T)
#
caledVis = np.array(caledVis)
#-------- SPW-specific phase using BP scan
GainP = np.array(GainP) # GainP[spw, pol, ant, time]
spwPhase = [0.0]* 2* spwNum
for ant_index in range(1,UseAntNum):
    for pol_index in range(2):
        spwPhase = spwPhase + [0.0]
        for spw_index in range(1,spwNum): spwPhase = spwPhase + [np.angle(GainP[spw_index, pol_index, ant_index].dot(GainP[0, pol_index, ant_index].conjugate()))]
    #
#
spwPhase = np.array(spwPhase).reshape([UseAntNum, 2, spwNum]); spwTwiddle = exp(1.0j *spwPhase)
QUsolution = XXYY2QU(PA, np.mean(caledVis[:,[0,3]], axis=0))
if 'QUMODEL' in locals():
    StokesBP = np.array(StokesDic[BPcal])
    if QUMODEL: QUsolution = np.array(StokesBP[[1,2]])/StokesBP[0]
else: QUsolution = XXYY2QU(PA, np.mean(caledVis[:,[0,3]], axis=0))
#-------- XY phase cal in Bandpass table
XYsign = np.ones(spwNum)
for spw_index in range(spwNum):
    XYphase = XY2Phase(QUsolution[1]* np.cos(2.0* PA) - QUsolution[0]* np.sin(2.0* PA), caledVis[spw_index][[1,2]])
    XYsign[spw_index] = np.sign(np.cos(XYphase))
    BPList[spw_index] = (BPList[spw_index].transpose(2,0,1)* spwTwiddle[:,:,spw_index]).transpose(1,2,0)
    BPList[spw_index][:,1] *= XYsign[spw_index]
    print 'SPW[%d] : XYphase = %6.1f [deg] sign = %3.0f' % (spwList[spw_index], 180.0*XYphase/pi, XYsign[spw_index])
#
#-------- Gain Equalization between X and Y
# PolEqualize.py is obsolete. The polarization gain equalization is already implemented in this script
#if 'PolEQ' in locals():
#    if PolEQ: execfile(SCR_DIR + 'PolEqualize.py')
#
#-------- Flux Density
ScanFlux, ScanSlope, ErrFlux, centerFreqList, FreqList = np.zeros([scanNum, spwNum, 4]), np.zeros([scanNum, spwNum, 4]), np.zeros([scanNum, spwNum, 4]), [], []
for spw_index in range(spwNum):
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index])
    centerFreqList.append( np.median(Freq)*1.0e-9 )
    FreqList = FreqList + [Freq]
#
print '---Flux densities of sources ---'
pp, polLabel, Pcolor = PdfPages('FL_' + prefix + '_' + UniqBands[band_index] + '.pdf'), ['I', 'Q', 'U', 'V'], ['black', 'blue', 'red', 'green']
XYD, XYC, scanTime = [], [], []     # XY delay and correlation
scanDic   = dict(zip(sourceList, [[]]*len(sourceList))) # Scan list index for each source
timeDic   = dict(zip(sourceList, [[]]*len(sourceList))) # Time index list for each source
PADic     = dict(zip(sourceList, [[]]*len(sourceList))) # PA list for each source
figFL = plt.figure(figsize = (11, 8))
figFL.suptitle(prefix + ' ' + UniqBands[band_index])
figFL.text(0.45, 0.05, 'Projected baseline [m]')
figFL.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90)
AmpCalChAvg = np.zeros([spwNum, 4, UseBlNum, timeSum], dtype=complex)
timePointer = 0
for scan_index in range(scanNum):
    scanFlag, DcalFlag = True, True
    sourceName = sourceList[sourceIDscan[scan_index]]
    scanDic[sourceName] = scanDic[sourceName] + [scan_index, 1]
    if scan_index > 0:
        for PL in IList: figFL.delaxes(PL)
        for PL in PList: figFL.delaxes(PL)
    SSO_flag = False
    if sourceName in SSOCatalog: SSO_flag = True; scanDic[sourceName][1] *= 0; DcalFlag = False
    #-------- UV distance
    timeStamp, UVW = GetUVW(msfile, spwList[0], scanList[scan_index]);  timeNum = len(timeStamp)
    scanTime = scanTime + [np.median(timeStamp)]
    timeNum = len(timeStamp)
    timeDic[sourceName] = range(timePointer, timePointer + timeNum)
    uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    #-------- Flagging
    if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
    else: flagIndex = range(timeNum)
    #-------- Parallactic Angle
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    if np.max(ElScan) == 0.0: AzScan, ElScan = AzElMatch(timeStamp, azelTime, 0, refantID, AZ, EL)
    PA = (AzEl2PA(AzScan, ElScan) + BandPA[band_index])[flagIndex]; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
    PADic[sourceName] = PA.tolist()
    #-------- Plot Frame
    IList, PList = [], []      # XY delay and correlation
    timeLabel = qa.time('%fs' % np.median(timeStamp), form='ymd')[0]
    text_src  = ' %02d %010s EL=%4.1f deg' % (scanList[scan_index], sourceName, 180.0* OnEL[scan_index]/pi)
    BPCaledXspec = []
    #-------- Subarray formation
    SAantMap, SAblMap, SAblInv, SAant0, SAant1 = antMap, blMap, blInv, ant0, ant1
    bpAntMap = indexList(antList[SAantMap],antList[antMap])
    Trx2antMap = indexList( antList[SAantMap], antList[TrxMap] )
    #-------- Baseline-based cross power spectra
    for spw_index in range(spwNum):
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], scanList[scan_index])
        timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); UseChNum = len(chRange)
        if np.max(abs(Xspec)) < 1.0e-9: continue
        #-------- Position offset phase correction
        if 'offAxis' in locals():
            lm = np.array(offAxis[sourceIDscan[scan_index]])
            Twiddle =  np.exp((0.0 + 1.0j)* np.outer(FreqList[spw_index], uvw[0:2].transpose(1,2,0).dot(lm)).reshape([chNum, UseBlNum, timeNum])* RADperHzMeterArcsec)
            tempSpec = CrossPolBL(Xspec[:,:,SAblMap]*Twiddle, SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        else:
            tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        #-------- Bandpass Calibration
        BPCaledXspec = BPCaledXspec + [(tempSpec / (BPList[spw_index][SAant0][:,polYindex]* BPList[spw_index][SAant1][:,polXindex].conjugate())).transpose(2,3,1,0)] 
    #
    #-------- Antenna-based Phase Solution
    BPCaledXspec = np.array(BPCaledXspec)   # BPCaledXspec[spw, pol, ch, bl, time]
    chAvgVis = np.mean(np.array(BPCaledXspec)[:,:,chRange], axis=(0,2)) # chAvgVis[pol, bl, time]
    if 'timeBunch' in locals():
        useTimeNum = timeNum / timeBunch * timeBunch
        leapNum = timeNum - useTimeNum
        timeAvgVis = np.mean(chAvgVis[:,:,range(useTimeNum)].reshape(polNum, UseBlNum, timeNum / timeBunch, timeBunch), axis=3)
        GainP = np.array([np.apply_along_axis(clphase_solve, 0, timeAvgVis[0]), np.apply_along_axis(clphase_solve, 0, timeAvgVis[3])]).repeat(timeBunch, axis=2)
        if leapNum > 0: GainP = np.append(GainP,  GainP[:,:,(useTimeNum-1):(useTimeNum)].repeat(leapNum, axis=2), axis=2)
    else:
        GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    #
    pCalVis = (BPCaledXspec.transpose(0,2,1,3,4) / (GainP[polYindex][:,SAant0]* GainP[polXindex][:,SAant1].conjugate()))[:,chRange]
    #-------- XY phase spectra
    for spw_index in range(spwNum):
        delayFact = (chNum + 0.0)/len(chRange)
        XYspec = np.mean(pCalVis[spw_index, :, 1:3, :], axis=(2,3)).T
        XYdelay, XYamp = delay_search(XYspec[:,0]); YXdelay, YXamp = delay_search(XYspec[:,1])
        XYD = XYD + [XYdelay* delayFact, YXdelay* delayFact]
        XYC = XYC + [np.mean(delay_cal(XYspec[:,0], XYdelay)), np.mean(delay_cal(XYspec[:,1], YXdelay))]
    #-------- Full-Stokes parameters
    text_Stokes = np.repeat('',spwNum).tolist()
    for spw_index in range(spwNum):
        StokesI_PL = figFL.add_subplot( 2, spwNum, spw_index + 1 )
        StokesP_PL = figFL.add_subplot( 2, spwNum, spwNum + spw_index + 1 )
        IList = IList + [StokesI_PL]
        PList = PList + [StokesP_PL]
        exp_Tau = np.exp(-(Tau0spec[spw_index] + scipy.interpolate.splev(np.median(timeStamp), exTauSP) ) / np.mean(np.sin(ElScan)))
        atmCorrect = 1.0 / exp_Tau
        TsysSPW = (TrxList[spw_index].transpose(2,0,1) + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))[Trx2antMap]    # [ant, pol, ch]
        SEFD = 2.0* kb* (TsysSPW * atmCorrect).transpose(2,1,0) / Ae[:,bpAntMap]   # SEFD[ch,pol,antMap]
        SAantNum = len(SAantMap); SAblNum = len(SAblMap)
        AmpCalVis = np.mean(np.mean(pCalVis[spw_index], axis=3)[:,[0,3]] * np.sqrt(SEFD[chRange][:,:,ant0[0:SAblNum]]* SEFD[chRange][:,:,ant1[0:SAblNum]]), axis=0)
        indivRelGain = abs(gainComplexVec(AmpCalVis.T)); indivRelGain /= np.percentile(indivRelGain, 75, axis=0)
        SEFD /= (indivRelGain**2).T
        AmpCalVis = (pCalVis[spw_index].transpose(3,0,1,2)* np.sqrt(SEFD[chRange][:,polYindex][:,:,SAant0]* SEFD[chRange][:,polXindex][:,:,SAant1])).transpose(3,2,1,0)
        text_Stokes[spw_index] = ' SPW%02d %5.1f GHz ' % (spwList[spw_index], centerFreqList[spw_index])
        Stokes = np.zeros([4,SAblNum], dtype=complex)
        for bl_index in range(SAblNum):
            Minv = InvMullerMatrix(Dcat[SAant1[bl_index], 0, spw_index], Dcat[SAant1[bl_index], 1, spw_index], Dcat[SAant0[bl_index], 0, spw_index], Dcat[SAant0[bl_index], 1, spw_index])
            Stokes[:,bl_index] = PS.reshape(4, 4*PAnum).dot(Minv.dot(np.mean(AmpCalVis[bl_index], axis=1)).reshape(4*PAnum)) / PAnum
        #
        StokesVis, StokesErr = Stokes.real, Stokes.imag
        #-------- Visibility slope vs uvdist using Stokes I
        percent75 = np.percentile(StokesVis[0], 75); sdvis = np.std(StokesVis[0])
        visFlag = np.where(abs(StokesVis[0] - percent75) < 2.0* sdvis )[0].tolist()
        weight = np.zeros(SAblNum); weight[visFlag] = 1.0/np.var(StokesVis[0][visFlag])
        P, W = np.c_[np.ones(SAblNum), uvDist[SAblMap]], np.diag(weight)
        PtWP_inv = scipy.linalg.inv(P.T.dot(W.dot(P)))
        solution, solerr = PtWP_inv.dot(P.T.dot(weight* StokesVis[0])),  np.sqrt(np.diag(PtWP_inv)) # solution[0]:intercept, solution[1]:slope
        slopeSNR = abs(solution[1]) / abs(solerr[1]) # ; print 'Slope SNR = ' + `slopeSNR`
        if slopeSNR < 3.0: solution[0], solution[1] = np.percentile(StokesVis[0][visFlag], 75),  0.0
        ScanFlux[scan_index, spw_index, 0], ScanSlope[scan_index, spw_index, 0], ErrFlux[scan_index, spw_index, 0] = solution[0], solution[1], solerr[0]
        for pol_index in range(1,4):
            ScanSlope[scan_index, spw_index, pol_index] = ScanSlope[scan_index, spw_index, 0] * np.median(StokesVis[pol_index])/ScanFlux[scan_index, spw_index, 0]
            solution[0] = (weight.dot(StokesVis[pol_index]) - ScanSlope[scan_index, spw_index, pol_index]* weight.dot(uvDist[SAblMap]))/(np.sum(weight))
            ScanFlux[scan_index, spw_index, pol_index] = solution[0]
            resid = StokesVis[pol_index] - ScanSlope[scan_index, spw_index, pol_index]* uvDist[SAblMap] - solution[0]; ErrFlux[scan_index, spw_index, pol_index] = np.sqrt(weight.dot(resid**2)/np.sum(weight))
        #
        #-------- Check SNR of Stokes I
        if np.max(ErrFlux[scan_index, spw_index]) > 0.2: DcalFlag = False; scanFlag = False
        if ScanFlux[scan_index, spw_index, 0] < 3.0* ErrFlux[scan_index, spw_index, 0]: DcalFlag = False; scanFlag = False
        if DcalFlag :
            AmpCalChAvg[spw_index][:,:,timePointer:timePointer+timeNum] = np.mean(AmpCalVis, axis=2).transpose(1,0,2)
        else:
            scanDic[sourceName][1] *= 0
        #
        #-------- Plot Stokes visibilities
        for pol_index in range(4): text_Stokes[spw_index] = text_Stokes[spw_index] + ' %7.4f (%.4f) ' % (ScanFlux[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index])
        StokesI_PL.plot( uvDist[SAblMap], StokesVis[0], '.', label=polLabel[0], color=Pcolor[0])
        StokesP_PL.plot( uvDist[SAblMap], StokesVis[1], '.', label=polLabel[1], color=Pcolor[1])
        StokesP_PL.plot( uvDist[SAblMap], StokesVis[2], '.', label=polLabel[2], color=Pcolor[2])
        StokesP_PL.plot( uvDist[SAblMap], StokesVis[3], '.', label=polLabel[3], color=Pcolor[3])
        text_Stokes[spw_index] = text_Stokes[spw_index] + '%6.3f   %6.1f ' % (100.0* np.sqrt(ScanFlux[scan_index, spw_index, 1]**2 + ScanFlux[scan_index, spw_index, 2]**2)/ScanFlux[scan_index, spw_index, 0], np.arctan2(ScanFlux[scan_index, spw_index, 2],ScanFlux[scan_index, spw_index, 1])*90.0/pi)
    #
    uvMin, uvMax, IMax = min(uvDist[SAblMap]), max(uvDist[SAblMap]), max(ScanFlux[scan_index,:,0])
    for spw_index in range(spwNum):
        StokesI_PL, StokesP_PL = IList[spw_index], PList[spw_index]
        if spw_index == 0: StokesI_PL.text(0.0, IMax*1.35, text_src)
        if spw_index == spwNum - 1: StokesI_PL.text(0.0, IMax*1.35, timeLabel)
        StokesI_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 0], ScanFlux[scan_index, spw_index, 0]+ uvMax* ScanSlope[scan_index, spw_index, 0]]), '-', color=Pcolor[0])
        StokesP_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 1], ScanFlux[scan_index, spw_index, 1]+ uvMax* ScanSlope[scan_index, spw_index, 1]]), '-', color=Pcolor[1])
        StokesP_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 2], ScanFlux[scan_index, spw_index, 2]+ uvMax* ScanSlope[scan_index, spw_index, 2]]), '-', color=Pcolor[2])
        StokesP_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 3], ScanFlux[scan_index, spw_index, 3]+ uvMax* ScanSlope[scan_index, spw_index, 3]]), '-', color=Pcolor[3])
        StokesI_PL.axis([0.0, uvMax, 0.0, 1.25*IMax]); StokesP_PL.axis([0.0, uvMax, -0.25*IMax, 0.25*IMax])
        StokesI_PL.text(0.0, 1.26*IMax, 'SPW%2d %5.1f GHz' % (spwList[spw_index], centerFreqList[spw_index]))
    #
    StokesI_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    StokesP_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    plt.show()
    figFL.savefig(pp, format='pdf')
    freqArray = np.array(centerFreqList)[range(spwNum)]; meanFreq = np.mean(freqArray); relFreq = freqArray - meanFreq
    if spwNum > 1:
        text_mean = ' mean  %5.1f GHz ' % (meanFreq)
        pflux, pfluxerr = np.zeros(4), np.zeros(4)
        spwStokesDic[sourceName] = []
        for pol_index in range(4):
            sol, solerr = linearRegression(relFreq, ScanFlux[scan_index, :, pol_index], ErrFlux[scan_index, :, pol_index] ); pflux[pol_index], pfluxerr[pol_index] = sol[0], solerr[0]
            text_mean = text_mean + ' %7.4f (%.4f) ' % (pflux[pol_index], pfluxerr[pol_index])
            if SSO_flag and (pol_index > 0):  sol = np.zeros(2)
            spwStokesDic[sourceName] = spwStokesDic[sourceName] + (sol[0] + sol[1]* relFreq).tolist()
        #
        text_mean = text_mean + '%6.3f   %6.1f' % (100.0* np.sqrt(pflux[1]**2 + pflux[2]**2)/pflux[0], np.arctan2(pflux[2],pflux[1])*90.0/pi)
        waveLength = 299.792458/meanFreq    # wavelength in mm
        text_ingest = '%10s, NE, NE, NE, NE, %.2fE+09, %6.3f, %6.4f, %6.3f, %6.4f, %5.1f, %5.1f, NE, NE, %s, %s, %s\n' % (sourceName, meanFreq, pflux[0], np.sqrt(0.0004*pflux[0]**2 + pfluxerr[0]**2), np.sqrt(pflux[1]**2 + pflux[2]**2)/pflux[0], np.sqrt(pfluxerr[1]**2 + pfluxerr[2]**2)/pflux[0], np.arctan2(pflux[2],pflux[1])*90.0/pi, np.sqrt(pfluxerr[1]**2 + pfluxerr[2]**2)/np.sqrt(pflux[1]**2 + pflux[2]**2)*90.0/pi, timeLabel.replace('/','-'), fluxCalText, prefix)
    #
    if scanFlag:
        logfile.write(text_src + ' ' + timeLabel + '\n'); print text_src + ' ' + timeLabel
        text_sd = ' SPW  Frequency     I                 Q                 U                 V                %Pol     EVPA '; logfile.write(text_sd + '\n'); print text_sd
        text_sd = ' --------------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
        for spw_index in range(spwNum): logfile.write(text_Stokes[spw_index] + '\n'); print text_Stokes[spw_index]
        text_sd = ' --------------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
        if spwNum > 1: logfile.write(text_mean + '\n'); print text_mean; ingestFile.write(text_ingest)
        text_sd = 'UV_min_max  %6.1f  %6.1f ' % (uvMin, uvMax); logfile.write(text_sd + '\n'); print text_sd
    else:
        '\033[91m %+.3f%+.3fi \033[0m'
        print '\033[91m' + text_src + ' ' + timeLabel + '\033[0m'
        text_sd = ' SPW  Frequency     I                 Q                 U                 V                %Pol     EVPA '; print '\033[91m' + text_sd + '\033[0m'
        text_sd = ' --------------------------------------------------------------------------------------------------------'; print '\033[91m' + text_sd + '\033[0m'
        for spw_index in range(spwNum): print '\033[91m' + text_Stokes[spw_index] + '\033[0m'
        text_sd = ' --------------------------------------------------------------------------------------------------------'; print '\033[91m' + text_sd + '\033[0m'
        if spwNum > 1: print '\033[91m' + text_mean + '\033[0m'
        text_sd = 'UV_min_max  %6.1f  %6.1f ' % (uvMin, uvMax); print '\033[91m' + text_sd + '\033[0m'
    #
    if COMPDB: 
        print ' -------- Comparison with ALMA Calibrator Catalog --------'
        au.searchFlux(sourcename='%s' % (sourceName), band=int(UniqBands[band_index][3:5]), date=timeText[0:10], maxrows=3)
    #
    if DcalFlag :
        timePointer += timeNum
    else:
        timeSum -= timeNum
    #
    logfile.write('\n'); print ''
#
AmpCalChAvg = AmpCalChAvg[:,:,:,range(timeSum)]
for PL in IList: figFL.delaxes(PL)
for PL in PList: figFL.delaxes(PL)
ingestFile.close()
plt.close('all')
pp.close()
np.save(prefix + '-' + UniqBands[band_index] + '.Flux.npy', ScanFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Ferr.npy', ErrFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Source.npy', np.array(sourceList)[sourceIDscan])
np.save(prefix + '-' + UniqBands[band_index] + '.AZEL.npy', np.array([scanTime, OnAZ, OnEL, OnPA]))
np.save(prefix + '-' + UniqBands[band_index] + '.XYC.npy', np.array(XYC).reshape([len(XYC)/spwNum/2, spwNum, 2]))
np.save(prefix + '-' + UniqBands[band_index] + '.XYD.npy', np.array(XYD).reshape([len(XYC)/spwNum/2, spwNum, 2]))
#---- D-term solution
StokesI, QCpUS, UCmQS = np.zeros(timeSum), np.zeros(timeSum), np.zeros(timeSum)
Dterm = np.zeros([UseAntNum, spwNum, 2], dtype=complex)
DtermDic['Freq'] = centerFreqList
StokesDic['Freq'] = np.mean(np.array(centerFreqList))
for sourceName in sourceList:
    if sourceName == '': continue
    if len(scanDic[sourceName]) == 0: del StokesDic[sourceName]
for spw_index in range(spwNum):
    for sourceName in sourceList:
        if len(scanDic[sourceName]) == 0: continue
        if scanDic[sourceName][1] != 1: continue
        chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = Freq * 1.0e-9
        timeIndex = timeDic[sourceName]
        PA = np.array(PADic[sourceName]); timeNum = len(PA)
        CS, SN = np.cos(2.0* PA), np.sin(2.0* PA)
        Isol, Qsol, Usol = spwStokesDic[sourceName][spw_index], spwStokesDic[sourceName][spwNum + spw_index], spwStokesDic[sourceName][2*spwNum + spw_index]
        QCpUS[timeIndex] = Qsol* CS + Usol* SN
        UCmQS[timeIndex] = Usol* CS - Qsol* SN
        StokesI[timeIndex] = Isol
    #
    Dterm[:,spw_index, 0], Dterm[:,spw_index,1] = VisMuiti_solveD(AmpCalChAvg[spw_index], QCpUS, UCmQS, [], [], StokesI)
#
for ant_index in range(UseAntNum):
    DtermDic[antList[antMap[ant_index]]] = Dterm[ant_index]
#
fileDic = open(prefix + '-B' + `int(UniqBands[band_index][3:5])` + '.D.dic', mode='wb'); pickle.dump(DtermDic, fileDic); fileDic.close()
fileDic = open(prefix + '-B' + `int(UniqBands[band_index][3:5])` + '.S.dic', mode='wb'); pickle.dump(StokesDic, fileDic); fileDic.close()
text_sd = 'D-term        '; logfile.write(text_sd); print text_sd,
for spw_index in range(spwNum):
    text_sd = '    SPW%02d                    ' % (spwList[spw_index]);  logfile.write(text_sd); print text_sd,
logfile.write('\n'); print ''
text_sd = 'Ant  '; logfile.write(text_sd); print text_sd,
for spw_index in range(spwNum):
    text_sd = '       Dx             Dy     ';  logfile.write(text_sd); print text_sd,
logfile.write('\n'); print ''
for ant_index in range(UseAntNum):
    text_sd = '%s  ' % (antList[antMap[ant_index]])
    text_fd = text_sd
    for spw_index in range(spwNum):
        for pol_index in range(2):
            if abs(Dterm[ant_index,spw_index,pol_index]) > 0.1:
                text_sd = text_sd + '\033[91m %+.3f%+.3fi \033[0m' % (Dterm[ant_index,spw_index,pol_index].real, Dterm[ant_index,spw_index,pol_index].imag)
            else:
                text_sd = text_sd + ' %+.3f%+.3fi ' % (Dterm[ant_index,spw_index,pol_index].real, Dterm[ant_index,spw_index,pol_index].imag)
            #
            text_fd = text_fd + ' %+.3f%+.3fi ' % (Dterm[ant_index,spw_index,pol_index].real, Dterm[ant_index,spw_index,pol_index].imag)
        #
    #
    logfile.write(text_fd + '\n'); print text_sd
#----
logfile.close()
msmd.close()
msmd.done()
del flagAnt, TrxFlag, gainFlag, refantID, AntID, Xspec, tempSpec, BPCaledXspec, BP_ant, Gain, GainP, Minv, SEFD, TrxList, TsysSPW, azelTime, azelTime_index, chAvgVis, W, refIndex
if 'refant' in locals(): del refant
