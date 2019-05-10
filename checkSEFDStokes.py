import sys
from scipy import stats
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
from matplotlib.backends.backend_pdf import PdfPages
RADperHzMeterArcsec = 2.0* pi / 299792458 / (180*3600/pi)
class END(Exception):
    pass
#
msmd.open(msfile)
#-------- Configure Array
print '---Checking array configulation'
antDia = np.ones(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
print '  -- usable antenna checking for BP scan : '
spwList = scnspw
gainFlag = np.ones([antNum])
for spw_index in range(spwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPScan)
    timeNum, chNum, blNum = Xspec.shape[3], Xspec.shape[1], Xspec.shape[2]; chRange, timeRange = range(int(0.05*chNum), int(0.95*chNum)), range(timeNum-4, timeNum-1)
    for polID in pPol:
        blD, blA = np.apply_along_axis(delay_search, 0, np.mean(Xspec[polID][chRange][:,:,timeRange], axis=2))
        blA = blA / (antDia[ANT0[0:blNum]]* antDia[ANT1[0:blNum]])
        errD, errA = np.where(abs(blD - np.median(blD)) > 4.0)[0].tolist(), np.where(abs(blA - np.median(blA)) > 0.5* np.median(blA))[0].tolist()
        errCount = np.zeros(antNum)
        for bl in set(errD) or set(errA): errCount[ list(Bl2Ant(bl)) ] += 1
        gainFlag[np.where(errCount > 2.5 )[0].tolist()] *= 0.0
    #
#
#-------- Check D-term files
Dloaded = False
if not 'DPATH' in locals(): DPATH = SCR_DIR
print '---Checking D-term files in ' + DPATH
DantList, noDlist, Dflag = [], [], np.ones([antNum])
Dpath = DPATH + 'DtermB' + `int(UniqBands[band_index][3:5])` + '/'
for ant_index in range(antNum):
    Dfile = Dpath + 'B' + `int(UniqBands[band_index][3:5])` + '-SPW0-' + antList[ant_index] + '.DSpec.npy'
    if os.path.exists(Dfile): DantList += [ant_index]
    else: noDlist += [ant_index]; Dflag[ant_index] *= 0.0
#
DantNum, noDantNum = len(DantList), len(noDlist)
print 'Antennas with D-term file (%d):' % (DantNum),
for ant_index in DantList: print '%s ' % antList[ant_index],
print ''
if noDantNum > 0:
    print 'Antennas without D-term file (%d) : ' % (noDantNum),
    for ant_index in noDlist: print '%s ' % antList[ant_index],
    sys.exit(' Run DtermTransfer first!!')
#
#-------- Load Tsys table
TrxList = []
for spw in spwList: TrxList = TrxList + [np.median(np.load(prefix +  '-' + UniqBands[band_index] + '-SPW' + `spw` + '.Trx.npy'), axis=3)]  # TrxList[spw][pol, ch, ant]
TrxAnts  = np.load(prefix +  '-' + UniqBands[band_index] + '.TrxAnt.npy') # TrxAnts[ant]
Tau0spec = np.load(prefix +  '-' + UniqBands[band_index] + '.Tau0.npy') # Tau0spec[spw][ch]
Tau0E    = np.load(prefix +  '-' + UniqBands[band_index] + '.TauE.npy') # Tau0E[spw, atmScan]
atmTimeRef = np.load(prefix +  '-' + UniqBands[band_index] + '.atmTime.npy') # atmTimeRef[atmScan]
TrxMap = indexList(TrxAnts, antList); TrxFlag = np.zeros([antNum]); TrxFlag[TrxMap] = 1.0
Tau0E = np.nanmedian(Tau0E, axis=0); Tau0E[np.isnan(Tau0E)] = np.nanmedian(Tau0E); Tau0E[np.isnan(Tau0E)] = 0.0
for spw_index in range(spwNum):
    TrxMed = np.median(TrxList[spw_index], axis=1)  # TrxMed[pol, ant]
    for pol_index in range(2): TrxFlag[np.where(abs(TrxMed[pol_index] - np.median(TrxMed[pol_index])) > 0.8* np.median(TrxMed[pol_index]))[0].tolist()] *= 0.0
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
print; print 'Dtrm',
for ant_index in range(antNum): print '   %.0f' % (Dflag[ant_index]),
print
flagAnt = flagAnt* TrxFlag* gainFlag* Dflag
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
if len(UseAnt) < 5: sys.exit('Too few usable antennas. Reduction failed.')
#-------- Check Scans for atmCal
ingestFile = open(prefix + '-' + UniqBands[band_index] + '-Ingest.log', 'w') 
text_sd = '#source,   RA,eRA,dec,edec,frequency,  flux,eflux,degree,edeg,EVPA,eEVPA,uvmin,uvmax,date\n'; ingestFile.write(text_sd)
logfile = open(prefix + '-' + UniqBands[band_index] + '-Flux.log', 'w') 
logfile.write(BPcalText + '\n'); logfile.write(EQcalText + '\n')
#-------- Tsys measurements
scanList = onsourceScans
msmd.close()
msmd.done()
#-------- Array Configuration
print '---Determining refant'
#flagList = np.where(np.median(chAvgTrx.reshape(UseAntNum, 2* spwNum), axis=1) > 2.0* np.percentile(chAvgTrx, 75))[0].tolist()  # Avoid too-high Trx
#flagList = unique(flagList + np.where(np.min(chAvgTrx.reshape(UseAntNum, 2* spwNum), axis=1) < 1.0 )[0].tolist()).tolist()     # Avoid too-low Trx
#if len(flagList) >0 :
#    for index in flagList: del UseAnt[index]
#
#text_sd = '  Usable antennas: '
#for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
#logfile.write(text_sd + '\n'); print text_sd
#text_sd = '  Flagged by Trx:  '
#for ants in antList[flagList].tolist(): text_sd = text_sd + ants + ' '
#logfile.write(text_sd + '\n'); print text_sd
#for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
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
AeNominal = 0.7* 0.25* np.pi* antDia**2      # Nominal Collecting Area
#-------- Load D-term file
DxList, DyList = [], []
for ant_index in range(antNum):
    Dpath = SCR_DIR + 'DtermB' + `int(UniqBands[band_index][3:5])` + '/'
    for spw_index in range(spwNum):
        Dfile = Dpath + 'B' + `int(UniqBands[band_index][3:5])` + '-SPW' + `spw_index` + '-' + antList[ant_index] + '.DSpec.npy'
        # print 'Loading %s' % (Dfile)
        Dterm = np.load(Dfile)
        DxList = DxList + [Dterm[1] + (0.0 + 1.0j)* Dterm[2]]
        DyList = DyList + [Dterm[3] + (0.0 + 1.0j)* Dterm[4]]
    #
#
chNum = np.array(DxList).shape[1]
DxSpec, DySpec = np.array(DxList).reshape([antNum, spwNum, chNum]), np.array(DyList).reshape([antNum, spwNum, chNum])
Dloaded = True
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
BPDone = False
print '---Generating antenna-based bandpass table'
BPList = []
for spw_index in range(spwNum):
    BP_ant, XY_BP, XYdelay, Gain = BPtable(msfile, spwList[spw_index], BPScan, blMap, blInv)     # BP_ant[antMap, pol, ch]
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
scanList = onsourceScans
GainP, AeSeqX, AeSeqY = [], [], []  # effective area x flux density of the equalizer
polXindex, polYindex, scan_index = (arange(4)//2).tolist(), (arange(4)%2).tolist(), onsourceScans.index(EQScan)
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PA = np.arctan2( np.sin(PA), np.cos(PA))
QUsolution = np.zeros(2)
if catalogStokesQ.get(EQcal) > 0.0 :
    QUsolution = np.array([catalogStokesQ.get(EQcal), catalogStokesU.get(EQcal)])
QCpUS = (QUsolution[0]* np.cos(2.0* PA) + QUsolution[1]* np.sin(2.0* PA)) / catalogStokesI.get(EQcal)
#
if len(atmTimeRef) > 5:
    #exTauSP  = UnivariateSpline(atmTimeRef, Tau0E, np.ones(len(atmTimeRef)), s=0.1*np.std(Tau0E), ext=3)
    exTauSP  = UnivariateSpline(np.append(np.append(atmTimeRef[0]-180.0,atmTimeRef), atmTimeRef[-1]+180.0), np.append(np.append(Tau0E[0], Tau0E), Tau0E[-1]), np.ones(len(atmTimeRef)+2), s=0.1*np.std(Tau0E))
else:
    tempTime = np.arange(np.min(atmTimeRef) - 3600.0,  np.max(atmTimeRef) + 3600.0, 300.0)
    tempTauE = np.repeat(np.median(Tau0E[spw_index]), len(tempTime))
    exTauSP = UnivariateSpline(tempTime, tempTauE, np.ones(len(tempTime)), s=0.1, ext=3)
#
for spw_index in range(spwNum):
    exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
    TsysEQScan = np.mean(TrxList[spw_index].transpose(2,0,1)[:,:,chRange] + Tcmb*exp_Tau[chRange] + tempAtm* (1.0 - exp_Tau[chRange]), axis=2)[Trx2antMap] # [antMap, pol]
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)[flagIndex]       # Cross Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)[pPol]  # chAvgVis[pPol,bl,time]
    GainP = GainP + [np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])]
    pCalVisX = np.mean(chAvgVis[0] / (1.0 + QCpUS) / (GainP[spw_index][0,ant0]* GainP[spw_index][0,ant1].conjugate()), axis=1)
    pCalVisY = np.mean(chAvgVis[1] / (1.0 - QCpUS) / (GainP[spw_index][1,ant0]* GainP[spw_index][1,ant1].conjugate()), axis=1)
    #-------- Antenna-based Gain
    AeSeqX = AeSeqX + [2.0* kb* TsysEQScan[:, 0]* abs(gainComplex(pCalVisX))**2] # Ae x Seq (X-pol)
    AeSeqY = AeSeqY + [2.0* kb* TsysEQScan[:, 1]* abs(gainComplex(pCalVisY))**2] # Ae x Seq (Y-pol)
#
GainP = np.array(GainP) # GainP[spw, pol, ant, time]
AeSeqX, AeSeqY = np.array(AeSeqX), np.array(AeSeqY) # AeSeqX[spw, antMap], AeSeqX[spw, antMap] : (relative) aperture efficiency, assuming 1 Jy
#
##-------- inter-SPW phasing using EQ scan
spwPhase = [0.0]* 2* spwNum
for ant_index in range(1,UseAntNum):
    for pol_index in range(2):
        spwPhase = spwPhase + [0.0]
        for spw_index in range(1,spwNum): spwPhase = spwPhase + [np.angle(GainP[spw_index, pol_index, ant_index].dot(GainP[0, pol_index, ant_index].conjugate()))]
    #
#
spwPhase = np.array(spwPhase).reshape([UseAntNum, 2, spwNum]); spwTwiddle = exp(1.0j *spwPhase)
for spw_index in range(1,spwNum): BPList[spw_index] = (BPList[spw_index].transpose(2,0,1)* spwTwiddle[:,:,spw_index]).transpose(1,2,0)
#-------- Flux models for solar system objects
msmd.done()
execfile(SCR_DIR + 'SSOflux.py'); logfile.write(FLScaleText + '\n')
######## Outputs from SSOflux.py :
#  SSOflux0[SSO, spw] : model flux density of solar system objects at zero spacing
#  SSOmodelVis[SSO, spw, bl] : predicted relative visibility (max=1.0)
#  uvFlag[SSO, spw, bl] : 0=resolved, 1=unresolved
########
flaggedBlList = list(set(range(blNum)) - set(blMap)); uvFlag[:,:,flaggedBlList] = 0.0
atmCorrect = np.exp(-outer( np.mean(Tau0spec, axis=1),  1.0/np.sin( np.array(OnEL)[indexList(np.array(SSOscanID), np.array(onsourceScans))]))).T
SSOflux = SSOflux0* atmCorrect  # SSOflux[SSO, spw] : attenuated SSO flux
uvFlag = np.min(uvFlag, axis=1) # all-SPW uv flag
##-------- Scaling with the flux calibrator
AeX, AeY = np.zeros([UseAntNum, spwNum, SSONum]), np.zeros([UseAntNum, spwNum, SSONum])
#-------- Sub-array with unflagged antennas (short baselines)
SSO_flag = np.ones(SSONum)
for sso_index in range(SSONum):
    scan_index = indexList(np.array(SSOscanID), np.array(onsourceScans))[sso_index]
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], SSOscanID[sso_index]); timeNum = len(timeStamp)
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    if np.max(ElScan) < 35.0/180.0*pi : continue    # Too low Elevation
    SAantMap, SAblMap, SAblInv = subArrayIndex(uvFlag[sso_index], refantID) # in Canonical ordering
    if len(SAantMap) < 4: continue #  Too few antennas
    print 'Subarray : ',; print antList[SAantMap]
    SAblIndex = indexList(np.array(SAblMap), np.array(blMap))
    SAant0, SAant1 = np.array(ant0)[SAblIndex].tolist(), np.array(ant1)[SAblIndex].tolist()
    bpAntMap = indexList(antList[SAantMap],antList[antMap])
    Trx2antMap = indexList( antList[SAantMap], antList[TrxMap] )
    for spw_index in range(spwNum):
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], SSOscanID[sso_index])
        chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        BPCaledXspec = (tempSpec / (BPList[spw_index][SAant0][:,polYindex]* BPList[spw_index][SAant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
        chAvgVis = (np.mean(BPCaledXspec[:, chRange], axis=1)[pPol].transpose(0,2,1)/SSOmodelVis[sso_index,spw_index,SAblMap]).transpose(0,2,1)
        GainX, GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[0]), np.apply_along_axis( gainComplex, 0, chAvgVis[1])
        #-------- Tsys
        Ta = SSOflux[sso_index, spw_index]* AeNominal[SAantMap] / (2.0* kb)
        exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
        TsysSPW = np.mean(TrxList[spw_index].transpose(2,0,1)[:,:,chRange] + Tcmb*exp_Tau[chRange] + tempAtm* (1.0 - exp_Tau[chRange]), axis=2)[Trx2antMap] # [antMap, pol]
        #-------- Aperture efficiency
        AeX[bpAntMap, spw_index, sso_index] = 2.0* kb* np.percentile(abs(GainX), 75, axis=1)**2 * (Ta + TsysSPW[:, 0]) / SSOflux[sso_index, spw_index]
        AeY[bpAntMap, spw_index, sso_index] = 2.0* kb* np.percentile(abs(GainY), 75, axis=1)**2 * (Ta + TsysSPW[:, 1]) / SSOflux[sso_index, spw_index]
    #
#
#-------- SSO Flagging
for sso_index in range(SSONum):
    for spw_index in range(spwNum):
        index = np.where(AeX[:, spw_index, sso_index] > 1.0)[0].tolist()
        if len(index) < 4: SSO_flag[sso_index] = 0.0; continue
        FLX_stat, FLY_stat = AeX[index, spw_index, sso_index]/AeNominal[np.array(antMap)[index].tolist()], AeY[index, spw_index, sso_index]/AeNominal[np.array(antMap)[index].tolist()]
        if np.median(FLX_stat) < 0.7: SSO_flag[sso_index] = 0.0 
        if np.median(FLY_stat) < 0.7: SSO_flag[sso_index] = 0.0 
        if np.median(FLX_stat) > 1.5: SSO_flag[sso_index] = 0.0 
        if np.median(FLY_stat) > 1.5: SSO_flag[sso_index] = 0.0 
        if np.percentile(FLX_stat, 75) / np.median(FLX_stat) > 1.5: SSO_flag[sso_index] = 0.0
        if np.percentile(FLY_stat, 75) / np.median(FLY_stat) > 1.5: SSO_flag[sso_index] = 0.0
    #
    try:
        if SSO_flag[sso_index] == 0.0: raise END
    except END:
        text_sd = '%s is not available as a flux calibrator.' % (sourceList[BandSSOList[sso_index]])
        logfile.write(text_sd + '\n'); print text_sd
    #
#
SSOUseList = np.where(SSO_flag == 1.0)[0].tolist()
if len(SSOUseList) > 0:
    execfile(SCR_DIR + 'SEFDStokes.py')
else:
    print 'No available Solar System Objects!! Try a-priori calibration.'
    msmd.close()
    msmd.done()
    #del flagAnt, TrxFlag, gainFlag, Dflag, AntID, BPCaledXspec, BP_ant, Gain, GainP, Minv, SEFD, TrxList, TsysSPW, TsysBL, azelTime, azelTime_index, chAvgVis, W
    execfile(SCR_DIR + 'aprioriStokes.py')
#
