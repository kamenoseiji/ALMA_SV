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
TrxAnts  = np.load(prefix +  '-' + UniqBands[band_index] + '.TrxAnt.npy') # TrxAnts[ant]
Tau0spec = np.load(prefix +  '-' + UniqBands[band_index] + '.Tau0.npy') # Tau0spec[spw][ch]
Trxspec  = np.load(prefix +  '-' + UniqBands[band_index] + '.Trx.npy')  # Trxspec[spw, ant, pol, ch]
Tau0E    = np.load(prefix +  '-' + UniqBands[band_index] + '.TauE.npy') # Tau0E[spw, atmScan]
atmTimeRef = np.load(prefix +  '-' + UniqBands[band_index] + '.atmTime.npy') # atmTimeRef[atmScan]
TrxMap = indexList(TrxAnts, antList); TrxFlag = np.zeros([antNum]); TrxFlag[TrxMap] = 1.0
Tau0E = np.nanmedian(Tau0E, axis=0); Tau0E[np.isnan(Tau0E)] = np.nanmedian(Tau0E); Tau0E[np.isnan(Tau0E)] = 0.0
TrxMed = np.median(Trxspec, axis=3)
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
BPDone = False
print '---Generating antenna-based bandpass table'
BPList = []
for spw_index in spwList:
    BP_ant, XY_BP, XYdelay, Gain = BPtable(msfile, spw_index, BPScan, blMap, blInv)     # BP_ant[antMap, pol, ch]
    BP_ant[:,1] *= XY_BP
    exp_Tau = np.exp(-Tau0spec[spw_index] / np.sin(BPEL))
    atmCorrect = 1.0 / exp_Tau
    TsysBPScan = atmCorrect* (Trxspec[spw_index][Trx2antMap] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau)) # [antMap, pol, ch]
    TsysBPShape = (TsysBPScan.transpose(2,0,1) / np.median(TsysBPScan, axis=2)).transpose(1,2,0)
    BPList = BPList + [BP_ant* np.sqrt(TsysBPShape)]
    #BPList = BPList + [BP_ant]
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
Trx2antMap = indexList( antList[antMap], antList[TrxMap] )
#
if len(atmTimeRef) > 5:
    exTauSP  = UnivariateSpline(atmTimeRef, Tau0E, np.ones(len(atmTimeRef)), s=0.1*np.std(Tau0E), ext=3)
else:
    tempTime = np.arange(np.min(atmTimeRef) - 3600.0,  np.max(atmTimeRef) + 3600.0, 300.0)
    tempTauE = np.repeat(np.median(Tau0E[spw_index]), len(tempTime))
    exTauSP = UnivariateSpline(tempTime, tempTauE, np.ones(len(tempTime)), s=0.1, ext=3)
#
for spw_index in range(spwNum):
    exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
    TsysEQScan = np.mean(Trxspec[spw_index,:,:,chRange].transpose(1,2,0) + Tcmb*exp_Tau[chRange] + tempAtm* (1.0 - exp_Tau[chRange]), axis=2)[Trx2antMap] # [antMap, pol]
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
        TsysSPW = np.mean(Trxspec[spw_index,:,:,chRange].transpose(1,2,0) + Tcmb*exp_Tau[chRange] + tempAtm* (1.0 - exp_Tau[chRange]), axis=2)[Trx2antMap] # [antMap, pol]
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
if len(SSOUseList) == 0: sys.exit('No available Solar System Objects!! Try a-priori calibration.')
EQflux = np.ones([2*spwNum])
#-------- Flux density of the equalizer
for spw_index in range(spwNum):
    FLX, FLY = [], []
    for sso_index in SSOUseList:
        index = np.where(AeX[:, spw_index, sso_index] > 1.0)[0].tolist()
        FLX = FLX + (AeSeqX[spw_index, index] / AeX[index, spw_index, sso_index]).tolist()
        FLY = FLY + (AeSeqY[spw_index, index] / AeY[index, spw_index, sso_index]).tolist()
    #
    EQflux[spw_index], EQflux[spw_index + spwNum] = np.median(np.array(FLX)), np.median(np.array(FLY))
#
#-------- Power-law fit
P = np.c_[ np.r_[np.log(centerFreqList),np.log(centerFreqList)], np.r_[np.ones(spwNum), np.zeros(spwNum)], np.r_[np.zeros(spwNum), np.ones(spwNum)]]
EQmodel = scipy.linalg.solve(np.dot(P.T, P), np.dot(P.T, np.log(EQflux)))   # alpha, logSx, logSy
EQflux = np.c_[np.exp(EQmodel[0]* np.log(centerFreqList) + EQmodel[1]), np.exp(EQmodel[0]* np.log(centerFreqList) + EQmodel[2])]
##-------- Aperture efficiencies
Ae     = np.c_[AeSeqX.T / EQflux[:,0], AeSeqY.T / EQflux[:,1]].reshape(UseAntNum, ppolNum, spwNum)
AEFF   = (Ae.transpose(1,2,0) / (0.25* pi*antDia[antMap]**2)).transpose(2,0,1)
np.save(prefix + '-' + UniqBands[band_index] + '.AntList.npy', antList[antMap]) 
np.save(prefix + '-' + UniqBands[band_index] + '.Aeff.npy', AEFF) 
text_sd =  ' Aeff:'; logfile.write(text_sd); print text_sd,
for spw_index in range(spwNum):
    for pol_index in range(2):
        text_sd = 'SPW%02d-%s ' % (spwList[spw_index], PolList[pol_index])
        logfile.write(text_sd); print text_sd,
    #
#
logfile.write('\n'); print ''
logjy = open(prefix + '-' + UniqBands[band_index] + '-JyK.log', 'w')
for ant_index in range(UseAntNum):
    text_sd = '%s :' % (antList[antMap[ant_index]]); logfile.write(text_sd); print text_sd,
    for spw_index in range(spwNum):
        for pol_index in range(2):
            #logjy.write( '%s %s %d %' % (prefix, antList[antMap[ant_index]], spwList[spw_index])
            text_sd = '  %4.1f%% ' % (100.0* AEFF[ant_index, pol_index, spw_index]); logfile.write(text_sd); print text_sd,
            text_jy = '%s %s %d %s %f %s %s %fe+9 %e %5.1f' % (prefix, antList[antMap[ant_index]], spwList[spw_index], PolList[pol_index], 2.0* kb / (0.25* AEFF[ant_index, pol_index, spw_index]* pi*antDia[antMap[ant_index]]**2), timeLabel, UniqBands[band_index], centerFreqList[spw_index], len(chRange)* chWid[0], tempAtm); logjy.write(text_jy + '\n')
        #
    #
    logfile.write('\n'); print ''; logjy.write('\n')
##
logfile.write('\n'); print ''
logjy.write('\n'); logjy.close()
#-------- XY phase using BP scan
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], BPScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PA = np.arctan2( np.sin(PA), np.cos(PA))
XYphase, caledVis = [], []
scan_index = onsourceScans.index(BPScan)
Trx2antMap = indexList( antList[antMap], antList[TrxMap] )
for spw_index in range(spwNum):
    exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
    atmCorrect = 1.0 / exp_Tau
    TsysSPW = (Trxspec[spw_index] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))[Trx2antMap] # [antMap, pol, ch]
    TsysBL  = np.sqrt( TsysSPW[ant0][:,polYindex]* TsysSPW[ant1][:,polXindex])
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPScan); timeNum = len(timeStamp)
    if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
    else : flagIndex = range(timeNum)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec * TsysBL/ (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    SEFD = 2.0* kb / Ae[:,:,spw_index].T
    caledVis.append(np.mean((chAvgVis / (GainP[polYindex][:,ant0]* GainP[polXindex][:,ant1].conjugate())).transpose(2, 0, 1)* np.sqrt(SEFD[polYindex][:,ant0]* SEFD[polXindex][:,ant1]), axis=2).T)
#
caledVis = np.array(caledVis)   # [spw, pol, time]
if 'QUMODEL' in locals():
    if QUMODEL: QUsolution = np.array([catalogStokesQ.get(BPcal), catalogStokesU.get(BPcal)])
else: QUsolution = XXYY2QU(PA, np.mean(caledVis[:,[0,3]], axis=0))
#-------- XY phase cal in Bandpass table
for spw_index in range(spwNum):
    XYphase = XY2Phase(PA, QUsolution[0], QUsolution[1], caledVis[spw_index][[1,2]])
    XYsign = np.sign(np.cos(XYphase))
    BPList[spw_index][:,1] *= XYsign
    print 'SPW[%d] : XY phase = %6.1f [deg] sign = %3.0f' % (spwList[spw_index], 180.0*XYphase/pi, XYsign)
#
#-------- Flux Density
ScanFlux = np.zeros([scanNum, spwNum, 4])
ScanSlope= np.zeros([scanNum, spwNum, 4])
ErrFlux  = np.zeros([scanNum, spwNum, 4])
print '---Flux densities of sources ---'
pp = PdfPages('FL_' + prefix + '_' + UniqBands[band_index] + '.pdf')
polLabel = ['I', 'Q', 'U', 'V']
Pcolor   = ['black', 'blue', 'red', 'green']
XYD, XYC, scanTime = [], [],[]      # XY delay and correlation
figFL = plt.figure(figsize = (11, 8))
figFL.suptitle(prefix + ' ' + UniqBands[band_index])
figFL.text(0.45, 0.05, 'Projected baseline [m]')
figFL.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90)
for scan_index in range(scanNum):
    if scan_index > 0:
        for PL in IList: figFL.delaxes(PL)
        for PL in PList: figFL.delaxes(PL)
    #
    #-------- UV distance
    timeStamp, UVW = GetUVW(msfile, spwList[spw_index], onsourceScans[scan_index])
    scanTime = scanTime + [np.median(timeStamp)]
    uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    if min(ElScan) < 20.0 / 180.0* pi: continue
    PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
    #-------- Prepare plots
    IList, PList = [], []      # XY delay and correlation
    text_time = qa.time('%fs' % np.median(timeStamp), form='ymd')[0]
    text_src  = ' %02d %010s EL=%4.1f deg' % (scanList[scan_index], sourceList[sourceIDscan[scan_index]], 180.0* OnEL[scan_index]/pi); logfile.write(text_src + ' ' + text_time + '\n'); print text_src + ' ' + text_time
    if(onsourceScans[scan_index] in SSOscanID):
        SSO_flag = True
        SSO_ID = SSOscanID.index(onsourceScans[scan_index])
        text_sd = ' SPW  Frequency    I               Q               U               V             | Model I'; logfile.write(text_sd + '\n'); print text_sd
    else:
        SSO_flag = False
        text_sd = ' SPW  Frequency    I                 Q                 U                 V                 %Pol     EVPA '; logfile.write(text_sd + '\n'); print text_sd
    text_sd = ' --------------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
    BPCaledXspec = []
    #-------- Sub-array formation
    SAantMap, SAblMap, SAblInv, SAant0, SAant1 = antMap, blMap, blInv, ant0, ant1
    if SSO_flag:
        SAantMap, SAblMap, SAblInv = subArrayIndex(uvFlag[SSO_ID], refantID) # antList[SAantMap] lists usable antennas
        SAblIndex = indexList(np.array(SAblMap), np.array(blMap))
        SAant0, SAant1 = np.array(ant0)[SAblIndex].tolist(), np.array(ant1)[SAblIndex].tolist()
    #
    bpAntMap = indexList(antList[SAantMap],antList[antMap])
    Trx2antMap = indexList( antList[SAantMap], antList[TrxMap] )
    SAantNum = len(SAantMap); SAblNum = len(SAblMap)
    if SAblNum < 6:
        text_sd = ' Only %d baselines for short enough sub-array. Skip!' % (SAblNum) ; logfile.write(text_sd + '\n'); print text_sd
        continue
    #
    #-------- Baseline-based cross power spectra
    for spw_index in range(spwNum):
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], onsourceScans[scan_index])
        timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); UseChNum = len(chRange)
        tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        #-------- Bandpass Calibration
        BPCaledXspec = BPCaledXspec + [(tempSpec / (BPList[spw_index][SAant0][:,polYindex]* BPList[spw_index][SAant1][:,polXindex].conjugate())).transpose(2,3,1,0)] # Bandpass Cal
    #
    #-------- Antenna-based Gain
    BPCaledXspec = np.array(BPCaledXspec)   # BPCaledXspec[spw, pol, ch, bl, time]
    chAvgVis = np.mean(BPCaledXspec[:, :, chRange], axis=(0,2))
    if(SSO_flag): chAvgVis =(np.mean(BPCaledXspec[:,:, chRange], axis=(0,2)).transpose(0,2,1) / SSOmodelVis[SSO_ID, spw_index, SAblMap]).transpose(0,2,1)
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    pCalVis = (BPCaledXspec.transpose(0,2,1,3,4) / (GainP[polYindex][:,ant0[0:SAblNum]]* GainP[polXindex][:,ant1[0:SAblNum]].conjugate()))[:,chRange]
    #-------- XY phase spectra
    for spw_index in range(spwNum):
        delayFact = (chNum + 0.0)/len(chRange)
        XYspec = np.mean(pCalVis[spw_index, :, 1:3, :], axis=(2,3)).T
        XYdelay, XYamp = delay_search(XYspec[:,0]); YXdelay, YXamp = delay_search(XYspec[:,1])
        XYD = XYD + [XYdelay* delayFact, YXdelay* delayFact]
        xyc, yxc = np.mean(delay_cal(XYspec[:,0], XYdelay)), np.mean(delay_cal(XYspec[:,1], YXdelay))
        XYC = XYC + [xyc, yxc]
        #text_sd = ' spw%d : Delay (samples) XY= %.3f YX=%.3f  phase (deg) XY= %.2f YX=%.2f' % (spwList[spw_index], XYdelay* delayFact, YXdelay* delayFact, 90.0* np.angle(np.exp((0.0 + 2.0j)* np.angle(xyc))) / pi, 90.0* np.angle(np.exp((0.0 + 2.0j)* np.angle(yxc))) / pi)
        #print text_sd
    #
    for spw_index in range(spwNum):
        StokesI_PL = figFL.add_subplot( 2, spwNum, spw_index + 1 )
        StokesP_PL = figFL.add_subplot( 2, spwNum, spwNum + spw_index + 1 )
        IList = IList + [StokesI_PL]
        PList = PList + [StokesP_PL]
        exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
        atmCorrect = 1.0 / exp_Tau
        TsysSPW = (Trxspec[spw_index] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))[Trx2antMap]    # [ant, pol, ch]
        if SSO_flag:
            Ta = SSOflux[sso_index, spw_index]* Ae[bpAntMap, :, spw_index]* np.mean(atmCorrect) / (2.0* kb)
            TsysSPW = (TsysSPW.transpose(2,0,1) + Ta).transpose(1,2,0)
        #
        SEFD = 2.0* kb* (TsysSPW * atmCorrect).transpose(2,1,0) / Ae[bpAntMap][:,:,spw_index].T   # SEFD[ch,pol,antMap]
        SAantNum = len(SAantMap); SAblNum = len(SAblMap)
        #-------- Additional equalizaiton
        if not SSO_flag:        # Additional equalization for point sources
            AmpCalVis = np.mean(np.mean(pCalVis[spw_index], axis=3)[:,[0,3]] * np.sqrt(SEFD[chRange][:,:,ant0[0:SAblNum]]* SEFD[chRange][:,:,ant1[0:SAblNum]]), axis=0)
            indivRelGain = abs(gainComplexVec(AmpCalVis.T)); indivRelGain /= np.percentile(indivRelGain, 75, axis=0)
            SEFD /= (indivRelGain**2).T
        #
        AmpCalVis = (pCalVis[spw_index].transpose(3,0,1,2)* np.sqrt(SEFD[chRange][:,polYindex][:,:,ant0[0:SAblNum]]* SEFD[chRange][:,polXindex][:,:,ant1[0:SAblNum]])).transpose(3,2,1,0)
        text_sd = ' SPW%02d %5.1f GHz ' % (spwList[spw_index], centerFreqList[spw_index]); logfile.write(text_sd); print text_sd,
        Stokes = np.zeros([4,SAblNum], dtype=complex)
        for bl_index in range(SAblNum):
            Minv = InvMullerVector(DxSpec[SAant1[bl_index], spw_index][chRange], DySpec[SAant1[bl_index], spw_index][chRange], DxSpec[SAant0[bl_index], spw_index][chRange], DySpec[SAant0[bl_index], spw_index][chRange], np.ones(UseChNum))
            Stokes[:,bl_index] = PS.reshape(4, 4*PAnum).dot(Minv.reshape(4, 4*UseChNum).dot(AmpCalVis[bl_index].reshape(4*UseChNum, PAnum)).reshape(4*PAnum)) / (PAnum* UseChNum)
        #
        StokesVis, StokesErr = Stokes.real, Stokes.imag
        if SSO_flag: StokesVis /= SSOmodelVis[SSO_ID, spw_index][SAblMap]
        percent75 = np.percentile(StokesVis[0], 75); sdvis = np.std(StokesVis[0])
        visFlag = np.where(abs(StokesVis[0] - percent75) < 2.0* sdvis )[0]
        if len(visFlag) < 2 : continue
        weight = np.zeros(SAblNum); weight[visFlag] = 1.0/np.var(StokesVis[0][visFlag])
        if SSO_flag: weight *= (SSOmodelVis[SSO_ID, spw_index][SAblMap]**2)
        P, W = np.c_[np.ones(len(weight)), uvDist[SAblMap]], np.diag(weight)
        PtWP_inv = scipy.linalg.inv(P.T.dot(W.dot(P)))
        solution, solerr = PtWP_inv.dot(P.T.dot(weight* StokesVis[0])),  np.sqrt(np.diag(PtWP_inv)) # solution[0]:intercept, solution[1]:slope
        if abs(solution[1]) < 5.0* solerr[1]: solution[0], solution[1] = np.median(StokesVis[0][visFlag]), 0.0
        ScanFlux[scan_index, spw_index, 0], ScanSlope[scan_index, spw_index, 0], ErrFlux[scan_index, spw_index, 0] = solution[0], solution[1], solerr[0]
        for pol_index in range(1,4):
            ScanSlope[scan_index, spw_index, pol_index] = ScanSlope[scan_index, spw_index, 0] * np.median(StokesVis[pol_index])/ScanFlux[scan_index, spw_index, 0]
            solution[0] = (weight.dot(StokesVis[pol_index]) - ScanSlope[scan_index, spw_index, pol_index]* weight.dot(uvDist[SAblMap]))/(np.sum(weight))
            ScanFlux[scan_index, spw_index, pol_index] = solution[0]
            resid = StokesVis[pol_index] - ScanSlope[scan_index, spw_index, pol_index]* uvDist[SAblMap] - solution[0]; ErrFlux[scan_index, spw_index, pol_index] = np.sqrt(weight.dot(resid**2)/np.sum(weight))
        #
        for pol_index in range(4):
            if len(visFlag) < 6:
                text_sd = 'Only %d vis.    ' % (len(visFlag)) ; logfile.write(text_sd + '\n'); print text_sd,
                continue
            #
            text_sd = '%7.4f (%.4f) ' % (ScanFlux[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index]); logfile.write(text_sd); print text_sd,
        #
        StokesI_PL.plot( uvDist[SAblMap], StokesVis[0], '.', label=polLabel[0], color=Pcolor[0])
        StokesP_PL.plot( uvDist[SAblMap], StokesVis[1], '.', label=polLabel[1], color=Pcolor[1])
        StokesP_PL.plot( uvDist[SAblMap], StokesVis[2], '.', label=polLabel[2], color=Pcolor[2])
        StokesP_PL.plot( uvDist[SAblMap], StokesVis[3], '.', label=polLabel[3], color=Pcolor[3])
        if(SSO_flag):
            text_sd = '| %6.3f ' % (SSOflux0[SSO_ID, spw_index]); logfile.write(text_sd); print text_sd,
            logfile.write('\n'); print ''
        else: 
            text_sd = '%6.3f   %6.1f ' % (100.0* np.sqrt(ScanFlux[scan_index, spw_index, 1]**2 + ScanFlux[scan_index, spw_index, 2]**2)/ScanFlux[scan_index, spw_index, 0], np.arctan2(ScanFlux[scan_index, spw_index, 2],ScanFlux[scan_index, spw_index, 1])*90.0/pi); logfile.write(text_sd); print text_sd,
            logfile.write('\n'); print ''
        #
        #
    #
    uvMin, uvMax, IMax = min(uvDist[blMap]), max(uvDist[blMap]), max(ScanFlux[scan_index,:,0])
    for spw_index in range(spwNum):
        StokesI_PL, StokesP_PL = IList[spw_index], PList[spw_index]
        if spw_index == 0: StokesI_PL.text(0.0, IMax*1.35, text_src)
        if spw_index == spwNum - 1: StokesI_PL.text(0.0, IMax*1.35, text_time)
        #for spw_index in range(spwNum):
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
    #
    freqArray = np.array(centerFreqList)[range(spwNum)]; meanFreq = np.mean(freqArray); relFreq = freqArray - meanFreq
    text_sd = ' --------------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
    pflux, pfluxerr = np.zeros(4), np.zeros(4)
    text_sd = ' mean  %5.1f GHz' % (meanFreq); logfile.write(text_sd); print text_sd,
    for pol_index in range(4):
        sol, solerr = linearRegression(relFreq, ScanFlux[scan_index, :, pol_index], ErrFlux[scan_index, :, pol_index] ); pflux[pol_index], pfluxerr[pol_index] = sol[0], solerr[0]
        text_sd = '%7.4f (%.4f) ' % (pflux[pol_index], pfluxerr[pol_index]) ; logfile.write(text_sd); print text_sd,
    #
    text_sd = '%6.3f   %6.1f ' % (100.0* np.sqrt(pflux[1]**2 + pflux[2]**2)/pflux[0], np.arctan2(pflux[2],pflux[1])*90.0/pi); logfile.write(text_sd); print text_sd,
    logfile.write('\n')
    text_sd = 'UV_min_max  %6.1f  %6.1f ' % (uvMin, uvMax); logfile.write(text_sd); print text_sd,
    logfile.write('\n \n'); print '\n'
    if not SSO_flag:
        waveLength = 299.792458/meanFreq    # wavelength in mm
        text_sd = '%s, NE, NE, NE, NE, %.2fE+09, %.3f, %.3f, %.3f, %.3f, %.2f, %.2f, %.2f, %.2f, %s\n' % (sourceList[sourceIDscan[scan_index]], meanFreq, pflux[0], pfluxerr[0], np.sqrt(pflux[1]**2 + pflux[2]**2)/pflux[0], np.sqrt(pfluxerr[1]**2 + pfluxerr[2]**2)/pflux[0], np.arctan2(pflux[2],pflux[1])*90.0/pi, np.sqrt(pfluxerr[1]**2 + pfluxerr[2]**2)/np.sqrt(pflux[1]**2 + pflux[2]**2)*90.0/pi, uvMin/waveLength, uvMax/waveLength, timeLabel[0:10].replace('/','-'))
        ingestFile.write(text_sd)
    #
    if COMPDB & (not SSO_flag):
        print ' -------- Comparison with ALMA Calibrator Catalog --------'
        au.searchFlux(sourcename='%s' % (sourceList[sourceIDscan[scan_index]]), band=int(UniqBands[band_index][3:5]), date=timeLabel[0:10], maxrows=3)
        print '\n'
    #
#
for PL in IList: figFL.delaxes(PL)
for PL in PList: figFL.delaxes(PL)
#
ingestFile.close()
logfile.close()
plt.close('all')
pp.close()
np.save(prefix + '-' + UniqBands[band_index] + '.Flux.npy', ScanFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Ferr.npy', ErrFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Source.npy', np.array(sourceList)[sourceIDscan])
np.save(prefix + '-' + UniqBands[band_index] + '.AZEL.npy', np.array([scanTime, OnAZ, OnEL, OnPA]))
np.save(prefix + '-' + UniqBands[band_index] + '.XYC.npy', np.array(XYC).reshape([len(XYC)/spwNum/2, spwNum, 2]))
np.save(prefix + '-' + UniqBands[band_index] + '.XYD.npy', np.array(XYD).reshape([len(XYC)/spwNum/2, spwNum, 2]))
msmd.close()
msmd.done()
del flagAnt, TrxFlag, gainFlag, Dflag, AntID, BPCaledXspec, BP_ant, Gain, GainP, Minv, SEFD, Trxspec, TsysSPW, TsysBL, azelTime, azelTime_index, chAvgVis, W
