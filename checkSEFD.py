import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
from matplotlib.backends.backend_pdf import PdfPages
#
msmd.open(msfile)
print '---Checking array configulation'
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
Tatm_OFS  = 15.0     # Ambient-load temperature - Atmosphere temperature
kb        = 1.38064852e3
#-------- Check Scans for atmCal
logfile = open(prefix + '-' + UniqBands[band_index] + '-Flux.log', 'w') 
logfile.write(FLScaleText + '\n'); logfile.write(BPcalText + '\n'); logfile.write(EQcalText + '\n')
print '---Checking time series in MS and atmCal scans'
#-------- Tsys measurements
execfile(SCR_DIR + 'TsysCal.py')
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
AeNominal = 0.6* 0.25* np.pi* antDia**2      # Nominal Collecting Area
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
##-------- Equalization using EQ scan
print '---Gain Equalization'
GainP, AeSeqX, AeSeqY = [], [],[]  # effective area x flux density of the equalizer
polXindex, polYindex, scan_index = (arange(4)//2).tolist(), (arange(4)%2).tolist(), onsourceScans.index(EQScan)
for spw_index in range(spwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = ParaPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)       # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0]* BPList[spw_index][ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    #-------- Antenna-based Gain
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainP = GainP + [np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])]
    #GainX, GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[0]), np.apply_along_axis( gainComplex, 0, chAvgVis[1])
    #BLphsX, BLphsY = GainX[ant0]* GainX[ant1].conjugate() / abs(GainX[ant0]* GainX[ant1]), GainY[ant0]* GainY[ant1].conjugate() / abs(GainY[ant0]* GainY[ant1])
    pCalVisX, pCalVisY = np.mean(chAvgVis[0] / (GainP[spw_index][0,ant0]* GainP[spw_index][0,ant1].conjugate()), axis=1), np.mean(chAvgVis[1] / (GainP[spw_index][1,ant0]* GainP[spw_index][1,ant1].conjugate()), axis=1)
    #-------- Antenna-based Gain
    AeSeqX = AeSeqX + [2.0* kb* chAvgTsys[antMap, spw_index, 0, scan_index]* abs(gainComplex(pCalVisX))**2] # Ae x Seq (X-pol)
    AeSeqY = AeSeqY + [2.0* kb* chAvgTsys[antMap, spw_index, 1, scan_index]* abs(gainComplex(pCalVisY))**2] # Ae x Seq (X-pol)
#
GainP = np.array(GainP) # GainP[spw, pol, ant, time]
AeSeqX, AeSeqY = np.array(AeSeqX), np.array(AeSeqY) # AeSeqX[spw, antMap], AeSeqX[spw, antMap] : (relative) aperture efficiency, assuming 1 Jy
#
##-------- Phase alignment between SPWs
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
execfile(SCR_DIR + 'SSOflux.py')
######## Outputs from SSOflux.py :
#  SSOflux0[SSO, spw] : model flux density of solar system objects at zero spacing
#  SSOmodelVis[SSO, spw, bl] : predicted relative visibility (max=1.0)
#  uvFlag[SSO, spw, bl] : 0=resolved, 1=unresolved; bl is orderd as blMap
########
SSOflux = SSOflux0* np.exp(-onTau.transpose(1,0)[indexList(np.array(SSOscanID), np.array(onsourceScans))])  # SSOflux[SSO, spw] : attenuated SSO flux
##-------- Scaling with the flux calibrator
AeX, AeY = np.zeros([UseAntNum, spwNum, SSONum]), np.zeros([UseAntNum, spwNum, SSONum])
#-------- Sub-array with unflagged antennas (short baselines)
SSO_flag = np.ones(SSONum)
for sso_index in range(SSONum):
    for spw_index in range(spwNum):
        SAantennas, SAbl, SAblFlag, SAant0, SAant1 = subArrayIndex(uvFlag[sso_index, spw_index]) # in Canonical ordering
        if len(SAantennas) < 3: continue #  Too few antennas
        SAantMap, SAblMap, SAblInv = np.array(antMap)[SAantennas].tolist(), np.array(blMap)[SAbl].tolist(), np.array(blInv)[SAbl].tolist()
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], SSOscanID[sso_index])
        chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        tempSpec = ParaPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        BPCaledXspec = (tempSpec / (BPList[spw_index][SAant0]* BPList[spw_index][SAant1].conjugate())).transpose(2,3,1,0) #  BPCaledXspe[pol, ch, bl, time]
        chAvgVis = (np.mean(BPCaledXspec[:, chRange], axis=1).transpose(0,2,1)/SSOmodelVis[sso_index, spw_index,SAbl]).transpose(0,2,1)
        GainX, GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[0]), np.apply_along_axis( gainComplex, 0, chAvgVis[1])
        Ta = SSOflux[sso_index, spw_index]* AeNominal[SAantMap] / (2.0* kb)
        AeX[SAantennas, spw_index, sso_index] = 2.0* kb* np.percentile(abs(GainX), 75, axis=1)**2 * (Ta + chAvgTsys[SAantMap, spw_index, 0, scan_index]) / SSOflux[sso_index, spw_index]
        AeY[SAantennas, spw_index, sso_index] = 2.0* kb* np.percentile(abs(GainY), 75, axis=1)**2 * (Ta + chAvgTsys[SAantMap, spw_index, 1, scan_index]) / SSOflux[sso_index, spw_index]
    #
#
#-------- SSO Flagging
for sso_index in range(SSONum):
    for spw_index in range(spwNum):
        index = np.where(AeX[:, spw_index, sso_index] > 1.0)[0].tolist()
        if len(index) < 3: SSO_flag[sso_index] = 0.0; continue
        FLX_stat, FLY_stat = AeX[index, spw_index, sso_index]/AeNominal[np.array(antMap)[index].tolist()], AeY[index, spw_index, sso_index]/AeNominal[np.array(antMap)[index].tolist()]
        if np.percentile(FLX_stat, 75) / np.median(FLX_stat) > 1.25: SSO_flag[sso_index] = 0.0
        if np.percentile(FLY_stat, 75) / np.median(FLY_stat) > 1.25: SSO_flag[sso_index] = 0.0
    #
#
SSOUseList = np.where(SSO_flag == 1.0)[0].tolist()
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
        text_sd = 'SPW%02d-%s ' % (spw[spw_index], PolList[pol_index])
        logfile.write(text_sd); print text_sd,
    #
#
logfile.write('\n'); print ''
for ant_index in range(UseAntNum):
    text_sd = '%s :' % (antList[antMap[ant_index]]); logfile.write(text_sd); print text_sd,
    for spw_index in range(spwNum):
        for pol_index in range(2): text_sd = '  %4.1f%% ' % (100.0* AEFF[ant_index, pol_index, spw_index]); logfile.write(text_sd); print text_sd,
    #
    logfile.write('\n'); print ''
#
logfile.write('\n'); print ''
#-------- Flux Density
ScanFlux = np.zeros([scanNum, spwNum])
ScanSlope= np.zeros([scanNum, spwNum])
ErrFlux  = np.zeros([scanNum, spwNum])
ScanEL     = np.zeros([scanNum])
print '---Flux densities of sources ---'
text_sd = 'Scan  Source     EL(deg) '
for spw_index in range(spwNum): text_sd = text_sd + ' SPW%02d %5.1f GHz ' %  (spw[spw_index], centerFreqList[spw_index])
logfile.write(text_sd + '\n'); print text_sd
for scan_index in range(scanNum):
    ScanEL[scan_index] = np.median(OnEL[:,scan_index])
    text_sd = ' --------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
    text_sd = ' %03d %010s  %4.1f    ' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0* ScanEL[scan_index]/pi )
    if(onsourceScans[scan_index] in SSOscanID):
        SSO_flag = T
        SSO_ID = SSOscanID.index(onsourceScans[scan_index])
    else:
        SSO_flag = F
    BPCaledXspec = []
    SEFD = np.ones([spwNum, 2, UseAntNum])
    #-------- Sub-array with unflagged antennas (short baselines)
    if SSO_flag:
        SAantennas, SAbl, SAblFlag, SAant0, SAant1 = subArrayIndex(uvFlag[SSO_ID, spwNum - 1])
        SAantMap, SAblMap, SAblInv = np.array(antMap)[SAantennas].tolist(), np.array(blMap)[SAbl].tolist(), np.array(blInv)[SAbl].tolist()
    else:
        SAantennas, SAbl, SAblFlag, SAant0, SAant1 = range(UseAntNum), range(UseBlNum), np.ones([blNum]), ant0, ant1
        SAantMap, SAblMap, SAblInv = antMap, blMap, blInv
        TA = 0.0
    #
    for spw_index in range(spwNum):
        atmCorrect = np.exp(-onTau[spw_index, scan_index])
        if SSO_flag: TA = Ae[SAantennas,:,spw_index]* SSOflux0[SSO_ID, spw_index]* atmCorrect  / (2.0* kb)
        SEFD[spw_index][:,SAantennas]  = (2.0* kb* (chAvgTsys[SAantMap,spw_index, :,scan_index] + TA) / (Ae[SAantennas,:,spw_index]* atmCorrect)).T
        SAantNum = len(SAantennas); SAblNum = len(SAblMap)
        if SAblNum < 3:
            text_sd = ' Only %d baselines for short enough sub-array. Skip!' % (SAblNum) ; logfile.write(text_sd + '\n'); print text_sd
            continue
        #
        #-------- UV distance
        timeStamp, UVW = GetUVW(msfile, spw[spw_index], onsourceScans[scan_index])
        uvw = np.mean(UVW[:,SAblMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], onsourceScans[scan_index])
        timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
        tempSpec = ParaPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # tempSpec[timeNum, blNum, polNum, chNum]
        #-------- Bandpass Calibration
        BPCaledXspec = BPCaledXspec + [(tempSpec / (BPList[spw_index][SAant0]* BPList[spw_index][SAant1].conjugate())).transpose(2,3,1,0)] # Bandpass Cal
        #BPCaledXspec = (tempSpec / (BPList[spw_index][SAant0]* BPList[spw_index][SAant1].conjugate())).transpose(2,3,1,0) # BPCaledXspec[polNum, chNum, blNum, timeNum]
    #
    if SAblNum < 3: continue
    BPCaledXspec = np.array(BPCaledXspec)   # BPCaledXspec[spw, pol, ch, bl, time]
    if(SSO_flag): chAvgVis =(np.mean(BPCaledXspec[:,:, chRange], axis=(0,2)).transpose(0,2,1) / SSOmodelVis[SSO_ID, spw_index][SAbl]).transpose(0,2,1)
    else: chAvgVis = np.mean(BPCaledXspec[:, :, chRange], axis=(0,2))
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
    pCalVis = (BPCaledXspec.transpose(0,2,1,3,4) / (GainP[:,ant0[0:SAblNum]]* GainP[:,ant1[0:SAblNum]].conjugate()))[:,chRange]
    for spw_index in range(spwNum):
        chAvgVis = np.mean(pCalVis[spw_index], axis=(0,3))* np.sqrt(SEFD[spw_index][:, ant0[0:SAblNum]]* SEFD[spw_index][:, ant1[0:SAblNum]])
        indivRelGain = abs(gainComplexVec(chAvgVis.T))
        SEFD[spw_index][:,SAantennas] *= ((np.percentile(indivRelGain, 75, axis=0) / indivRelGain)**2).T
    #
    AmpCalVis = (pCalVis.transpose(1,4,0,2,3)* np.sqrt(SEFD[:,:,ant0[0:SAblNum]]* SEFD[:,:,ant1[0:SAblNum]])).transpose(2,4,3,0,1) # AmpCalVis[spw,bl,pol,ch,time]
    StokesI = np.mean(AmpCalVis, axis=(2, 3, 4)).real
    for spw_index in range(spwNum):
        visFlag = np.where(abs(StokesI[spw_index] - np.percentile(StokesI[spw_index], 75))/np.percentile(StokesI[spw_index], 75) < 0.2 )[0]
        if len(visFlag) < 3:
            text_sd = 'Only %d vis.    ' % (len(visFlag)) ; logfile.write(text_sd + '\n'); print text_sd,
            continue
        #
        weight = np.zeros(SAblNum); weight[visFlag] = 1.0/np.var(StokesI[spw_index][visFlag])
        P, W = np.c_[np.ones(SAblNum), uvDist], np.diag(weight)
        PtWP_inv = scipy.linalg.inv(np.dot(P.T, np.dot(W, P)))
        solution, solerr = PtWP_inv.dot(P.T.dot(weight* StokesI[spw_index])),  np.sqrt(np.diag(PtWP_inv))
        if solution[1] > -2.0* solerr[1]: solution[0] = np.median(StokesI[spw_index]); solution[1] = 0.0
        ScanFlux[scan_index, spw_index], ScanSlope[scan_index, spw_index], ErrFlux[scan_index, spw_index] = solution[0], solution[1], solerr[0]
        text_sd = text_sd + '  %6.3f (%.3f) ' % (ScanFlux[scan_index, spw_index], ErrFlux[scan_index, spw_index])
        #
    #
    logfile.write(text_sd); print text_sd
    if(SSO_flag):
        text_sd = '       SSO model         '
        for spw_index in range(spwNum): text_sd = text_sd + '  %6.3f         ' % (SSOflux0[SSO_ID, spw_index])
        logfile.write(text_sd + '\n'); print text_sd
    #
    if COMPDB:
        if not SSO_flag:
            print ' -------- Comparison with ALMA Calibrator Catalog --------'
            au.searchFlux(sourcename='%s' % (sourceList[sourceIDscan[scan_index]]), band=int(UniqBands[band_index][3:5]), date=timeLabel[0:10], maxrows=3)
            print '\n'
        #
    #
#
logfile.close()
np.save(prefix + '-' + UniqBands[band_index] + '.Flux.npy', ScanFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Ferr.npy', ErrFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Source.npy', np.array(sourceList)[sourceIDscan])
np.save(prefix + '-' + UniqBands[band_index] + '.EL.npy', ScanEL)
msmd.done()
