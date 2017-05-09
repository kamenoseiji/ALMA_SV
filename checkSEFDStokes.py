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
#-------- Definitions
def residTskyTransfer( param, Tamb, secz, Tsky, weight ):
    exp_Tau = np.exp( -param[1]* secz )
    return weight* (Tsky - (param[0] + 2.718* exp_Tau  + Tamb* (1.0 - exp_Tau)))
#
def residTskyTransfer2( param, Tamb, Tau0, secz, Tsky, weight ):
    exp_Tau = np.exp( -Tau0* secz )
    return weight* (Tsky - (param[0] + 2.718* exp_Tau  + Tamb* (1.0 - exp_Tau)))
#
def get_progressbar_str(progress):
    MAX_LEN = 48
    BAR_LEN = int(MAX_LEN * progress)
    return ('[' + '=' * BAR_LEN + ('>' if BAR_LEN < MAX_LEN else '') + ' ' * (MAX_LEN - BAR_LEN) + '] %.1f%%' % (progress * 100.))
#
Tatm_OFS  = 15.0     # Ambient-load temperature - Atmosphere temperature
AeNominal = 0.6* 0.25* np.pi* antDia**2      # Nominal Collecting Area
kb        = 1.38064852e3
#-------- Check Scans for atmCal
logfile = open(prefix + '-' + UniqBands[band_index] + '-Flux.log', 'w') 
logfile.write(FLScaleText + '\n')
logfile.write(BPcalText + '\n')
logfile.write(EQcalText + '\n')
print '---Checking time series in MS and atmCal scans'
tb.open(msfile); timeXY = tb.query('ANTENNA1 == 0 && ANTENNA2 == 0 && DATA_DESC_ID == '+`spw[0]`).getcol('TIME'); tb.close()
OnTimeIndex = []
for scan_index in range(scanNum):
    OnTimeIndex.append( indexList(msmd.timesforscan(onsourceScans[scan_index]), timeXY) )
offTime = sort( list(set(timeXY) & set(timeOFF)) )
ambTime = sort( list(set(timeXY) & set(timeAMB)) )
hotTime = sort( list(set(timeXY) & set(timeHOT)) )
#
offTimeIndex = indexList( offTime, timeXY)
ambTimeIndex = indexList( ambTime, timeXY)
hotTimeIndex = indexList( hotTime, timeXY)
#-------- Load D-term file
DxList, DyList = [], []
for ant_index in range(UseAntNum):
    Dpath = SCR_DIR + 'DtermB' + `int(UniqBands[band_index][3:5])` + '/'
    for spw_index in range(4):
        Dfile = Dpath + 'B' + `int(UniqBands[band_index][3:5])` + '-SPW' + `spw_index` + '-' + antList[antMap[ant_index]] + '.DSpec.npy'
        Dterm = np.load(Dfile)
        DxList = DxList + [Dterm[1] + (0.0 + 1.0j)* Dterm[2]]
        DyList = DyList + [Dterm[3] + (0.0 + 1.0j)* Dterm[4]]
    #
#
chNum = np.array(DxList).shape[1]
DxSpec, DySpec = np.array(DxList).reshape([UseAntNum, spwNum, chNum]), np.array(DyList).reshape([UseAntNum, spwNum, chNum])
#-------- Bandpass Table
print '---Generating antenna-based bandpass table'
XYdelayList, BPList = [], []
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
########
##-------- Equalization using EQ scan
GainP, AeSeqX, AeSeqY = [], [], []  # effective area x flux density of the equalizer
polXindex, polYindex, scan_index = (arange(4)//2).tolist(), (arange(4)%2).tolist(), onsourceScans.index(EQScan)
for spw_index in range(spwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)       # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)[pPol]
    GainP = GainP + [np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])]
    pCalVisX, pCalVisY = np.mean(chAvgVis[0] / (GainP[spw_index][0,ant0]* GainP[spw_index][0,ant1].conjugate()), axis=1), np.mean(chAvgVis[1] / (GainP[spw_index][1,ant0]* GainP[spw_index][1,ant1].conjugate()), axis=1)
    #-------- Antenna-based Gain
    AeSeqX = AeSeqX + [2.0* kb* chAvgTsys[antMap, spw_index, 0, scan_index]* abs(gainComplex(pCalVisX))**2] # Ae x Seq (X-pol)
    AeSeqY = AeSeqY + [2.0* kb* chAvgTsys[antMap, spw_index, 1, scan_index]* abs(gainComplex(pCalVisY))**2] # Ae x Seq (Y-pol)
#
#
GainP = np.array(GainP) # GainP[spw, pol, ant, time]
AeSeqX, AeSeqY = np.array(AeSeqX), np.array(AeSeqY)
######## Outputs from SSOflux.py :
# AeSeqX[spw, antMap], AeSeqX[spw, antMap] : (relative) aperture efficiency, assuming 1 Jy
#
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
for spw_index in range(1,spwNum):
    BPList[spw_index] = (BPList[spw_index].transpose(2,0,1)* spwTwiddle[:,:,spw_index]).transpose(1,2,0)
#
#-------- Flux models for solar system objects
msmd.done()
execfile(SCR_DIR + 'SSOflux.py'); logfile.write(FLScaleText + '\n')
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
        if len(SAantennas) < 4: continue #  Too few antennas
        SAantMap, SAblMap, SAblInv = np.array(antMap)[SAantennas].tolist(), np.array(blMap)[SAbl].tolist(), np.array(blInv)[SAbl].tolist()
        #print 'Subarray : ',; print antList[SAantMap]
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], SSOscanID[sso_index])
        chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        BPCaledXspec = (tempSpec / (BPList[spw_index][SAant0][:,polYindex]* BPList[spw_index][SAant1][:,polXindex].conjugate())).transpose(2,3,1,0) #  BPCaledXspe[pol, ch, bl, time]
        chAvgVis = (np.mean(BPCaledXspec[:, chRange], axis=1)[pPol].transpose(0,2,1)/SSOmodelVis[sso_index, spw_index,SAbl]).transpose(0,2,1)
        GainX, GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[0]), np.apply_along_axis( gainComplex, 0, chAvgVis[1])
        Ta = SSOflux[sso_index, spw_index]* AeNominal[SAantMap] / (2.0* kb)
        AeX[SAantennas, spw_index, sso_index] = 2.0* kb* np.percentile(abs(GainX), 75, axis=1)**2 * (Ta + chAvgTsys[SAantMap, spw_index, 0, scan_index]) / SSOflux[sso_index, spw_index]
        AeY[SAantennas, spw_index, sso_index] = 2.0* kb* np.percentile(abs(GainY), 75, axis=1)**2 * (Ta + chAvgTsys[SAantMap, spw_index, 1, scan_index]) / SSOflux[sso_index, spw_index]
        #
        #index = np.where(AeX[:, spw_index, sso_index] > 1.0)[0].tolist()
        #if np.std(np.append(AeX[index, spw_index, sso_index]/AeNominal[SAantMap], AeY[index, spw_index, sso_index]/AeNominal[SAantMap])) > 0.2:
        #    SSO_flag[sso_index] *= 0.0
        #
    #
#
#-------- SSO Flagging
for sso_index in range(SSONum):
    for spw_index in range(spwNum):
        index = np.where(AeX[:, spw_index, sso_index] > 1.0)[0].tolist()
        FLX_stat, FLY_stat = AeX[index, spw_index, sso_index]/AeNominal[np.array(antMap)[index].tolist()], AeY[index, spw_index, sso_index]/AeNominal[np.array(antMap)[index].tolist()]
        if np.std(FLX_stat) / np.median(FLX_stat) > 0.3: SSO_flag[sso_index] = 0.0
        if np.std(FLY_stat) / np.median(FLY_stat) > 0.3: SSO_flag[sso_index] = 0.0
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
    EQflux[spw_index], EQflux[spw_index + spwNum] = np.percentile(np.array(FLX), 25), np.percentile(np.array(FLY), 25)
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
        for pol_index in range(2):
            text_sd = '  %4.1f%% ' % (100.0* AEFF[ant_index, pol_index, spw_index]); logfile.write(text_sd); print text_sd,
        #
    #
    logfile.write('\n'); print ''
##
logfile.write('\n'); print ''
#-------- XY phase using BP scan
interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spw[0], BPScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PA = np.arctan2( np.sin(PA), np.cos(PA))
XYphase, caledVis = [], []
scan_index = onsourceScans.index(BPScan)
for spw_index in range(spwNum):
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], BPScan)
    timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    SEFD = 2.0* kb* chAvgTsys[antMap,spw_index, :,scan_index] / Ae[:,:,spw_index]
    caledVis.append(np.mean((chAvgVis / (GainP[polYindex][:,ant0]* GainP[polXindex][:,ant1].conjugate())).transpose(2, 0, 1)* np.sqrt(SEFD[ant0][:,polYindex].T* SEFD[ant1][:,polXindex].T), axis=2).T)
#
caledVis = np.array(caledVis)   # [spw, pol, time]
#QUsolution = XXYY2QU(PA, np.mean(caledVis[:,[0,3]], axis=0))
QUsolution = np.array([catalogStokesQ.get(BPcal), catalogStokesU.get(BPcal)])
#-------- XY phase cal in Bandpass table
for spw_index in range(spwNum):
    XYphase = XY2Phase(PA, QUsolution[0], QUsolution[1], caledVis[spw_index][[1,2]])
    XYsign = np.sign(np.cos(XYphase))
    BPList[spw_index][:,1] *= XYsign
    print 'SPW[%d] : XY phase = %6.1f [deg] sign = %3.0f' % (spw[spw_index], 180.0*XYphase/pi, XYsign)
#
#-------- Flux Density
ScanFlux = np.zeros([scanNum, spwNum, 4])
ScanSlope= np.zeros([scanNum, spwNum, 4])
ErrFlux  = np.zeros([scanNum, spwNum, 4])
ScanEL     = np.zeros([scanNum])
print '---Flux densities of sources ---'
pp = PdfPages('FL_' + prefix + '_' + UniqBands[band_index] + '.pdf')
polLabel = ['I', 'Q', 'U', 'V']
Pcolor   = ['black', 'blue', 'red', 'green']
for scan_index in range(scanNum):
    figScan = plt.figure(scan_index, figsize = (11, 8))
    figScan.suptitle(prefix + ' ' + UniqBands[band_index])
    figScan.text(0.75, 0.95, qa.time('%fs' % timeStamp[0], form='ymd')[0])
    figScan.text(0.45, 0.05, 'Projected baseline [m]')
    figScan.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90)
    ScanEL[scan_index] = np.median(OnEL[:,scan_index])
    text_sd = ' %02d %010s EL=%4.1f deg' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0* ScanEL[scan_index]/pi ); logfile.write(text_sd + '\n'); print text_sd
    figScan.text(0.05, 0.95, text_sd)
    if(onsourceScans[scan_index] in SSOscanID):
        SSO_flag = T
        SSO_ID = SSOscanID.index(onsourceScans[scan_index])
        text_sd = ' SPW  Frequency    I               Q               U               V             | Model I'; logfile.write(text_sd + '\n'); print text_sd
    else:
        SSO_flag = F
        text_sd = ' SPW  Frequency    I               Q               U               V               %Pol     EVPA '; logfile.write(text_sd + '\n'); print text_sd
    #
    text_sd = ' ------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
    BPCaledXspec = []
    SEFD = np.ones([spwNum, 2, UseAntNum])
    for spw_index in range(spwNum):
        #text_sd = 'SPW%02d %5.1f GHz ' % (spw[spw_index], centerFreqList[spw_index]); logfile.write(text_sd); print text_sd,
        atmCorrect = np.exp(-onTau[spw_index, scan_index])
        #-------- Sub-array with unflagged antennas (short baselines)
        if SSO_flag:
            SAantennas, SAbl, SAblFlag, SAant0, SAant1 = subArrayIndex(uvFlag[SSO_ID, spw_index])
            SAantMap, SAblMap, SAblInv = np.array(antMap)[SAantennas].tolist(), np.array(blMap)[SAbl].tolist(), np.array(blInv)[SAbl].tolist()
            TA = Ae[SAantennas,:,spw_index]* SSOflux0[SSO_ID, spw_index]* atmCorrect  / (2.0* kb)
        else:
            SAantennas, SAbl, SAblFlag, SAant0, SAant1 = range(UseAntNum), range(UseBlNum), np.ones([blNum]), ant0, ant1
            SAantMap, SAblMap, SAblInv = antMap, blMap, blInv
            TA = 0.0
        #
        SEFD[spw_index][:,SAantennas]  = (2.0* kb* (chAvgTsys[SAantMap,spw_index, :,scan_index] + TA) / (Ae[SAantennas,:,spw_index]* atmCorrect)).T
        SAantNum = len(SAantennas); SAblNum = len(SAblMap)
        if SAblNum < 6:
            text_sd = ' Only %d baselines for short enough sub-array. Skip!' % (SAblNum) ; logfile.write(text_sd + '\n'); print text_sd
            continue
        #
        #
        #-------- UV distance
        timeStamp, UVW = GetUVW(msfile, spw[spw_index], onsourceScans[scan_index])
        uvw = np.mean(UVW[:,SAblMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], onsourceScans[scan_index])
        timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); UseChNum = len(chRange)
        #timeThresh = np.median(diff(timeStamp))
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
        tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        #-------- Bandpass Calibration
        BPCaledXspec = BPCaledXspec + [(tempSpec / (BPList[spw_index][SAant0][:,polYindex]* BPList[spw_index][SAant1][:,polXindex].conjugate())).transpose(2,3,1,0)] # Bandpass Cal
    #
    #-------- Antenna-based Gain
    if SAblNum < 6: continue
    BPCaledXspec = np.array(BPCaledXspec)   # BPCaledXspec[spw, pol, ch, bl, time]
    if(SSO_flag): chAvgVis =(np.mean(BPCaledXspec[:,:, chRange], axis=(0,2)).transpose(0,2,1) / SSOmodelVis[SSO_ID, spw_index][SAbl]).transpose(0,2,1)
    else: chAvgVis = np.mean(BPCaledXspec[:, :, chRange], axis=(0,2))
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    pCalVis = (BPCaledXspec.transpose(0,2,1,3,4) / (GainP[polYindex][:,ant0[0:SAblNum]]* GainP[polXindex][:,ant1[0:SAblNum]].conjugate()))[:,chRange]
    for spw_index in range(spwNum):
        chAvgVis = np.mean(pCalVis[spw_index], axis=(0,3))[[0,3]]* np.sqrt(SEFD[spw_index][:, ant0[0:SAblNum]]* SEFD[spw_index][:, ant1[0:SAblNum]])
        indivRelGain = abs(gainComplexVec(chAvgVis.T))
        SEFD[spw_index][:,SAantennas] *= ((np.percentile(indivRelGain, 75, axis=0) / indivRelGain)**2).T
    #
    AmpCalVis = (pCalVis.transpose(1,4,0,2,3)* np.sqrt(SEFD[:,polYindex][:,:,ant0[0:SAblNum]]* SEFD[:,polXindex][:,:,ant1[0:SAblNum]])).transpose(2,4,3,0,1) # AmpCalVis[spw,bl,pol,ch,time]
    for spw_index in range(spwNum):
        StokesI_PL = figScan.add_subplot( 2, spwNum, spw_index + 1 )
        StokesP_PL = figScan.add_subplot( 2, spwNum, spwNum + spw_index + 1 )
        text_sd = ' SPW%02d %5.1f GHz' % (spw[spw_index], centerFreqList[spw_index]); logfile.write(text_sd); print text_sd,
        #pCalVis = (pCalVis.transpose(0,3,1,2)* np.sqrt(SEFD[ant0[0:SAblNum]][:,polYindex].T* SEFD[ant0[0:SAblNum]][:,polXindex].T))[chRange].transpose(3, 2, 0, 1)
        Stokes = np.zeros([4,SAblNum], dtype=complex)
        for bl_index in range(SAblNum):
            Minv = InvMullerVector(DxSpec[SAant1[bl_index], spw_index][chRange], DySpec[SAant1[bl_index], spw_index][chRange], DxSpec[SAant0[bl_index], spw_index][chRange], DySpec[SAant0[bl_index], spw_index][chRange], np.ones(UseChNum))
            Stokes[:,bl_index] = PS.reshape(4, 4*PAnum).dot(Minv.reshape(4, 4*UseChNum).dot(AmpCalVis[spw_index][bl_index].reshape(4*UseChNum, PAnum)).reshape(4*PAnum)) / (PAnum* UseChNum)
        #
        StokesVis, StokesErr = Stokes.real, Stokes.imag
        for pol_index in range(4):
            visFlag = np.where(abs(StokesVis[pol_index] - np.percentile(StokesVis[pol_index], 75))/np.percentile(StokesVis[0], 75) < 0.2 )[0]
            if len(visFlag) < 6:
                text_sd = 'Only %d vis.    ' % (len(visFlag)) ; logfile.write(text_sd + '\n'); print text_sd,
                continue
            #
            weight = np.zeros(SAblNum); weight[visFlag] = 1.0/np.var(StokesVis[pol_index][visFlag])
            P, W = np.c_[np.ones(SAblNum), uvDist], np.diag(weight)
            PtWP_inv = scipy.linalg.inv(np.dot(P.T, np.dot(W, P))) 
            solution, solerr = PtWP_inv.dot(P.T.dot(weight* StokesVis[pol_index])),  np.sqrt(np.diag(PtWP_inv))
            if solution[1] > -2.0* solerr[1]: solution[0] = np.median(StokesVis[pol_index]); solution[1] = 0.0
            ScanFlux[scan_index, spw_index, pol_index], ScanSlope[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index] = solution[0], solution[1], solerr[0] 
            text_sd = '%6.3f (%.3f) ' % (ScanFlux[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index]); logfile.write(text_sd); print text_sd,
        #
        StokesI_PL.plot( uvDist, StokesVis[0], '.', label=polLabel[0], color=Pcolor[0])
        StokesP_PL.plot( uvDist, StokesVis[1], '.', label=polLabel[1], color=Pcolor[1])
        StokesP_PL.plot( uvDist, StokesVis[2], '.', label=polLabel[2], color=Pcolor[2])
        StokesP_PL.plot( uvDist, StokesVis[3], '.', label=polLabel[3], color=Pcolor[3])
        if(SSO_flag):
            text_sd = '| %6.3f ' % (SSOflux0[SSO_ID, spw_index]); logfile.write(text_sd); print text_sd,
            logfile.write('\n'); print ''
        else: 
            text_sd = '%6.3f   %6.1f ' % (100.0* np.sqrt(ScanFlux[scan_index, spw_index, 1]**2 + ScanFlux[scan_index, spw_index, 2]**2)/ScanFlux[scan_index, spw_index, 0], np.arctan2(ScanFlux[scan_index, spw_index, 2],ScanFlux[scan_index, spw_index, 1])*90.0/pi); logfile.write(text_sd); print text_sd,
            logfile.write('\n'); print ''
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
    #
    freqArray = np.array(centerFreqList)[range(spwNum)]; meanFreq = np.mean(freqArray); relFreq = freqArray - meanFreq
    text_sd = ' ------------------------------------------------------------------------------------------------'; logfile.write(text_sd); print text_sd,
    logfile.write('\n'); print ''
    pflux, pfluxerr = np.zeros(4), np.zeros(4)
    text_sd = ' mean  %5.1f GHz' % (meanFreq); logfile.write(text_sd); print text_sd,
    for pol_index in range(4):
        sol, solerr = linearRegression(relFreq, ScanFlux[scan_index, :, pol_index], ErrFlux[scan_index, :, pol_index] ); pflux[pol_index], pfluxerr[pol_index] = sol[0], solerr[0]
        text_sd = '%6.3f (%.3f) ' % (pflux[pol_index], pfluxerr[pol_index]) ; logfile.write(text_sd); print text_sd,
    #
    text_sd = '%6.3f   %6.1f ' % (100.0* np.sqrt(pflux[1]**2 + pflux[2]**2)/pflux[0], np.arctan2(pflux[2],pflux[1])*90.0/pi); logfile.write(text_sd); print text_sd,
    logfile.write('\n'); print '\n'
    if COMPDB:
        if not SSO_flag:
            print ' -------- Comparison with ALMA Calibrator Catalog --------'
            au.searchFlux(sourcename='%s' % (sourceList[sourceIDscan[scan_index]]), band=int(UniqBands[band_index][3:5]), date=timeLabel[0:10], maxrows=3)
            print '\n'
        #
    #
#
logfile.close()
plt.close('all')
pp.close()
np.save(prefix + '-' + UniqBands[band_index] + '.Flux.npy', ScanFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Ferr.npy', ErrFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Source.npy', np.array(sourceList)[sourceIDscan])
np.save(prefix + '-' + UniqBands[band_index] + '.EL.npy', ScanEL)
msmd.done()
