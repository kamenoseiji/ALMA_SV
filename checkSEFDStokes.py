import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
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
##-------- Equalization using EQ scan
AeSeqX, AeSeqY = [], []  # effective area x flux density of the equalizer
polXindex, polYindex, scan_index = (arange(4)//2).tolist(), (arange(4)%2).tolist(), onsourceScans.index(EQScan)
for spw_index in range(spwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)       # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)[pPol]
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
    pCalVisX, pCalVisY = np.mean(chAvgVis[0] / (GainP[0,ant0]* GainP[0,ant1].conjugate()), axis=1), np.mean(chAvgVis[1] / (GainP[1,ant0]* GainP[1,ant1].conjugate()), axis=1)
    #-------- Antenna-based Gain
    AeSeqX = AeSeqX + [2.0* kb* chAvgTsys[:, spw_index, 0, scan_index]* abs(gainComplex(pCalVisX))**2] # Ae x Seq (X-pol)
    AeSeqY = AeSeqY + [2.0* kb* chAvgTsys[:, spw_index, 1, scan_index]* abs(gainComplex(pCalVisY))**2] # Ae x Seq (Y-pol)
#
AeSeqX, AeSeqY = np.array(AeSeqX), np.array(AeSeqY)
#-------- Flux models for solar system objects
msmd.done()
execfile(SCR_DIR + 'SSOflux.py'); logfile.write(FLScaleText + '\n')
SSOflux = SSOflux0* np.exp(-onTau.transpose(1,0)[indexList(np.array(SSOscanID), np.array(onsourceScans))])
##-------- Scaling with the flux calibrator
AeX, AeY = [], []
scan_index = onsourceScans.index(FCScan)
#-------- Sub-array with unflagged antennas (short baselines)
SAantennas, SAbl, SAblFlag, SAant0, SAant1 = subArrayIndex(FCSFlag[0])
SAantList = antList[np.array(antMap)[SAantennas]] 
SAantMap = np.array(antMap)[SAantennas].tolist()
SAantList = antList[SAantMap]
SAblMap = np.array(blMap)[SAbl].tolist()
SAblInv = np.array(blInv)[SAbl].tolist()
SAantNum = len(SAantMap); SAblNum = len(SAblMap)
for spw_index in range(spwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], FCScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); UseChNum = len(chRange)
    #-------- Baseline-based cross power spectra
    tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)       # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec / (BPList[spw_index][np.array(ant0)[SAbl].tolist()][:,polYindex]* BPList[spw_index][np.array(ant1)[SAbl].tolist()][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    chAvgVis = (np.mean(BPCaledXspec[:, chRange], axis=1)[pPol].transpose(0,2,1)/FCSmodelVis[spw_index,SAbl]).transpose(0,2,1)
    GainX, GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[0]), np.apply_along_axis( gainComplex, 0, chAvgVis[1])
    Ta = SSOflux[FCS_ID, spw_index]* AeNominal[SAantennas] / (2.0* kb)
    AeX = AeX + (2.0* kb* np.median(abs(GainX), axis=1)**2 * (Ta + chAvgTsys[SAantennas, spw_index, 0, scan_index]) / SSOflux[FCS_ID, spw_index]).tolist()
    AeY = AeY + (2.0* kb* np.median(abs(GainY), axis=1)**2 * (Ta + chAvgTsys[SAantennas, spw_index, 1, scan_index]) / SSOflux[FCS_ID, spw_index]).tolist()
    #
#
AeX, AeY = np.array(AeX).reshape(spwNum, SAantNum), np.array(AeY).reshape(spwNum, SAantNum)
##-------- Flux density of the equalizer, aligned in the power law
EQflux = np.r_[np.median((AeSeqX[:,SAantennas]/AeX), axis=1), np.median((AeSeqY[:,SAantennas]/AeY), axis=1)]
P = np.c_[ np.r_[np.log(centerFreqList),np.log(centerFreqList)], np.r_[np.ones(spwNum), np.zeros(spwNum)], np.r_[np.zeros(spwNum), np.ones(spwNum)]]
EQmodel = scipy.linalg.solve(np.dot(P.T, P), np.dot(P.T, np.log(EQflux)))   # alpha, logSx, logSy
EQflux = np.c_[np.exp(EQmodel[0]* np.log(centerFreqList) + EQmodel[1]), np.exp(EQmodel[0]* np.log(centerFreqList) + EQmodel[2])]
Ae     = np.c_[AeSeqX.T / EQflux[:,0], AeSeqY.T / EQflux[:,1]].reshape(UseAntNum, ppolNum, spwNum)
AEFF   = (Ae.transpose(1,2,0) / (0.25* pi*antDia**2)).transpose(2,0,1)
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
    SEFD = 2.0* kb* chAvgTsys[:,spw_index, :,scan_index] / Ae[:,:,spw_index]
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
    print 'SPW[%d] : XY phase = %6.1f sign = %3.0f' % (spw[spw_index], XYphase, XYsign)
#
#-------- Flux Density
ScanFlux = np.zeros([scanNum, spwNum, 4])
ErrFlux  = np.zeros([scanNum, spwNum, 4])
ScanEL     = np.zeros([scanNum])
print '---Flux densities of sources ---'
for scan_index in range(scanNum):
    ScanEL[scan_index] = np.median(OnEL[:,scan_index])
    text_sd = ' %02d %010s EL=%4.1f deg' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0* ScanEL[scan_index]/pi ); logfile.write(text_sd + '\n'); print text_sd
    text_sd = ' SPW  Frequency    I               Q               U               V             | Model I'; logfile.write(text_sd + '\n'); print text_sd
    text_sd = ' -----------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
    if(onsourceScans[scan_index] in SSOscanID):
        SSO_flag = T
        SSO_ID = SSOscanID.index(onsourceScans[scan_index])
    else:
        SSO_flag = F
    for spw_index in range(spwNum):
        text_sd = 'SPW%02d %5.1f GHz ' % (spw[spw_index], centerFreqList[spw_index]); logfile.write(text_sd); print text_sd,
        atmCorrect = np.exp(-onTau[spw_index, scan_index])
        #-------- Sub-array with unflagged antennas (short baselines)
        if SSO_flag:
            SAantennas, SAbl, SAblFlag, SAant0, SAant1 = subArrayIndex(uvFlag[SSO_ID, spw_index])
            SAantMap, SAblMap, SAblInv = np.array(antMap)[SAantennas].tolist(), np.array(blMap)[SAbl].tolist(), np.array(blInv)[SAbl].tolist()
            TA = Ae[:,:,spw_index]* SSOflux0[SSO_ID, spw_index]* atmCorrect  / (2.0* kb)
        else:
            SAantennas, SAbl, SAblFlag, SAant0, SAant1 = range(UseAntNum), range(UseBlNum), np.ones([blNum]), ant0, ant1
            SAantMap, SAblMap, SAblInv = antMap, blMap, blInv
            TA = 0.0
        #
        SEFD = 2.0* kb* (chAvgTsys[:,spw_index, :,scan_index] + TA) / (Ae[:,:,spw_index]* atmCorrect)
        SAantNum = len(SAantennas); SAblNum = len(SAblMap)
        if SAblNum < 3:
            text_sd = ' Only %d baselines for short enough sub-array. Skip!' % (SAblNum) ; logfile.write(text_sd + '\n'); print text_sd
            continue
        #
        #
        #-------- UV distance
        timeStamp, UVW = GetUVW(msfile, spw[spw_index], onsourceScans[scan_index])
        uvw = np.mean(UVW[:,SAblMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], onsourceScans[scan_index])
        timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        timeThresh = np.median(diff(timeStamp))
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
        tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        #-------- Bandpass Calibration
        BPCaledXspec = (tempSpec / (BPList[spw_index][np.array(ant0)[SAbl].tolist()][:,polYindex]* BPList[spw_index][np.array(ant1)[SAbl].tolist()][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
        #-------- Antenna-based Gain
        if(SSO_flag): chAvgVis =(np.mean(BPCaledXspec[:, chRange], axis=1).transpose(0,2,1) / SSOmodelVis[SSO_ID, spw_index][SAbl]).transpose(0,2,1)
        else: chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
        GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
        pCalVis = BPCaledXspec.transpose(1,0,2,3) / (GainP[polYindex][:,ant0[0:SAblNum]]* GainP[polXindex][:,ant1[0:SAblNum]].conjugate())
        pCalVis = (pCalVis.transpose(0,3,1,2)* np.sqrt(SEFD[np.array(ant0)[SAbl].tolist()][:,polYindex].T* SEFD[np.array(ant1)[SAbl].tolist()][:,polXindex].T))[chRange].transpose(3, 2, 0, 1)
        Stokes = np.zeros([4,SAblNum], dtype=complex)
        for bl_index in range(SAblNum):
            Minv = InvMullerVector(DxSpec[ant1[bl_index], spw_index][chRange], DySpec[ant1[bl_index], spw_index][chRange], DxSpec[ant0[bl_index], spw_index][chRange], DySpec[ant0[bl_index], spw_index][chRange], np.ones(UseChNum))
            Stokes[:,bl_index] = PS.reshape(4, 4*PAnum).dot(Minv.reshape(4, 4*UseChNum).dot(pCalVis[bl_index].reshape(4*UseChNum, PAnum)).reshape(4*PAnum)) / (PAnum* UseChNum)
        #
        StokesVis, StokesErr = Stokes.real, Stokes.imag
        for pol_index in range(4):
            visFlag = np.where(abs(StokesVis[pol_index] - np.median(StokesVis[pol_index]))/np.median(StokesVis[0]) < 0.2 )[0]
            if len(visFlag) < 4:
                text_sd = 'Only %d vis.    ' % (len(visFlag)) ; logfile.write(text_sd + '\n'); print text_sd,
                continue
            #
            weight = np.zeros(SAblNum); weight[visFlag] = 1.0/np.var(StokesVis[pol_index][visFlag])
            P, W = np.c_[np.ones(SAblNum), uvDist], np.diag(weight)
            PtWP_inv = scipy.linalg.inv(np.dot(P.T, np.dot(W, P))) 
            ScanFlux[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index] = np.dot(PtWP_inv, np.dot(P.T, np.dot(W, StokesVis[pol_index])))[0], np.sqrt(PtWP_inv[0,0])
            text_sd = '%6.3f (%.3f) ' % (ScanFlux[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index]); logfile.write(text_sd); print text_sd,
        #
        if(SSO_flag): text_sd = '| %6.3f ' % (SSOflux0[SSO_ID, spw_index]); logfile.write(text_sd); print text_sd,
        logfile.write('\n'); print ''
    #
    logfile.write('\n'); print ''
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
