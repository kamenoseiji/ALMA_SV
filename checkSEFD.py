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
Tatm_OFS  = 15.0     # Ambient-load temperature - Atmosphere temperature
AeNominal = 0.6* 0.25* np.pi* antDia**2      # Nominal Collecting Area
kb        = 1.38064852e3
#-------- Check Scans for atmCal
logfile = open(prefix + '-' + UniqBands[band_index] + '-Flux.log', 'w') 
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
#-------- Bandpass Table
print '---Generating antenna-based bandpass table'
BPList = []
for spw_index in spw:
    BP_ant, XY_BP, XYdelay, Gain = BPtable(msfile, spw_index, BPScan, blMap, blInv)
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
print '---Gain Equalization'
AeSeqX, AeSeqY = [], []  # effective area x flux density of the equalizer
polXindex, polYindex, scan_index = (arange(4)//2).tolist(), (arange(4)%2).tolist(), onsourceScans.index(EQScan)
for spw_index in range(spwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = ParaPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)       # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0]* BPList[spw_index][ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    #-------- Antenna-based Gain
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainX, GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[0]), np.apply_along_axis( gainComplex, 0, chAvgVis[1])
    BLphsX, BLphsY = GainX[ant0]* GainX[ant1].conjugate() / abs(GainX[ant0]* GainX[ant1]), GainY[ant0]* GainY[ant1].conjugate() / abs(GainY[ant0]* GainY[ant1])
    pCalVisX, pCalVisY = np.mean(chAvgVis[0] / BLphsX, axis=1), np.mean(chAvgVis[1] / BLphsY, axis=1)
    #-------- Antenna-based Gain
    AeSeqX = AeSeqX + [2.0* kb* chAvgTsys[:, spw_index, 0, scan_index]* abs(gainComplex(pCalVisX))**2] # Ae x Seq (X-pol)
    AeSeqY = AeSeqY + [2.0* kb* chAvgTsys[:, spw_index, 1, scan_index]* abs(gainComplex(pCalVisY))**2] # Ae x Seq (X-pol)
#
AeSeqX, AeSeqY = np.array(AeSeqX), np.array(AeSeqY)
#-------- Flux models for solar system objects
msmd.done()
execfile(SCR_DIR + 'SSOflux.py')
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
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    #-------- Baseline-based cross power spectra
    tempSpec = ParaPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)       # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec / (BPList[spw_index][np.array(ant0)[SAbl].tolist()]* BPList[spw_index][np.array(ant1)[SAbl].tolist()].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    chAvgVis = (np.mean(BPCaledXspec[:, chRange], axis=1).transpose(0,2,1)/FCSmodelVis[spw_index,SAbl]).transpose(0,2,1)
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

#-------- Flux Density
ScanFlux = np.zeros([scanNum, spwNum])
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
    for spw_index in range(spwNum):
        # atmCorrect = np.exp(-Tau0med[spw_index]/ np.sin(np.median(OnEL[:, scan_index])))
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
        tempSpec = ParaPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # tempSpec[timeNum, blNum, polNum, chNum]
        #-------- Bandpass Calibration
        BPCaledXspec = (tempSpec / (BPList[spw_index][np.array(ant0)[SAbl].tolist()]* BPList[spw_index][np.array(ant1)[SAbl].tolist()].conjugate())).transpose(2,3,1,0) # BPCaledXspec[polNum, chNum, blNum, timeNum]
        #-------- Antenna-based Gain
        if(SSO_flag): chAvgVis =(np.mean(BPCaledXspec[:, chRange], axis=1).transpose(0,2,1) / SSOmodelVis[SSO_ID, spw_index][SAbl]).transpose(0,2,1)
        else: chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)  # chAvgVis[polNum, blNum, timeNum]
        GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])  # GainP[polNum, antNum, timeNum]
        pCalVis = np.mean(chAvgVis / (GainP[:,ant0[0:SAblNum]]* GainP[:,ant1[0:SAblNum]].conjugate()), axis=2)    #  pCalVi[polNum, blNum]
        pCalVis = pCalVis* np.sqrt(SEFD[np.array(ant0)[SAbl].tolist()].T* SEFD[np.array(ant1)[SAbl].tolist()].T) # pCalVis[polNum,blNum]
        IVis = np.mean(pCalVis, axis=0).real
        visFlag = np.where(abs(IVis - np.median(IVis))/np.median(IVis) < 0.2)[0]
        if len(visFlag) < 4: text_sd = text_sd + ' Too few vis '; continue
        weight = np.zeros(SAblNum); weight[visFlag] = 1.0/np.var(IVis[visFlag])
        P, W = np.c_[np.ones(SAblNum), uvDist], np.diag(weight)
        PtWP_inv = scipy.linalg.inv(np.dot(P.T, np.dot(W, P)))
        ScanFlux[scan_index, spw_index], ErrFlux[scan_index, spw_index] = np.dot(PtWP_inv, np.dot(P.T, np.dot(W, IVis)))[0], np.sqrt(PtWP_inv[0,0])
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
