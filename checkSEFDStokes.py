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
Dfile = open(SCR_DIR + 'Dterm.B' + `int(UniqBands[band_index][3:5])` + '.data')
Dlines = Dfile.readlines()
Dfile.close()
Dx, Dy = np.zeros([antNum, spwNum], dtype=complex), np.zeros([antNum, spwNum], dtype=complex)
for ant_index in range(antNum):
    for Dline in Dlines:
        if antList[ant_index] in Dline:
            Dx[ant_index, int(Dline.split(' ')[1])] = float(Dline.split(' ')[2]) + (0.0+1.0j)* float(Dline.split(' ')[3])
            Dy[ant_index, int(Dline.split(' ')[1])] = float(Dline.split(' ')[4]) + (0.0+1.0j)* float(Dline.split(' ')[5])
        #
    #
#
#-------- Bandpass Table
print '---Generating antenna-based bandpass table'
XYdelayList, XYList, BPList = [], [], []
for spw_index in spw:
    BP_ant, XY_BP, XYdelay = BPtable(msfile, spw_index, BPScan, blMap, blInv)
    BPList = BPList + [BP_ant]
    XYList = XYList + [XY_BP]
    XYdelayList = XYdelayList + [XYdelay]
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
#-------- Load autocorrelation power spectra
print '---Loading autocorr power spectra'
OnSpecList, OffSpecList, AmbSpecList, HotSpecList = [], [], [], []
for ant_index in range(UseAntNum):
    for spw_index in range(spwNum):
        progress = (1.0* ant_index* spwNum + spw_index + 1.0) / (UseAntNum* spwNum)
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
        timeXY, Pspec = GetPSpec(msfile, antMap[ant_index], spw[spw_index])
        OffSpecList.append(Pspec[pPol][:,:,offTimeIndex])
        AmbSpecList.append(Pspec[pPol][:,:,ambTimeIndex])
        HotSpecList.append(Pspec[pPol][:,:,hotTimeIndex])
        for scan_index in range(scanNum):
            OnSpecList.append(np.mean( Pspec[pPol][:,:,OnTimeIndex[scan_index]], axis=2 ))
        #
    #
#
sys.stderr.write('\n')
sys.stderr.flush()
#-------- Load Az El position
azelTime, AntID, AZ, EL = GetAzEl(msfile)
OnAZ, OnEL, OffEL = np.ones([UseAntNum, scanNum]), np.ones([UseAntNum, scanNum]), np.ones([UseAntNum, len(offTimeIndex)])
for ant_index in range(UseAntNum):
    azelTime_index = np.where( AntID == antMap[ant_index] )[0].tolist()
    for scan_index in range(scanNum):
        refTime = np.median(msmd.timesforscan(onsourceScans[scan_index]))
        OnAZ[ant_index, scan_index] = AZ[azelTime_index[argmin(abs(azelTime[azelTime_index] - refTime))]]
        OnEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - refTime))]]
    #
    for time_index in range(len(offTimeIndex)):
        OffEL[ant_index, time_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - offTime[time_index]))]]
    #
#
secZ = 1.0 / np.sin( OffEL )
#-------- Time-interpolation of ambient and hot
print '---Analyzing Trec and Tsky using atmCal scans'
chAvgTrx, chAvgTsky, chAvgTsys = np.zeros([UseAntNum, spwNum, 2, len(offTime)]), np.zeros([UseAntNum, spwNum, 2, len(offTime)]), np.zeros([UseAntNum, spwNum, 2, scanNum])
TrxFlag, TsysFlag, onTau = np.ones([UseAntNum, spwNum, 2, len(offTime)]), np.ones([UseAntNum, spwNum, 2, scanNum]), np.zeros([spwNum, scanNum])
TrxList, TskyList = [], []
tempAmb, tempHot  = np.zeros([UseAntNum]), np.zeros([UseAntNum])
for ant_index in range(UseAntNum):
    tempAmb[ant_index], tempHot[ant_index] = GetLoadTemp(msfile, antMap[ant_index], spw[0])
    for spw_index in range(spwNum):
        AntSpwIndex = ant_index* spwNum + spw_index
        chNum = AmbSpecList[AntSpwIndex].shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        Trx, Tsky = np.zeros([2, chNum, len(offTime)]), np.zeros([2, chNum, len(offTime)])
        for pol_index in range(ppolNum):
            ambSpec, hotSpec = AmbSpecList[AntSpwIndex][pol_index], HotSpecList[AntSpwIndex][pol_index]
            SPL_amb = UnivariateSpline(ambTime, np.mean(ambSpec[chRange], axis=0), s=0.001)
            SPL_hot = UnivariateSpline(hotTime, np.mean(hotSpec[chRange], axis=0), s=0.001)
            for time_index in range(len(offTime)):
                Psamb = np.mean(ambSpec, axis=1)*SPL_amb(offTime[time_index])/np.mean(ambSpec[chRange])
                Pshot = np.mean(hotSpec, axis=1)*SPL_hot(offTime[time_index])/np.mean(hotSpec[chRange])
                Psoff = OffSpecList[AntSpwIndex][pol_index][:,time_index]
                Trx[pol_index, :, time_index] = (tempHot[ant_index]* Psamb - Pshot* tempAmb[ant_index]) / (Pshot - Psamb)
                Tsky[pol_index, :, time_index]= (Psoff* (tempHot[ant_index] - tempAmb[ant_index]) + tempAmb[ant_index]* Pshot - tempHot[ant_index]* Psamb) / (Pshot - Psamb)
                Phot, Pamb, Poff = np.mean(Pshot[chRange]), np.mean(Psamb[chRange]), np.mean(Psoff[chRange])
                chAvgTrx[ant_index, spw_index, pol_index, time_index] = (tempHot[ant_index]* Pamb - Phot* tempAmb[ant_index]) / (Phot - Pamb)
                chAvgTsky[ant_index, spw_index, pol_index, time_index]= (Poff* (tempHot[ant_index] - tempAmb[ant_index]) + tempAmb[ant_index]* Phot - tempHot[ant_index]* Pamb) / (Phot - Pamb)
            #-------- Trx flagging
            TrxTemp = chAvgTrx[ant_index, spw_index, pol_index]
            TrxFlag[ant_index, spw_index, pol_index, np.where( abs(TrxTemp - np.median(TrxTemp)) > 0.2* np.median(TrxTemp))[0].tolist()] = 0.0
            TrxFlag[ant_index, spw_index, pol_index, np.where( TrxTemp < 1.0)[0].tolist()] = 0.0
        #
        TrxList.append(Trx)
        TskyList.append(Tsky)
    #
#
#-------- Tau and TantN fitting
param = [0.0, 0.05]     # Initial Parameter
Tau0, TantN, Trms = np.zeros([UseAntNum, spwNum, 2]), np.zeros([UseAntNum, spwNum, 2]), np.zeros([UseAntNum, spwNum, 2])
for ant_index in range(UseAntNum):
    for spw_index in range(spwNum):
        for pol_index in range(2):
            fit = scipy.optimize.leastsq(residTskyTransfer, param, args=(tempAmb[ant_index]-Tatm_OFS, secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index]))
            TantN[ant_index, spw_index, pol_index] = fit[0][0]
            Tau0[ant_index, spw_index, pol_index]  = fit[0][1]
        #
    #
    TantN[ant_index] = np.median(TantN[ant_index])* np.ones([spwNum, 2])
#
Tau0err, Tau0med = np.std(Tau0, axis=(0,2)), np.mean(Tau0, axis=(0,2))
#-------- Trx Transfer
for ant_index in range(UseAntNum):
    for spw_index in range(spwNum):
        for pol_index in range(ppolNum):
            flgIndex = np.where( TrxFlag[ant_index, spw_index, pol_index] == 0 )[0].tolist()
            vldIndex = np.where( TrxFlag[ant_index, spw_index, pol_index] == 1 )[0].tolist()
            if(len(vldIndex) > 0):
                chAvgTrx[ant_index, spw_index, pol_index, flgIndex] = np.median( chAvgTrx[ant_index, spw_index, pol_index, vldIndex] )
                chAvgTsky[ant_index, spw_index, pol_index, flgIndex]= (tempAmb[ant_index]-Tatm_OFS)* (1.0 - np.exp(-Tau0[ant_index, spw_index, pol_index]* secZ[ant_index, flgIndex])) + TantN[ant_index, spw_index, pol_index]
            #
        #
    #
#
#-------- Tsys for scans
for scan_index in range(scanNum):
    OnTimeRange = timeXY[OnTimeIndex[scan_index]]
    ambTimeIndex = argmin(abs(ambTime - OnTimeRange[0]))
    offTimeIndex = argmin(abs(offTime - OnTimeRange[0]))
    chAvgTsys[:,:,:,scan_index] = chAvgTrx[:,:,:,ambTimeIndex] + chAvgTsky[:,:,:,offTimeIndex] 
    TsysFlag[:,:,:,scan_index] = TrxFlag[:,:,:,ambTimeIndex]
#
#-------- Tau on source
for scan_index in range(scanNum):
    OnTimeRange = timeXY[OnTimeIndex[scan_index]]
    offTimeIndex = argmin(abs(offTime - OnTimeRange[0]))
    tempTau = -np.log((chAvgTsky[:,:,:,offTimeIndex] - TantN - np.median(tempAmb) + Tatm_OFS) / (2.718 - np.median(tempAmb) + Tatm_OFS))
    onTau[:,scan_index] = np.median(tempTau.transpose(1,2,0).reshape(spwNum, -1), axis=1)
#
for spw_index in range(spwNum): text_sd = 'SPW=%d : Tau(zenith) = %6.4f +- %6.4f' % (spw[spw_index], Tau0med[spw_index], Tau0err[spw_index]); logfile.write(text_sd + '\n'); print text_sd
#
np.save(prefix +  '-' + UniqBands[band_index] + '.Trx.npy', TrxList) 
np.save(prefix +  '-' + UniqBands[band_index] + '.Tsky.npy', TskyList) 
np.save(prefix +  '-' + UniqBands[band_index] + '.TrxFlag.npy', TrxFlag) 
np.save(prefix +  '-' + UniqBands[band_index] + '.Tau0.npy', Tau0) 
#-------- Antenna-dependent leakage noise
param = [0.0]
text_sd = ' TantN:'; logfile.write(text_sd); print text_sd,
for spw_index in range(spwNum): text_sd = ' PolX   SPW%02d  PolY           |' % (spw[spw_index]); logfile.write(text_sd); print text_sd,
#
logfile.write('\n'); print ' '
text_sd = ' ----:--------------------------------+-------------------------------+-------------------------------+-------------------------------+'; logfile.write(text_sd + '\n'); print text_sd
for ant_index in range(UseAntNum):
    text_sd = antList[antMap[ant_index]] + ' : '; logfile.write(text_sd); print text_sd,
    for spw_index in range(spwNum):
        for pol_index in range(2):
            fit = scipy.optimize.leastsq(residTskyTransfer2, param, args=(tempAmb[ant_index]-Tatm_OFS, Tau0med[spw_index], secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index]))
            TantN[ant_index, spw_index, pol_index] = fit[0][0]
            resid = residTskyTransfer([fit[0][0], Tau0med[spw_index]], tempAmb[ant_index]-Tatm_OFS, secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index])
            Trms[ant_index, spw_index, pol_index]  = sqrt(np.dot(resid, resid) / len(np.where(TrxFlag[ant_index, spw_index, pol_index] > 0.0)[0]))
            text_sd = '%4.1f (%4.1f) K ' % (TantN[ant_index, spw_index, pol_index], Trms[ant_index, spw_index, pol_index])
            logfile.write(text_sd); print text_sd,
        text_sd = '|'; logfile.write(text_sd); print text_sd,
    logfile.write('\n'); print ''
logfile.write('\n'); print ''
np.save(prefix +  '-' + UniqBands[band_index] + '.TantN.npy', TantN) 
#-------- Trx
text_sd = ' Trec: '; logfile.write(text_sd); print text_sd,
for spw_index in range(spwNum): text_sd = ' SPW%02d X        Y |' % (spw[spw_index]); logfile.write(text_sd); print text_sd,
logfile.write('\n'); print ' '
text_sd =  ' ----:--------------------+-------------------+-------------------+-------------------+'; logfile.write(text_sd + '\n'); print text_sd
for ant_index in range(UseAntNum):
    text_sd =  antList[antMap[ant_index]] + ' : '; logfile.write(text_sd); print text_sd,
    for spw_index in range(spwNum):
        for pol_index in range(2):
            text_sd = '%6.1f K' % (np.median(chAvgTrx[ant_index, spw_index, pol_index]))
            logfile.write(text_sd); print text_sd,
        text_sd = '|'; logfile.write(text_sd); print text_sd,
    logfile.write('\n'); print ' '
print ' '
#-------- Tsys
text_sd = ' Tsys: '; logfile.write(text_sd); print text_sd,
for spw_index in range(spwNum): text_sd = ' SPW%02d X        Y |' % (spw[spw_index]); logfile.write(text_sd); print text_sd,
logfile.write('\n'); print ' '
text_sd =  ' ----:--------------------+-------------------+-------------------+-------------------+'; logfile.write(text_sd + '\n'); print text_sd
for ant_index in range(UseAntNum):
    text_sd =  antList[antMap[ant_index]] + ' : '; logfile.write(text_sd); print text_sd,
    for spw_index in range(spwNum):
        for pol_index in range(2):
            text_sd = '%6.1f K' % (np.median(chAvgTrx[ant_index, spw_index, pol_index]) + np.median(chAvgTsky[ant_index, spw_index, pol_index]))
            logfile.write(text_sd); print text_sd,
        text_sd = '|'; logfile.write(text_sd); print text_sd,
    logfile.write('\n'); print ' '
#-------- Plot optical depth
if PLOTTAU: plotTau(prefix + '_' + UniqBands[band_index], antList[antMap], spw, secZ, (chAvgTsky.transpose(3,0,1,2) - TantN).transpose(1,2,3,0), np.median(tempAmb) - Tatm_OFS, Tau0med, TrxFlag, 2.0*np.median(chAvgTsky), PLOTFMT) 
if PLOTTSYS: plotTsys(prefix + '_' + UniqBands[band_index], antList[antMap], ambTime, spw, TrxList, TskyList, PLOTFMT)
##-------- Equalization using EQ scan
AeSeqX, AeSeqY = [], []  # effective area x flux density of the equalizer
polXindex, polYindex, scan_index = (arange(4)//2).tolist(), (arange(4)%2).tolist(), onsourceScans.index(EQScan)
for spw_index in range(spwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)       # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    #-------- Antenna-based Gain
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)[pPol]
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
#-------- XY phase using EQ scan
XYphase, caledVis = [], []
scan_index = onsourceScans.index(EQScan)
for spw_index in range(spwNum):
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], EQScan)
    timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    timeThresh = np.median(diff(timeStamp))
    AzScan, ElScan = np.zeros(timeNum), np.zeros(timeNum)
    for time_index in range(timeNum):
        AzScan[time_index], ElScan[time_index] = AzElMatch(timeStamp[time_index], azelTime, AntID, UseAnt[refantID], timeThresh, AZ, EL)
    #
    PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    #XYdlSpec = delay_cal(np.ones([chNum], dtype=complex), XYdelayList[spw_index])
    XYdlSpec = XYList[spw_index].conjugate()
    BPCaledXspec[1] = (BPCaledXspec[1].transpose(1,2,0)* XYdlSpec).transpose(2,0,1)
    BPCaledXspec[2] = (BPCaledXspec[2].transpose(1,2,0)* XYdlSpec.conjugate()).transpose(2,0,1)
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    Gain = np.r_[np.apply_along_axis(gainComplex, 0, chAvgVis[0]), np.apply_along_axis(gainComplex, 0, chAvgVis[3])].reshape(2, UseAntNum, timeNum)
    Gain = Gain / abs(Gain)
    SEFD = 2.0* kb* chAvgTsys[:,spw_index, :,scan_index] / Ae[:,:,spw_index]
    caledVis.append(np.mean((chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())).transpose(2, 0, 1)* np.sqrt(SEFD[ant0][:,polYindex].T* SEFD[ant1][:,polXindex].T), axis=2).T)
#
caledVis = np.array(caledVis)
QUsolution = XXYY2QU(PA, np.mean(caledVis[:,[0,3]], axis=0))
for spw_index in range(spwNum):
    XYphase.append(XY2Phase(PA, QUsolution[0], QUsolution[1], caledVis[spw_index][[1,2]]))
    print 'SPW[%d] : XY phase = %5.1f deg' % (spw[spw_index], 180.0* XYphase[spw_index]/pi)
#
XYphase = np.array(XYphase)
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
        AzScan, ElScan = np.zeros(timeNum), np.zeros(timeNum)
        for time_index in range(timeNum):
            AzScan[time_index], ElScan[time_index] = AzElMatch(timeStamp[time_index], azelTime, AntID, UseAnt[refantID], timeThresh, AZ, EL)
        #
        PA = np.mean(AzEl2PA(AzScan, ElScan)) + BandPA[band_index]; PS = InvPAMatrix(PA)
        tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        #-------- Bandpass Calibration
        BPCaledXspec = (tempSpec / (BPList[spw_index][np.array(ant0)[SAbl].tolist()][:,polYindex]* BPList[spw_index][np.array(ant1)[SAbl].tolist()][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
        #-------- XY delay cal
        #XYdlSpec = delay_cal(np.ones([chNum], dtype=complex), XYdelayList[spw_index])
        XYdlSpec = XYList[spw_index].conjugate()
        BPCaledXspec[1] = (BPCaledXspec[1].transpose(1,2,0)* XYdlSpec).transpose(2,0,1)
        BPCaledXspec[2] = (BPCaledXspec[2].transpose(1,2,0)* XYdlSpec.conjugate()).transpose(2,0,1)
        #-------- Antenna-based Gain
        if(SSO_flag): chAvgVis =(np.mean(BPCaledXspec[:, chRange], axis=1).transpose(0,2,1) / SSOmodelVis[SSO_ID, spw_index][SAbl]).transpose(0,2,1)
        else: chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
        Gain = np.r_[np.apply_along_axis(gainComplex, 0, chAvgVis[0]), np.apply_along_axis(gainComplex, 0, chAvgVis[3])* np.exp((0.0 + 1.0j)* XYphase[spw_index])].reshape(2, SAantNum, timeNum)
        Gain = Gain / abs(Gain)
        pCalVis = np.mean(chAvgVis / (Gain[polYindex][:,ant0[0:SAblNum]]* Gain[polXindex][:,ant1[0:SAblNum]].conjugate()), axis=2)* np.sqrt(SEFD[np.array(ant0)[SAbl].tolist()][:,polYindex].T* SEFD[np.array(ant1)[SAbl].tolist()][:,polXindex].T)
        StokesVis, StokesErr = np.zeros([4,SAblNum]), np.zeros([4,SAblNum])
        for bl_index in range(SAblNum):
            Minv = InvMullerMatrix(Dx[SAantMap[ant1[bl_index]], spw_index], Dy[SAantMap[ant1[bl_index]], spw_index], Dx[SAantMap[ant0[bl_index]], spw_index], Dy[SAantMap[ant0[bl_index]], spw_index])
            Stokes = np.dot(PS, np.dot(Minv, pCalVis[:,bl_index]))
            StokesVis[:,bl_index], StokesErr[:,bl_index] = Stokes.real, Stokes.imag
        #
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
            # ScanFlux[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index] = np.median(StokesVis[pol_index]), np.std(StokesVis[pol_index])/np.sqrt(SAblNum - 1.0)
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
