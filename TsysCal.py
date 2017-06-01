def residTskyTransfer( param, Tamb, secz, Tsky, weight ):
    exp_Tau = np.exp( -param[1]* secz )
    return weight* (Tsky - (param[0] + 2.718* exp_Tau  + Tamb* (1.0 - exp_Tau)))
#
def residTskyTransfer0( param, Tamb, secz, Tsky, weight ):
    exp_Tau = np.exp( -param[0]* secz )
    return weight* (Tsky - (2.718* exp_Tau  + Tamb* (1.0 - exp_Tau)))
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
kb        = 1.38064852e3
#-------- Check Ambient Load Timing
print '---Checking time for ambient and hot load'
timeOFF, timeAMB, timeHOT = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT")
if len(timeAMB) == 0:
    #timeON  = msmd.timesforintent("CALIBRATE_ATMOSPHERE#ON_SOURCE")
    timeXY, Pspec = GetPSpec(msfile, 0, spwList[0])
    timeNum, chNum = Pspec.shape[2], Pspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    chAvgPower = np.mean(Pspec[0][chRange], axis=0)
    offTimeIndex = indexList(timeOFF, timeXY)
    hotTimeIndex = (np.array(offTimeIndex) - 1).tolist()
    ambTimeIndex = (np.array(offTimeIndex) - 2).tolist()
    #edge = np.where( abs(np.diff( chAvgPower )) > 0.1* np.median(chAvgPower))[0].tolist()
    #atmList = np.array(edge)[np.where( abs( np.diff(timeXY[edge]) - np.median(np.diff(timeXY[edge]))) < 0.5* np.median(np.diff(timeXY[edge])))[0].tolist()]
    #atmList = atmList.reshape(len(atmList)/2, 2)
    #ambTimeIndex, hotTimeIndex, offTimeIndex = atmList[:,0].tolist(), atmList[:,1].tolist(), (atmList[:,1] + 1).tolist()
    ambTime, hotTime, offTime = timeXY[ambTimeIndex], timeXY[hotTimeIndex], timeXY[offTimeIndex]
else:
    tb.open(msfile); timeXY = tb.query('ANTENNA1 == 0 && ANTENNA2 == 0 && DATA_DESC_ID == '+`spw[0]`).getcol('TIME'); tb.close()
    offTime, ambTime, hotTime = sort( list(set(timeXY) & set(timeOFF)) ), sort( list(set(timeXY) & set(timeAMB)) ), sort( list(set(timeXY) & set(timeHOT)) )
    offTimeIndex, ambTimeIndex, hotTimeIndex = indexList(offTime, timeXY),  indexList(ambTime, timeXY),  indexList(hotTime, timeXY)
#
OnTimeIndex = []
for scan_index in range(scanNum): OnTimeIndex.append( indexList(msmd.timesforscan(onsourceScans[scan_index]), timeXY) )
#-------- Load autocorrelation power spectra
print '---Loading autocorr power spectra'
OnSpecList, OffSpecList, AmbSpecList, HotSpecList = [], [], [], []
#maxTimeIndex = max(offTimeIndex + ambTimeIndex + hotTimeIndex + max(OnTimeIndex))
for ant_index in range(antNum):
    for spw_index in range(spwNum):
        progress = (1.0* ant_index* spwNum + spw_index + 1.0) / (antNum* spwNum)
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
        timeXY, Pspec = GetPSpec(msfile, ant_index, spw[spw_index])
        #if len(timeXY) < maxTimeIndex:
        #    chNum, addTimeNum = Pspec.shape[1], maxTimeIndex - len(timeXY) + 1
        #    Pspec = np.concatenate((Pspec, np.zeros([polNum, chNum, addTimeNum])), axis=2)
        #
        OffSpecList.append(Pspec[pPol][:,:,offTimeIndex])
        AmbSpecList.append(Pspec[pPol][:,:,ambTimeIndex])
        HotSpecList.append(Pspec[pPol][:,:,hotTimeIndex])
        for scan_index in range(scanNum):
            #OnSpecList.append(np.mean( Pspec[pPol][:,:,OnTimeIndex[scan_index]], axis=2 ))
            OnSpecList.append(np.median( Pspec[pPol][:,:,OnTimeIndex[scan_index]], axis=2 ))
        #
    #
#
sys.stderr.write('\n')
sys.stderr.flush()
#-------- Load Az El position
azelTime, AntID, AZ, EL = GetAzEl(msfile)
OnAZ, OnEL, OffEL = np.ones([antNum, scanNum]), np.ones([antNum, scanNum]), np.ones([antNum, len(offTimeIndex)])
for ant_index in range(antNum):
    azelTime_index = np.where( AntID == ant_index )[0].tolist()
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
chAvgTrx, chAvgTsky, chAvgTsys = np.zeros([antNum, spwNum, 2, len(offTime)]), np.zeros([antNum, spwNum, 2, len(offTime)]), np.zeros([antNum, spwNum, 2, scanNum])
TrxFlag, TsysFlag, onTau = np.ones([antNum, spwNum, 2, len(offTime)]), np.ones([antNum, spwNum, 2, scanNum]), np.zeros([spwNum, scanNum])
TrxList, TskyList = [], []
tempAmb, tempHot  = np.zeros([antNum]), np.zeros([antNum])
for ant_index in range(antNum):
    tempAmb[ant_index], tempHot[ant_index] = GetLoadTemp(msfile, ant_index, spw[0])
    if tempAmb[ant_index] < 250: tempAmb[ant_index] += 273.15
    if tempHot[ant_index] < 300: tempHot[ant_index] += 273.15
    for spw_index in range(spwNum):
        AntSpwIndex = ant_index* spwNum + spw_index
        chNum = AmbSpecList[AntSpwIndex].shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        Trx, Tsky = np.zeros([2, chNum, len(offTime)]), np.zeros([2, chNum, len(offTime)])
        for pol_index in range(ppolNum):
            ambSpec, hotSpec = AmbSpecList[AntSpwIndex][pol_index], HotSpecList[AntSpwIndex][pol_index]
            for time_index in range(len(offTime)):
                #Psamb, Pshot = np.mean(ambSpec, axis=1), np.mean(hotSpec, axis=1)
                Psamb, Pshot = np.median(ambSpec, axis=1), np.median(hotSpec, axis=1)
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
            if 'SSOScanIndex' in locals(): TrxFlag[ant_index, spw_index, pol_index, SSOScanIndex] = 0.0
        #
        TrxList.append(Trx)
        TskyList.append(Tsky)
    #
#
#-------- Tau and TantN fitting
param = [0.0, 0.05]     # Initial Parameter
Tau0, TantN, Trms = np.zeros([antNum, spwNum, 2]), np.zeros([antNum, spwNum, 2]), np.zeros([antNum, spwNum, 2])
for ant_index in range(antNum):
    for spw_index in range(spwNum):
        for pol_index in range(2):
            if max(secZ[0]) - min(secZ[0]) > 0.5:
                fit = scipy.optimize.leastsq(residTskyTransfer, param, args=(tempAmb[ant_index]-Tatm_OFS, secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index]))
                TantN[ant_index, spw_index, pol_index] = fit[0][0]
                Tau0[ant_index, spw_index, pol_index]  = fit[0][1]
            else:
                param = [0.05]
                fit = scipy.optimize.leastsq(residTskyTransfer0, param, args=(tempAmb[ant_index]-Tatm_OFS, secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index]))
                Tau0[ant_index, spw_index, pol_index]  = fit[0][0]
                TantN[ant_index, spw_index, pol_index] = Tatm_OFS
        #
    #
    TantN[ant_index] = np.median(TantN[ant_index])* np.ones([spwNum, 2])
#
Tau0err, Tau0med = np.std(Tau0, axis=(0,2)), np.mean(Tau0, axis=(0,2))
#-------- Trx Transfer
for ant_index in range(antNum):
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
chAvgTrx = np.median(chAvgTrx, axis=3)
#-------- Tsys for scans
for scan_index in range(scanNum):
    OnTimeRange = timeXY[OnTimeIndex[scan_index]]
    ambTimeIndex = argmin(abs(ambTime - OnTimeRange[0]))
    offTimeIndex = argmin(abs(offTime - OnTimeRange[0]))
    chAvgTsys[:,:,:,scan_index] = chAvgTrx + chAvgTsky[:,:,:,offTimeIndex] 
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
np.save(prefix +  '-' + bandName + '.Trx.npy', TrxList) 
np.save(prefix +  '-' + bandName + '.Tsky.npy', TskyList) 
np.save(prefix +  '-' + bandName + '.TrxFlag.npy', TrxFlag) 
np.save(prefix +  '-' + bandName + '.Tau0.npy', Tau0) 
#-------- Antenna-dependent leakage noise
param = [0.0]
text_sd = ' TantN:'; logfile.write(text_sd); print text_sd,
for spw_index in range(spwNum): text_sd = ' PolX   SPW%02d  PolY           |' % (spw[spw_index]); logfile.write(text_sd); print text_sd,
#
logfile.write('\n'); print ' '
text_sd = ' ----:--------------------------------+-------------------------------+-------------------------------+-------------------------------+'; logfile.write(text_sd + '\n'); print text_sd
for ant_index in range(antNum):
    text_sd = antList[ant_index] + ' : '; logfile.write(text_sd); print text_sd,
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
np.save(prefix +  '-' + bandName + '.TantN.npy', TantN) 
#-------- Trx
text_sd = ' Trec: '; logfile.write(text_sd); print text_sd,
for spw_index in range(spwNum): text_sd = ' SPW%02d X        Y |' % (spw[spw_index]); logfile.write(text_sd); print text_sd,
logfile.write('\n'); print ' '
text_sd =  ' ----:--------------------+-------------------+-------------------+-------------------+'; logfile.write(text_sd + '\n'); print text_sd
for ant_index in range(antNum):
    text_sd =  antList[ant_index] + ' : '; logfile.write(text_sd); print text_sd,
    for spw_index in range(spwNum):
        for pol_index in range(2):
            text_sd = '%6.1f K' % (chAvgTrx[ant_index, spw_index, pol_index])
            logfile.write(text_sd); print text_sd,
        text_sd = '|'; logfile.write(text_sd); print text_sd,
    logfile.write('\n'); print ' '
print ' '
#-------- Tsys
text_sd = ' Tsys: '; logfile.write(text_sd); print text_sd,
for spw_index in range(spwNum): text_sd = ' SPW%02d X        Y |' % (spw[spw_index]); logfile.write(text_sd); print text_sd,
logfile.write('\n'); print ' '
text_sd =  ' ----:--------------------+-------------------+-------------------+-------------------+'; logfile.write(text_sd + '\n'); print text_sd
for ant_index in range(antNum):
    text_sd =  antList[ant_index] + ' : '; logfile.write(text_sd); print text_sd,
    for spw_index in range(spwNum):
        for pol_index in range(2):
            text_sd = '%6.1f K' % (chAvgTrx[ant_index, spw_index, pol_index] + np.median(chAvgTsky[ant_index, spw_index, pol_index]))
            logfile.write(text_sd); print text_sd,
        text_sd = '|'; logfile.write(text_sd); print text_sd,
    logfile.write('\n'); print ' '
#-------- Plot optical depth
if PLOTTAU: plotTau(prefix + '_' + bandName, antList, spw, secZ, (chAvgTsky.transpose(3,0,1,2) - TantN).transpose(1,2,3,0), np.median(tempAmb) - Tatm_OFS, Tau0med, TrxFlag, 2.0*np.median(chAvgTsky), PLOTFMT) 
if PLOTTSYS: plotTsys(prefix + '_' + bandName, antList, ambTime, spw, TrxList, TskyList, PLOTFMT)
