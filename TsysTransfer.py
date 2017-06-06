#-------- Load autocorrelation power spectra
print '---Loading autocorr power spectra'
timeXY = np.load(prefix +  '-' + '.autocorrTime.npy') 
OffSpecList = np.load(prefix +  '-' + '.OffSpecList.npy') 
AmbSpecList = np.load(prefix +  '-' + '.AmbSpecList.npy') 
HotSpecList = np.load(prefix +  '-' + '.HotSpecList.npy') 
OnSpecList = np.load(prefix +  '-' + '.OnSpecList.npy') 
#-------- Load Az El position
azelTime, AntID, AZ, EL = GetAzEl(msfile)
OnAZ, OnEL, OffEL = np.ones([UseAntNum, scanNum]), np.ones([UseAntNum, scanNum]), np.ones(UseAntNum)
for ant_index in range(UseAntNum):
    azelTime_index = np.where( AntID == antMap[ant_index] )[0].tolist()
    for scan_index in range(scanNum):
        refTime = np.median(msmd.timesforscan(onsourceScans[scan_index]))
        OnAZ[ant_index, scan_index] = AZ[azelTime_index[argmin(abs(azelTime[azelTime_index] - refTime))]]
        OnEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - refTime))]]
    #
    #for time_index in range(len(offTimeIndex)):
    OffEL[ant_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - np.median(offTime[offTimeIndex])))]]
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
            for time_index in range(len(offTime)):
                Psamb = np.mean(ambSpec, axis=1)
                Pshot = np.mean(hotSpec, axis=1)
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
            #fit = scipy.optimize.leastsq(residTskyTransfer, param, args=(tempAmb[ant_index]-Tatm_OFS, secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index]))
            #TantN[ant_index, spw_index, pol_index] = fit[0][0]
            #Tau0[ant_index, spw_index, pol_index]  = fit[0][1]
            TantN[ant_index, spw_index, pol_index] = Tatm_OFS
            Tau0[ant_index, spw_index, pol_index]  = -np.log(1.0 - np.median(chAvgTsky[ant_index, spw_index, pol_index]) / (tempAmb[ant_index]-Tatm_OFS))* np.sin( OnEL[ant_index,0] )
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
            if(len(flgIndex) > 0):
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
    chAvgTsys[:,:,:,scan_index] = chAvgTrx[:,:,:,ambTimeIndex] + (np.median(tempAmb) - Tatm_OFS)* (1.0 - np.exp(-Tau0 / np.sin(np.median(OnEL[:,scan_index])))) + TantN
    TsysFlag[:,:,:,scan_index] = TrxFlag[:,:,:,ambTimeIndex]
#
#-------- Tau on source
for scan_index in range(scanNum):
    #OnTimeRange = timeXY[OnTimeIndex[scan_index]]
    #offTimeIndex = argmin(abs(offTime - OnTimeRange[0]))
    #tempTau = -np.log((chAvgTsky[:,:,:,offTimeIndex] - TantN - np.median(tempAmb) + Tatm_OFS) / (2.718 - np.median(tempAmb) + Tatm_OFS))
    onTau[:,scan_index] = np.median(np.median(Tau0, axis=2), axis=0) / np.sin( np.median(OnEL[:,scan_index]) )
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
#if PLOTTAU: plotTau(prefix + '_' + UniqBands[band_index], antList[antMap], spw, secZ, (chAvgTsky.transpose(3,0,1,2) - TantN).transpose(1,2,3,0), np.median(tempAmb) - Tatm_OFS, Tau0med, TrxFlag, 2.0*np.median(chAvgTsky), PLOTFMT) 
#if PLOTTSYS: plotTsys(prefix + '_' + UniqBands[band_index], antList[antMap], ambTime, spw, TrxList, TskyList, PLOTFMT)
