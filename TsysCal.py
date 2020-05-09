# TsysCal.py 
# inputs
#    SCR_DIR    : (character) path to the script (e.g. '/home/skameno/ALMA_SV/' )
#    PLOTTAU    : (boolean)   plot optical depths (True or False)
#    PLOTTSYS   : (boolean)   plot Trx and Tsys spectra (True or False)
#    prefix     : (character) UID name (e.g. 'uid___A002_Xc02418_X3f67' )
#
# outputs
#  chAvgTsys[band* scan* ant* spw* pol] : List of channel-averaged system noise temperature
#  TsysSpec[band* scan* ant* spw* pol][ch] : List of Tsys spectrum
#  TsysFlag[ant, spw, pol, scan] : 0 = invalid, 1=valid
#  Tau0med[spw] : median-value of the zenith optical depth
#  onTau[spw, scan] : on-source optical depth
#
#  They include all of antennas (even if flagged) in MS order
#
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
#-------- Get atmCal scans
def scanAtmSpec(msfile, useAnt, scanList, spwList, timeOFF=0, timeON=0, timeAMB=0, timeHOT=0):
    timeList, offSpecList, ambSpecList, hotSpecList = [], [], [], []
    antNum, scanNum, spwNum = len(useAnt), len(scanList), len(spwList)
    scanTimeList = []
    scanFlag = range(scanNum)
    for scan_index in range(scanNum):
        scanID = scanList[scan_index]
        interval, scanTimeRec = GetTimerecord(msfile, 0, 0, spwList[0], scanID)
        offTime, ambTime, hotTime = sort(list(set(scanTimeRec) & set(timeOFF))), sort(list(set(scanTimeRec) & set(timeAMB))), sort(list(set(scanTimeRec) & set(timeHOT)))
        if (len(offTime)* len(ambTime)* len(hotTime) == 0):
            scanFlag.remove(scan_index)
        else:
            scanTimeList = scanTimeList + [scanTimeRec]
        #
    #
    scanList = np.array(scanList)[scanFlag].tolist()
    scanNum = len(scanTimeList)
    for ant_index in range(antNum):
        progress = (1.0* ant_index + 1.0) / antNum
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
        for spwID in spwList:
            timeXY, Pspec = GetPSpec(msfile, useAnt[ant_index], spwID)
            chNum = Pspec.shape[1]
            pPol = [0,1]
            if Pspec.shape[0] == 4: pPol = [0,3]
            for scan_index in range(scanNum):
                scanID = scanList[scan_index]
                scanTimeRec = scanTimeList[scan_index]
                offTime, ambTime, hotTime = sort(list(set(scanTimeRec) & set(timeOFF))), sort(list(set(scanTimeRec) & set(timeAMB))), sort(list(set(scanTimeRec) & set(timeHOT)))
                offTimeIndex, ambTimeIndex, hotTimeIndex = indexList(offTime, timeXY),  indexList(ambTime, timeXY),  indexList(hotTime, timeXY)
                if len(offTimeIndex) * len(ambTimeIndex) * len(hotTimeIndex) == 0: continue   # Unuseable atmCal scan
                offTimeIndex, ambTimeIndex, hotTimeIndex = indexList(offTime, timeXY)[-1],  indexList(ambTime, timeXY)[-1],  indexList(hotTime, timeXY)[-1]
                if ((ant_index == 0) & (spwID == spwList[0]) & (len(ambTime) > 0)): timeList = timeList + [offTime[-1]] # Record off-position time
                #
                if ((ant_index == 0) & (spwID == spwList[0]) & (len(ambTime) == 0)):    # No available ambient load data
                    chRange = range(int(0.05*chNum), int(0.95*chNum)); chAvgPower = np.mean(Pspec[0][chRange], axis=0)
                    offTimeIndex = indexList(timeOFF, timeXY)
                    ambhotTimeIndex = indexList(timeON, timeXY)
                    ambhotThresh = 0.5*(max(chAvgPower[ambhotTimeIndex]) + min(chAvgPower[ambhotTimeIndex]))
                    hotTimeIndex = np.array(ambhotTimeIndex)[np.where( chAvgPower[ambhotTimeIndex] > ambhotThresh )[0].tolist()].tolist()
                    ambTimeIndex = np.array(ambhotTimeIndex)[np.where( chAvgPower[ambhotTimeIndex] < ambhotThresh )[0].tolist()].tolist()
                    timeList = timeList + [np.median(timeXY[offTimeIndex])]
                #
                if len(ambTime) > 0 :
                    offSpecList = offSpecList + [Pspec[pPol][:,:,offTimeIndex]]
                    ambSpecList = ambSpecList + [Pspec[pPol][:,:,ambTimeIndex]]
                    hotSpecList = hotSpecList + [Pspec[pPol][:,:,hotTimeIndex]]
                    #index += 1
                else:
                    offSpecList = offSpecList + [np.median(Pspec[pPol][:,:,offTimeIndex], axis=2)]
                    ambSpecList = ambSpecList + [np.median(Pspec[pPol][:,:,ambTimeIndex], axis=2)]
                    hotSpecList = hotSpecList + [np.median(Pspec[pPol][:,:,hotTimeIndex], axis=2)]
                #
            #
        #
    #            
    sys.stderr.write('\n'); sys.stderr.flush()
    return np.array(timeList), offSpecList, ambSpecList, hotSpecList, scanList
#
#-------- Log Trx
def LogTrx(antList, spwList, freqList, scanList, timeRef, Trx, logFile):
    antNum, spwNum, scanNum = len(antList), len(spwList), Trx[0].shape[3]
    for scan_index in range(scanNum):
        text_sd =  'Scan %d : %s' % (scanList[scan_index], qa.time('%fs' % (timeRef[scan_index]), form='fits')[0]); logFile.write(text_sd + '\n'); print text_sd
        text_sd = 'Trx  : '
        for spw_index in range(spwNum): text_sd = text_sd + ' SPW%02d  %6.1f GHz |' % (spwList[spw_index], freqList[spw_index])
        logFile.write(text_sd + '\n'); print text_sd
        text_sd = ' Pol : '
        for spw_index in range(spwNum): text_sd = text_sd + '     X        Y    |'
        logFile.write(text_sd + '\n'); print text_sd
        text_sd =  ' ----:--------------------+-------------------+-------------------+-------------------+'; logFile.write(text_sd + '\n'); print text_sd
        for ant_index in range(antNum):
            text_sd =  antList[ant_index] + ' : '; logFile.write(text_sd); print text_sd,
            for spw_index in range(spwNum):
                for pol_index in range(2):
                    text_sd = '%6.1f K' % (np.median(Trx[spw_index], axis=1)[pol_index, ant_index, scan_index])
                    logFile.write(text_sd); print text_sd,
                text_sd = '|'; logFile.write(text_sd); print text_sd,
            logFile.write('\n'); print ' '
        print ' '
    return
#
#-------- Trx and Tsky
def TrxTskySpec(useAnt, tempAmb, tempHot, spwList, scanList, ambSpec, hotSpec, offSpec):
    TrxList, TskyList = [], []
    useAntNum, spwNum, scanNum =  len(useAnt), len(spwList), len(scanList)
    outLierFlag = np.zeros([spwNum, 2, useAntNum]) 
    for spw_index in range(len(spwList)):
        chNum = ambSpec[spw_index* scanNum].shape[1]
        chRange = range(int(0.02*chNum), int(0.99*chNum)); chOut = sort(list(set(range(chNum)) - set(chRange))).tolist()
        TrxSpec = np.zeros([2, chNum, useAntNum, scanNum])
        TskySpec = np.zeros([2, chNum, useAntNum, scanNum])
        for pol_index in range(2):
            for ant_index in range(useAntNum):
                for scan_index in range(scanNum):
                    index = (ant_index* spwNum + spw_index)* scanNum + scan_index
                    Psamb, Pshot, Psoff = ambSpec[index][pol_index], hotSpec[index][pol_index], offSpec[index][pol_index]
                    TrxSpec[pol_index, :, ant_index, scan_index] = (tempHot[ant_index]* Psamb - Pshot* tempAmb[ant_index]) / (Pshot - Psamb)
                    TskySpec[pol_index, :, ant_index, scan_index] = (Psoff* (tempHot[ant_index] - tempAmb[ant_index]) + tempAmb[ant_index]* Pshot - tempHot[ant_index]* Psamb) / (Pshot - Psamb)
                    TrxSpec[pol_index, chOut, ant_index, scan_index] = np.median(TrxSpec[pol_index, chRange, ant_index, scan_index])       # extraporate band-edge
                    TskySpec[pol_index, chOut, ant_index, scan_index] = np.median(TskySpec[pol_index, chRange, ant_index, scan_index])     # extraporate band-edge
                #
            #
        #
        # TrxSpec = np.median(TrxSpec, axis=3).transpose(2,0,1)   # Trx[ant, pol, ch] : Trx is time independent
        chAvgTrx = np.median(TrxSpec[:, chRange], axis=(1, 3)).T    # chAvgTrx[ant, pol]
        #-------- Trx transfer to outlier antennas
        for pol_index in range(2):
            medTrx = np.median(chAvgTrx[:,pol_index])
            flagIndex = np.where(abs(chAvgTrx[:,pol_index] - medTrx) > 0.8* medTrx )[0].tolist()
            if len(flagIndex) > 0:
                outLierFlag[spw_index, pol_index, flagIndex] = 1.0
                #TrxSpec[flagIndex, pol_index] = np.median(TrxSpec[:,pol_index], axis=0)
                for scan_index in range(scanNum):
                    TskySpec[pol_index, :, flagIndex, scan_index] = np.median(TskySpec[pol_index, :, :, scan_index], axis=1)
                #
            #
        #
        TrxList = TrxList + [TrxSpec]
        TskyList = TskyList + [np.mean(TskySpec, axis=0)]
    #
    return  TrxList, TskyList, outLierFlag
#
#-------- Zenith opacity fitting
def tau0SpecFit(tempAtm, secZ, useAnt, spwList, TskyList):
    Tau0List, TantNList = [], []
    scanNum, useAntNum, spwNum = len(secZ), len(useAnt), len(spwList)
    Tau0Excess = np.zeros([spwNum, scanNum])
    if (np.max(secZ) - np.min(secZ)) < 0.5:
        for spw_index in range(spwNum):
            chNum = TskyList[spw_index].shape[0]
            TantNList = TantNList + [np.zeros([useAntNum, chNum])]
            Tau0Med = np.zeros(chNum)
            for ch_index in range(chNum):
                param = [0.05]
                #-------- Fit for Tau0 (fixed TantN)
                fit = scipy.optimize.leastsq(residTskyTransfer0, param, args=(tempAtm, secZ, np.nanmedian(TskyList[spw_index][ch_index], axis=0), np.ones(scanNum)))
                Tau0Med[ch_index]  = fit[0][0]
            #
            Tau0List  = Tau0List  + [Tau0Med]
            Tau0Excess[spw_index] = residTskyTransfer0([np.median(Tau0Med)], tempAtm, secZ, np.median(TskyList[spw_index], axis=(0,1)), np.ones(scanNum) ) / (tempAtm - Tcmb)* np.exp(-np.median(Tau0Med)* secZ) / secZ
        #
        return Tau0List, Tau0Excess, TantNList
    #
    for spw_index in range(spwNum):
        param = [0.05] # Initial parameter [Tau0]
        chNum = TskyList[spw_index].shape[0]
        Tau0, TantN = np.zeros([useAntNum, chNum]), np.zeros([useAntNum, chNum])
        for ant_index in range(useAntNum):
            for ch_index in range(chNum):
                #-------- Fit for Tau0 (without TantN)
                fit = scipy.optimize.leastsq(residTskyTransfer0, param, args=(tempAtm, secZ, TskyList[spw_index][ch_index, ant_index], np.ones(scanNum)))
                Tau0[ant_index, ch_index]  = fit[0][0]
            #
        #
        Tau0Med = np.median(Tau0, axis=0)   # Tau0 is independent on antenna
        param = [0.0]
        for ant_index in range(useAntNum):
            for ch_index in range(chNum):
                #-------- Fit for TantN (fixed Tau0)
                fit = scipy.optimize.leastsq(residTskyTransfer2, param, args=(tempAtm, Tau0Med[ch_index], secZ, TskyList[spw_index][ch_index, ant_index], np.ones(scanNum)))
                TantN[ant_index, ch_index]  = fit[0][0]
            #
        #
        TskyResid = np.median((TskyList[spw_index].transpose(2,1,0) - TantN), axis=1)
        for ch_index in range(chNum):
            param = [Tau0Med[ch_index]]
            #-------- Fit for Tau0 (fixed TantN)
            fit = scipy.optimize.leastsq(residTskyTransfer0, param, args=(tempAtm, secZ, TskyResid[:,ch_index], np.ones(scanNum)))
            Tau0Med[ch_index]  = fit[0][0]
        #
        Tau0Excess[spw_index] = residTskyTransfer0([np.median(Tau0Med)], tempAtm, secZ, np.median(TskyResid, axis=1), np.ones(scanNum) ) / (tempAtm - Tcmb)* np.exp(-np.median(Tau0Med)* secZ) / secZ
        Tau0List  = Tau0List  + [Tau0Med]
        TantNList = TantNList + [TantN]
        #
    #
    #
    return Tau0List, Tau0Excess, TantNList
#
#-------- Check MS file
msfile = wd + prefix + '.ms'
tempAtm = GetTemp(msfile)
if tempAtm != tempAtm: tempAtm = 270.0; print 'Cannot get ambient-load temperature ... employ 270.0 K, instead.'
antList = GetAntName(msfile)
antNum = len(antList)
if 'flagAnt' not in locals(): flagAnt = np.ones(antNum)
if 'antFlag' in locals():
    index =  indexList(np.array(antFlag), antList)
    if len(index) > 0: flagAnt[index] = 0.0
useAnt = np.where(flagAnt == 1.0)[0].tolist(); useAntNum = len(useAnt)
#-------- Check SPWs
print '---Checking spectral windows and scans with atmCal for ' + prefix
if 'atmSPWs' not in locals():
    atmSPWs = GetAtmSPWs(msfile)
atmBandNames = GetBandNames(msfile); UniqBands = unique(atmBandNames).tolist(); NumBands = len(UniqBands)
msmd.open(msfile)
atmspwLists, atmscanLists = [], []
for band_index in range(NumBands):
    bandAtmSPWs = np.array(atmSPWs)[indexList(np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()
    #bandAtmSPWs = np.array(atmSPWs)
    atmspwLists = atmspwLists + [bandAtmSPWs]
    if 'spwFlag' in locals():
        flagIndex = indexList(np.array(spwFlag), np.array(atmspwLists[band_index]))
        for index in flagIndex: del atmspwLists[band_index][index]
    #
    atmscanList = list(set(msmd.scansforspw(atmspwLists[band_index][0]))& set(msmd.scansforintent("CALIBRATE_ATMOSPHERE*"))); atmscanList.sort(); atmscanLists = atmscanLists + [atmscanList]
    print ' ',
    print UniqBands[band_index] + ': atmSPW=' + `atmspwLists[band_index]`
#
# atmSPWs[band] : SPWs used in atmCal scans
# bpSPWs[band]  : SPWs used in bandpass scan (i.e. SPWs for OBS_TARGET)
# atmscanLists[band] : atmCal scans
# bpscanLists[band]  : scans on source
#
print '---Checking time for ambient and hot load'
timeOFF, timeON, timeAMB, timeHOT = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#ON_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT")
if len(timeAMB) == 0:
    timeXY, Pspec = GetPSpec(msfile, 0, atmSPWs[0])
    timeNum, chNum = Pspec.shape[2], Pspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    chAvgPower = np.mean(Pspec[0][chRange], axis=0)
    onTimeIndex  = indexList(timeON, timeXY)
    onTime, onPower = timeXY[onTimeIndex], chAvgPower[onTimeIndex]
    hot_index, amb_index = np.where(onPower >  np.median(onPower))[0].tolist(), np.where(onPower <  np.median(onPower))[0].tolist()
    hotStart = [hot_index[0]] + np.array(hot_index)[np.where(np.diff(onTime[hot_index]) > 60.0)[0] + 1].tolist()
    ambStart = [amb_index[0]] + np.array(amb_index)[np.where(np.diff(onTime[amb_index]) > 60.0)[0] + 1].tolist()
    hotTimeIndex, ambTimeIndex = np.array(onTimeIndex)[list(set(hot_index) - set(hotStart))], np.array(onTimeIndex)[list(set(amb_index) - set(ambStart))]
    timeAMB, timeHOT = timeXY[ambTimeIndex], timeXY[hotTimeIndex]
#
azelTime, AntID, AZ, EL = GetAzEl(msfile)
# timeOFF : mjd of CALIBRATE_ATMOSPHERE#OFF_SOURCE
# timeON  : mjd of CALIBRATE_ATMOSPHERE#ON_SOURCE (becore Cycle 3, ambient + hot loads
# timeAMB : mjd of CALIBRATE_ATMOSPHERE#AMBIENT (after Cycle 3)
# timeHOT : mjd of CALIBRATE_ATMOSPHERE#HOT (after Cycle 3)
#
#-------- Get Load Temperature
tempAmb, tempHot  = np.zeros([useAntNum]), np.zeros([useAntNum])
for ant_index in range(useAntNum):
    tempAmb[ant_index], tempHot[ant_index] = GetLoadTemp(msfile, useAnt[ant_index], atmspwLists[0][0])
    if tempAmb[ant_index] < 240: tempAmb[ant_index] += 273.15       # Old MS describes the load temperature in Celsius
    if tempHot[ant_index] < 240: tempHot[ant_index] += 273.15       #
#

#-------- Trx, TantN, and Tau0
for band_index in range(NumBands):
    tsysLog = open(prefix + '-' + UniqBands[band_index] + '-Tsys.log', 'w')
    #-------- Trx
    atmTimeRef, offSpec, ambSpec, hotSpec, scanList = scanAtmSpec(msfile, useAnt, atmscanLists[band_index], atmspwLists[band_index], timeOFF, timeON, timeAMB, timeHOT)
    atmscanLists[band_index] = scanList
    atmscanNum, spwNum = len(atmscanLists[band_index]), len(atmspwLists[band_index])
    TrxList, TskyList, outLierFlag = TrxTskySpec(useAnt, tempAmb, tempHot, atmspwLists[band_index], atmscanLists[band_index], ambSpec, hotSpec, offSpec)
    for spw_index in range(spwNum):
        np.save(prefix +  '-' + UniqBands[band_index] + '-SPW' + `atmspwLists[band_index][spw_index]` + '.Trx.npy', TrxList[spw_index])    # TxList[spw][pol,ch,ant,scan]
    #-------- Az and El position at atmCal and onsource scans
    AtmEL = np.ones([useAntNum, atmscanNum])
    for ant_index in range(useAntNum):
        azelTime_index = np.where( AntID == useAnt[ant_index] )[0].tolist()
        if len(azelTime_index) == 0: azelTime_index = np.where( AntID == 0 )[0].tolist()
        for scan_index in range(atmscanNum): AtmEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - atmTimeRef[scan_index]))]]
    #
    atmsecZ  = 1.0 / np.sin( np.median(AtmEL, axis=0) )
    #-------- Tsky and TantN
    Tau0, Tau0Excess, TantN = tau0SpecFit(tempAtm - 5.0, atmsecZ, useAnt, atmspwLists[band_index], TskyList)
    for spw_index in range(spwNum):
        TrxList[spw_index] = (TrxList[spw_index].transpose(0,3,2,1) + TantN[spw_index]).transpose(0,3,2,1)
    #
    SPWfreqList, freqList = [], []
    for spw_index in range(spwNum):
        chNum, chWid, freq = GetChNum(msfile, atmspwLists[band_index][spw_index]); freqList = freqList + [freq*1.0e-9]; SPWfreqList = SPWfreqList + [np.median(freq)*1.0e-9]
    #
    LogTrx(antList[useAnt], atmspwLists[band_index], SPWfreqList, scanList, atmTimeRef, TrxList, tsysLog)
    for spw_index in range(spwNum):
        text_sd = 'SPW=%d : Tau(zenith) = %6.4f' % (atmspwLists[band_index][spw_index], np.median(Tau0[spw_index]))
        tsysLog.write(text_sd + '\n'); print text_sd
    #
    tsysLog.close()
    #-------- Save to npy files
    np.save(prefix +  '-' + UniqBands[band_index] + '.TrxAnt.npy', antList[useAnt])     # antList[ant]
    np.save(prefix +  '-' + UniqBands[band_index] + '.atmTime.npy', atmTimeRef)         # [scan]
    np.save(prefix +  '-' + UniqBands[band_index] + '.TauE.npy', Tau0Excess)            # [spw][scan]
    for spw_index in range(spwNum):
        np.save(prefix +  '-' + UniqBands[band_index] + '-SPW' + `atmspwLists[band_index][spw_index]` + '.TrxFreq.npy', freqList[spw_index])    # freqList[spw]
        np.save(prefix +  '-' + UniqBands[band_index] + '-SPW' + `atmspwLists[band_index][spw_index]` + '.Trx.npy', TrxList[spw_index])    # [spw][ant, pol, ch]
        np.save(prefix +  '-' + UniqBands[band_index] + '-SPW' + `atmspwLists[band_index][spw_index]` + '.Tau0.npy', Tau0[spw_index])      # [spw][ch]
    #
    #---- Plots
    if not 'PLOTFMT' in locals():   PLOTFMT = 'pdf'
    if PLOTTAU: plotTau(prefix + '_' + UniqBands[band_index], atmspwLists[band_index], freqList, Tau0) 
    if PLOTTSYS: plotTsys(prefix + '_' + UniqBands[band_index], antList[useAnt], atmspwLists[band_index], freqList, atmTimeRef, TrxList, TskyList)
#
#-------- Plot optical depth
msmd.close()
msmd.done()
