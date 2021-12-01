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
SunAngleTsysLimit = 5.0 # [deg] 
if 'PLOTTAU'  not in locals(): PLOTTAU  = False
if 'PLOTTSYS' not in locals(): PLOTTSYS = False
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
            if Pspec.shape[0] == 1: pPol = [0]
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
def LogTrx(antList, spwList, freqList, scanList, timeRef, Trx, TantN, logFile):
    antNum, spwNum, scanNum, polNum = len(antList), len(spwList), Trx[0].shape[3], Trx[0].shape[0]
    for scan_index in range(scanNum):
        text_sd =  'Scan %d : %s' % (scanList[scan_index], qa.time('%fs' % (timeRef[scan_index]), form='fits')[0]); logFile.write(text_sd + '\n'); print text_sd
        text_sd = 'Trx  : '
        for spw_index in range(spwNum): text_sd = text_sd + ' SPW%03d  %6.1f GHz |' % (spwList[spw_index], freqList[spw_index])
        logFile.write(text_sd + '\n'); print text_sd
        text_sd = ' Pol : '
        for spw_index in range(spwNum): text_sd = text_sd + '     X        Y     |'
        logFile.write(text_sd + '\n'); print text_sd
        text_sd =  ' ----:-'
        for spw_index in range(spwNum): text_sd = text_sd + '--------------------+'
        logFile.write(text_sd + '\n'); print text_sd
        for ant_index in range(antNum):
            text_sd =  antList[ant_index] + ' : '
            for spw_index in range(spwNum):
                for pol_index in range(polNum):
                    text_sd = text_sd + '%7.1f K ' % (np.median(Trx[spw_index], axis=1)[pol_index, ant_index, scan_index])
                text_sd = text_sd + '|'
            logFile.write(text_sd + '\n'); print text_sd
        #
    #
    #-------- Log mean and SD in Trx
    text_sd = 'mean : '
    for spw_index in range(spwNum): text_sd = text_sd + '                 SPW%03d  %6.1f GHz |' % (spwList[spw_index], freqList[spw_index])
    logFile.write(text_sd + '\n'); print text_sd
    text_sd = ' Pol : '
    for spw_index in range(spwNum): text_sd = text_sd + ' X mean (   sd)    Y mean (   sd)   |'
    logFile.write(text_sd + '\n'); print text_sd
    text_sd =  '-----:-'
    for spw_index in range(spwNum): text_sd = text_sd + '------------------------------------+'
    logFile.write(text_sd + '\n'); print text_sd
    for ant_index in range(antNum):
        text_sd =  antList[ant_index] + ' : '
        for spw_index in range(spwNum):
            for pol_index in range(polNum):
                text_sd = text_sd + '%7.1f (%5.1f) K ' % (np.median(Trx[spw_index], axis=(1,3))[pol_index, ant_index], np.std(np.median(Trx[spw_index], axis=1)[pol_index, ant_index]) )
            text_sd = text_sd + '|'
        logFile.write(text_sd + '\n'); print text_sd
    #
    #-------- Log TantN
    text_sd = 'TantN: '
    for spw_index in range(spwNum): text_sd = text_sd + ' SPW%03d   |' % (spwList[spw_index])
    logFile.write(text_sd + '\n'); print text_sd
    text_sd =  ' ----:-'
    for spw_index in range(spwNum): text_sd = text_sd + '----------+'
    logFile.write(text_sd + '\n'); print text_sd
    for ant_index in range(antNum):
        text_sd =  antList[ant_index] + ' : '
        for spw_index in range(spwNum):
            text_sd = text_sd + '%7.1f K ' % (np.median(TantN[spw_index][ant_index]))
            text_sd = text_sd + '|'
        logFile.write(text_sd + '\n'); print text_sd
    #
    text_sd =  ' ----:-'
    for spw_index in range(spwNum): text_sd = text_sd + '----------+'
    logFile.write(text_sd + '\n'); print text_sd
    return
#
#-------- Trx and Tsky
def TrxTskySpec(useAnt, tempAmb, tempHot, spwList, scanList, ambSpec, hotSpec, offSpec):
    TrxList, TskyList = [], []
    useAntNum, spwNum, scanNum, polNum =  len(useAnt), len(spwList), len(scanList), ambSpec[0].shape[0]
    scanFlag = np.ones([spwNum, polNum, useAntNum, scanNum]) 
    for spw_index in range(len(spwList)):
        chNum = ambSpec[spw_index* scanNum].shape[1]
        chRange = range(int(0.05*chNum), int(0.95*chNum)); chOut = sort(list(set(range(chNum)) - set(chRange))).tolist()
        #chRange = range(int(0.02*chNum), int(0.99*chNum)); chOut = sort(list(set(range(chNum)) - set(chRange))).tolist()
        TrxSpec = -np.ones([polNum, chNum, useAntNum, scanNum])
        TskySpec = -np.ones([polNum, chNum, useAntNum, scanNum])
        for pol_index in range(polNum):
            for ant_index in range(useAntNum):
                for scan_index in range(scanNum):
                    index = (ant_index* spwNum + spw_index)* scanNum + scan_index
                    Psamb, Pshot, Psoff = ambSpec[index][pol_index], hotSpec[index][pol_index], offSpec[index][pol_index]
                    if np.max(Pshot[chRange]/Psamb[chRange]) > tempHot[ant_index] / tempAmb[ant_index] : continue       # negative Trx
                    if np.min(Pshot[chRange]/Psamb[chRange]) < 1.01 : continue                                          # infinite Trx
                    TrxSpec[pol_index, chRange, ant_index, scan_index] = (tempHot[ant_index]* Psamb[chRange] - Pshot[chRange]* tempAmb[ant_index]) / (Pshot - Psamb)[chRange]
                    TskySpec[pol_index, chRange, ant_index, scan_index]= (Psoff[chRange]* (tempHot[ant_index] - tempAmb[ant_index]) + tempAmb[ant_index]* Pshot[chRange] - tempHot[ant_index]* Psamb[chRange]) / (Pshot - Psamb)[chRange]
                    TrxSpec[pol_index, chOut, ant_index, scan_index] = np.median(TrxSpec[pol_index, chRange, ant_index, scan_index])       # extraporate band-edge
                    TskySpec[pol_index, chOut, ant_index, scan_index] = np.median(TskySpec[pol_index, chRange, ant_index, scan_index])     # extraporate band-edge
                #
            #
        #
        #-------- Flag negative Trx
        chAvgTrx = np.mean(TrxSpec[:,chRange], axis=1)
        scanFlag[spw_index] = (np.sign(chAvgTrx - 10.0) + 1.0 )/2       # Flag (Trx < 10.0 K) out, scanFlag[spw, pol, ant, scan]
        chAvgTrx = (np.sum(scanFlag[spw_index]* chAvgTrx, axis=2) / np.sum(scanFlag[spw_index] + 1.0e-9, axis=2)).T # chAvgTrx[ant, pol]
        TskySpec = np.sum(TskySpec.transpose(1,0,2,3)* scanFlag[spw_index], axis=1) / np.sum(scanFlag[spw_index] + 1.0e-9, axis=0)
        TrxList = TrxList + [TrxSpec]
        TskyList = TskyList + [TskySpec]
    #
    return  TrxList, TskyList, scanFlag
#
#-------- Smooth time-variable Tau
def tauSMTH( timeSample, TauE ):
    if len(timeSample) > 5:
        SplineWeight = np.ones(len(timeSample) + 4)
        flagIndex = (np.where(abs(TauE - np.median(TauE))/np.std(TauE) > 3.0)[0] + 2).tolist()
        SplineWeight[flagIndex] = 0.01
        tempTime = np.append([timeSample[0]-500.0, timeSample[0]-300.0], np.append(timeSample, [timeSample[-1]+300.0, timeSample[-1]+500.0]))
        tempTauE = np.append([TauE[0], TauE[0]], np.append(TauE, [TauE[-1], TauE[-1]]))
        #smthTau = scipy.interpolate.splrep(tempTime, tempTauE, k=3, w=SplineWeight, t=tempTime[range(2, len(tempTime)-2, 3)] - 60.0 )
        smthTau = scipy.interpolate.splrep(tempTime, tempTauE, k=3, w=SplineWeight, t=tempTime[range(1, len(tempTime), 2)] - 60.0 )
    else:
        tempTime = np.arange(np.min(timeSample) - 3600.0,  np.max(timeSample) + 3600.0, 300.0)
        tempTauE = np.repeat(np.median(TauE), len(tempTime))
        smthTau = scipy.interpolate.splrep(tempTime, tempTauE, k=3)
    #
    return smthTau
#
#-------- Zenith opacity fitting
#Tau0, Tau0Excess, Tau0Coef, TantN = tau0SpecFit(tempAtm, atmsecZ, useAnt, atmspwLists[band_index], TskyList, scanFlag)
def tau0SpecFit(tempAtm, secZ, useAnt, spwList, TskyList, scanFlag):
    Tau0List, TantNList, Tau0Coef = [], [], []
    scanNum, useAntNum, spwNum = len(secZ), len(useAnt), len(spwList)
    Tau0Excess = np.zeros([spwNum, scanNum])
    #-------- Case1: Single atmCal scan
    if scanNum < 2:
        for spw_index in range(spwNum):
            chNum = TskyList[spw_index].shape[0]
            TantNList = TantNList + [np.zeros([useAntNum, chNum])]
            Tau0List  = Tau0List  + [ -np.log( (np.median(TskyList[spw_index], axis=(1,2)) - tempAtm) / (Tcmb - tempAtm) ) / secZ ]
            Tau0Coef = Tau0Coef + [np.zeros(2)]
        return Tau0List, Tau0Excess, Tau0Coef, TantNList
    #    
    #-------- Case2: Multiple atmCal scans, but insuffcient SecZ coverage
    print 'SD(secZ) = %f' % (np.std(secZ))
    if np.std(secZ) < 0.25:
        for spw_index in range(spwNum):
            scanWeight = np.sum(scanFlag[spw_index], axis=(0,1))
            chNum = TskyList[spw_index].shape[0]
            TantNList = TantNList + [np.zeros([useAntNum, chNum])]
            Tau0Med = np.zeros(chNum)
            for ch_index in range(chNum):
                param = [0.05]
                #-------- Fit for Tau0 (fixed TantN)
                fit = scipy.optimize.leastsq(residTskyTransfer0, param, args=(tempAtm, secZ, np.nanmedian(TskyList[spw_index][ch_index], axis=0), scanWeight / (np.var(TskyList[spw_index][ch_index], axis=0) + 1e-3) ))
                Tau0Med[ch_index]  = fit[0][0]
            #
            Tau0List  = Tau0List  + [Tau0Med]
            Tau0Excess[spw_index] = residTskyTransfer0([np.median(Tau0Med)], tempAtm, secZ, np.median(TskyList[spw_index], axis=(0,1)), scanWeight ) / (tempAtm - Tcmb)* np.exp(-np.median(Tau0Med)* secZ) / secZ / (scanWeight + 1e-3)
        #
        #-------- Tau0Excess dependent on secZ 
        for spw_index in range(spwNum):
            seczSum, seczVar, tauRes = np.sum(secZ),  secZ.dot(secZ), secZ.dot(Tau0Excess[spw_index])
            detTau = scanNum* seczVar - seczSum**2
            coef = np.array([[seczVar, -seczSum],[-seczSum, scanNum]]).dot( np.array([np.sum(Tau0Excess[spw_index]), secZ.dot(Tau0Excess[spw_index])])) / detTau
            Tau0Excess[spw_index] = Tau0Excess[spw_index] - coef[0] - coef[1]* secZ
            Tau0Coef = Tau0Coef + [coef]
        #
        return Tau0List, Tau0Excess, Tau0Coef, TantNList
    #
    #-------- Case3: Suffcient SecZ coverage
    for spw_index in range(spwNum):
        param = [0.0, 0.05] # Initial parameter [TantN, Tau0]
        chNum = TskyList[spw_index].shape[0]
        Tau0, TantN = param[1]* np.ones([useAntNum, chNum]), np.zeros([useAntNum, chNum])
        #-------- Fit for Tau0 (without TantN)
        for ant_index in range(useAntNum):
            scanWeight = scanFlag[spw_index, 0, ant_index] * scanFlag[spw_index, 1, ant_index]
            if len(np.where(scanWeight > 0)[0]) > 6:    # at least 6 points to fit
                for ch_index in range(chNum):
                    fit = scipy.optimize.leastsq(residTskyTransfer, param, args=(tempAtm, secZ, TskyList[spw_index][ch_index, ant_index], scanWeight))
                    TantN[ant_index, ch_index] = fit[0][0]
                    Tau0[ant_index, ch_index]  = fit[0][1]
                #
            #
        #
        Tau0Med = np.median(Tau0, axis=0)   # Tau0 is independent on antenna
        #-------- Fit for TantN (fixed Tau0)
        for ant_index in range(useAntNum):
            scanWeight = scanFlag[spw_index, 0, ant_index] * scanFlag[spw_index, 1, ant_index]
            if len(np.where(scanWeight > 0)[0]) > 1:
                for ch_index in range(chNum):
                    param = Tau0[ant_index, ch_index] 
                    fit = scipy.optimize.leastsq(residTskyTransfer2, param, args=(tempAtm, Tau0Med[ch_index], secZ, TskyList[spw_index][ch_index, ant_index], scanWeight / (np.var(TskyList[spw_index][ch_index], axis=0) + 1e-3) ))
                    TantN[ant_index, ch_index]  = fit[0][0]
                #
            #
        #
        TskyResid = np.median((TskyList[spw_index].transpose(2,1,0) - TantN), axis=1)
        #-------- Fit for Tau0 (fixed TantN)
        scanWeight = np.sum(scanFlag[spw_index], axis=(0,1))
        for ch_index in range(chNum):
            param = [Tau0Med[ch_index]]
            fit = scipy.optimize.leastsq(residTskyTransfer0, param, args=(tempAtm, secZ, TskyResid[:,ch_index], scanWeight / (np.var(TskyList[spw_index][ch_index], axis=0) + 1e-2)))
            Tau0Med[ch_index]  = fit[0][0]
        #
        Tau0Excess[spw_index] = residTskyTransfer0([np.median(Tau0Med)], tempAtm, secZ, np.median(TskyResid, axis=1), scanWeight ) / (tempAtm - Tcmb)* np.exp(-np.median(Tau0Med)* secZ) / secZ / (scanWeight + 1e-3)
        Tau0List  = Tau0List  + [Tau0Med]
        TantNList = TantNList + [TantN]
        #
    #
    #-------- Tau0Excess dependent on secZ 
    for spw_index in range(spwNum):
        seczSum, seczVar, tauRes = np.sum(secZ),  secZ.dot(secZ), secZ.dot(Tau0Excess[spw_index])
        detTau = scanNum* seczVar - seczSum**2
        coef = np.array([[seczVar, -seczSum],[-seczSum, scanNum]]).dot( np.array([np.sum(Tau0Excess[spw_index]), secZ.dot(Tau0Excess[spw_index])])) / detTau
        Tau0Excess[spw_index] = Tau0Excess[spw_index] - coef[0] - coef[1]* secZ
        Tau0Coef = Tau0Coef + [coef]
    #
    return Tau0List, Tau0Excess, Tau0Coef, TantNList
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
    #atmSPWs = list( set(GetBPcalSPWs(msfile)) & set(GetAtmSPWs(msfile)) ); atmSPWs.sort()
    atmSPWs = GetAtmSPWs(msfile); atmSPWs.sort()
atmBandNames = GetBandNames(msfile, atmSPWs); UniqBands = unique(atmBandNames).tolist()
if UniqBands == []: UniqBands = BandList
NumBands = len(UniqBands)
msmd.open(msfile)
atmspwLists, atmscanLists = [], []
for band_index in range(NumBands):
    bandAtmSPWs = np.array(atmSPWs)[indexList(np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()
    if atmBandNames == []: bandAtmSPWs = np.array(atmSPWs)
    atmspwLists = atmspwLists + [bandAtmSPWs]
    if 'spwFlag' in locals():
        flagIndex = indexList(np.array(spwFlag), np.array(atmspwLists[band_index]))
        for index in flagIndex: del atmspwLists[band_index][index]
    #
    atmscanList = list(set(msmd.scansforspw(atmspwLists[band_index][0]))& set(msmd.scansforintent("CALIBRATE_ATMOSPHERE*")))
    atmscanList.sort()
    atmscanLists = atmscanLists + [atmscanList]
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
    for band_index in range(NumBands):
        atmSPW = atmspwLists[band_index]
        timeXY, Pspec = GetPSpec(msfile, 0, atmSPW[0])
        timeNum, chNum = Pspec.shape[2], Pspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        chAvgPower = np.mean(Pspec[0][chRange], axis=0)
        onTimeIndex  = indexList(timeON, timeXY)
        onTime, onPower = timeXY[onTimeIndex], chAvgPower[onTimeIndex]
        hot_index, amb_index = np.where(onPower >  np.median(onPower))[0].tolist(), np.where(onPower <  np.median(onPower))[0].tolist()
        hotStart = [hot_index[0]] + np.array(hot_index)[np.where(np.diff(onTime[hot_index]) > 60.0)[0] + 1].tolist()
        ambStart = [amb_index[0]] + np.array(amb_index)[np.where(np.diff(onTime[amb_index]) > 60.0)[0] + 1].tolist()
        if len(hot_index) > len(hotStart):
            hotTimeIndex, ambTimeIndex = np.array(onTimeIndex)[list(set(hot_index) - set(hotStart))], np.array(onTimeIndex)[list(set(amb_index) - set(ambStart))]
        else:
            hotTimeIndex, ambTimeIndex = np.array(onTimeIndex)[hot_index], np.array(onTimeIndex)[amb_index]
        timeAMB = np.append(timeAMB, timeXY[ambTimeIndex])
        timeHOT = np.append(timeHOT, timeXY[hotTimeIndex])
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
Tau0Max = np.zeros(NumBands)
sourceList, polList = GetSourceList(msfile)
SunAngleSourceList = GetSunAngle(msfile)
for band_index in range(NumBands):
    tsysLog = open(prefix + '-' + UniqBands[band_index] + '-Tsys.log', 'w')
    #-------- Trx
    atmTimeRef, offSpec, ambSpec, hotSpec, scanList = scanAtmSpec(msfile, useAnt, atmscanLists[band_index], atmspwLists[band_index], timeOFF, timeON, timeAMB, timeHOT)
    atmscanLists[band_index] = scanList
    print 'atmCal scans = ' + `atmscanLists[band_index]`
    atmscanNum, spwNum = len(atmscanLists[band_index]), len(atmspwLists[band_index])
    TrxList, TskyList, scanFlag = TrxTskySpec(useAnt, tempAmb, tempHot, atmspwLists[band_index], atmscanLists[band_index], ambSpec, hotSpec, offSpec)
    #-------- Check sun angle 
    for scan_index in range(len(scanList)):
        sourceID = msmd.sourceidforfield(msmd.fieldsforscan(scanList[scan_index])[0])
        SunAngle = SunAngleSourceList[sourceID]
        if SunAngle < SunAngleTsysLimit: scanFlag[:,:,:,scan_index] *= 0.01
    #
    for spw_index in range(spwNum):
        np.save(prefix +  '-' + UniqBands[band_index] + '-SPW' + `atmspwLists[band_index][spw_index]` + '.Trx.npy', TrxList[spw_index])    # TxList[spw][pol,ch,ant,scan]
    #-------- Az and El position at atmCal and onsource scans
    AtmEL = np.ones([useAntNum, atmscanNum])
    for ant_index in range(useAntNum):
        azelTime_index = np.where( AntID == useAnt[ant_index] )[0].tolist()
        if len(azelTime_index) == 0: azelTime_index = np.where( AntID == useAnt[ant_index + 1] )[0].tolist()
        for scan_index in range(atmscanNum): AtmEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - atmTimeRef[scan_index]))]]
    #
    atmsecZ  = 1.0 / np.sin( np.median(AtmEL, axis=0) )
    #-------- Tsky and TantN
    Tau0, Tau0Excess, Tau0Coef, TantN = tau0SpecFit(tempAtm, atmsecZ, useAnt, atmspwLists[band_index], TskyList, scanFlag)
    SPWfreqList, freqList = [], []
    Tau0med = np.zeros(spwNum)
    for spw_index in range(spwNum):
        chNum, chWid, freq = GetChNum(msfile, atmspwLists[band_index][spw_index]); freqList = freqList + [freq*1.0e-9]; SPWfreqList = SPWfreqList + [np.median(freq)*1.0e-9]
        Tau0med[spw_index] = np.median(Tau0[spw_index])
    #
    Tau0Max[band_index] = np.max(Tau0med)
    #-------- Log Tau0 (mean opacity at zenith)
    LogTrx(antList[useAnt], atmspwLists[band_index], SPWfreqList, scanList, atmTimeRef, TrxList, TantN, tsysLog)
    text_sd = 'Tau0 :  '
    for spw_index in range(spwNum): text_sd = text_sd + ' %7.5f | ' % (Tau0med[spw_index])
    tsysLog.write(text_sd + '\n'); print text_sd
    tsysLog.close()
    #-------- Save to npy files
    np.save(prefix +  '-' + UniqBands[band_index] + '.TrxAnt.npy', antList[useAnt])     # antList[ant]
    np.save(prefix +  '-' + UniqBands[band_index] + '.atmTime.npy', atmTimeRef)         # [scan]
    np.save(prefix +  '-' + UniqBands[band_index] + '.TauE.npy', Tau0Excess)            # [spw][scan]
    for spw_index in range(spwNum):
        np.save(prefix +  '-' + UniqBands[band_index] + '-SPW' + `atmspwLists[band_index][spw_index]` + '.TrxFreq.npy', freqList[spw_index])    # freqList[spw]
        np.save(prefix +  '-' + UniqBands[band_index] + '-SPW' + `atmspwLists[band_index][spw_index]` + '.Trx.npy', TrxList[spw_index])         # [spw][ant, pol, ch]
        np.save(prefix +  '-' + UniqBands[band_index] + '-SPW' + `atmspwLists[band_index][spw_index]` + '.TantN.npy', TantN[spw_index])         # [spw][ant, ch]
        np.save(prefix +  '-' + UniqBands[band_index] + '-SPW' + `atmspwLists[band_index][spw_index]` + '.Tau0.npy', Tau0[spw_index])           # [spw][ch]
        np.save(prefix +  '-' + UniqBands[band_index] + '-SPW' + `atmspwLists[band_index][spw_index]` + '.Tau0C.npy', Tau0Coef[spw_index])      # [spw][2]
    #
    #---- Plots
    if not 'PLOTFMT' in locals():   PLOTFMT = 'pdf'
    if PLOTTAU:
        plotTauSpec(prefix + '_' + UniqBands[band_index], atmspwLists[band_index], freqList, Tau0) 
        plotTauFit(prefix + '_' + UniqBands[band_index], antList[useAnt], atmspwLists[band_index], atmsecZ, tempAtm, Tau0, TantN, TskyList, np.min(scanFlag, axis=1)) 
        if len(atmscanLists[band_index]) > 5: plotTau0E(prefix + '_' + UniqBands[band_index], atmTimeRef, atmspwLists[band_index], Tau0, Tau0Excess, np.min(scanFlag, axis=(1,2))) 
    if PLOTTSYS: plotTsys(prefix + '_' + UniqBands[band_index], antList[useAnt], atmspwLists[band_index], freqList, atmTimeRef, TrxList, TskyList)
#
#-------- Plot optical depth
msmd.close()
msmd.done()
