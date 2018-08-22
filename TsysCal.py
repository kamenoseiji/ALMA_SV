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
def LogTrx(antList, spwList, Trx, text_sd, logFile):
    antNum, spwNum = len(antList), len(spwList)
    for spw_index in range(spwNum): text_sd = text_sd + '  SPW%02d X        Y |' % (spwList[spw_index])
    logFile.write(text_sd + '\n'); print text_sd
    text_sd =  ' ----:--------------------+-------------------+-------------------+-------------------+'; logFile.write(text_sd + '\n'); print text_sd
    for ant_index in range(antNum):
        text_sd =  antList[ant_index] + ' : '; logFile.write(text_sd); print text_sd,
        for spw_index in range(spwNum):
            for pol_index in range(2):
                text_sd = '%6.1f K' % (Trx[ant_index, spw_index, pol_index])
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
        chRange = range(int(0.05*chNum), int(0.95*chNum)); chOut = sort(list(set(range(chNum)) - set(chRange))).tolist()
        TrxSpec = np.zeros([2, chNum, useAntNum, scanNum])
        TskySpec = np.zeros([2, chNum, useAntNum, scanNum])
        for pol_index in range(2):
            for ant_index in range(useAntNum):
                for scan_index in range(scanNum):
                    index = (ant_index* spwNum + spw_index)* scanNum + scan_index
                    Psamb, Pshot, Psoff = ambSpec[index][pol_index], hotSpec[index][pol_index], offSpec[index][pol_index]
                    TrxSpec[pol_index, :, ant_index, scan_index] = (tempHot[ant_index]* Psamb - Pshot* tempAmb[ant_index]) / (Pshot - Psamb)
                    TskySpec[pol_index, :, ant_index, scan_index] = (Psoff* (tempHot[ant_index] - tempAmb[ant_index]) + tempAmb[ant_index]* Pshot - tempHot[ant_index]* Psamb) / (Pshot - Psamb)
                    TrxSpec[pol_index, chOut, ant_index, scan_index] = np.median(TrxSpec[pol_index, chRange, ant_index, scan_index])
                    TskySpec[pol_index, chOut, ant_index, scan_index] = np.median(TskySpec[pol_index, chRange, ant_index, scan_index])
                #
            #
        #
        TrxSpec = np.median(TrxSpec, axis=3).transpose(2,0,1)   # Trx[ant, pol, ch] : Trx is time independent
        chAvgTrx = np.median(TrxSpec[:,:,chRange], axis=2)        # chAvgTrx[ant, pol]
        #-------- Trx transfer to outlier antennas
        for pol_index in range(2):
            medTrx = np.median(chAvgTrx[:,pol_index])
            flagIndex = np.where(abs(chAvgTrx[:,pol_index] - medTrx) > 0.8* medTrx )[0].tolist()
            if len(flagIndex) > 0:
                outLierFlag[spw_index, pol_index, flagIndex] = 1.0
                TrxSpec[flagIndex, pol_index] = np.median(TrxSpec[:,pol_index], axis=0)
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
                fit = scipy.optimize.leastsq(residTskyTransfer0, param, args=(tempAtm, secZ, np.median(TskyList[spw_index], axis=1)[ch_index], np.ones(scanNum)))
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
msfile = prefix + '.ms'
tempAtm = GetTemp(msfile)
antList = GetAntName(msfile)
antNum = len(antList)
if 'flagAnt' not in locals(): flagAnt = np.ones(antNum)
if 'antFlag' in locals(): flagAnt[indexList(antFlag, antList)] = 0.0
useAnt = np.where(flagAnt == 1.0)[0].tolist(); useAntNum = len(useAnt)
#-------- Check SPWs
msmd.open(msfile)
print '---Checking spectral windows and scans with atmCal for ' + prefix
atmSPWs = list(set(msmd.tdmspws()) & set(msmd.spwsforintent("CALIBRATE_ATMOSPHERE*"))); atmSPWs.sort()
bpSPWs  = msmd.spwsforintent("CALIBRATE_BANDPASS*").tolist(); bpSPWs.sort()
if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_DELAY*").tolist(); bpSPWs.sort()
atmspwNames, bpspwNames = msmd.namesforspws(atmSPWs), msmd.namesforspws(bpSPWs)
bpSPWs = np.array(bpSPWs)[indexList(np.array(atmspwNames), np.array(bpspwNames))].tolist(); bpspwNames = msmd.namesforspws(bpSPWs)
atmBandNames, atmPattern = [], r'RB_..'
for spwName in atmspwNames : atmBandNames = atmBandNames + re.findall(atmPattern, spwName)
UniqBands = unique(atmBandNames).tolist(); NumBands = len(UniqBands)
atmspwLists, bpspwLists, atmscanLists, bpscanLists = [], [], [], []
for band_index in range(NumBands):
    bandAtmSPWs = np.array(atmSPWs)[indexList(np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()
    bandBpSPWs  = np.array(bpSPWs)[indexList( np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()
    atmspwLists = atmspwLists + [bandAtmSPWs]
    bpspwLists  = bpspwLists  + [bandBpSPWs]
    atmscanList = list(set(msmd.scansforspw(atmspwLists[band_index][0]))& set(msmd.scansforintent("CALIBRATE_ATMOSPHERE*"))); atmscanList.sort(); atmscanLists = atmscanLists + [atmscanList]
    bpscanList  = list(set(msmd.scansforspw(bpspwLists[band_index][0])) & set(msmd.scansforintent("*#ON_SOURCE"))); bpscanList.sort(); bpscanLists = bpscanLists + [bpscanList]
    print ' ',
    print UniqBands[band_index] + ': atmSPW=' + `atmspwLists[band_index]` + ' bpSPW=' + `bpspwLists[band_index]`
    #
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
    offTimeIndex = indexList(timeOFF, timeXY)
    gapList = np.where(np.diff(chAvgPower) > 5.0* np.std(chAvgPower[offTimeIndex]))[0]
    gapList = np.append(gapList, [max(gapList) + 6])
    hotTimeIndex, ambTimeIndex = [], []
    for gap_index in range(0, len(gapList)-1, 2):
        ambTimeIndex = ambTimeIndex + range((gapList[gap_index] + 2), (gapList[gap_index + 1]-1))
        hotTimeIndex = hotTimeIndex + range((gapList[gap_index+1] + 2), (gapList[gap_index + 2]-1))
    #
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
    atmscanNum, scanNum, spwNum = len(atmscanLists[band_index]), len(bpscanLists[band_index]), len(atmspwLists[band_index])
    TrxList, TskyList, outLierFlag = TrxTskySpec(useAnt, tempAmb, tempHot, atmspwLists[band_index], atmscanLists[band_index], ambSpec, hotSpec, offSpec)
    np.save(prefix +  '-' + UniqBands[band_index] + '.Trx.npy', TrxList)    # TxList[spw][ant,pol,ch]
    #-------- Az and El position at atmCal and onsource scans
    AtmEL, OnEL = np.ones([useAntNum, atmscanNum]), np.ones([useAntNum, scanNum])
    scanTimeRef = np.zeros(scanNum)
    for scan_index in range(scanNum): scanTimeRef[scan_index] = np.median(msmd.timesforscan(bpscanLists[band_index][scan_index]))
    for ant_index in range(useAntNum):
        azelTime_index = np.where( AntID == useAnt[ant_index] )[0].tolist()
        if len(azelTime_index) == 0: azelTime_index = np.where( AntID == 0 )[0].tolist()
        for scan_index in range(atmscanNum): AtmEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - atmTimeRef[scan_index]))]]
        for scan_index in range(scanNum):     OnEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - scanTimeRef[scan_index]))]]
    #
    OnsecZ, atmsecZ  = 1.0 / np.sin( np.median(OnEL, axis=0) ), 1.0 / np.sin( np.median(AtmEL, axis=0) )
    #-------- Tsky and TantN
    Tau0, Tau0Excess, TantN = tau0SpecFit(tempAtm - 5.0, atmsecZ, useAnt, atmspwLists[band_index], TskyList)
    for spw_index in range(spwNum): TrxList[spw_index] = (TrxList[spw_index].transpose(1,0,2) - TantN[spw_index]).transpose(1,0,2)
    LogTrx(antList[useAnt], atmspwLists[band_index], np.mean(np.array(TrxList), axis=3).transpose(1,0,2), 'Trec : ', tsysLog)
    LogTrx(antList[useAnt], atmspwLists[band_index], np.mean(np.array([TantN,TantN]), axis=3).transpose(2,1,0), 'TantN: ', tsysLog)
    freqList = []
    for spw_index in range(spwNum):
        chNum, chWid, freq = GetChNum(msfile, atmspwLists[band_index][spw_index]); freqList = freqList + [freq*1.0e-9]
        text_sd = 'SPW=%d : Tau(zenith) = %6.4f' % (atmspwLists[band_index][spw_index], np.median(Tau0[spw_index]))
        tsysLog.write(text_sd + '\n'); print text_sd
    #
    tsysLog.close()
    #-------- Save to npy files
    np.save(prefix +  '-' + UniqBands[band_index] + '.Trx.npy', TrxList)    # [spw][ant, pol, ch]
    np.save(prefix +  '-' + UniqBands[band_index] + '.Tau0.npy', Tau0)      # [spw][ch]
    np.save(prefix +  '-' + UniqBands[band_index] + '.TauE.npy', Tau0Excess)# [spw][scan]
#
    #---- Plots
    if not 'PLOTFMT' in locals():   PLOTFMT = 'pdf'
    if PLOTTAU: plotTau(prefix + '_' + UniqBands[band_index], atmspwLists[band_index], freqList, Tau0) 
    if PLOTTSYS: plotTsys(prefix + '_' + UniqBands[band_index], antList[useAnt], atmspwLists[band_index], freqList, atmTimeRef, TrxList, TskyList)
#
#-------- Plot optical depth
msmd.close()
msmd.done()
