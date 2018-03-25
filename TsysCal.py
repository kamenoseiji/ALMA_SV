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
def scanAtmSpec(msfile, antNum, scanList, spwList, timeOFF=0, timeON=0, timeAMB=0, timeHOT=0):
    timeList, offSpecList, ambSpecList, hotSpecList = [], [], [], []
    scanNum, spwNum = len(scanList), len(spwList)
    pPol = [0,1]
    index = 0
    for ant_index in range(antNum):
        progress = (1.0* ant_index + 1.0) / antNum
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
        for spwID in spwList:
            for scanID in scanList:
                # print '----%s SPW%02d Scan%02d Index%02d' % (antList[ant_index], spwID, scanID, index)
                timeXY, Pspec = GetPSpecScan(msfile, ant_index, spwID, scanID)
                chNum = Pspec.shape[1]
                if Pspec.shape[0] == 4: pPol = [0,3]
                offTime, ambTime, hotTime = sort( list(set(timeXY) & set(timeOFF)) ), sort( list(set(timeXY) & set(timeAMB)) ), sort( list(set(timeXY) & set(timeHOT)) )
                offTimeIndex, ambTimeIndex, hotTimeIndex = indexList(offTime, timeXY)[-1],  indexList(ambTime, timeXY)[-1],  indexList(hotTime, timeXY)[-1]
                if ((ant_index == 0) & (spwID == spwList[0]) & (len(timeAMB) > 0)):
                    timeList = timeList + [offTime[-1]]
                # 
                if ((ant_index == 0) & (spwID == spwList[0]) & (len(timeAMB) == 0)):
                    chRange = range(int(0.05*chNum), int(0.95*chNum)); chAvgPower = np.mean(Pspec[0][chRange], axis=0)
                    offTimeIndex = indexList(timeOFF, timeXY)
                    ambhotTimeIndex = indexList(timeON, timeXY)
                    ambhotThresh = 0.5*(max(chAvgPower[ambhotTimeIndex]) + min(chAvgPower[ambhotTimeIndex]))
                    hotTimeIndex = np.array(ambhotTimeIndex)[np.where( chAvgPower[ambhotTimeIndex] > ambhotThresh )[0].tolist()].tolist()
                    ambTimeIndex = np.array(ambhotTimeIndex)[np.where( chAvgPower[ambhotTimeIndex] < ambhotThresh )[0].tolist()].tolist()
                    timeList = timeList + [np.median(timeXY[offTimeIndex])]
                #
                if len(timeAMB) > 0 :
                    offSpecList = offSpecList + [Pspec[pPol][:,:,offTimeIndex]]
                    ambSpecList = ambSpecList + [Pspec[pPol][:,:,ambTimeIndex]]
                    hotSpecList = hotSpecList + [Pspec[pPol][:,:,hotTimeIndex]]
                    index += 1
                else:
                    offSpecList = offSpecList + [np.median(Pspec[pPol][:,:,offTimeIndex], axis=2)]
                    ambSpecList = ambSpecList + [np.median(Pspec[pPol][:,:,ambTimeIndex], axis=2)]
                    hotSpecList = hotSpecList + [np.median(Pspec[pPol][:,:,hotTimeIndex], axis=2)]
                #
            #
        #
    #            
    sys.stderr.write('\n'); sys.stderr.flush()
    return np.array(timeList), offSpecList, ambSpecList, hotSpecList
#
Tatm_OFS  = 15.0     # Ambient-load temperature - Atmosphere temperature
Tcmb      = 2.725    # CMB temperature
#-------- Check MS file
msfile = prefix + '.ms'
tempAtm = GetTemp(msfile)
antList = GetAntName(msfile)
antNum = len(antList)
if 'flagAnt' not in locals(): flagAnt = np.ones(antNum)
#-------- Check SPWs
msmd.open(msfile)
print '---Checking spectral windows and scans with atmCal for ' + prefix
atmSPWs = list(set(msmd.tdmspws()) & set(msmd.spwsforintent("CALIBRATE_ATMOSPHERE*"))); atmSPWs.sort()
bpSPWs  = msmd.spwsforintent("CALIBRATE_BANDPASS*").tolist(); bpSPWs.sort()
atmspwNames, bpspwNames = msmd.namesforspws(atmSPWs), msmd.namesforspws(bpSPWs)
bpSPWs = np.array(bpSPWs)[indexList(np.array(atmspwNames), np.array(bpspwNames))].tolist(); bpspwNames = msmd.namesforspws(bpSPWs)
atmBandNames, atmPattern = [], r'RB_..'
for spwName in atmspwNames : atmBandNames = atmBandNames + re.findall(atmPattern, spwName)
UniqBands = unique(atmBandNames).tolist(); NumBands = len(UniqBands)
atmspwLists, bpspwLists, atmscanLists, bpscanLists = [], [], [], []
for band_index in range(NumBands):
    atmspwLists = atmspwLists + [np.array(atmSPWs)[indexList(np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()]
    bpspwLists  = bpspwLists  + [np.array(bpSPWs)[indexList( np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()]
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
azelTime, AntID, AZ, EL = GetAzEl(msfile)
# timeOFF : mjd of CALIBRATE_ATMOSPHERE#OFF_SOURCE
# timeON  : mjd of CALIBRATE_ATMOSPHERE#ON_SOURCE (becore Cycle 3, ambient + hot loads
# timeAMB : mjd of CALIBRATE_ATMOSPHERE#AMBIENT (after Cycle 3)
# timeHOT : mjd of CALIBRATE_ATMOSPHERE#HOT (after Cycle 3)
#
OnElList = []
for band_index in range(NumBands):
    TrxList, TskyList, chAvgTsys = [], [], []
    tsysLog = open(prefix + '-' + UniqBands[band_index] + '-Tsys.log', 'w')
    atmTimeRef, offSpec, ambSpec, hotSpec = scanAtmSpec(msfile, antNum, atmscanLists[band_index], atmspwLists[band_index], timeOFF, timeON, timeAMB, timeHOT)
    atmscanNum, scanNum = len(atmscanLists[band_index]), len(bpscanLists[band_index])
    #-------- Az and El position at atmCal and onsource scans
    AtmEL, OnEL = np.ones([antNum, atmscanNum]), np.ones([antNum, scanNum])
    scanTimeRef = np.zeros(scanNum)
    for scan_index in range(scanNum): scanTimeRef[scan_index] = np.median(msmd.timesforscan(bpscanLists[band_index][scan_index]))
    for ant_index in range(antNum):
        if flagAnt[ant_index] < 1.0: continue
        azelTime_index = np.where( AntID == ant_index )[0].tolist()
        for scan_index in range(atmscanNum): AtmEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - atmTimeRef[scan_index]))]]
        for scan_index in range(scanNum):     OnEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - scanTimeRef[scan_index]))]]
    #
    OnsecZ, atmsecZ  = 1.0 / np.sin( OnEL ), 1.0 / np.sin( AtmEL )
    OnElList = OnElList + [OnEL]
    #-------- Time-interpolation of ambient and hot
    print '---Analyzing ' + UniqBands[band_index] + ' Trec and Tsky using atmCal scans'
    spwNum = len(atmspwLists[band_index])
    chAvgTrx, chAvgTsky = np.zeros([antNum, spwNum, 2, atmscanNum]), np.zeros([antNum, spwNum, 2, atmscanNum])
    TrxFlag = np.ones([antNum, spwNum, 2, atmscanNum])
    tempAmb, tempHot  = np.zeros([antNum]), np.zeros([antNum])
    for ant_index in range(antNum):
        # if flagAnt[ant_index] < 1.0: continue
        tempAmb[ant_index], tempHot[ant_index] = GetLoadTemp(msfile, ant_index, atmspwLists[band_index][0])
        if tempAmb[ant_index] < 240: tempAmb[ant_index] += 273.15       # Old MS describes the load temperature in Celsius
        if tempHot[ant_index] < 240: tempHot[ant_index] += 273.15       #
        for spw_index in range(spwNum):
            chNum = offSpec[spw_index* atmscanNum].shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); chOut = sort(list(set(range(chNum)) - set(chRange)))
            TrxSpec, TskySpec = np.zeros([2, chNum, atmscanNum]), np.zeros([2, chNum, atmscanNum])
            for pol_index in range(2):
                for scan_index in range(atmscanNum):
                    index = (ant_index* spwNum + spw_index)* atmscanNum + scan_index
                    Psoff, Psamb, Pshot = offSpec[index][pol_index], ambSpec[index][pol_index], hotSpec[index][pol_index]
                    TrxSpec[pol_index, :, scan_index] = (tempHot[ant_index]* Psamb - Pshot* tempAmb[ant_index]) / (Pshot - Psamb)
                    TskySpec[pol_index, :, scan_index] = (Psoff* (tempHot[ant_index] - tempAmb[ant_index]) + tempAmb[ant_index]* Pshot - tempHot[ant_index]* Psamb) / (Pshot - Psamb)
                    TrxSpec[pol_index, chOut, scan_index] = np.median(TrxSpec[pol_index, chRange, scan_index])
                    TskySpec[pol_index, chOut, scan_index] = np.median(TskySpec[pol_index, chRange, scan_index])
                    Phot, Pamb, Poff = np.mean(Pshot[chRange]), np.mean(Psamb[chRange]), np.mean(Psoff[chRange])
                    chAvgTrx[ant_index, spw_index, pol_index, scan_index] = (tempHot[ant_index]* Pamb - Phot* tempAmb[ant_index]) / (Phot - Pamb)
                    chAvgTsky[ant_index, spw_index, pol_index, scan_index]= (Poff* (tempHot[ant_index] - tempAmb[ant_index]) + tempAmb[ant_index]* Phot - tempHot[ant_index]* Pamb) / (Phot - Pamb)
                #-------- Trx flagging
                TrxTemp = chAvgTrx[ant_index, spw_index, pol_index]
                TrxFlag[ant_index, spw_index, pol_index, np.where( abs(TrxTemp - np.median(TrxTemp)) > 0.2* np.median(TrxTemp))[0].tolist()] = 0.0
                TrxFlag[ant_index, spw_index, pol_index, np.where( TrxTemp < 1.0)[0].tolist()] = 0.0
                if 'SSOScanIndex' in locals(): TrxFlag[ant_index, spw_index, pol_index, SSOScanIndex] = 0.0
            #
            TrxList = TrxList + [np.median(TrxSpec, axis=2)] # Trx is time invariable
            TskyList= TskyList+ [np.mean(TskySpec, axis=0)]  # Sky is unpolarized
        #
    #
    chAvgTsky = np.mean(chAvgTsky, axis=2)  # Sky is unpolarized
    #-------- Tau and TantN fitting
    param = [0.0, 0.05]     # Initial Parameter
    Tau0, TantN, Tau0err, Tau0med = np.zeros([antNum, spwNum]), np.zeros([antNum, spwNum]), np.zeros(spwNum), np.zeros(spwNum)
    Tau0spec = []
    if max(atmsecZ[0]) - min(atmsecZ[0]) > 0.5:
        for ant_index in range(antNum):
            if flagAnt[ant_index] < 1.0: continue
            for spw_index in range(spwNum):
                fit = scipy.optimize.leastsq(residTskyTransfer, param, args=(tempAtm-Tatm_OFS, atmsecZ[ant_index], chAvgTsky[ant_index, spw_index], np.min(TrxFlag, axis=2)[ant_index, spw_index]))
                TantN[ant_index, spw_index] = fit[0][0]
                Tau0[ant_index, spw_index]  = fit[0][1]
    #
    else: Tau0 = 0.05*np.ones([antNum, spwNum]); TantN = Tatm_OFS* np.ones([antNum, spwNum])
    for spw_index in range(spwNum):
        chNum = TskyList[spw_index].shape[0]; chTau = np.zeros([antNum, chNum])
        for ant_index in range(antNum):
            if flagAnt[ant_index] < 1.0: continue
            index = ant_index* spwNum + spw_index
            for ch_index in range(chNum):
                fit = scipy.optimize.leastsq(residTskyTransfer0, [Tau0[ant_index, spw_index]], args=(tempAtm-Tatm_OFS, atmsecZ[ant_index], TskyList[index][ch_index]-TantN[ant_index,spw_index], np.min(TrxFlag, axis=2)[ant_index, spw_index]))
                chTau[ant_index, ch_index] = fit[0][0]
            #
        #   
        Tau0spec = Tau0spec + [np.median(chTau, axis=0)]   # Tau is antnna independent 
        Tau0err[spw_index] = np.median(np.std(chTau, axis=0))/sqrt(antNum - 1.0)
        Tau0med[spw_index] = np.median(chTau)
    #
    #-------- Trx Transfer
    chAvgTrx = np.median(chAvgTrx, axis=3)
    bpSPWs = bpspwLists[band_index]; bpspwNum = len(bpSPWs)
    #-------- Tsys for onsource scans
    for scan_index in range(scanNum):
        for ant_index in range(antNum):
            for spw_index in range(bpspwNum):
                index = ant_index* spwNum + spw_index
                exp_Tau = np.exp(-OnsecZ[ant_index, scan_index]* Tau0spec[spw_index] )
                for pol_index in range(2):
                    Tsys = TantN[ant_index, spw_index] + TrxList[index][pol_index] + Tcmb* exp_Tau + (tempAtm-Tatm_OFS)* (1.0 - exp_Tau)
                    # TsysList = TsysList + [Tsys]
                    chAvgTsys = chAvgTsys + [np.median(Tsys)]
                #
            #
        #
    #
    #-------- Log TantN
    text_sd = 'TantN:----------+---------+---------+---------+\n'; tsysLog.write(text_sd); print text_sd,
    for ant_index in range(antNum):
        if flagAnt[ant_index] < 1.0: continue
        text_sd = antList[ant_index] + ' : '; tsysLog.write(text_sd); print text_sd,
        for spw_index in range(spwNum):
            text_sd = '%5.1f K' % (TantN[ant_index, spw_index])
            tsysLog.write(text_sd); print text_sd,
            text_sd = '|'; tsysLog.write(text_sd); print text_sd,
        tsysLog.write('\n'); print ''
    tsysLog.write('\n'); print ''
    #-------- Log Trx
    text_sd = 'Trec : '
    for spw_index in range(spwNum): text_sd = text_sd + '  SPW%02d X        Y |' % (atmspwLists[band_index][spw_index])
    tsysLog.write(text_sd + '\n'); print text_sd
    text_sd =  ' ----:--------------------+-------------------+-------------------+-------------------+'; tsysLog.write(text_sd + '\n'); print text_sd
    for ant_index in range(antNum):
        if flagAnt[ant_index] < 1.0: continue
        text_sd =  antList[ant_index] + ' : '; tsysLog.write(text_sd); print text_sd,
        for spw_index in range(spwNum):
            for pol_index in range(2):
                text_sd = '%6.1f K' % (chAvgTrx[ant_index, spw_index, pol_index])
                tsysLog.write(text_sd); print text_sd,
            text_sd = '|'; tsysLog.write(text_sd); print text_sd,
        tsysLog.write('\n'); print ' '
    print ' '
    #-------- Log Tsys
    text_sd = 'Tsys : '
    for spw_index in range(spwNum): text_sd = text_sd + '  SPW%02d X        Y |' % (atmspwLists[band_index][spw_index])
    tsysLog.write(text_sd + '\n'); print text_sd
    text_sd =  ' ----:--------------------+-------------------+-------------------+-------------------+'; tsysLog.write(text_sd + '\n'); print text_sd
    for ant_index in range(antNum):
        if flagAnt[ant_index] < 1.0: continue
        text_sd =  antList[ant_index] + ' : '; tsysLog.write(text_sd); print text_sd,
        for spw_index in range(spwNum):
            for pol_index in range(2):
                index = (((arange(scanNum)* antNum + ant_index)* spwNum + spw_index)* 2 + pol_index).tolist()
                text_sd = '%6.1f K' % (np.median(np.array(chAvgTsys)[index]))
                tsysLog.write(text_sd); print text_sd,
            text_sd = '|'; tsysLog.write(text_sd); print text_sd,
        tsysLog.write('\n'); print ' '
    print ' '
    for spw_index in range(spwNum): text_sd = 'SPW=%d : Tau(zenith) = %6.4f +- %6.4f' % (atmspwLists[band_index][spw_index], Tau0med[spw_index], Tau0err[spw_index]); tsysLog.write(text_sd + '\n'); print text_sd
    tsysLog.close()
    np.save(prefix +  '-' + UniqBands[band_index] + '.Trx.npy', TrxList) 
    np.save(prefix +  '-' + UniqBands[band_index] + '.Tsky.npy', TskyList) 
    np.save(prefix +  '-' + UniqBands[band_index] + '.TrxFlag.npy', TrxFlag) 
    np.save(prefix +  '-' + UniqBands[band_index] + '.Tau0.npy', Tau0spec) 
    #---- Plots
    if not 'PLOTFMT' in locals():   PLOTFMT = 'pdf'
    if PLOTTAU: plotTau(prefix + '_' + UniqBands[band_index], antList, atmspwLists[band_index], atmsecZ, (chAvgTsky.transpose(2,0,1) - TantN).transpose(1,2,0), tempAtm - Tatm_OFS, Tau0med, np.min(TrxFlag, axis=2), 2.0*np.median(chAvgTsky), PLOTFMT) 
    if PLOTTSYS: plotTsys(prefix + '_' + UniqBands[band_index], antList, atmTimeRef, atmspwLists[band_index], TrxList, TskyList, PLOTFMT)
#
#-------- Plot optical depth
msmd.close()
msmd.done()
