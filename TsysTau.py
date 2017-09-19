execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
def residTskyTransfer( param, Tamb, secz, Tsky, weight=1.0 ):
    exp_Tau = np.exp( -param[1]* secz )
    return weight* (Tsky - (param[0] + 2.718* exp_Tau  + Tamb* (1.0 - exp_Tau)))
#
def residTskyTransfer0( param, Tamb, secz, Tsky, weight=1.0 ):
    exp_Tau = np.exp( -param[0]* secz )
    return weight* (Tsky - (2.718* exp_Tau  + Tamb* (1.0 - exp_Tau)))
#
def residTskyTransfer2( param, Tamb, Tau0, secz, Tsky, weight=1.0 ):
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
msfile = prefix + '.ms'
antList = GetAntName(msfile)
antNum = len(antList)
msmd.open(msfile)
print '---Checking time for ambient and hot load'
timeOFF, timeAMB, timeHOT, atmScanList = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT"), msmd.scansforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE").tolist()
scanNum = len(atmScanList)
if len(timeAMB) == 0:
    ambTimeIndex, hotTimeIndex = [], []
    timeXY, Pspec = GetPSpec(msfile, 0, spwList[0])
    timeNum, chNum = Pspec.shape[2], Pspec.shape[1]
    if 'chRange' not in locals(): chRange = range(int(0.05*chNum), int(0.95*chNum))
    chAvgPower = np.mean(Pspec[0][chRange], axis=0)
    offTimeIndex = indexList(timeOFF, timeXY)
    calTimeIndex = indexList(msmd.timesforintent("CALIBRATE_ATMOSPHERE#ON_SOURCE"), timeXY)
    scanBreakIndex = np.where( np.diff(timeXY[calTimeIndex]) > 30.0)[0].tolist() + [len(calTimeIndex)-1]
    startBreak = 0
    for break_index in scanBreakIndex:
        scanSeq = np.array(calTimeIndex)[startBreak:break_index]
        ambHotBreak = np.where( np.diff(timeXY[scanSeq]) > 2.0* np.median(np.diff(timeXY[scanSeq])))[0][0]
        ambTimeIndex = ambTimeIndex + scanSeq[1:(ambHotBreak-1)].tolist()
        hotTimeIndex = hotTimeIndex + scanSeq[(ambHotBreak+1):(len(scanSeq)-1)].tolist()
        startBreak = break_index + 1
    #
    ambTime, hotTime, offTime = timeXY[ambTimeIndex], timeXY[hotTimeIndex], timeXY[offTimeIndex]
else:
    tb.open(msfile); timeXY = tb.query('ANTENNA1 == 0 && ANTENNA2 == 0 && DATA_DESC_ID == '+`spwList[0]`).getcol('TIME'); tb.close()
    offTime, ambTime, hotTime = sort( list(set(timeXY) & set(timeOFF)) ), sort( list(set(timeXY) & set(timeAMB)) ), sort( list(set(timeXY) & set(timeHOT)) )
    offTimeIndex, ambTimeIndex, hotTimeIndex = indexList(offTime, timeXY),  indexList(ambTime, timeXY),  indexList(hotTime, timeXY)
#
#-------- Check scan timings
OffTimeIndex = []
for scanID in atmScanList: OffTimeIndex.append( indexList(msmd.timesforscan(scanID), timeXY) )
#-------- Load autocorrelation power spectra
print '---Loading autocorr power spectra'
OffSpecList, AmbSpecList, HotSpecList = [], [], []
spwNum, ppolNum = len(spwList), len(pPol)
for ant_index in range(antNum):
    for spw_index in range(spwNum):
        progress = (1.0* ant_index* spwNum + spw_index + 1.0) / (antNum* spwNum)
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
        timeXY, Pspec = GetPSpec(msfile, ant_index, spwList[spw_index])
        OffSpecList.append(Pspec[pPol][:,:,offTimeIndex])
        AmbSpecList.append(Pspec[pPol][:,:,ambTimeIndex])
        HotSpecList.append(Pspec[pPol][:,:,hotTimeIndex])
    #
#
sys.stderr.write('\n')
sys.stderr.flush()
#-------- Load Az El position
azelTime, AntID, AZ, EL = GetAzEl(msfile)
OffEL = np.ones([antNum, len(offTimeIndex)])
for ant_index in range(antNum):
    azelTime_index = np.where( AntID == ant_index )[0].tolist()
    for time_index in range(len(offTimeIndex)):
        OffEL[ant_index, time_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - offTime[time_index]))]]
    #
#
secZ = 1.0 / np.sin( OffEL )
#-------- Time-interpolation of ambient and hot
print '---Analyzing Trec and Tsky using atmCal scans'
chAvgTrx, chAvgTsky, chAvgTsys = np.zeros([antNum, spwNum, 2, len(offTime)]), np.zeros([antNum, spwNum, 2, len(offTime)]), np.zeros([antNum, spwNum, 2, scanNum])
TrxList, TskyList = [], []
tempAmb, tempHot  = np.zeros([antNum]), np.zeros([antNum])
for ant_index in range(antNum):
    tempAmb[ant_index], tempHot[ant_index] = GetLoadTemp(msfile, ant_index, spwList[0])
    if tempAmb[ant_index] < 250: tempAmb[ant_index] += 273.15
    if tempHot[ant_index] < 300: tempHot[ant_index] += 273.15
    for spw_index in range(spwNum):
        AntSpwIndex = ant_index* spwNum + spw_index
        chNum = AmbSpecList[AntSpwIndex].shape[1]
        Trx, Tsky = np.zeros([2, chNum, len(offTime)]), np.zeros([2, chNum, len(offTime)])
        for pol_index in range(ppolNum):
            ambSpec, hotSpec = AmbSpecList[AntSpwIndex][pol_index], HotSpecList[AntSpwIndex][pol_index]
            for time_index in range(len(offTime)):
                Psamb, Pshot = np.median(ambSpec, axis=1), np.median(hotSpec, axis=1)
                Psoff = OffSpecList[AntSpwIndex][pol_index][:,time_index]
                Trx[pol_index, :, time_index] = (tempHot[ant_index]* Psamb - Pshot* tempAmb[ant_index]) / (Pshot - Psamb)
                Tsky[pol_index, :, time_index]= (Psoff* (tempHot[ant_index] - tempAmb[ant_index]) + tempAmb[ant_index]* Pshot - tempHot[ant_index]* Psamb) / (Pshot - Psamb)
            #
        #
        TrxList = TrxList + [Trx]
        TskyList= TskyList + [Tsky]
    #
#
TrxSpec  = np.array(TrxList).reshape(antNum, spwNum, 2, chNum, len(offTime))
TskySpec = np.array(TskyList).reshape(antNum, spwNum, 2, chNum, len(offTime))
#-------- Tau and TantN fitting
param = [0.0, 0.05]     # Initial Parameter
Tau0, TantN = np.zeros([spwNum, chNum]), np.zeros([spwNum, chNum])
for spw_index in range(spwNum):
    for ch_index in chRange:
        if 'Tant' in locals():
            param = [0.05]
            fit = scipy.optimize.leastsq(residTskyTransfer0, param, args=(np.median(tempAmb)-Tatm_OFS, np.median(secZ,axis=0), np.median(TskySpec, axis=(0,2))[spw_index, ch_index]-Tant[spw_index]))
            Tau0[spw_index, ch_index]  = fit[0][0]
            TantN[spw_index, ch_index] = Tant[spw_index]
        else:
            fit = scipy.optimize.leastsq(residTskyTransfer, param, args=(np.median(tempAmb)-Tatm_OFS, np.median(secZ,axis=0), np.median(TskySpec, axis=(0,2))[spw_index, ch_index]))
            TantN[spw_index, ch_index] = fit[0][0]
            Tau0[spw_index, ch_index]  = fit[0][1]
        #
    #
#
np.save(prefix +  '.Trx.npy', TrxSpec) 
np.save(prefix +  '.Tau0.npy', Tau0) 
np.save(prefix +  '.TantN.npy', TantN) 
#-------- Antenna-dependent leakage noise
msmd.close()
msmd.done()
