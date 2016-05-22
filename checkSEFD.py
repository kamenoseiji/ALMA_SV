import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
#
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
#-------- Check Scans for atmCal
print '---Checking time series in MS and atmCal scans'
tb.open(msfile); timeXY = tb.query('ANTENNA1 == 0 && ANTENNA2 == 0 && DATA_DESC_ID == '+`spw[0]`).getcol('TIME'); tb.close()
OnTimeIndex = []
sourceIDscan = []
for scan_index in range(scanNum):
    OnTimeIndex.append( indexList(msmd.timesforscan(onsourceScans[scan_index]), timeXY) )
    sourceIDscan.append( msmd.fieldsforscan(onsourceScans[scan_index])[0])
#
offTime = sort( list(set(timeXY) & set(timeOFF)) )
ambTime = sort( list(set(timeXY) & set(timeAMB)) )
hotTime = sort( list(set(timeXY) & set(timeHOT)) )
#
offTimeIndex = indexList( offTime, timeXY)
ambTimeIndex = indexList( ambTime, timeXY)
hotTimeIndex = indexList( hotTime, timeXY)
#-------- Bandpass Table
print '---Generating antenna-based bandpass table'
XYdelayList = []
BPList = []
for spw_index in spw:
    BP_ant, XYdelay = BPtable(msfile, spw_index, BPScan, blMap, blInv)
    BPList = BPList + [BP_ant]
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
sys.stderr.write('\n')
sys.stderr.flush()
secZ = 1.0 / np.sin( OffEL )
#-------- Time-interpolation of ambient and hot
print '---Analyzing Trec and Tsky using atmCal scans'
chAvgTrx = np.zeros([UseAntNum, spwNum, 2, len(offTime)])
chAvgTsky= np.zeros([UseAntNum, spwNum, 2, len(offTime)])
chAvgTsys= np.zeros([UseAntNum, spwNum, 2, scanNum])
TrxFlag  = np.ones([UseAntNum, spwNum, 2, len(offTime)])
TsysFlag = np.ones([UseAntNum, spwNum, 2, scanNum])
onTau = np.zeros([spwNum, scanNum])
TrxList, TskyList = [], []
tempAmb, tempHot  = np.zeros([UseAntNum]), np.zeros([UseAntNum])
for ant_index in range(UseAntNum):
    tempAmb[ant_index], tempHot[ant_index] = GetLoadTemp(msfile, antMap[ant_index], spw[0])
    for spw_index in range(spwNum):
        AntSpwIndex = ant_index* spwNum + spw_index
        chNum = AmbSpecList[AntSpwIndex].shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        Trx = np.zeros([2, chNum, len(offTime)])
        Tsky= np.zeros([2, chNum, len(offTime)])
        for pol_index in range(ppolNum):
            ambSpec = AmbSpecList[AntSpwIndex][pol_index]
            hotSpec = HotSpecList[AntSpwIndex][pol_index]
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
#-------- Tsys for scans
for scan_index in range(scanNum):
    OnTimeRange = timeXY[OnTimeIndex[scan_index]]
    ambTimeIndex = argmin(abs(ambTime - OnTimeRange[0]))
    offTimeIndex = argmin(abs(offTime - OnTimeRange[0]))
    chAvgTsys[:,:,:,scan_index] = chAvgTrx[:,:,:,ambTimeIndex] + chAvgTsky[:,:,:,offTimeIndex] 
    TsysFlag[:,:,:,scan_index] = TrxFlag[:,:,:,ambTimeIndex]
#
#-------- Tau and TantN fitting
param = [3.0, 0.05]     # Initial Parameter
Tau0, TantN, Trms = np.zeros([UseAntNum, spwNum, 2]), np.zeros([UseAntNum, spwNum, 2]), np.zeros([UseAntNum, spwNum, 2])
for ant_index in range(UseAntNum):
    for spw_index in range(spwNum):
        for pol_index in range(2):
            fit = scipy.optimize.leastsq(residTskyTransfer, param, args=(tempAmb[ant_index]-20.0, secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index]))
            TantN[ant_index, spw_index, pol_index] = fit[0][0]
            Tau0[ant_index, spw_index, pol_index]  = fit[0][1]
        #
    #
#
Tau0err = np.std( Tau0, axis=(0,2) )
Tau0med = np.mean( Tau0, axis=(0,2) )
#-------- Tau on source
for scan_index in range(scanNum):
    OnTimeRange = timeXY[OnTimeIndex[scan_index]]
    offTimeIndex = argmin(abs(offTime - OnTimeRange[0]))
    onTau[:,scan_index] = -np.median(np.log((chAvgTsky[:,:,:,offTimeIndex] - TantN - np.median(tempAmb) + 20.0) / (2.718 - np.median(tempAmb) + 20.0)), axis=(0, 2))
#
for spw_index in range(spwNum):
    print 'SPW=%d : Tau(zenith) = %6.4f +- %6.4f' % (spw[spw_index], Tau0med[spw_index], Tau0err[spw_index])
#
np.save(prefix + '.Trx.npy', TrxList) 
np.save(prefix + '.Tsky.npy', TskyList) 
np.save(prefix + '.TrxFlag.npy', TrxFlag) 
np.save(prefix + '.Tau0.npy', Tau0) 
#-------- Antenna-dependent leakage noise
param = [5.0]
print 'TantN: ',
for spw_index in range(spwNum):
    print 'PolX    SPW%02d  PolY           |' % (spw[spw_index]),
#
print ' '
print '-----:--------------------------------+-------------------------------+-------------------------------+-------------------------------+'
for ant_index in range(UseAntNum):
    print antList[antMap[ant_index]] + ' : ',
    for spw_index in range(spwNum):
        for pol_index in range(2):
            fit = scipy.optimize.leastsq(residTskyTransfer2, param, args=(tempAmb[ant_index]-20.0, Tau0med[spw_index], secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index]))
            TantN[ant_index, spw_index, pol_index] = fit[0][0]
            resid = residTskyTransfer([fit[0][0], Tau0med[spw_index]], tempAmb[ant_index]-20.0, secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index])
            Trms[ant_index, spw_index, pol_index]  = sqrt(np.dot(resid, resid) / len(np.where(TrxFlag[ant_index, spw_index, pol_index] > 0.0)[0]))
            print '%4.1f (%4.1f) K ' % (TantN[ant_index, spw_index, pol_index], Trms[ant_index, spw_index, pol_index]),
        print '|',
    print ' '
print ' '
np.save(prefix + '.TantN.npy', TantN) 
#-------- Trx
print 'Trec : ',
for spw_index in range(spwNum): print 'SPW%02d  X        Y |' % (spw[spw_index]),
print ' '
print '-----:--------------------+-------------------+-------------------+-------------------+'
for ant_index in range(UseAntNum):
    print antList[antMap[ant_index]] + ' : ',
    for spw_index in range(spwNum):
        for pol_index in range(2):
            print '%6.1f K' % (np.median(chAvgTrx[ant_index, spw_index, pol_index])),
        print '|',
    print ' '
print ' '
#-------- Tsys
print 'Tsys : ',
for spw_index in range(spwNum): print 'SPW%02d  X        Y |' % (spw[spw_index]),
print ' '
print '-----:--------------------+-------------------+-------------------+-------------------+'
for ant_index in range(UseAntNum):
    print antList[antMap[ant_index]] + ' : ',
    for spw_index in range(spwNum):
        for pol_index in range(2):
            print '%6.1f K' % (np.median(chAvgTrx[ant_index, spw_index, pol_index]) + np.median(chAvgTsky[ant_index, spw_index, pol_index])),
        print '|',
    print ' '
#-------- Plot optical depth
if PLOTTAU: plotTau(prefix + ' ' + UniqBands[band_index], antList[antMap], spw, secZ, (chAvgTsky.transpose(3,0,1,2) - TantN).transpose(1,2,3,0), np.median(tempAmb) - 20.0, Tau0med, TrxFlag, 2.0*np.median(chAvgTsky), PLOTFMT) 
if PLOTTSYS: plotTsys(prefix + ' ' + UniqBands[band_index], antList[antMap], ambTime, spw, TrxList, TskyList, PLOTFMT)
##-------- Equalization using Bandpass scan
GainAnt = []
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
for spw_index in range(spwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    PA = AzEl2PA(np.median(OnAZ[:,onsourceScans.index(EQScan)]), np.median(OnEL[:,onsourceScans.index(EQScan)]))
    PS = InvPAMatrix(PA)
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)       # Cross Polarization Baseline Mapping
    Xspec = (tempSpec / (BPList[spw_index][ant0][:,polXindex]* BPList[spw_index][ant1][:,polYindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    #-------- XY delay cal
    XYdlSpec = delay_cal( np.ones([chNum], dtype=complex), XYdelayList[spw_index] )
    Xspec[1] = (Xspec[1].transpose(1,2,0)* XYdlSpec).transpose(2,0,1)
    Xspec[2] = (Xspec[2].transpose(1,2,0)* XYdlSpec.conjugate()).transpose(2,0,1)
    #-------- Antenna-based Gain
    chAvgVis = np.mean(Xspec[:, chRange], axis=1)
    GainX = np.apply_along_axis( gainComplex, 0, chAvgVis[0])
    GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[3])
    VisXY = np.array([np.median(gainCalVis( chAvgVis[0], GainX, GainX)), np.median(gainCalVis( chAvgVis[1], GainX, GainY)), np.median(gainCalVis( chAvgVis[2], GainY, GainX)), np.median(gainCalVis( chAvgVis[3], GainY, GainY))])
    StokesVis = np.dot(PS, VisXY)
    #print '%f %f %f %f' % (StokesVis[0], StokesVis[1], StokesVis[2], StokesVis[3])
    #
    #-------- Phase Cal and channel average
    BLphsX, BLphsY = GainX[ant0]* GainX[ant1].conjugate() / abs(GainX[ant0]* GainX[ant1]), GainY[ant0]* GainY[ant1].conjugate() / abs(GainY[ant0]* GainY[ant1])
    pCalVisX, pCalVisY = np.mean(chAvgVis[0] / BLphsX, axis=1), np.mean(chAvgVis[3] / BLphsY, axis=1)
    #-------- Antenna-based Gain
    GainAnt = GainAnt + [gainComplex(pCalVisX)]
    GainAnt = GainAnt + [gainComplex(pCalVisY)]
#
GainEq = abs(np.array(GainAnt)).reshape([spwNum, 2, UseAntNum]).transpose(2,0,1)    # Ant, SPW, Pol
#-------- Flux models for solar system objects
SSONum = len(SSOList)
timeLabel = qa.time('%fs' % (timeXY[0]), form='ymd')[0]
SSOflux0= []
SSOsize = []
centerFreqList = []
primaryBeam = np.ones([UseBlNum])
for bl_index in range(UseBlNum):
    beam0, beam1 = 1.0/antDia[ant0[bl_index]], 1.0/antDia[ant1[bl_index]] 
    primaryBeam[bl_index] = np.sqrt(2.0/ ((beam0)**2 + (beam1)**2 )) * beam0* beam1
#
for spw_index in range(spwNum): 
    chNum, chWid, Freq = GetChNum(msfile, spw[spw_index])
    centerFreqList.append( np.median(Freq)*1.0e-9 )
#
for ssoIndex in range(SSONum):
    for spw_index in range(spwNum): 
        text_Freq = '%6.2fGHz' % (centerFreqList[spw_index])
        SSOmodel = predictcomp(objname=sourceList[SSOList[ssoIndex]], standard="Butler-JPL-Horizons 2012", minfreq=text_Freq, maxfreq=text_Freq, nfreqs=1, prefix="", antennalist="aca.cycle3.cfg", epoch=timeLabel, showplot=T)
        SSOflux0.append(SSOmodel['spectrum']['bl0flux']['value'])
    #
    SSOsize.append(SSOmodel['shape']['majoraxis']['value']* pi / 21600.0)   # arcmin -> rad, diameter -> radius
#
plt.close('all')
SSOflux0= np.array(SSOflux0).reshape(SSONum, spwNum)     # [SSO, spw]
uvFlag = np.ones([SSONum, spwNum, UseBlNum])
SSOmodelVis = []
SSOscanID   = []
for ssoIndex in range(SSONum):
    UVlimit = 0.32 / SSOsize[ssoIndex]                              # Maximum uv distane(lambda) available for the SSO size
    scanID = list(set( msmd.scansforfield(SSOList[ssoIndex]).tolist()) & set(onsourceScans))[0]; SSOscanID.append(scanID)
    if( scanID == FCScan):
        FCS_ID = ssoIndex
        print 'Flux Calibrator is %s at %s' % (sourceList[SSOList[ssoIndex]], timeLabel)
    #
    timeStamp, UVW = GetUVW(msfile, spw[spw_index], scanID)
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    for spw_index in range(spwNum):
        uvWave = uvDist* centerFreqList[spw_index] / 0.299792458    # UV distance in wavelength
        uvFlag[ssoIndex, spw_index, np.where( uvWave > UVlimit )[0].tolist()] = 0.0
        SSOmodelVis.append(diskVisBeam(SSOsize[ssoIndex], uvWave, 1.13* 0.299792458* primaryBeam/centerFreqList[spw_index]))
    #
#
SSOflux = SSOflux0* np.exp(-onTau.transpose(1,0)[indexList(np.array(SSOscanID), np.array(onsourceScans))])
SSOmodelVis = np.array(SSOmodelVis).reshape(SSONum, spwNum, UseBlNum)
FCSmodelVis = SSOmodelVis[FCS_ID]
FCSFlag     = uvFlag[FCS_ID]
#
##-------- Scaling with the flux calibrator
medSF, sdSF = [], []; del GainAnt
for spw_index in range(spwNum):
    #-------- Sub-array with unflagged antennas (short baselines)
    SAantennas, SAblMap, SAblFlag, SAant0, SAant1 = subArranIndex(uvFlag[FCS_ID, spw_index])
    SAantNum = len(SAantennas); SAblNum = SAantNum* (SAantNum - 1)/2
    if SAantNum < 4:
        print 'Too few antennas for %s. Change REFANT or Flux Calibrator.'
        sys.exit(0)
    #
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], FCScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    PA = AzEl2PA(np.median(OnAZ[:,onsourceScans.index(FCScan)]), np.median(OnEL[:,onsourceScans.index(FCScan)]))
    PS = InvPAMatrix(PA)
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1) 
    Xspec = (tempSpec / (BPList[spw_index][ant0][:,polXindex]* BPList[spw_index][ant1][:,polYindex].conjugate())).transpose(2,3,1,0)[:,:,SAblMap]
    Xspec = (Xspec.transpose(1,3,2,0)/(GainEq[SAant0, spw_index][:,polXindex]* GainEq[SAant1, spw_index][:,polYindex].conjugate())).transpose(3,0,2,1)
    #-------- XY delay cal
    XYdlSpec = delay_cal( np.ones([chNum], dtype=complex), XYdelayList[spw_index] )
    Xspec[1] = (Xspec[1].transpose(1,2,0)* XYdlSpec).transpose(2,0,1)
    Xspec[2] = (Xspec[2].transpose(1,2,0)* XYdlSpec.conjugate()).transpose(2,0,1)
    #-------- Antenna-based Gain
    chAvgVis =(np.mean(Xspec[:, chRange], axis=1).transpose(0,2,1) / FCSmodelVis[spw_index, SAblMap]* SAblFlag).transpose(0,2,1)
    GainX, GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[0]), np.apply_along_axis( gainComplex, 0, chAvgVis[3])
    #
    #-------- Phase Cal and channel average
    BLphsX = GainX[ant0[0:SAblNum]]* GainX[ant1[0:SAblNum]].conjugate() / abs(GainX[ant0[0:SAblNum]]* GainX[ant1[0:SAblNum]])
    BLphsY = GainY[ant0[0:SAblNum]]* GainY[ant1[0:SAblNum]].conjugate() / abs(GainY[ant0[0:SAblNum]]* GainY[ant1[0:SAblNum]])
    pCalVisX, pCalVisY = np.mean(chAvgVis[0] / BLphsX, axis=1), np.mean(chAvgVis[3] / BLphsY, axis=1)
    #-------- Antenna-based Gain
    GainAnt = abs(gainComplex(pCalVisX))**2; medSF = medSF + [np.median(GainAnt)]; sdSF = sdSF + [np.std(GainAnt)/np.sqrt(SAantNum-1)]
    GainAnt = abs(gainComplex(pCalVisY))**2; medSF = medSF + [np.median(GainAnt)]; sdSF = sdSF + [np.std(GainAnt)/np.sqrt(SAantNum-1)]
#
medSF, sdSF = np.array(medSF).reshape([spwNum, 2]), np.array(sdSF).reshape([spwNum, 2])
scaleFact = (SSOflux[FCS_ID] / medSF.transpose(1,0)).transpose(1,0)
FCS_Eq = GainEq/ np.sqrt(scaleFact)
SEFD   = 1.0 / FCS_Eq**2
AEFF   = ((2761.297* chAvgTsys[:,:,:,onsourceScans.index(FCScan)]/SEFD).transpose(1,2,0)/(0.25* pi*antDia**2)).transpose(2,0,1)
print 'Aeff :',
for spw_index in range(spwNum):
    for pol_index in range(2):
        print 'SPW%02d-%s ' % (spw[spw_index], PolList[pol_index]),
    #
#
print ''
for ant_index in range(UseAntNum):
    print '%s :' % (antList[antMap[ant_index]]),
    for spw_index in range(spwNum):
        for pol_index in range(2):
            print '  %4.1f%% ' % (100.0* AEFF[ant_index, spw_index, pol_index]),
        #
    #
    print ''
#
#-------- Antenna-based Gain
print '---Flux densities of sources ---'
print 'Scan   Source  EL   ',
for spw_index in range(spwNum): print 'SPW%02d         ' % (spw[spw_index]),
print ' '
for scan_index in range(scanNum):
    print '%02d %010s %4.1f ' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0* np.median(OnEL[:,scan_index])/pi ),
    PA = AzEl2PA(np.median(OnAZ[:,scan_index]), np.median(OnEL[:,scan_index]))
    PS = InvPAMatrix(PA)
    if(onsourceScans[scan_index] in SSOscanID):
        SSO_flag = T
        SSO_ID = SSOscanID.index(onsourceScans[scan_index])
    else:
        SSO_flag = F
    for spw_index in range(spwNum):
        #atmCorrect = np.exp(-Tau0med[spw_index]/ np.sin(np.median(OnEL[:, scan_index])))
        atmCorrect = np.exp(-onTau[spw_index, scan_index])
        #-------- Sub-array with unflagged antennas (short baselines)
        if SSO_flag:
            SAantennas, SAblMap, SAblFlag, SAant0, SAant1 = subArranIndex(uvFlag[SSO_ID, spw_index])
        else:
            SAantennas, SAblMap, SAblFlag, SAant0, SAant1 = range(UseAntNum), range(UseBlNum), np.ones([blNum]), ant0, ant1
        #
        SAantNum = len(SAantennas); SAblNum = SAantNum* (SAantNum - 1)/2
        if SAantNum < 4:
            print ' too few ants ',
            continue
        #
        #
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], onsourceScans[scan_index])
        chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)[:,SAblMap]       # Cross Polarization Baseline Mapping
        #-------- Bandpass Calibration
        GainBL = FCS_Eq[SAant0,spw_index][:,polXindex]* FCS_Eq[SAant1,spw_index][:,polYindex]* atmCorrect
        BP_bl = ((BPList[spw_index][SAant0][:,polXindex]* BPList[spw_index][SAant1][:,polYindex].conjugate()).transpose(2,0,1)* GainBL).transpose(1,2,0)
        Xspec = (tempSpec / BP_bl).transpose(2,3,1,0) # Bandpass Cal
        #-------- XY delay cal
        XYdlSpec = delay_cal(np.ones([chNum], dtype=complex), XYdelayList[spw_index])
        Xspec[1] = (Xspec[1].transpose(1,2,0)* XYdlSpec).transpose(2,0,1)
        Xspec[2] = (Xspec[2].transpose(1,2,0)* XYdlSpec.conjugate()).transpose(2,0,1)
        #-------- Antenna-based Gain
        if(SSO_flag):
            chAvgVis =(np.mean(Xspec[:, chRange], axis=1).transpose(0,2,1) / SSOmodelVis[SSO_ID, spw_index, SAblMap]).transpose(0,2,1)
        else:
            chAvgVis = np.mean(Xspec[:, chRange], axis=1)
        #
        GainX = np.apply_along_axis( gainComplex, 0, chAvgVis[0]) #; Xflag = np.where( np.std(abs(GainX), axis=1) < np.median(abs(GainX), axis=1))[0]
        GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[3]) #; Yflag = np.where( np.std(abs(GainY), axis=1) < np.median(abs(GainY), axis=1))[0]
        # Tflag = list( set(Xflag) & set(Yflag) )
        # print len(Tflag)
        # if len(Tflag) < 4: print ' low SNR   ',; continue
        #VisXY = np.array([np.median(gainCalVis( chAvgVis[0], GainX, GainX)), np.median(gainCalVis( chAvgVis[1], GainX, GainY)), np.median(gainCalVis( chAvgVis[2], GainY, GainX)), np.median(gainCalVis( chAvgVis[3], GainY, GainY))])
        #StokesVis = np.dot(PS, VisXY)
        #print '%f %f %f %f' % (StokesVis[0], StokesVis[1], StokesVis[2], StokesVis[3])
        #
        #-------- Phase Cal and channel average
        BLphsX = GainX[ant0[0:SAblNum]]* GainX[ant1[0:SAblNum]].conjugate() / abs(GainX[ant0[0:SAblNum]]* GainX[ant1[0:SAblNum]])
        BLphsY = GainY[ant0[0:SAblNum]]* GainY[ant1[0:SAblNum]].conjugate() / abs(GainY[ant0[0:SAblNum]]* GainY[ant1[0:SAblNum]])
        pCalVisX = np.mean((chAvgVis[0]/BLphsX), axis=1)
        pCalVisY = np.mean((chAvgVis[3]/BLphsY), axis=1)
        #-------- Antenna-based Gain
        GainAntX, GainAntY = abs(gainComplex(pCalVisX))**2, abs(gainComplex(pCalVisY))**2
        #-------- Check 
        gainXflag = np.where(abs(GainAntX - np.median(GainAntX)) / np.median(GainAntX) > 0.16)[0].tolist()
        gainYflag = np.where(abs(GainAntY - np.median(GainAntY)) / np.median(GainAntY) > 0.16)[0].tolist()
        Xants, Yants = list(set(range(SAantNum)) - set(gainXflag)), list(set(range(SAantNum)) - set(gainYflag))
        Iants = list( set(Xants) & set(Yants) )
        if len(Iants) < 2: print ' too low SNR  ',; continue
        GainI = (GainAntX[Iants] + GainAntY[Iants])/2.0
        meanI, sdI = np.mean(GainI), np.std(GainI)/sqrt(len(GainI) - 1.0)
        print '%6.3f (%3.1f%%) ' % (meanI, 100.0* sdI/meanI),
    #
    if(SSO_flag):
        for spw_index in range(spwNum):
            print ' %6.3f ' % (SSOflux0[SSO_ID, spw_index]),
        #
    #
    print ' '
#