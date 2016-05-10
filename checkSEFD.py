import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import matplotlib.cm as cm
execfile(SCR_DIR + 'interferometry.py')
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
#-------- Procedures
msfile = wd + prefix + '.ms'
#-------- Check Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
blNum = antNum* (antNum - 1) / 2
if refant:
    try:
        refantID = np.where( antList == refant )[0][0]
    except:
        refantID = 0
        print 'Warning: refant ' + refant + ' is not found in the antenna list.'
    #
else:
    refantID = 0
#
print 'Use ' + antList[refantID] + ' as the refant.'
#-------- Baseline Mapping
Refant = [refantID] + range(refantID) + range(refantID+1, antNum)
blMap = range(blNum)
blInv = [False]* blNum
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
for bl_index in range(blNum):
    #ants = Bl2Ant(bl_index)
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(Refant[ant0[bl_index]], Refant[ant1[bl_index]])
#
print `len(np.where( blInv )[0])` + ' baselines are inverted.'

#-------- Check SPWs of atmCal
msmd.open(msfile)
print '---Checking spectral windows'
spw = list(set(msmd.tdmspws()) & set(msmd.spwsforintent("CALIBRATE_ATMOSPHERE*")))
spwNum = len(spw)
#-------- Check source list
print '---Checking source list'
# sourceList = msmd.fieldnames()
sourceList, posList = GetSourceList(msfile) 
numSource = len(sourceList)
SSOList   = np.where( (np.array(posList)[:,0] == 0.0) & (np.array(posList)[:,1] == 0.0) )[0].tolist()   # Solar System Objects
#-------- Check MJD for Ambient Load
print '---Checking time for ambient and hot load'
timeOFF = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE")
timeAMB = msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT")
timeHOT = msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT")
#-------- Check bandpass and on-source scans
print '---Checking bandpass and on-source time'
BPScan        = msmd.scansforintent("CALIBRATE_BANDPASS#ON_SOURCE")[0]
FCScan        = msmd.scansforintent("CALIBRATE_FLUX#ON_SOURCE")[0]
onsourceScans = [BPScan] + [FCScan] + msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE").tolist()
scanNum = len(onsourceScans)
#-------- Check Scans for atmCal
print '---Checking time series in MS and atmCal scans'
tb.open(msfile); timeXY = tb.query('ANTENNA1 == 0 && ANTENNA2 == 0 && DATA_DESC_ID == '+`spw[0]`).getcol('TIME'); tb.close()
OnTimeIndex = []
for scan_index in range(scanNum):
    OnTimeIndex.append( indexList(msmd.timesforscan(onsourceScans[scan_index]), timeXY) )
#
polNum = msmd.ncorrforpol(msmd.polidfordatadesc(spw[0]))
pPol, cPol = [0,1], []  # parallel and cross pol
PolList = ['X', 'Y']
if polNum == 4:
    pPol, cPol = [0,3], [1,2]  # parallel and cross pol
#
ppolNum, cpolNum = len(pPol), len(cPol)
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
    BP_ant, XYdelay = BPtable(msfile, spw_index, BPScan)
    BPList = BPList + [BP_ant]
    XYdelayList = XYdelayList + [XYdelay]
#
if PLOTBP:
    figAnt = PlotBP(msfile, antList, spw, BPList)
    fileExt = '.pdf'
    if PLOTFMT == 'png': fileExt = '.png'
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index)
        plotFigFileName = 'BP_' + prefix + '_' + antList[ant_index] + '_REF' + antList[0] + '_Scan' + `BPScan` + fileExt
        figAnt.savefig(plotFigFileName)
    #
    plt.close('all')
#
#
#-------- Load autocorrelation power spectra
print '---Loading autocorr power spectra'
OnSpecList, OffSpecList, AmbSpecList, HotSpecList = [], [], [], []
antDia = np.ones([antNum])
for ant_index in range(antNum):
    antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
    for spw_index in range(spwNum):
        progress = (1.0* ant_index* spwNum + spw_index) / (antNum* spwNum)
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
        timeXY, Pspec = GetPSpec(msfile, ant_index, spw[spw_index])
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
OnEL, OffEL = np.ones([antNum, scanNum]), np.ones([antNum, len(offTimeIndex)])
for ant_index in range(antNum):
    azelTime_index = np.where( AntID == ant_index )[0].tolist()
    for scan_index in range(scanNum):
        refTime = np.median(msmd.timesforscan(onsourceScans[scan_index]))
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
chAvgTrx = np.zeros([antNum, spwNum, 2, len(offTime)])
chAvgTsky= np.zeros([antNum, spwNum, 2, len(offTime)])
chAvgTsys= np.zeros([antNum, spwNum, 2, scanNum])
TrxFlag  = np.ones([antNum, spwNum, 2, len(offTime)])
TsysFlag = np.ones([antNum, spwNum, 2, scanNum])
TrxList, TskyList = [], []
tempAmb, tempHot  = np.zeros([antNum]), np.zeros([antNum])
for ant_index in range(antNum):
    tempAmb[ant_index], tempHot[ant_index] = GetLoadTemp(msfile, ant_index, spw[0])
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
            #
            #-------- Trx flagging
            TrxTemp = chAvgTrx[ant_index, spw_index, pol_index]
            TrxFlag[ant_index, spw_index, pol_index, np.where( abs(TrxTemp - np.median(TrxTemp)) > 0.2* np.median(TrxTemp))[0].tolist()] = 0.0
            TrxFlag[ant_index, spw_index, pol_index, np.where( TrxTemp < 1.0)[0].tolist()] = 0.0
            #-------- Tsys for scans
            for scan_index in range(scanNum):
                OnTimeRange = timeXY[OnTimeIndex[scan_index]]
                chAvgTsys[ant_index, spw_index, pol_index, scan_index] = chAvgTrx[ant_index, spw_index, pol_index, argmin(abs(ambTime - OnTimeRange[0]))] + chAvgTsky[ant_index, spw_index, pol_index, argmin(abs(offTime - OnTimeRange[0]))] 
                TsysFlag[ant_index, spw_index, pol_index, scan_index] = TrxFlag[ant_index, spw_index, pol_index, argmin(abs(ambTime - OnTimeRange[0]))]
            #
        #
        TrxList.append(Trx)
        TskyList.append(Tsky)
    #
#
#-------- Tau and TantN fitting
param = [5.0, 0.05]     # Initial Parameter
Tau0 = np.zeros([antNum, spwNum, 2])
TantN= np.zeros([antNum, spwNum, 2])
Trms = np.zeros([antNum, spwNum, 2])
for ant_index in range(antNum):
    for spw_index in range(spwNum):
        for pol_index in range(2):
            fit = scipy.optimize.leastsq(residTskyTransfer, param, args=(tempAmb[ant_index], secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index]))
            TantN[ant_index, spw_index, pol_index] = fit[0][0]
            Tau0[ant_index, spw_index, pol_index]  = fit[0][1]
        #
    #
#
Tau0err = np.std( Tau0, axis=(0,2) )
Tau0med = np.mean( Tau0, axis=(0,2) )
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
for ant_index in range(antNum):
    print antList[ant_index] + ' : ',
    for spw_index in range(spwNum):
        for pol_index in range(2):
            fit = scipy.optimize.leastsq(residTskyTransfer2, param, args=(tempAmb[ant_index], Tau0med[spw_index], secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index]))
            TantN[ant_index, spw_index, pol_index] = fit[0][0]
            resid = residTskyTransfer([fit[0][0], Tau0med[spw_index]], tempAmb[ant_index], secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index])
            Trms[ant_index, spw_index, pol_index]  = sqrt(np.dot(resid, resid) / len(np.where(TrxFlag[ant_index, spw_index, pol_index] > 0.0)[0]))
            print '%4.1f (%4.1f) K ' % (TantN[ant_index, spw_index, pol_index], Trms[ant_index, spw_index, pol_index]),
        #
        print '|',
    #
    print ' '
#
print ' '
np.save(prefix + '.TantN.npy', TantN) 
#-------- Trx
print 'Trec : ',
for spw_index in range(spwNum):
    print 'SPW%02d  X        Y |' % (spw[spw_index]),
#
print ' '
print '-----:--------------------+-------------------+-------------------+-------------------+'
for ant_index in range(antNum):
    print antList[ant_index] + ' : ',
    for spw_index in range(spwNum):
        for pol_index in range(2):
            print '%6.1f K' % (np.median(chAvgTrx[ant_index, spw_index, pol_index])),
        #
        print '|',
    #
    print ' '
#
print ' '
#-------- Tsys
print 'Tsys : ',
for spw_index in range(spwNum):
    print 'SPW%02d  X        Y |' % (spw[spw_index]),
#
print ' '
print '-----:--------------------+-------------------+-------------------+-------------------+'
for ant_index in range(antNum):
    print antList[ant_index] + ' : ',
    for spw_index in range(spwNum):
        for pol_index in range(2):
            print '%6.1f K' % (np.median(chAvgTrx[ant_index, spw_index, pol_index]) + np.median(chAvgTsky[ant_index, spw_index, pol_index])),
        #
        print '|',
    #
    print ' '
#
#-------- Plot optical depth
if PLOTTAU:
    TrmThresh = 2.0* np.median(Trms)
    plotMax = 2.0 * np.median(chAvgTsky)
    airmass = np.arange( 1.0, 1.25*np.max(secZ), 0.01)
    figTau = plt.figure(0, figsize = (11,8))
    figTau.suptitle(prefix + ' Optical Depth')
    figTau.text(0.45, 0.05, 'Airmass')
    figTau.text(0.03, 0.45, 'Sky Temperature [K]', rotation=90)
    for spw_index in range(spwNum):
        for pol_index in range(2):
            TskyPL = figTau.add_subplot(2, spwNum, spwNum* pol_index + spw_index + 1 )
            TskyPL.axis([1.0, 1.25*np.max(secZ), 0.0, plotMax])
            TskyPL.plot( airmass, 2.713* np.exp(-Tau0med[spw_index]* airmass) + tempAmb[ant_index]* (1.0 - np.exp(-Tau0med[spw_index]* airmass)), '-')
            for ant_index in range(antNum):
                plotTsky = chAvgTsky[ant_index, spw_index, pol_index] - TantN[ant_index, spw_index, pol_index]
                TskyPL.scatter( secZ[ant_index], plotTsky, s=15* TrxFlag[ant_index, spw_index, pol_index] + 1, color=cm.gist_ncar( float(ant_index) / antNum ), alpha=0.25, label = antList[ant_index])
            #
            text_sd = 'Pol %s Tau(zenith)=%6.4f' % (PolList[pol_index], Tau0med[spw_index])
            TskyPL.text(1.01, 0.95* plotMax, text_sd, fontsize='9')
            if pol_index == 0:
                TskyPL.set_title('SPW ' + `spw[spw_index]`)
            #
        #
    #
    TskyPL.legend(loc = 'lower right', prop={'size' :7}, numpoints = 1)
    if PLOTFMT == 'png':
        figTau.savefig('TAU_' + prefix + '.png')
    else :
        figTau.savefig('TAU_' + prefix + '.pdf')
    plt.close('all')
#
if PLOTTSYS:
    #-------- Plots for Tsys spectra
    timeThresh = 30.0
    TimePlot = ambTime[np.where( diff(ambTime) < timeThresh)[0].tolist()]
    numTimePlot = len(TimePlot)
    plotMax = 1.5 * np.median(Trx + Tsky)
    #-------- Prepare Plots
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index, figsize = (8, 11))
        figAnt.suptitle(prefix + ' ' + antList[ant_index])
        figAnt.text(0.45, 0.05, 'Frequency [GHz]')
        figAnt.text(0.03, 0.45, 'Tsys (solid) and Trec (dotted) [K]', rotation=90)
    #
    #-------- Plot BP
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index)
        for spw_index in range(spwNum):
            AntSpwIndex = ant_index* spwNum + spw_index
            chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq = 1.0e-9* Freq  # GHz
            for scan_index in range(numTimePlot):
                TsysPL = figAnt.add_subplot(numTimePlot, spwNum, spwNum* scan_index + spw_index + 1 )
                timeRange = np.where( abs( ambTime - TimePlot[scan_index]) < timeThresh )[0].tolist()
                timeLabel = qa.time('%fs' % (TimePlot[scan_index]), form='fits')[0]
                for pol_index in range(2):
                    plotTrx  = np.mean(TrxList[AntSpwIndex][pol_index][:, timeRange], axis=1)
                    plotTsys = np.mean(TskyList[AntSpwIndex][pol_index][:, timeRange], axis=1) + plotTrx
                    TsysPL.plot( Freq, plotTsys, ls='steps-mid', label = 'Tsys_Pol=' + PolList[pol_index])
                    TsysPL.plot( Freq, plotTrx,  ls=':', label = 'Trec_Pol=' + PolList[pol_index])
                #
                TsysPL.axis([np.min(Freq), np.max(Freq), 0.0, plotMax])
                if scan_index == 0:
                    TsysPL.set_title('SPW ' + `spw[spw_index]`)
                #
                if scan_index < numTimePlot - 1:
                    TsysPL.set_xticklabels([])
                #
                if spw_index == 0:
                    TsysPL.text(np.min(Freq), 0.8* plotMax, timeLabel, fontsize='8')
                else:
                    TsysPL.set_yticklabels([])
                #
            #
        #
        if PLOTFMT == 'png':
            figAnt.savefig('TSYS_' + prefix + '_' + antList[ant_index] + '.png')
        else :
            figAnt.savefig('TSYS_' + prefix + '_' + antList[ant_index] + '.pdf')
        #
    #
    plt.close('all')
#
#



#-------- Flux models for solar system objects
SSONum = len(SSOList)
timeLabel = qa.time('%fs' % (timeXY[0]), form='ymd')[0]
SSOflux = []
SSOsize = []
centerFreqList = []
primaryBeam = np.ones([blNum])
for bl_index in range(blNum):
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
        SSOflux.append(SSOmodel['spectrum']['bl0flux']['value']* np.exp(Tau0med[spw_index] / np.sin(SSOmodel['azel']['m1']['value'])))
    #
    SSOsize.append(SSOmodel['shape']['majoraxis']['value']* pi / 21600.0)   # arcmin -> rad, diameter -> radius
#
plt.close('all')
SSOflux = np.array(SSOflux).reshape(SSONum, spwNum)     # [SSO, spw]
uvFlag = np.ones([SSONum, spwNum, blNum])
SSOmodelVis = []
SSOscanID   = []
for ssoIndex in range(SSONum):
    UVlimit = 0.32 / SSOsize[ssoIndex]                              # Maximum uv distane(lambda) available for the SSO size
    scanID = list(set( msmd.scansforfield(SSOList[ssoIndex]).tolist()) & set(onsourceScans))[0]; SSOscanID.append(scanID)
    timeStamp, UVW = GetUVW(msfile, spw[spw_index], scanID)
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    for spw_index in range(spwNum):
        uvWave = uvDist* centerFreqList[spw_index] / 0.299792458    # UV distance in wavelength
        uvFlag[ssoIndex, spw_index, np.where( uvWave > UVlimit )[0].tolist()] = 0.0
        SSOmodelVis.append(SSOflux[ssoIndex, spw_index]* diskVisBeam(SSOsize[ssoIndex], uvWave, 1.13* 0.299792458* primaryBeam/centerFreqList[spw_index]))
    #
#
SSOmodelVis = np.array(SSOmodelVis).reshape(SSONum, spwNum, blNum)

#-------- Antenna-based Gain
GainAnt = []
print '---Equalization based on SEFD'
for scan_index in range(scanNum):
    if(onsourceScans[scan_index] in SSOscanID):
        SSO_flag = T
        SSO_ID = SSOscanID.index(onsourceScans[scan_index])
    else:
        SSO_flag = F
    for spw_index in range(spwNum):
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], onsourceScans[scan_index])
        chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        Xspec = Xspec[pPol]; Xspec = Xspec[:,:,blMap]
        Ximag = Xspec.transpose(0,1,3,2).imag * (-2.0* np.array(blInv) + 1.0)
        Xspec.imag = Ximag.transpose(0,1,3,2)
        #-------- Bandpass Calibration
        BP_bl = BPList[spw_index][ant0]* BPList[spw_index][ant1].conjugate()
        Xspec = (Xspec.transpose(3,2,0,1) / BP_bl).transpose(2,3,1,0)
        #-------- Antenna-based Gain
        if(SSO_flag):
            chAvgVis =(np.mean(Xspec[:, chRange], axis=1).transpose(0,2,1) / SSOmodelVis[SSO_ID, spw_index]).transpose(0,2,1)
        else:
            chAvgVis = np.mean(Xspec[:, chRange], axis=1)
        # 
        GainX = np.apply_along_axis( gainComplex, 0, chAvgVis[0])
        GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[1])
        #
        #-------- Phase Cal and channel average
        BLphsX = GainX[ant0]* GainX[ant1].conjugate() / abs(GainX[ant0]* GainX[ant1])
        BLphsY = GainY[ant0]* GainY[ant1].conjugate() / abs(GainY[ant0]* GainY[ant1])
        pCalVisX = np.mean(chAvgVis[0] / BLphsX, axis=1)
        pCalVisY = np.mean(chAvgVis[1] / BLphsY, axis=1)
        #-------- Antenna-based Gain
        GainAnt = GainAnt + [gainComplex(pCalVisX)]
        GainAnt = GainAnt + [gainComplex(pCalVisY)]
    #
#
SSOscanList = indexList(np.array(SSOscanID), np.array(onsourceScans))
AE = 2761.297*((abs(np.array(GainAnt))**2).reshape(scanNum, spwNum, ppolNum, antNum).transpose(3,1,2,0)*chAvgTsys)[:,:,:,SSOscanList]
AEFF = (AE.transpose(1,2,3,0) / (0.25* pi*antDia**2)).transpose(3,0,1,2)

##-------- Equalization using Bandpass scan
#scan_index = onsourceScans.index(BPScan)
#BP_SEFD = 1.0 /abs(GainAnt[:,:,:,scan_index])**2    # SEFD assuming Flux = 1 Jy
#BP_AEFF = 2761.297 * ((chAvgTsys[:,:,:,scan_index] / BP_SEFD).transpose(1,2,0) / (0.25* pi*antDia**2)).transpose(2,0,1)
#tempFlux = abs(GainAnt).transpose(3,0,1,2)**2 * BP_SEFD
#-------- Scaling
msmd.done()
