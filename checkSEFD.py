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
sourceList = msmd.fieldnames()
numSource = len(sourceList)
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
print '---Checking atmCal scans'
spw_index = 0; ant_index = 0
timeXY, Pspec = GetPSpec(msfile, ant_index, spw[spw_index])
polNum, chNum = Pspec.shape[0], Pspec.shape[1]
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
#-------- Prepare BP and Delay to store
print '---Generating antenna-based bandpass table'
XYdelayList = []
BPList = []
for spw_index in spw:
    BP_ant, XYdelay = BPtable(msfile, spw_index, BPScan)
    BPList = BPList + [BP_ant]
    XYdelayList = XYdelayList + [XYdelay]
#
plotMax = 1.25
if PLOTBP:
    #-------- Prepare Plots
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index, figsize = (11, 8))
        figAnt.suptitle(prefix + ' ' + antList[ant_index] + ' Scan = ' + `BPScan`)
        figAnt.text(0.45, 0.05, 'Frequency [GHz]')
        figAnt.text(0.03, 0.45, 'Bandpass Amplitude and Phase', rotation=90)
    #
    #-------- Plot BP
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index)
        for spw_index in range(spwNum):
            chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq = 1.0e-9* Freq  # GHz
            BPampPL = figAnt.add_subplot( 2, spwNum, spw_index + 1 )
            BPphsPL = figAnt.add_subplot( 2, spwNum, spw_index + spwNum + 1 )
            for pol_index in range(ppolNum):
                plotBP = BPList[spw_index][ant_index, pol_index]
                BPampPL.plot( Freq, abs(plotBP), ls='steps-mid', label = 'Pol=' + PolList[pol_index])
                BPampPL.axis([np.min(Freq), np.max(Freq), 0.0, 1.25* plotMax])
                BPampPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
                BPampPL.yaxis.offsetText.set_fontsize(10)
                BPampPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
                BPphsPL.plot( Freq, np.angle(plotBP), '.', label = 'Pol=' + PolList[pol_index])
                BPphsPL.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
            #
            BPampPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            BPphsPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
            BPampPL.text( np.min(Freq), 1.1* plotMax, 'SPW=' + `spw[spw_index]` + ' Amp')
            BPphsPL.text( np.min(Freq), 2.5, 'SPW=' + `spw[spw_index]` + ' Phase')
        #
        if PLOTFMT == 'png':
            figAnt.savefig('BP_' + prefix + '_' + antList[ant_index] + '-REF' + antList[0] + '_Scan' + `BPScan` + '.png')
        else :
            figAnt.savefig('BP_' + prefix + '_' + antList[ant_index] + '-REF' + antList[0] + '_Scan' + `BPScan` + '.pdf')
        #
    #
    plt.close('all')
#
"""
#
#-------- Load Az El position
azelTime, AntID, AZ, EL = GetAzEl(msfile)
#-------- Load autocorrelation power spectra
print '---Loading autocorr power spectra'
OffSpec = np.zeros([antNum, spwNum, 2, chNum, len(offTime)])    # 2 stand for num of pPol
AmbSpec = np.zeros([antNum, spwNum, 2, chNum, len(ambTime)])
HotSpec = np.zeros([antNum, spwNum, 2, chNum, len(hotTime)])
OnSpec  = np.zeros([antNum, spwNum, 2, chNum, scanNum])
OffEL   = np.zeros([antNum, len(offTime)])
scanEL  = np.zeros([antNum, scanNum])
for ant_index in range(antNum):
    for spw_index in range(spwNum):
        progress = (1.0* ant_index* spwNum + spw_index) / (antNum* spwNum)
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress))
        sys.stderr.flush()
        timeXY, Pspec = GetPSpec(msfile, ant_index, spw[spw_index])
        OffSpec[ant_index, spw_index] = Pspec[pPol][:,:,offTimeIndex]
        AmbSpec[ant_index, spw_index] = Pspec[pPol][:,:,ambTimeIndex]
        HotSpec[ant_index, spw_index] = Pspec[pPol][:,:,hotTimeIndex]
        for scan_index in range(scanNum):
            OnTimeIndex = indexList(msmd.timesforscan(onsourceScans[scan_index]), timeXY)
            OnSpec[ant_index, spw_index, :, :, scan_index] = np.median( Pspec[pPol][:,:,OnTimeIndex], axis=2 )
        #
    #
    azelTime_index = np.where( AntID == ant_index )[0].tolist()
    for scan_index in range(scanNum):
        refTime = np.median(msmd.timesforscan(onsourceScans[scan_index]))
        scanEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - refTime))]]
    #
    for time_index in range(len(offTimeIndex)):
        OffEL[ant_index, time_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - offTime[time_index]))]]
    #
#
sys.stderr.write('\n')
sys.stderr.flush()
secZ = 1.0 / np.sin( OffEL )
#-------- Time-interpolation of ambient and hot
print '---Analyzing time variability of autocorr power'
timeAvgOffSpec = np.mean( OffSpec, axis=4 )
timeAvgAmbSpec = np.mean( AmbSpec, axis=4 )
timeAvgHotSpec = np.mean( HotSpec, axis=4 )
timeAvgOnSpec  = np.mean( OnSpec, axis=4 )
chAvgOff = (np.mean(OffSpec[:,:,:,chRange], axis=3).transpose(3,0,1,2) / np.mean(timeAvgOffSpec[:,:,:,chRange],axis=3)).transpose(1,2,3,0)
chAvgAmb = (np.mean(AmbSpec[:,:,:,chRange], axis=3).transpose(3,0,1,2) / np.mean(timeAvgAmbSpec[:,:,:,chRange],axis=3)).transpose(1,2,3,0)
chAvgHot = (np.mean(HotSpec[:,:,:,chRange], axis=3).transpose(3,0,1,2) / np.mean(timeAvgHotSpec[:,:,:,chRange],axis=3)).transpose(1,2,3,0)
chAvgOn  = (np.mean(OnSpec[:,:,:,chRange], axis=3).transpose(3,0,1,2)  / np.mean(timeAvgOnSpec[:,:,:,chRange],axis=3)).transpose(1,2,3,0)
#-------- Off-source Tsys
Trx = np.zeros([antNum, spwNum, 2, chNum, len(offTime)])
Tsky= np.zeros([antNum, spwNum, 2, chNum, len(offTime)])
chAvgTrx = np.zeros([antNum, spwNum, 2, len(offTime)])
chAvgTsky= np.zeros([antNum, spwNum, 2, len(offTime)])
TrxFlag  = np.ones([antNum, spwNum, 2, len(offTime)])
tempAmb  = np.zeros([antNum])
tempHot  = np.zeros([antNum])
for ant_index in range(antNum):
    tempAmb[ant_index], tempHot[ant_index] = GetLoadTemp(msfile, ant_index, spw[0])
    for spw_index in range(spwNum):
        for pol_index in range(2):
            SPL_amb = UnivariateSpline(ambTime, chAvgAmb[ant_index, spw_index, pol_index], s=0.001)
            SPL_hot = UnivariateSpline(hotTime, chAvgHot[ant_index, spw_index, pol_index], s=0.001)
            for time_index in range(len(offTime)):
                Psamb = timeAvgAmbSpec[ant_index, spw_index, pol_index]* SPL_amb(offTime[time_index])
                Pshot = timeAvgHotSpec[ant_index, spw_index, pol_index]* SPL_hot(offTime[time_index])
                Psoff = OffSpec[ant_index, spw_index, pol_index, :, time_index]
                Trx[ant_index, spw_index, pol_index, :, time_index] = (tempHot[ant_index]* Psamb - Pshot* tempAmb[ant_index]) / (Pshot - Psamb) 
                Tsky[ant_index, spw_index, pol_index, :, time_index]= (Psoff* (tempHot[ant_index] - tempAmb[ant_index]) + tempAmb[ant_index]* Pshot - tempHot[ant_index]* Psamb) / (Pshot - Psamb)
                Phot, Pamb, Poff = np.mean(Pshot[chRange]), np.mean(Psamb[chRange]), np.mean(Psoff[chRange])
                chAvgTrx[ant_index, spw_index, pol_index, time_index] = (tempHot[ant_index]* Pamb - Phot* tempAmb[ant_index]) / (Phot - Pamb)
                chAvgTsky[ant_index, spw_index, pol_index, time_index]= (Poff* (tempHot[ant_index] - tempAmb[ant_index]) + tempAmb[ant_index]* Phot - tempHot[ant_index]* Pamb) / (Phot - Pamb)
            #
            #-------- Trx flagging
            TrxTemp = chAvgTrx[ant_index, spw_index, pol_index]
            TrxFlag[ant_index, spw_index, pol_index, np.where( abs(TrxTemp - np.median(TrxTemp)) > 0.1* np.mean(TrxTemp))[0].tolist()] = 0.0
            TrxFlag[ant_index, spw_index, pol_index, np.where( TrxTemp < 1.0)[0].tolist()] = 0.0
        #
    #
#
param = [5.0, 0.05]
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
Tau0 = np.median(np.median(Tau0, axis=0), axis=1)
for spw_index in range(spwNum):
    print 'SPW=%d : Tau(zenith) = %6.4f' % (spw[spw_index], Tau0[spw_index])
#
np.save(prefix + '.Trx.npy', Trx) 
np.save(prefix + '.Tsky.npy', Tsky) 
np.save(prefix + '.TrxFlag.npy', TrxFlag) 
np.save(prefix + '.Tau0.npy', Tau0) 
np.save(prefix + '.TantN.npy', TantN) 
msmd.close()
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
            fit = scipy.optimize.leastsq(residTskyTransfer2, param, args=(tempAmb[ant_index], Tau0[spw_index], secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index]))
            TantN[ant_index, spw_index, pol_index] = fit[0][0]
            resid = residTskyTransfer([fit[0][0], Tau0[spw_index]], tempAmb[ant_index], secZ[ant_index], chAvgTsky[ant_index, spw_index, pol_index], TrxFlag[ant_index, spw_index, pol_index])
            Trms[ant_index, spw_index, pol_index]  = sqrt(np.dot(resid, resid) / len(np.where(TrxFlag[ant_index, spw_index, pol_index] > 0.0)[0]))
            #print '%s SPW=%d %s : TantN=%6.3f K  Trms=%6.3f K' % (antList[ant_index], spw[spw_index], PolList[pol_index], fit[0][0], Trms[ant_index, spw_index, pol_index])
            print '%4.1f (%4.1f) K ' % (TantN[ant_index, spw_index, pol_index], Trms[ant_index, spw_index, pol_index]),
        #
        print '|',
    #
    print ' '
#
print ' '
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
            TskyPL.plot( airmass, 2.713* np.exp(-Tau0[spw_index]* airmass) + tempAmb[ant_index]* (1.0 - np.exp(-Tau0[spw_index]* airmass)), '-')
            for ant_index in range(antNum):
                plotTsky = chAvgTsky[ant_index, spw_index, pol_index] - TantN[ant_index, spw_index, pol_index]
                TskyPL.scatter( secZ[ant_index], plotTsky, s=15* TrxFlag[ant_index, spw_index, pol_index] + 1, color=cm.gist_ncar( float(ant_index) / antNum ), alpha=0.25, label = antList[ant_index])
            #
            text_sd = 'Pol %s Tau(zenith)=%6.4f' % (PolList[pol_index], Tau0[spw_index])
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
            chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq = 1.0e-9* Freq  # GHz
            for scan_index in range(numTimePlot):
                TsysPL = figAnt.add_subplot(numTimePlot, spwNum, spwNum* scan_index + spw_index + 1 )
                timeRange = np.where( abs( ambTime - TimePlot[scan_index]) < timeThresh )[0].tolist()
                timeLabel = qa.time('%fs' % (TimePlot[scan_index]), form='fits')[0]
                for pol_index in range(2):
                    plotTsys = np.mean(Tsky[ant_index, spw_index, pol_index][:, timeRange] + Trx[ant_index, spw_index, pol_index][:, timeRange], axis=1)
                    plotTrx  = np.mean(Trx[ant_index, spw_index, pol_index][:, timeRange], axis=1)
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
"""
#-------- Antenna-based Gain
GainAnt = []
for scan_index in range(scanNum):
    for spw_index in range(spwNum):
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], onsourceScans[scan_index])
        Xspec = Xspec[pPol]; Xspec = Xspec[:,:,blMap]
        Ximag = Xspec.transpose(0,1,3,2).imag * (-2.0* np.array(blInv) + 1.0)
        Xspec.imag = Ximag.transpose(0,1,3,2)
        #-------- Bandpass Calibration
        BP_bl = BPList[spw_index][ant0]* BPList[spw_index][ant1].conjugate()
        Xspec = (Xspec.transpose(3,2,0,1) / BP_bl).transpose(2,3,1,0)
        #-------- Antenna-based Gain
        chAvgVis = np.mean( Xspec[:, chRange], axis=1 )
        GainX = np.apply_along_axis( gainComplex, 0, chAvgVis[0])
        GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[1])
        #-------- Phase Cal and channel average
        BLphsX = GainX[ant0]* GainX[ant1].conjugate() / abs(GainX[ant0]* GainX[ant1])
        BLphsY = GainY[ant0]* GainY[ant1].conjugate() / abs(GainY[ant0]* GainY[ant1])
        pCalVisX = np.mean((Xspec[0] / BLphsX)[chRange], axis=(0,2))
        pCalVisY = np.mean((Xspec[1] / BLphsY)[chRange], axis=(0,2))
        #-------- Antenna-based Gain
        GainAnt = GainAnt + [gainComplex(pCalVisX)]
        GainAnt = GainAnt + [gainComplex(pCalVisY)]
    #
#

    #
#
"""
#-------- Baseline-based bandpass
BP_bl = (np.ones([blNum, spwNum, 2, chNum], dtype=complex)* BP_ant[ant0]* BP_ant[ant1].conjugate()).transpose(1, 2, 3, 0)
#-------- Plot bandpass
"""
