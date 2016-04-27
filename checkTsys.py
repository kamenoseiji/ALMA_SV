import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
#
#-------- Definitions
#-------- Procedures
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile)
antNum = len(antList)
#-------- Check SPWs of atmCal
msmd.open(msfile)
print '---Checking spectral windows'
TDMspw_atmCal = list(set(msmd.tdmspws()) & set(msmd.spwsforintent("CALIBRATE_ATMOSPHERE*")))
spwNum = len(TDMspw_atmCal)
#-------- Check source list
print '---Checking source list'
sourceList = msmd.fieldnames()
numSource = len(sourceList)
#-------- Check MJD for Ambient Load
print '---Checking time for ambient and hot load'
timeOFF = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE")
timeAMB = msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT")
timeHOT = msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT")
#-------- Check Scans on-source
print '---Checking on-source time'
BPScans       = msmd.scansforintent("CALIBRATE_BANDPASS#ON_SOURCE")
onsourceScans = msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE")
scanNum = len(onsourceScans)
#
spw_index = 0; ant_index = 0
timeXY, Pspec = GetPSpec(msfile, ant_index, TDMspw_atmCal[spw_index])
polNum, chNum = Pspec.shape[0], Pspec.shape[1]
pPol, cPol = [0,1], []  # parallel and cross pol
if polNum == 4:
    pPol, cPol = [0,3], [1,2]  # parallel and cross pol
#
offTime = sort( list(set(timeXY) & set(timeOFF)) ).tolist()
ambTime = sort( list(set(timeXY) & set(timeAMB)) ).tolist()
hotTime = sort( list(set(timeXY) & set(timeHOT)) ).tolist()
#
offTimeIndex = indexList( offTime, timeXY)
ambTimeIndex = indexList( ambTime, timeXY)
hotTimeIndex = indexList( hotTime, timeXY)
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
        timeXY, Pspec = GetPSpec(msfile, ant_index, TDMspw_atmCal[spw_index])
        OffSpec[ant_index, spw_index] = Pspec[pPol][:,:,offTimeIndex]
        AmbSpec[ant_index, spw_index] = Pspec[pPol][:,:,ambTimeIndex]
        HotSpec[ant_index, spw_index] = Pspec[pPol][:,:,hotTimeIndex]
        for scan_index in range(len(onsourceScans)):
            OnTimeIndex = indexList(msmd.timesforscan(onsourceScans[scan_index]), timeXY)
            OnSpec[ant_index, spw_index, :, :, scan_index] = np.median( Pspec[pPol][:,:,OnTimeIndex], axis=2 )
        #
    #
    azelTime_index = np.where( AntID == ant_index )[0].tolist()
    for scan_index in range(len(onsourceScans)):
        refTime = np.median(msmd.timesforscan(onsourceScans[scan_index]))
        scanEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - refTime))]]
    #
    for time_index in range(len(offTimeIndex)):
        OffEL[ant_index, time_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - offTime[time_index]))]]
    #
#
#-------- Time-interpolation of ambient and hot
print '---Analyzing time variability of autocorr power'
chRange = range(int(0.05*chNum), int(0.95*chNum))
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
Tsys= np.zeros([antNum, spwNum, 2, chNum, len(offTime)])
chAvgTrx = np.zeros([antNum, spwNum, 2, len(offTime)])
chAvgTsys= np.zeros([antNum, spwNum, 2, len(offTime)])
for ant_index in range(antNum):
    tempAmb, tempHot = GetLoadTemp(msfile, ant_index, TDMspw_atmCal[0])
    for spw_index in range(spwNum):
        for pol_index in range(2):
            SPL_amb = UnivariateSpline(ambTime, chAvgAmb[ant_index, spw_index, pol_index], s=0.001)
            SPL_hot = UnivariateSpline(hotTime, chAvgHot[ant_index, spw_index, pol_index], s=0.001)
            for time_index in range(len(offTime)):
                Psamb = timeAvgAmbSpec[ant_index, spw_index, pol_index]* SPL_amb(offTime[time_index])
                Pshot = timeAvgHotSpec[ant_index, spw_index, pol_index]* SPL_hot(offTime[time_index])
                Psoff = OffSpec[ant_index, spw_index, pol_index, :, time_index]
                Trx[ant_index, spw_index, pol_index, :, time_index] = (tempHot* Psamb - Pshot* tempAmb) / (Pshot - Psamb) 
                Tsys[ant_index, spw_index, pol_index, :, time_index] = (Psoff* tempAmb) / (Psamb - Psoff)
                Phot, Pamb, Poff = np.mean(Pshot[chRange]), np.mean(Psamb[chRange]), np.mean(Psoff[chRange])
                chAvgTrx[ant_index, spw_index, pol_index, time_index] = (tempHot* Pamb - Phot* tempAmb) / (Phot - Pamb)
                chAvgTsys[ant_index, spw_index, pol_index, time_index]= (Poff* tempAmb) / (Pamb - Poff)
            #
        #
    #
#
def residTskyTransfer( param, Tamb, secz, Tsky ):
    exp_Tau = np.exp( -param[1]* secz )
    return Tsky - (param[0] + 2.718* exp_Tau  + Tamb* (1.0 - exp_Tau))
#
def residTskyTransfer2( param, Tamb, Tau0, secz, Tsky ):
    exp_Tau = np.exp( -Tau0* secz )
    return Tsky - (param[0] + 2.718* exp_Tau  + Tamb* (1.0 - exp_Tau))
#
param = [5.0, 0.05]
Tau0 = np.zeros([antNum, spwNum, 2])
TantN= np.zeros([antNum, spwNum, 2])
for ant_index in range(antNum):
    for spw_index in range(spwNum):
        for pol_index in range(2):
            Tsky = chAvgTsys[ant_index, spw_index, pol_index] - chAvgTrx[ant_index, spw_index, pol_index]
            fit = scipy.optimize.leastsq(residTskyTransfer, param, args=(tempAmb, 1.0/np.sin(OffEL[ant_index]), Tsky))
            TantN[ant_index, spw_index, pol_index] = fit[0][0]
            Tau0[ant_index, spw_index, pol_index]  = fit[0][1]
        #
    #
#
Tau0 = np.median(Tau0, axis=(0,2))
for spw_index in range(spwNum):
    print 'SPW=%d : Tau(zenith) = %6.4f' % (TDMspw_atmCal[spw_index], Tau0[spw_index])
#
#-------- Antenna-dependent leakage noise
PolList = ['X', 'Y']
param = [5.0]
for ant_index in range(antNum):
    for spw_index in range(spwNum):
        for pol_index in range(2):
            Tsky = chAvgTsys[ant_index, spw_index, pol_index] - chAvgTrx[ant_index, spw_index, pol_index]
            fit = scipy.optimize.leastsq(residTskyTransfer2, param, args=(tempAmb, Tau0[spw_index], 1.0/np.sin(OffEL[ant_index]), Tsky))
            TantN[ant_index, spw_index, pol_index] = fit[0][0]
            print '%s SPW=%d %s : %6.3f K' % (antList[ant_index], TDMspw_atmCal[spw_index], PolList[pol_index], fit[0][0])
        #
    #
#

msmd.close()
"""
#-------- Plots
if BPPLOT:
    #-------- Prepare Plots
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index, figsize = (11, 8))
        figAnt.suptitle(prefix + ' ' + antList[ant_index] + ' Scan = ' + `BPscan[0]`)
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
                plotBP = BP_ant[ant_index, spw_index, pol_index]
                BPampPL.plot( Freq, abs(plotBP), ls='steps-mid', label = 'Pol=' + polName[ppol[pol_index]])
                BPampPL.axis([np.min(Freq), np.max(Freq), 0.0, 1.25* plotMax])
                BPampPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
                BPampPL.yaxis.offsetText.set_fontsize(10)
                BPampPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
                BPphsPL.plot( Freq, np.angle(plotBP), '.', label = 'Pol=' + polName[ppol[pol_index]])
                BPphsPL.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
            #
            BPampPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            BPphsPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
            BPampPL.text( np.min(Freq), 1.1* plotMax, 'SPW=' + `spw[spw_index]` + ' Amp')
            BPphsPL.text( np.min(Freq), 2.5, 'SPW=' + `spw[spw_index]` + ' Phase')
        #
        if PLOTFMT == 'png':
            figAnt.savefig('BP_' + prefix + '_' + antList[ant_index] + '-REF' + antList[0] + '_Scan' + `BPscan[0]` + '.png')
        else :
            figAnt.savefig('BP_' + prefix + '_' + antList[ant_index] + '-REF' + antList[0] + '_Scan' + `BPscan[0]` + '.pdf')
        #
    #
    plt.close('all')
#
np.save(prefix + '-REF' + antList[0] + '.Delay.npy', Delay_ant) 
"""
