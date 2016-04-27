import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
def indexList( refArray, motherArray ):
    IL = []
    for currentItem in refArray:
        IL = IL + np.where( motherArray == currentItem )[0].tolist()
    #
    return IL
#
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
timeAMB = msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT")
timeHOT = msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT")
#-------- Check Scans on-source
print '---Checking on-source time'
onsourceScans = msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE")
scanNum = len(onsourceScans)
#
spw_index = 0; ant_index = 0
timeXY, Pspec = GetPSpec(msfile, ant_index, TDMspw_atmCal[spw_index])
polNum, chNum = Pspec.shape[0], Pspec.shape[1]
ambTime = sort( list(set(timeXY) & set(timeAMB)) ).tolist()
hotTime = sort( list(set(timeXY) & set(timeHOT)) ).tolist()
#
ambTimeIndex = indexList( ambTime, timeXY)
hotTimeIndex = indexList( hotTime, timeXY)
#
AmbSpec = np.zeros([antNum, spwNum, polNum, chNum, len(ambTime)])
HotSpec = np.zeros([antNum, spwNum, polNum, chNum, len(hotTime)])
OnSpec  = np.zeros([antNum, spwNum, polNum, chNum, scanNum])

for ant_index in range(antNum):
    tempAmb, tempHot = GetLoadTemp(msfile, ant_index, TDMspw_atmCal[0])
    for spw_index in range(spwNum):
        timeXY, Pspec = GetPSpec(msfile, ant_index, TDMspw_atmCal[spw_index])
        AmbSpec[ant_index, spw_index] = Pspec[:,:,ambTimeIndex]
        HotSpec[ant_index, spw_index] = Pspec[:,:,hotTimeIndex]
        for scan_index in range(len(onsourceScans)):
            OnTimeIndex = indexList(msmd.timesforscan(onsourceScans[scan_index]), timeXY)
            OnSpec[ant_index, spw_index, :, :, scan_index] = np.median( Pspec[:,:,OnTimeIndex], axis=2 )
        #
    #
#
SPL_amb = UnivariateSpline(ambTime, AmbSpec[5,0,0,20], s=0.01*np.median(AmbSpec[5,0,0,20]))


for sourceIndex in range(numSource):
    scanList = list( set(msmd.scansforfield(sourceList[sourceIndex])) & set(onsourceScans) )
    print sourceList[sourceIndex], scanList
#
msmd.close()
"""
for scan_index in range(len(atmCalScan)):
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, TDMspw_atmCal[spw_index], atmCalScan[scan_index])   # Xspec[pol, ch, bl, time]
    currentAmbTime = list(set(timeStamp) & set(timeSKY))
    ambTime = ambTime + currentAmbTime
    for time_index in range(len(currentAmbTime)):
        ambSpec = Pspec[:,:,:,amb_index]
#
#ambPspec = np.median(Pspec[:,:,:,ambTimeIndex], axis=3)
"""
"""
blMap = range(blNum)
blInv = [False]* blNum		# True -> inverted baseline
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(rllefanteants[0]], refant[ants[1]])
#
#-------- Procedures
#
#-------- Prepare BP and Delay to store
chNum, chWid, Freq = GetChNum(msfile, spw[0])
BP_ant    = np.ones([antNum, spwNum, ppolNum, chNum], dtype=complex)
Delay_ant = np.zeros([antNum, spwNum, (ppolNum + cpolNum)])
XYdelay = np.zeros(spwNum)
BPXY = np.ones([chNum, blNum], dtype=complex)
BPYX = np.ones([chNum, blNum], dtype=complex)
#-------- Loop for SPW
for spw_index in range(spwNum):
    chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq = 1.0e-9* Freq  # GHz
    chRange = range(int(round(chNum/chBunch * 0.05)), int(round(chNum/chBunch * 0.95)))
    XPspec = np.zeros([scanNum, (ppolNum + cpolNum), chNum, blNum], dtype=complex)
    scanWeight = np.ones([scanNum, ppolNum])
    for scan_index in range(scanNum):
        interval, timeStamp = GetTimerecord(msfile, 0, 0, ppol[0], spw[spw_index], BPscan[scan_index])
        timeNum = len(timeStamp)
        print ' Loading SPW=' + `spw[spw_index]` + ' Scan=' + `BPscan[scan_index]`,
        #
        #------- Load Cross-power spectrum
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], BPscan[scan_index])   # Xspec[pol, ch, bl, time]
        Xspec = Xspec[:,:,blMap]
        Ximag = Xspec.transpose(0,1,3,2).imag* (-2.0* np.array(blInv) + 1.0)
        Xreal = Xspec.transpose(0,1,3,2).real
        if cpolNum > 0: # Full polarization pairs
            Xspec[0].imag = Ximag[0].transpose(0,2,1)
            Xspec[1].real = (Xreal[1]*(1.0 - np.array(blInv)) + Xreal[2]* np.array(blInv)).transpose(0,2,1)
            Xspec[1].imag = (Ximag[1]*(1.0 - np.array(blInv)) + Ximag[2]* np.array(blInv)).transpose(0,2,1)
            Xspec[2].real = (Xreal[2]*(1.0 - np.array(blInv)) + Xreal[1]* np.array(blInv)).transpose(0,2,1)
            Xspec[2].imag = (Ximag[2]*(1.0 - np.array(blInv)) + Ximag[1]* np.array(blInv)).transpose(0,2,1)
            Xspec[3].imag = Ximag[3].transpose(0,2,1)
            chAvgXX = np.mean(Xspec[0,chRange], axis=0 )
            chAvgYY = np.mean(Xspec[3,chRange], axis=0 )
        #
        else:   # parallel polarization only
            Xspec[0].imag = Ximag[0].transpose(0,2,1)
            Xspec[1].imag = Ximag[1].transpose(0,2,1)
            chAvgXX = np.mean(Xspec[0,chRange], axis=0 )
            chAvgYY = np.mean(Xspec[1,chRange], axis=0 )
        #
        #-------- Antenna-based Gain Cal
        GainX = np.apply_along_axis( gainComplex, 0, chAvgXX); scanWeight[scan_index, 0] = np.sum(abs(GainX)**2)
        GainY = np.apply_along_axis( gainComplex, 0, chAvgYY); scanWeight[scan_index, 1] = np.sum(abs(GainY)**2)
        if cpolNum > 0: # Full polarization pairs
            for ch_index in range(chNum):
                Xspec[0, ch_index] = gainCalVis( Xspec[0,ch_index], GainX, GainX)
                Xspec[1, ch_index] = gainCalVis( Xspec[1,ch_index], GainX, GainY)
                Xspec[2, ch_index] = gainCalVis( Xspec[2,ch_index], GainY, GainX)
                Xspec[3, ch_index] = gainCalVis( Xspec[3,ch_index], GainY, GainY)
            #
            XCspec = np.mean(Xspec, axis=3)[cpol]                         # Time Average and Select Pol
        else:
            for ch_index in range(chNum):
                Xspec[0, ch_index] = gainCalVis( Xspec[0,ch_index], GainX, GainX)
                Xspec[1, ch_index] = gainCalVis( Xspec[1,ch_index], GainY, GainY)
            #
        #
        #-------- Time Average
        XPspec[scan_index,ppol[0]] = np.mean(Xspec, axis=3)[ppol[0]] * scanWeight[scan_index, 0]  # Time Average and Select Pol
        XPspec[scan_index,ppol[1]] = np.mean(Xspec, axis=3)[ppol[1]] * scanWeight[scan_index, 1]  # Time Average and Select Pol
        if cpolNum > 0: # Full polarization pairs
            XPspec[scan_index,cpol[0]] = np.mean(Xspec, axis=3)[cpol[0]] * scanWeight[scan_index, 0]  # Time Average and Select Pol
            XPspec[scan_index,cpol[1]] = np.mean(Xspec, axis=3)[cpol[1]] * scanWeight[scan_index, 1]  # Time Average and Select Pol
        #
        print 'Weight = %6.1e %6.1e' % (scanWeight[scan_index, 0], scanWeight[scan_index, 1])
    #
    #-------- Antenna-based bandpass spectra
    for pol_index in range(ppolNum):
        #-------- Solution (BL -> Ant)
        BP_ant[:,spw_index, pol_index] = np.apply_along_axis(gainComplex, 0, np.mean(XPspec, axis=0)[ppol[pol_index]].T)
    #
    #-------- Bandpass Correction for Cross-pol
    if cpolNum > 0: # Full polarization pairs
        for bl_index in range(blNum):
            ants = Bl2Ant(bl_index)
            BPXY[:,bl_index] = BP_ant[ants[0], spw_index, 0]* BP_ant[ants[1], spw_index, 1].conjugate()
            BPYX[:,bl_index] = BP_ant[ants[0], spw_index, 1]* BP_ant[ants[1], spw_index, 0].conjugate()
        #
        XC = np.mean( (XCspec[0] / BPXY), axis=1 ) + np.mean( (XCspec[1] / BPYX), axis=1 ).conjugate()
        XYdelay[spw_index], amp = delay_search( XC[chRange] )
        XYdelay[spw_index] *= (float(chNum) / float(len(chRange)))
    #
#
print 'XY delay [sample] = ' + `XYdelay`
if plotMax == 0.0:
    plotMax = 1.5* np.median(abs(BP_ant))
#
#-------- Save CalTables
for spw_index in range(spwNum):
    np.save(prefix + '-REF' + antList[0] + '-SPW' + `spw[spw_index]` + '-BPant.npy', BP_ant[:,spw_index]) 
    if cpolNum > 0: # Full polarization pairs
        np.save(prefix + '-REF' + antList[0] + '-SPW' + `spw[spw_index]` + '-XYdelay.npy', XYdelay[spw_index]) 
    #
#
np.save(prefix + '.Ant.npy', antList) 
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
