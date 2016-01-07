#---- Script for Band-3 Astroholograpy Data
import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
#
#-------- Definitions
antNum = len(refant)
blNum = antNum* (antNum - 1) / 2 
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
spwNum  = len(spw)
scanNum  = len(scan)
ppolNum  = len(ppol)
cpolNum  = len(cpol)
#
#-------- Procedures
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile)[refant]
blMap = range(blNum)
blInv = [False]* blNum		# True -> inverted baseline
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(refant[ants[0]], refant[ants[1]])
#
#-------- Prepare BP and Delay to store
chNum, chWid, Freq = GetChNum(msfile, spw[0])
AvgSpec   = np.ones([spwNum, chNum], dtype=complex)
Delay_ant = np.zeros([antNum, spwNum, (ppolNum + cpolNum)])
XYdelay = np.zeros(spwNum)
BPXY = np.ones([chNum, blNum], dtype=complex)
BPYX = np.ones([chNum, blNum], dtype=complex)
#-------- Loop for SPW
for spw_index in range(spwNum):
    chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq = 1.0e-9* Freq  # GHz
    chRange = range(int(round(chNum/chBunch * 0.05)), int(round(chNum/chBunch * 0.95)))
    XPspec = np.zeros([scanNum, ppolNum, chNum], dtype=complex)
    scanWeight = np.ones([scanNum, ppolNum])
    BP_ant = load( wd + BPtable[spw_index])
    for scan_index in range(scanNum):
        interval, timeStamp = GetTimerecord(msfile, 0, 0, ppol[0], spw[spw_index], scan[scan_index])
        timeNum = len(timeStamp)
        print ' Loading SPW=' + `spw[spw_index]` + ' Scan=' + `scan[scan_index]`,
        #
        #------- Load Cross-power spectrum
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], scan[scan_index])   # Xspec[pol, ch, bl, time]
        Xspec = Xspec[:,:,blMap]
        Xspec = (Xspec.transpose(3,2,0,1) / (BP_ant[ant1].conjugate()* BP_ant[ant0])).transpose(2, 3, 1, 0)
        Ximag = Xspec.transpose(0,1,3,2).imag* (-2.0* np.array(blInv) + 1.0)
        Xspec.imag = Ximag.transpose(0,1,3,2)
        #
        chAvgXX = np.mean(Xspec[0,chRange], axis=0 )
        chAvgYY = np.mean(Xspec[1,chRange], axis=0 )
        #
        #-------- Antenna-based Gain Cal
        GainX = np.apply_along_axis( gainComplex, 0, chAvgXX); scanWeight[scan_index, 0] = np.sum(abs(GainX)**2)
        GainY = np.apply_along_axis( gainComplex, 0, chAvgYY); scanWeight[scan_index, 1] = np.sum(abs(GainY)**2)
        WeightBLX = abs(np.mean(GainX, axis=1))[ant0]* abs(np.mean(GainX, axis=1))[ant1]
        WeightBLY = abs(np.mean(GainY, axis=1))[ant0]* abs(np.mean(GainY, axis=1))[ant1]
        for ch_index in range(chNum):
            Xspec[0, ch_index] = gainCalVis( Xspec[0,ch_index], GainX, GainX)
            Xspec[1, ch_index] = gainCalVis( Xspec[1,ch_index], GainY, GainY)
        #
        #-------- Time Average
        XPspec[scan_index,ppol[0]] = np.mean(np.mean(Xspec, axis=3)[ppol[0]] * WeightBLX, axis=1)  # Time Average and Select Pol
        XPspec[scan_index,ppol[1]] = np.mean(np.mean(Xspec, axis=3)[ppol[1]] * WeightBLY, axis=1)  # Time Average and Select Pol
        print 'Weight = %6.1e %6.1e' % (scanWeight[scan_index, 0], scanWeight[scan_index, 1])
    #
    AvgSpec[spw_index] = np.mean( XPspec.transpose(2,0,1)* scanWeight, axis=(1,2) )
    plt.plot( Freq[chRange], abs(AvgSpec[spw_index, chRange]), ls='steps-mid')
#
np.save(prefix + '-Scan' + `scan[0]` + '.AvSpec.npy', AvgSpec) 
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
