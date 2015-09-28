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
spwNum  = len(spw)
polNum  = len(pol)
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
#-------- Procedures
interval, timeStamp = GetTimerecord(msfile, 0, 0, pol[0], spw[0], BPscan)
timeNum = len(timeStamp)
#
#-------- Prepare BP and Delay to store
chNum, chWid, Freq = GetChNum(msfile, spw[0])
BP_ant    = np.ones([antNum, spwNum, polNum, chNum], dtype=complex)
Delay_ant = np.zeros([antNum, spwNum, polNum])
#-------- Loop for SPW
for spw_index in range(spwNum):
    print ' Loading SPW = ' + `spw[spw_index]`
    chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq = 1.0e-9* Freq  # GHz
    chRange = range(int(round(chNum/chBunch * 0.05)), int(round(chNum/chBunch * 0.95)))
    #------- Load Cross-power spectrum
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], BPscan)   # Xspec[pol, ch, bl, time]
    XPspec = np.mean(Xspec, axis=3)[ppol][:,:,blMap]                         # Time Average and Select Pol
    XCspec = np.mean(Xspec, axis=3)[cpol][:,:,blMap]                         # Time Average and Select Pol
    #-------- Baseline-based cross power spectra
    Ximag = XPspec.imag * (-2.0* np.array(blInv) + 1.0)   # Complex conjugate for
    XPspec.imag = Ximag                                   # inversed baselines
    XYreal = XCspec[0].real* (1.0 - np.array(blInv)) + XCspec[1].real* np.array(blInv)  # Complex conjugate and swap(XY, YX)
    YXreal = XCspec[1].real* (1.0 - np.array(blInv)) + XCspec[0].real* np.array(blInv)  # for inverted baselines
    XYimag = XCspec[0].imag* (1.0 - np.array(blInv)) - XCspec[1].imag* np.array(blInv)  #
    YXimag = XCspec[1].imag* (1.0 - np.array(blInv)) - XCspec[0].imag* np.array(blInv)  #
    XCspec[0] = XYreal + 1.0j* XYimag
    XCspec[1] = YXreal + 1.0j* YXimag
    #-------- Antenna-based bandpass spectra
    for pol_index in range(polNum):
        #-------- Solution (BL -> Ant)
        BP_ant[:,spw_index, pol_index] = np.apply_along_axis(gainComplex, 0, XPspec[pol_index].T)
    #
    #-------- BP Cal for cross-pol
#
if plotMax == 0.0:
    plotMax = 1.5* np.median(abs(BP_ant))
#
#-------- Save CalTables
np.save(prefix + '.BPant.npy', BP_ant) 
np.save(prefix + '.Ant.npy', antList) 
#-------- Plots
if BPPLOT:
    #-------- Prepare Plots
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index, figsize = (11, 8))
        figAnt.suptitle(prefix + ' ' + antList[ant_index] + ' Scan = ' + `BPscan`)
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
            for pol_index in range(polNum):
                plotBP = BP_ant[ant_index, spw_index, pol_index]
                BPampPL.plot( Freq, abs(plotBP), ls='steps-mid', label = 'Pol=' + polName[pol_index])
                BPampPL.axis([np.min(Freq), np.max(Freq), 0.0, 1.25* plotMax])
                BPampPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
                BPampPL.yaxis.offsetText.set_fontsize(10)
                BPampPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
                BPphsPL.plot( Freq, np.angle(plotBP), '.', label = 'Pol=' + polName[pol_index])
                BPphsPL.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
            #
            BPampPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            BPphsPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
            BPampPL.text( np.min(Freq), 1.1* plotMax, 'SPW=' + `spw[spw_index]` + ' Amp')
            BPphsPL.text( np.min(Freq), 2.5, 'SPW=' + `spw[spw_index]` + ' Phase')
        #
        if PLOTFMT == 'png':
            figAnt.savefig('BP_' + prefix + '_' + antList[ant_index] + '_Scan' + `BPscan` + '.png')
        else :
            figAnt.savefig('BP_' + prefix + '_' + antList[ant_index] + '_Scan' + `BPscan` + '.pdf')
        #
    #
    np.save(prefix + '.Delay.npy', Delay_ant) 
    plt.close('all')
#
