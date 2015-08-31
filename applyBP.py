#---- Script for Band-3 Astroholograpy Data
import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
#
#-------- Load BP table
BP_ant = np.load(prefix + '.BPant.npy')     # BP_ant[ANT, SPW, POL, CH]
AntListBP = np.load(prefix + '.Ant.npy').tolist()
#-------- Definitions
antNum = len(refant)
blNum = antNum* (antNum - 1) / 2 
spwNum  = len(spw)
polNum  = len(pol)
#
#-------- Procedures
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile)[refant].tolist()
BPantMap  = antIndex(AntListBP, antList)
blMap = range(blNum)
blInv = [False]* blNum		# True -> inverted baseline
basAntBL = []
endAntBL = []
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(refant[ants[0]], refant[ants[1]])
    basAntBL.append(ants[1])
    endAntBL.append(ants[0])
#
#-------- Each SPW
for ant_index in range(antNum - 1):
    BLwithBaseAnt = np.where( array(basAntBL) == ant_index)[0].tolist() # BL index that have ant_index as the base antenna
#
for ant_index in range(1,antNum):
    BLwithEndAnt = np.where( array(endAntBL) == ant_index)[0].tolist() # BL index that have ant_index as the end antenna
#
for spw_index in range(spwNum):
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], scan)   # Xspec[pol, ch, bl, time]
    Xspec = Xspec[:,:,blMap]                                              # select BL
    #-------- Baseline-based cross power spectra
    Ximag = Xspec.transpose(0,1,3,2).imag * (-2.0* np.array(blInv) + 1.0)   # Complex conjugate for
    Xspec.imag = Ximag.transpose(0,1,3,2)                                   # inversed baselines
    tempSpec = Xspec.transpose(0,3,2,1)
#
    #-------- BP calibration for base antenna
    for ant_index in range(antNum - 1):
        antID = BPantMap[ant_index]
        BLwithBaseAnt = np.where( array(basAntBL) == ant_index)[0].tolist() # BL index that have ant_index as the base antenna
        tempSpec[0,:,BLwithBaseAnt] /= BP_ant[antID, spw_index, 0].conjugate()      # Pol XX  / BP X
        tempSpec[1,:,BLwithBaseAnt] /= BP_ant[antID, spw_index, 0].conjugate()      # Pol XY  / BP X
        tempSpec[2,:,BLwithBaseAnt] /= BP_ant[antID, spw_index, 1].conjugate()      # Pol YX  / BP Y
        tempSpec[3,:,BLwithBaseAnt] /= BP_ant[antID, spw_index, 1].conjugate()      # Pol YY  / BP Y
    #
    #-------- BP calibration for end antenna
    for ant_index in range(1,antNum):
        antID = BPantMap[ant_index]
        BLwithEndAnt = np.where( array(endAntBL) == ant_index)[0].tolist() # BL index that have ant_index as the end antenna
        tempSpec[0,:,BLwithEndAnt] /= BP_ant[antID, spw_index, 0]      # Pol XX  / BP X
        tempSpec[1,:,BLwithEndAnt] /= BP_ant[antID, spw_index, 1]      # Pol XY  / BP Y
        tempSpec[2,:,BLwithEndAnt] /= BP_ant[antID, spw_index, 0]      # Pol YX  / BP X
        tempSpec[3,:,BLwithEndAnt] /= BP_ant[antID, spw_index, 1]      # Pol YY  / BP Y
    #
#
"""
#
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
    Xspec = Xspec[pol]; Xspec = Xspec[:,:,blMap]                            # select polarization and BL
    #-------- Baseline-based cross power spectra
    Ximag = Xspec.transpose(0,1,3,2).imag * (-2.0* np.array(blInv) + 1.0)   # Complex conjugate for
    Xspec.imag = Ximag.transpose(0,1,3,2)                                   # inversed baselines
    tempVis = np.mean(Xspec, axis=3)                                        # tempVis[pol, ch, bl] (time-averaged)
    #-------- Antenna-based bandpass spectra
    for pol_index in range(polNum):
        #-------- Delay Determination and calibration
        if DELAYCAL :
            Delay_ant[:, spw_index, pol_index], delayCalXspec = delayCalSpec(tempVis[pol_index].T, chRange )
        else :
            delayCalXspec = tempVis[pol_index].T
        #
        #-------- Solution (BL -> Ant)
        BP_ant[:,spw_index, pol_index] = np.apply_along_axis(gainComplex, 0, delayCalXspec)
    #
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
"""
