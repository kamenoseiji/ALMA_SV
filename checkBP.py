#---- Script for Band-3 Astroholograpy Data
from scipy import stats
import matplotlib.pyplot as plt
execfile(SCR_DIR + 'interferometry.py')

#Real and Imaginary Solution
def antBP(Xspec):
	blNum, chNum = Xspec.shape[0], Xspec.shape[1]
	antNum =  Bl2Ant(blNum)[0]
	#
	BP_ant = np.zeros([antNum, chNum], dtype=complex)
	bl_err = np.std(np.angle(Xspec[:, 1:(chNum-1)]), axis=1) 
	for ch_index in range(chNum):
		BP_ant[:, ch_index] = clcomplex_solve(Xspec[:, ch_index], bl_err)
	#
	return BP_ant
#-------- Definitions
antNum = len(refant)
blNum = antNum* (antNum - 1) / 2 
spwNum  = len(spw)
polNum  = len(pol)
#
#-------- Procedures
timeRange = np.zeros([2])	# Start and End time period 
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
#-------- Prepare Plots
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index, figsize = (8, 11))
    figAnt.suptitle(prefix + ' ' + antList[ant_index] + ' Scan = ' + `BPscan`)
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'Bandpass Amplitude and Phase', rotation=90)
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
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], BPscan)
    Xspec = Xspec[pol]; Xspec = Xspec[:,:,blMap]
    #-------- Baseline-based cross power spectra
    for bl_index in range(blNum):
        if blInv[bl_index]:
            Xspec[:,:,bl_index,:] = Xspec[:,:,bl_index,:].conjugate()
        #
    #
    #-------- Delay Determination and calibration
    tempVis = np.mean(Xspec, axis=3)    # tempVis[pol, ch, bl]
    for pol_index in range(polNum):
        if DELAYCAL :
            Delay_ant[:, spw_index, pol_index], delayCalXspec = delayCalSpec(tempVis[pol_index].T, chRange )
        else :
            delayCalXspec = tempVis[pol_index].T
        #
        BP_ant[:,spw_index, pol_index] = antBP(delayCalXspec)
    #
#
plotMax = np.max(abs(BP_ant))
#-------- Plot BP
for spw_index in range(spwNum):
    for pol_index in range(polNum):
        for ant_index in range(antNum):
            plotBP = BP_ant[ant_index, spw_index, pol_index]
            figAnt = plt.figure(ant_index)
            BPampPL = figAnt.add_subplot( 4, spwNum, spwNum* (2* pol_index    ) + spw_index + 1 )
            BPampPL.plot( Freq, abs(plotBP), ls='steps-mid', label = 'Amp: SPW=' + `spw[spw_index]`)
            BPampPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            BPampPL.set_ylim(0.0, 1.5* plotMax )
            BPphsPL = figAnt.add_subplot( 4, spwNum, spwNum* (2* pol_index + 1) + spw_index + 1 )
            BPphsPL.plot( Freq, np.angle(plotBP), '.', label = 'Phase: SPW=' + `spw[spw_index]`)
            BPphsPL.set_ylim(-math.pi, math.pi)
            BPphsPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        #
    #
#
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index)
    figAnt.savefig(prefix + '_' + antList[ant_index] + '.pdf')
#
plt.close('all')
np.save(prefix + '.BPant.npy', BP_ant) 
np.save(prefix + '.Ant.npy', antList) 
np.save(prefix + '.Delay.npy', Delay_ant) 
