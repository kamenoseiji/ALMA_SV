#---- Script for Band-3 Astroholograpy Data
from scipy import stats
import matplotlib.pyplot as plt
execfile(SCR_DIR + 'interferometry.py')
#-------- Load BP and Delay
if BPCAL:
    try: 
        BP_ant = np.load( wd + BPprefix + '.BPant.npy' )
    except:
        BPCAL = False
    #
#
def gainComplex( vis ):
    return(clcomplex_solve(vis, 1.0e-8/abs(vis)))
#
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
interval, timeStamp = GetTimerecord(msfile, 0, 0, pol[0], spw[0], TGscan)
timeNum = len(timeStamp)
#
#-------- Prepare Plots
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index, figsize = (8, 11))
    figAnt.suptitle(prefix + ' ' + antList[ant_index] + ' Scan = ' + `TGscan`)
    figAnt.text(0.45, 0.05, 'MJD [sec]')
    figAnt.text(0.03, 0.45, 'Gain Amplitude and Phase', rotation=90)
#
Gain_ant = np.ones([antNum, spwNum, polNum, timeNum], dtype=complex)
#-------- Baseline-based bandpass
if BPCAL:
    BP_bl = np.ones([spwNum, polNum, chNum, blNum], dtype=complex)
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        BP_bl[:,:,:,bl_index] = BP_ant[ants[0]]* BP_ant[ants[1]].conjugate()
    #
#
#-------- Loop for SPW
for spw_index in range(spwNum):
    print ' Loading SPW = ' + `spw[spw_index]`
    chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq = 1.0e-9* Freq  # GHz
    chRange = range(int(round(chNum/chBunch * 0.05)), int(round(chNum/chBunch * 0.95)))
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], TGscan)
    Xspec = Xspec[pol]; Xspec = Xspec[:,:,blMap]
    #-------- Baseline-based cross power spectra
    Ximag = Xspec.transpose(0,1,3,2).imag * (-2.0* np.array(blInv) + 1.0)
    Xspec.imag = Ximag.transpose(0,1,3,2)
    #-------- Delay Determination and calibration
    if BPCAL:
        tempVis = np.mean( Xspec.transpose(3,0,1,2) / BP_bl[spw_index], axis=2).transpose(1,2,0)
    else:
        tempVis = np.mean(Xspec, axis=1)    # tempVis[pol, ch, bl]
    #
    for pol_index in range(polNum):
        vis_bl  = tempVis[pol_index]
        Gain_ant[:, spw_index, pol_index] = np.apply_along_axis(gainComplex, 0, vis_bl )
    #
#
plotMax = np.max(abs(Gain_ant))
#-------- Plot Gain
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index)
    for pol_index in range(polNum):
        GAampPL = figAnt.add_subplot( 4, 1, 2* pol_index + 1 )
        GAphsPL = figAnt.add_subplot( 4, 1, 2* pol_index + 2 )
        for spw_index in range(spwNum):
            plotGA = Gain_ant[ant_index, spw_index, pol_index]
            GAampPL.plot( timeStamp, abs(plotGA), ls='steps-mid', label = 'SPW=' + `spw[spw_index]`)
            GAphsPL.plot( timeStamp, np.angle(plotGA), '.', label = 'SPW=' + `spw[spw_index]`)
        #
        GAampPL.set_ylim(0.0, 1.25* plotMax )
        GAphsPL.set_ylim(-math.pi, math.pi)
        GAampPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
        GAphsPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
        GAampPL.text( np.min(timeStamp), 1.1* plotMax, polName[pol_index] + ' Amp')
        GAphsPL.text( np.min(timeStamp), 2.5, polName[pol_index] + ' Phase')
    #
    figAnt.savefig('GA_' + prefix + '_' + antList[ant_index] + '.pdf')
#
np.save(prefix + '.GA.npy', Gain_ant) 
plt.close('all')
