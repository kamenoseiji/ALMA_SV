#---- Script for Band-3 Astroholograpy Data
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
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
interval, timeStamp = GetTimerecord(msfile, 0, 0, spw[0], BPscan)
timeNum = len(timeStamp)
#
#-------- Prepare Plots
for spw_index in range(spwNum):
    figSPW = plt.figure(spw_index, figsize = (64, 64))
    #figSPW = plt.figure(spw_index, figsize = (11, 11))
    figSPW.suptitle(prefix + ' SPW = ' + `spw[spw_index]` + ' Scan = ' + `BPscan`)
    figSPW.text(0.45, 0.05, 'Frequency [GHz]')
    figSPW.text(0.05, 0.5, 'Phase', rotation=90)
    figSPW.text(0.95, 0.5, 'Amplitude', rotation=-90)
#
#-------- Plot BP
for spw_index in range(spwNum):
    print ' Loading SPW = ' + `spw[spw_index]`
    figSPW = plt.figure(spw_index)
    chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq = 1.0e-9* Freq  # GHz
    chRange = range(int(round(chNum/chBunch * 0.05)), int(round(chNum/chBunch * 0.95)))
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], BPscan)
    Xspec = Xspec[pol]; Xspec = Xspec[:,:,blMap]
    tempVis = np.mean(Xspec, axis=3)    # tempVis[pol, ch, bl]
    tempAC  = np.mean(Pspec, axis=3)    # tempVis[pol, ch, bl]
    if plotMax == 0.0:
        pMax = 2.0* np.median(abs(tempVis))
        #pMax = np.max(abs(tempVis))
    else:
        pMax = plotMax
    #
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        BLampPL = figSPW.add_subplot( antNum, antNum, antNum* ants[1] + ants[0] + 1 )
        BLphsPL = figSPW.add_subplot( antNum, antNum, antNum* ants[0] + ants[1] + 1 )
        for pol_index in range(polNum):
            plotVis = tempVis[pol_index, :, bl_index]
            BLampPL.plot( Freq, abs(plotVis), ls='steps-mid', label = 'Pol=' + polName[pol_index])
            #BLampPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
            #BLampPL.yaxis.offsetText.set_fontsize(10)
            #BLampPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
            BLphsPL.plot( Freq, np.angle(plotVis), '.', label = 'Pol=' + polName[pol_index])
        #
        #BLampPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
        #BLphsPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        #BLampPL.text( np.min(Freq), 1.1* pMax, 'SPW=' + `spw[spw_index]` + ' Amp')
        #BLphsPL.text( np.min(Freq), 2.5, 'SPW=' + `spw[spw_index]` + ' Phase')
        BLampPL.axis([np.min(Freq), np.max(Freq), 0.0, 1.25* pMax])
        BLampPL.xaxis.set_major_locator(plt.NullLocator())
        BLampPL.yaxis.set_major_locator(plt.NullLocator())
        BLphsPL.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
        if ants[1] == 0:
            BLampPL.set_title( antList[ants[0]] )
            BLphsPL.set_ylabel( antList[ants[0]] )
        #
        if ants[0] < antNum - 1:
            BLphsPL.xaxis.set_major_locator(plt.NullLocator())
        #
        if ants[1] > 0:
            BLphsPL.yaxis.set_major_locator(plt.NullLocator())
        #
    #
    for ant_index in range(antNum):
        ACampPL = figSPW.add_subplot( antNum, antNum, (antNum + 1)* ant_index + 1 )
        ACampPL.patch.set_facecolor('pink')
        for pol_index in range(polNum):
            ACampPL.plot( Freq, abs(tempAC[pol[pol_index], :, ant_index]), ls='steps-mid', label = 'Pol=' + polName[pol_index])
            ACampPL.xaxis.set_major_locator(plt.NullLocator())
            ACampPL.yaxis.set_major_locator(plt.NullLocator())
        #
        if ant_index == 0:
            ACampPL.set_title( antList[0] )
            ACampPL.set_ylabel( antList[0] )
            ACampPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
        #
    #
    figSPW.savefig('BP_' + prefix + '_Scan' + `BPscan` + '_SPW' + `spw[spw_index]` + '.png')
#
#-------- Save CalTables
#np.save(prefix + '.BPant.npy', BP_ant) 
#np.save(prefix + '.Ant.npy', antList) 
plt.close('all')
