#---- Script for Band-3 Astroholograpy Data
import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
#
#-------- Procedures
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile); antNum = len(antList); blNum = antNum* (antNum - 1)/2
#-------- Configure Array
print '---Checking array configulation'
flagAnt = indexList(antFlag, antList)
UseAnt = list(set(range(antNum)) - set(flagAnt)); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
blMap, blInv= range(UseBlNum), [False]* UseBlNum
try:
    refantID = np.where(antList[UseAnt] == refant )[0][0]
except:
    ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
    for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
    timeStamp, UVW = GetUVW(msfile, spw[0], BPscan)
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = bestRefant(uvDist)
#
print '  Use ' + antList[UseAnt[refantID]] + ' as the refant.'
#-------- Baseline Mapping
print '---Baseline Mapping'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
for bl_index in range(UseBlNum):
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
#-------- Bandpass Table
print '---Generating antenna-based bandpass table'
BPList, XYdelayList = [], []
spwNum = len(spw)
for spw_index in spw:
    BP_ant, XYdelay = BPtable(msfile, spw_index, BPscan, blMap, blInv)
    BPList = BPList + [BP_ant]
    XYdelayList = XYdelayList + [XYdelay]
    print 'SPW%d: XY delay [sample] = %f' % (spw_index, XYdelay)
#
BP_ant = np.array(BPList).transpose(1,0,2,3)
XYdelay = np.array(XYdelayList)
ppolNum = BP_ant.shape[2]
PolList = ['X', 'Y']
#
#-------- Save CalTables
np.save(prefix + '.Ant.npy', antList[antMap]) 
for spw_index in range(spwNum):
    np.save(prefix + '-REF' + antList[UseAnt[refantID]] + '-SPW' + `spw[spw_index]` + '-BPant.npy', BP_ant[:,spw_index]) 
    np.save(prefix + '-REF' + antList[UseAnt[refantID]] + '-SPW' + `spw[spw_index]` + '-XYdelay.npy', XYdelay[spw_index]) 
#
#-------- Plots
if BPPLOT:
    plotMax = 1.5* np.median(abs(BP_ant))
    #-------- Prepare Plots
    for ant_index in range(UseAntNum):
        figAnt = plt.figure(ant_index, figsize = (11, 8))
        figAnt.suptitle(prefix + ' ' + antList[antMap[ant_index]] + ' Scan = ' + `BPscan`)
        figAnt.text(0.45, 0.05, 'Frequency [GHz]')
        figAnt.text(0.03, 0.45, 'Bandpass Amplitude and Phase', rotation=90)
    #
    #-------- Plot BP
    for ant_index in range(UseAntNum):
        figAnt = plt.figure(ant_index)
        for spw_index in range(spwNum):
            chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq = 1.0e-9* Freq  # GHz
            BPampPL = figAnt.add_subplot( 2, spwNum, spw_index + 1 )
            BPphsPL = figAnt.add_subplot( 2, spwNum, spw_index + spwNum + 1 )
            for pol_index in range(ppolNum):
                plotBP = BP_ant[ant_index, spw_index, pol_index]
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
            figAnt.savefig('BP_' + prefix + '_' + antList[antMap[ant_index]] + '-REF' + antList[UseAnt[refantID]] + '_Scan' + `BPscan` + '.png')
        else :
            figAnt.savefig('BP_' + prefix + '_' + antList[antMap[ant_index]] + '-REF' + antList[UseAnt[refantID]] + '_Scan' + `BPscan` + '.pdf')
        #
    #
    plt.close('all')
#
