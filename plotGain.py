#---- Script for Band-3 Astroholograpy Data
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
#-------- Load tables
antList, timeStamp, Gain = np.load(antFile), np.load(timeFile), np.load(GainFile)
#-------- Plots
pp = PdfPages('GA_' + prefix + '.pdf')
plotMax = 1.5* np.median(abs(Gain))
antNum = len(antList)
#-------- Prepare Plots
for ant_index in range(antNum):
figAnt = plt.figure(ant_index, figsize = (11, 8))
figAnt.suptitle(prefix + ' ' + antList[ant_index])
figAnt.text(0.45, 0.05, 'Time')
figAnt.text(0.03, 0.45, 'Gain Amplitude and Phase', rotation=90)
#
#-------- Plot Gain
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
            #figAnt.savefig('BP_' + prefix + '_' + antList[antMap[ant_index]] + '-REF' + antList[UseAnt[refantID]] + '_Scan' + `BPscan` + '.pdf')
            figAnt.savefig(pp, format='pdf')
        #
    #
    plt.close('all')
    pp.close()
#
