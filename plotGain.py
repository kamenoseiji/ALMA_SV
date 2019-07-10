#---- Script for Band-3 Astroholograpy Data
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
#-------- Load tables
antList, timeStamp, Gain = np.load(antFile), np.load(timeFile), np.load(GainFile)
#-------- Plots
pp = PdfPages('GA_' + GainFile + '.pdf')
antNum = len(antList)
#-------- Prepare Plots
figAmp, figPhs = plt.figure(figsize = (8, 11)), plt.figure(figsize = (8, 11))
figAmp.suptitle(GainFile + ' Gain Amplitude'); figPhs.suptitle(GainFile + ' Gain Phase')
figAmp.text(0.45, 0.05, 'Time'); figPhs.text(0.45, 0.05, 'Time')
figAmp.text(0.03, 0.45, 'Gain Amplitude', rotation=90); figPhs.text(0.03, 0.45, 'Gain Phase [rad]', rotation=90)
#
#-------- Plot Gain
for ant_index in range(antNum):
    plotMax = 2.1* np.median(abs(Gain[ant_index]))
    plotMin = 0.0* np.min(abs(Gain[ant_index]))
    AmpPL = figAmp.add_subplot( int(np.ceil(antNum/2.0)), 2, ant_index + 1 )
    PhsPL = figPhs.add_subplot( int(np.ceil(antNum/2.0)), 2, ant_index + 1 )
    #PhsPL = figAnt.add_subplot( antNum, 1, ant_index + 1 )
    #plotX, plotY = Gain[ant_index, 0], Gain[ant_index, 0]
    #for pol_index in range(2): AmpPL.plot( timeStamp, abs(Gain[ant_index, pol_index]), ls='steps-mid')
    for pol_index in range(2): AmpPL.plot( timeStamp, abs(Gain[ant_index, pol_index]), '.')
    for pol_index in range(2): PhsPL.plot( timeStamp, np.angle(Gain[ant_index, pol_index]), '.')
    ampMed = np.median(abs(Gain[ant_index]), axis=1)
    AmpPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
    AmpPL.yaxis.offsetText.set_fontsize(3)
    PhsPL.yaxis.offsetText.set_fontsize(3)
    AmpPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    AmpPL.tick_params(labelsize=4)
    PhsPL.tick_params(labelsize=4)
    AmpPL.axis([np.min(timeStamp), np.max(timeStamp), plotMin, plotMax], fontsize=3)
    PhsPL.axis([np.min(timeStamp), np.max(timeStamp), -pi, pi], fontsize=3)
    #
    #AmpPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
    #PhsPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    text_sd = '%s : Gain(median) = (%.2f%% %.2f%%)' % (antList[ant_index], 100.0* ampMed[0], 100.0* ampMed[1])
    AmpPL.text( np.min(timeStamp), 0.8* plotMin + 0.2* plotMax, text_sd, fontsize=5)
    PhsPL.text( np.min(timeStamp), -0.5*pi, antList[ant_index], fontsize=5)
#
figAmp.savefig(pp, format='pdf')
figPhs.savefig(pp, format='pdf')
plt.close('all')
pp.close()
