#---- Script for Band-3 Astroholograpy Data
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
import datetime
#-------- Load tables
DT = []
antList, timeStamp, Gain = np.load(antFile), np.load(timeFile), np.load(GainFile)
for mjdSec in timeStamp.tolist(): DT.append(datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f'))
#-------- Plots
pp = PdfPages('GA_' + GainFile + '.pdf')
#pp = PdfPages('DP_' + GainFile + '.pdf')
antNum = len(antList)
#-------- Prepare Plots
figAmp, figPhs = plt.figure(figsize = (8, 11)), plt.figure(figsize = (8, 11))
figAmp.suptitle(GainFile + ' Gain Amplitude'); figPhs.suptitle(GainFile + ' Gain Phase')
figAmp.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d'))); figPhs.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')))
figAmp.text(0.03, 0.7, 'Gain Amplitude = sqrt(correlated flux / SEFD)', rotation=90); figPhs.text(0.03, 0.55, 'Gain Phase [deg]', rotation=90)
#
#plotMax = 1.5* np.median(abs(Gain))
plotMax = 1.1* np.max(abs(Gain))
plotMin = 0.0
#-------- Plot Gain
for ant_index in range(antNum):
    # plotMax = 1.1* np.max(abs(Gain[ant_index]))
    # plotMin = 0.0* np.min(abs(Gain[ant_index]))
    AmpPL = figAmp.add_subplot( int(np.ceil(antNum/2.0)), 2, ant_index + 1 )
    PhsPL = figPhs.add_subplot( int(np.ceil(antNum/2.0)), 2, ant_index + 1 )
    #PhsPL = figAnt.add_subplot( antNum, 1, ant_index + 1 )
    #plotX, plotY = Gain[ant_index, 0], Gain[ant_index, 0]
    #for pol_index in range(2): AmpPL.plot( timeStamp, abs(Gain[ant_index, pol_index]), ls='steps-mid')
    for pol_index in range(polNum): AmpPL.plot( DT, abs(Gain[ant_index, pol_index]), '.')
    for pol_index in range(polNum): PhsPL.plot( DT, np.angle(Gain[ant_index, pol_index])*180.0/pi, '.')
    #PhsPL.plot( DT, np.angle(Gain[ant_index, 1]* Gain[ant_index, 0].conjugate())*180.0/pi, '.')
    ampMed = np.median(abs(Gain[ant_index]), axis=1)
    AmpPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
    AmpPL.yaxis.offsetText.set_fontsize(3)
    PhsPL.yaxis.offsetText.set_fontsize(3)
    AmpPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    AmpPL.tick_params(labelsize=4)
    PhsPL.tick_params(labelsize=4)
    AmpPL.axis([np.min(DT), np.max(DT), plotMin, plotMax], fontsize=3)
    PhsPL.axis([np.min(DT), np.max(DT), -180.0, 180.0], fontsize=3)
    #
    #AmpPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
    #PhsPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    if polNum == 2: text_sd = '%s : Gain(median) = (%.2f%% %.2f%%)' % (antList[ant_index], 100.0* ampMed[0], 100.0* ampMed[1])
    else: text_sd = '%s : Gain(median) = (%.2f%%)' % (antList[ant_index], 100.0* ampMed[0])
    AmpPL.text( 0.05, 1.02, text_sd, transform=AmpPL.transAxes, fontsize=5)
    PhsPL.text( 0.05, 1.02, antList[ant_index], transform=PhsPL.transAxes, fontsize=5)
#
figAmp.savefig(pp, format='pdf')
figPhs.savefig(pp, format='pdf')
plt.close('all')
pp.close()
