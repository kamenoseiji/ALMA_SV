#---- Script for Band-3 Astroholograpy Data
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
import datetime
#-------- Load tables
DT = []
antList, timeStamp, Delay = np.load(antFile), np.load(timeFile), np.load(DelayFile)
Delay *= 1.0e9
for mjdSec in timeStamp.tolist(): DT.append(datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f'))
DTmin = datetime.datetime.strptime(qa.time('%fs' % (np.min(timeStamp) - 5.0), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f')
DTmax = datetime.datetime.strptime(qa.time('%fs' % (np.max(timeStamp) + 5.0), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f')
#-------- Plots
pp = PdfPages('DL_' + DelayFile + '.pdf')
antNum = len(antList)
polNum = Delay.shape[1]
#-------- Prepare Plots
figDL = plt.figure(figsize = (8, 11))
figDL.suptitle(DelayFile + ' Residual Delay')
figDL.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')))
figDL.text(0.03, 0.7, 'Residual Delay [ns]', rotation=90)
plotMax = 1.1* np.max(abs(Delay))
plotMin = -plotMax
#-------- Plot Gain
for ant_index in range(antNum):
    # plotMax = 1.1* np.max(abs(Gain[ant_index]))
    # plotMin = 0.0* np.min(abs(Gain[ant_index]))
    DLPL = figDL.add_subplot( int(np.ceil(antNum/2.0)), 2, ant_index + 1 )
    for pol_index in range(polNum): DLPL.plot( DT, Delay[ant_index, pol_index], '.')
    DLMed = np.median(Delay[ant_index], axis=1)
    DLPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
    DLPL.yaxis.offsetText.set_fontsize(3)
    DLPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    DLPL.tick_params(labelsize=4)
    DLPL.axis([DTmin, DTmax, plotMin, plotMax], fontsize=3)
    #
    text_sd = '%s : Delay(median) = (%.3f %.3f)' % (antList[ant_index], DLMed[0], DLMed[1])
    DLPL.text( 0.05, 1.02, text_sd, transform=DLPL.transAxes, fontsize=4)
#
figDL.savefig(pp, format='pdf')
plt.close('all')
pp.close()
