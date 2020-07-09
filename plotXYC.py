import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
import datetime
#-------- Load tables
DT = []
timeStamp, XYC, XYV, XYPH, AZEL = np.load(wd + timeFile), np.load(wd + xycFile), np.load(wd + xyvFile), np.load(wd + xypFile), np.load(wd + azelFile)
for mjdSec in timeStamp.tolist(): DT.append(datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f'))
DT = np.array(DT)
PA = np.angle(np.exp((0.0 + 1.0j)* AZEL[3]))
XYV = 0.5*(XYV[0] + XYV[1].conjugate()); XYV *= np.sign(XYV.real)
XYC = 0.5*(XYC[0] + XYC[1].conjugate()); XYC *= np.sign(XYC.real)
#-------- Plots
pp = PdfPages('XYP_' + xycFile + '.pdf')
#-------- Prepare Plots
figXYP = plt.figure(figsize = (8, 11))
figXYP.suptitle(xyvFile + ' XY Phase')
figXYP.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')))
plotMax, plotMin = 180.0, -180.0
#-------- Plot XY Phase
DTPL = figXYP.add_subplot( 3, 1, 1)
CTPL = figXYP.add_subplot( 3, 1, 2)
PAPL = figXYP.add_subplot( 3, 1, 3)
DTPL.grid()
DTPL.plot(DT, XYPH*180.0/pi, '-', label='XY phase correction')
DTPL.plot(DT, np.angle(XYV)*180.0/pi, '.', label='Raw XY correlation')
CTPL.grid()
CTPL.plot(DT, np.arctan2(XYC.imag, XYC.real)*180.0/pi, '.', label='residual after correction')
DTPL.set_ylabel('XY Phase [deg]')
DTPL.legend(loc='best')
CTPL.set_ylabel('XY Phase residual [deg]')
CTPL.legend(loc='best')
#-------- Plot X-feed orientation
PAPL.grid()
PAPL.plot(DT, PA* 180.0/pi, 'o')
PAPL.set_ylabel('X-feed PA [deg]')
figXYP.savefig(pp, format='pdf')
plt.close('all')
pp.close()
