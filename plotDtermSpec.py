import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
RADDEG = 180.0/ pi
#-------- Load tables
fileNum = len(DSpecList)
DList = []
pdfFileName = 'D_'
for DSpecFile in DSpecList:
    DList = DList + [np.load(DSpecFile)]
    pdfFileName = pdfFileName + DSpecFile
#
#-------- PDF file
pdfFileName = pdfFileName + '.pdf'
pp = PdfPages(pdfFileName)
#-------- Prepare Plots
#-------- Plot D-term
figD = plt.figure(figsize = (11, 8))
figD.suptitle('D-term spectra')
for file_index in range(fileNum):
    DampPL = figD.add_subplot(2, fileNum, file_index + 1)
    DphsPL = figD.add_subplot(2, fileNum, fileNum + file_index + 1)
    Dspec = DList[file_index]
    Freq = Dspec[0]; chNum = len(Freq); chRange = range(int(0.05*chNum), int(0.95*chNum))
    Dx, Dy = Dspec[1] + (0.0 + 1.0j)* Dspec[2], Dspec[3] + (0.0 + 1.0j)* Dspec[4]
    DampPL.grid()
    DampPL.plot( Freq, abs(Dx), ls='steps-mid', label = 'Dx' )
    DampPL.plot( Freq, abs(Dy), ls='steps-mid', label = 'Dy' )
    DampPL.set_title(DSpecList[file_index])
    DampPL.set_ylim(-0.15, 0.15)
    DphsPL.grid()
    DphsPL.plot( Freq, RADDEG* np.angle(Dx), '.', label = 'Dx' )
    DphsPL.plot( Freq, RADDEG* np.angle(Dy), '.', label = 'Dy' )
    DphsPL.set_xlabel('Frequency [GHz]')
    DphsPL.set_ylim(-180.0, 180.0)
    if file_index == 0:
        DampPL.set_ylabel('D-term amp')
        DphsPL.set_ylabel('D-term Phase [deg]')
        DampPL.legend(loc='best')
        DphsPL.legend(loc='best')
#
figD.savefig(pp, format='pdf')
#-------- Plot D-term comparison
if fileNum == 2:
    figD = plt.figure(figsize = (11, 8))
    figD.suptitle('D-term comparison ' + DSpecList[1] + ' / ' + DSpecList[0])
    DampPL = figD.add_subplot(2, 2, 1)
    DphsPL = figD.add_subplot(2, 2, 3)
    DampHS = figD.add_subplot(2, 2, 2)
    DphsHS = figD.add_subplot(2, 2, 4)
    Dx0, Dy0 = DList[0][1] + (0.0 + 1.0j)* DList[0][2], DList[0][3] + (0.0 + 1.0j)* DList[0][4]
    Dx1, Dy1 = DList[1][1] + (0.0 + 1.0j)* DList[1][2], DList[1][3] + (0.0 + 1.0j)* DList[1][4]
    DxR, DyR = Dx1 / Dx0, Dy1 / Dy0
    # Correlation between D-term spectra
    DcorrX = Dx1[chRange].dot(Dx0[chRange].conjugate()) / np.sqrt( Dx0[chRange].dot(Dx0[chRange].conjugate()) * Dx1[chRange].dot(Dx1[chRange].conjugate()) )
    DcorrY = Dy1[chRange].dot(Dy0[chRange].conjugate()) / np.sqrt( Dy0[chRange].dot(Dy0[chRange].conjugate()) * Dy1[chRange].dot(Dy1[chRange].conjugate()) )
    ampRatioDx, ampRatioDy = abs(np.mean(Dx1[chRange])/np.mean(Dx0[chRange])), abs(np.mean(Dy1[chRange])/np.mean(Dy0[chRange]))
    text_sd = '%s %d AmpRatio/PhaseDiff/Corr:  %.2f  %.1f  %.2f  %.2f  %.1f  %.2f' % (antName, BB, ampRatioDx, RADDEG* np.angle(DcorrX), abs(DcorrX), ampRatioDy, RADDEG* np.angle(DcorrY), abs(DcorrY))
    print text_sd
    # plot D-term ratio
    DampPL.grid()
    DampPL.plot( Freq, abs(DxR), ls='steps-mid', label = 'Dx ratio' )
    DampPL.plot( Freq, abs(DyR), ls='steps-mid', label = 'Dy ratio' )
    DampPL.set_ylabel('Amplitude ratio')
    DampPL.legend(loc='best')
    #
    DphsPL.grid()
    DphsPL.plot( Freq, RADDEG* np.angle(DxR), '.', label = 'Dx phase diff' )
    DphsPL.plot( Freq, RADDEG* np.angle(DyR), '.', label = 'Dy phase diff' )
    DphsPL.set_xlabel('Frequency [GHz]')
    DphsPL.set_ylim(-180.0, 180.0)
    DphsPL.set_ylabel('Phase diff [deg]')
    DphsPL.legend(loc='best')
    #
    DampHS.hist( abs(DxR[chRange]), bins=12, alpha=0.25, histtype='stepfilled', label='Dx ratio')
    DampHS.hist( abs(DyR[chRange]), bins=12, alpha=0.25, histtype='stepfilled', label='Dy ratio')
    DampHS.set_xlabel('Amplitude ratio')
    DampHS.legend(loc='best')
    DphsHS.hist( RADDEG* np.angle(DxR[chRange]), bins=12, alpha=0.25, histtype='stepfilled', label='Dx phase diff')
    DphsHS.hist( RADDEG* np.angle(DyR[chRange]), bins=12, alpha=0.25, histtype='stepfilled', label='Dy phase diff')
    DphsHS.set_xlabel('Phase diff [deg]')
    DphsHS.legend(loc='best')
    #
    figD.savefig(pp, format='pdf')
#
plt.close('all')
pp.close()
