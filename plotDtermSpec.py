import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
RADDEG = 180.0/ pi
#-------- Load tables
fileNum = len(DSpecList)
DList = []
pdfFileName = 'D_comp_' + antName
for DSpecFile in DSpecList:
    DList = DList + [np.load(DSpecFile)]
    # pdfFileName = pdfFileName + DSpecFile
#
Freq = DList[0][0]; chNum = len(Freq); chRange = range(int(0.05*chNum), int(0.95*chNum))
#-------- PDF file
pdfFileName = pdfFileName + '.pdf'
pp = PdfPages(pdfFileName)
#-------- Prepare Plots
#-------- Plot D-term
figD = plt.figure(figsize = (8, 11))
figD.suptitle('D-term spectra ' + antName)
'''
for file_index in range(fileNum):
    DrePL = figD.add_subplot(4, fileNum, file_index + 1)
    DimPL = figD.add_subplot(4, fileNum, fileNum + file_index + 1)
    Dspec = DList[file_index]
    Freq = Dspec[0]; chNum = len(Freq); chRange = range(int(0.05*chNum), int(0.95*chNum))
    Dx, Dy = Dspec[1] + (0.0 + 1.0j)* Dspec[2], Dspec[3] + (0.0 + 1.0j)* Dspec[4]
    #---- Real part
    DrePL.grid()
    DrePL.plot( Freq, Dx.real, ls='steps-mid', label = 'Re(Dx)' )
    DrePL.plot( Freq, Dy.real, ls='steps-mid', label = 'Re(Dy)' )
    DrePL.set_title(DSpecList[file_index], fontsize=7)
    DrePL.set_ylim(-0.05, 0.05)
    #---- Imaginary part
    DimPL.grid()
    DimPL.plot( Freq, Dx.imag, ls='steps-mid', label = 'Im(Dx)' )
    DimPL.plot( Freq, Dy.imag, ls='steps-mid', label = 'Im(Dy)' )
    DimPL.set_ylim(-0.05, 0.05)
    if file_index == 0:
        DrePL.set_ylabel('D-term real part')
        DimPL.set_ylabel('D-term imag parg')
        DrePL.legend(loc='best', prop={"size":7})
        DimPL.legend(loc='best', prop={"size":7})
    #
#
'''
#figD.savefig(pp, format='pdf')
#-------- Plot D-term comparison
#figD = plt.figure(figsize = (11, 8))
#figD.suptitle('D-term comparison ' + DSpecList[1] + ' / ' + DSpecList[0])
Dx0, Dy0 = DList[0][1] + (0.0 + 1.0j)* DList[0][2], DList[0][3] + (0.0 + 1.0j)* DList[0][4]
Dx1, Dy1 = DList[1][1] + (0.0 + 1.0j)* DList[1][2], DList[1][3] + (0.0 + 1.0j)* DList[1][4]
DxPL = figD.add_subplot(4, 2, 1)
DyPL = figD.add_subplot(4, 2, 2)
DxPL.grid()
DxPL.plot( Freq, Dx0.real, ls='steps-mid', color='black', label = 'ReDx' + labelList[0] )
DxPL.plot( Freq, Dx0.imag, linestyle=':', drawstyle='steps-mid', color='black', label = 'ImDx' + labelList[0] )
DxPL.plot( Freq, Dx1.real, ls='steps-mid', color='red', label = 'ReDx' + labelList[1] )
DxPL.plot( Freq, Dx1.imag, linestyle=':', drawstyle='steps-mid', color='red', label = 'ImDx' + labelList[1] )
DxPL.legend(loc='best', prop={"size":7})
DyPL.grid()
DyPL.plot( Freq, Dy0.real, ls='steps-mid', color='black', label = 'ReDy' + labelList[0] )
DyPL.plot( Freq, Dy0.imag, linestyle=':', drawstyle='steps-mid', color='black', label = 'ImDy' + labelList[0] )
DyPL.plot( Freq, Dy1.real, ls='steps-mid', color='red', label = 'ReDy' + labelList[1] )
DyPL.plot( Freq, Dy1.imag, linestyle=':', drawstyle='steps-mid', color='red', label = 'ImDy' + labelList[1] )
DyPL.legend(loc='best', prop={"size":7})
#-------- comparison
CxPL = figD.add_subplot(4, 2, 3)
CyPL = figD.add_subplot(4, 2, 4)
CxPL.grid()
CxPL.plot( Dx0.real, Dx1.real, '.', color='black', label = 'ReDx')
CxPL.plot( Dx0.imag, Dx1.imag, '.', color='red', label = 'ImDx')
CxPL.legend(loc='best', prop={"size":7})
CyPL.grid()
CyPL.plot( Dy0.real, Dy1.real, '.', color='black', label = 'ReDy')
CyPL.plot( Dy0.imag, Dy1.imag, '.', color='red', label = 'ImDy')
CyPL.legend(loc='best', prop={"size":7})
#
#---- Difference
diffDx, diffDy = Dx1 - Dx0, Dy1 - Dy0
# plot D-term diff spectra
DrePL = figD.add_subplot(4, 2, 5)
DimPL = figD.add_subplot(4, 2, 6)
DrePL.grid()
DrePL.plot( Freq, diffDx.real, ls='steps-mid', label = 'Re(diff Dx)' )
DrePL.plot( Freq, diffDy.real, ls='steps-mid', label = 'Re(diff Dy)' )
DrePL.set_ylabel('D-term diff.')
DrePL.legend(loc='best', prop={"size":7})
#
DimPL.grid()
DimPL.plot( Freq, diffDx.imag, ls='steps-mid', label = 'Im(diff Dx)' )
DimPL.plot( Freq, diffDy.imag, ls='steps-mid', label = 'Im(diff Dy)' )
DimPL.legend(loc='best', prop={"size":7})
DimPL.set_xlabel('Frequency [GHz]', fontsize=9)
#---- plot Histogram
DreHS = figD.add_subplot(4, 2, 7)
DimHS = figD.add_subplot(4, 2, 8)
DreHS.hist( diffDx[chRange].real, bins=12, alpha=0.25, histtype='stepfilled', label='Re(diff Dx)')
DreHS.hist( diffDy[chRange].real, bins=12, alpha=0.25, histtype='stepfilled', label='Re(diff Dy)')
DreHS.set_ylabel('Number')
DreHS.set_xlabel('Real-part diff')
DreHS.set_xlim(-0.05, 0.05)
DreHS.legend(loc='best', prop={"size":7})
#
DimHS.hist( diffDx[chRange].imag, bins=12, alpha=0.25, histtype='stepfilled', label='Im(diff Dx)')
DimHS.hist( diffDy[chRange].imag, bins=12, alpha=0.25, histtype='stepfilled', label='Im(diff Dy)')
DimHS.set_xlabel('Imag-part diff')
DimHS.set_xlim(-0.05, 0.05)
DimHS.legend(loc='best', prop={"size":7})
#
figD.savefig(pp, format='pdf')
plt.close('all')
pp.close()
