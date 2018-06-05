import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
#
#-------- Plot optical depth
def plotTau(prefix, spwList, freqList, Tau0spec):
    figTau = plt.figure(0, figsize = (11,8))
    figTau.suptitle(prefix + ' Zenith Opacity')
    figTau.text(0.45, 0.05, 'Frequency [GHz]')
    figTau.text(0.03, 0.45, 'Optical Depth', rotation=90)
    plotMax = 1.05* np.max(Tau0spec)
    spwNum = len(spwList)
    for spw_index in range(spwNum):
        chNum = len(freqList[spw_index]); chRange = range(int(0.05*chNum), int(0.95*chNum))
        TauPL = figTau.add_subplot(1, spwNum, spw_index + 1 )
        TauPL.axis([np.min(freqList[spw_index]), np.max(freqList[spw_index]), 0.0, plotMax])
        TauPL.plot(freqList[spw_index][chRange], Tau0spec[spw_index][chRange], ls='steps-mid')
        text_sd = 'SPW = %d' % (spwList[spw_index]); TauPL.text(np.min(freqList[spw_index]), 1.01* plotMax, text_sd, fontsize='8')
    #
    figTau.savefig('TAU_' + prefix + '.pdf')
    plt.close('all')
    return
#
#-------- Plot Tsys and Trx Spectrum
def plotTsys(prefix, antList, spwList, freqList, atmTime, TrxList, TskyList):
    pp = PdfPages('TSYS_' + prefix + '.pdf')
    #-------- Plots for Tsys spectra
    antNum, spwNum, scanNum  = len(antList), len(spwList), len(atmTime)
    plotMax = np.max(TrxList) + np.max(TskyList)
    PolList = ['X', 'Y']
    #-------- Prepare Plots
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index, figsize = (8, 11))
        figAnt.suptitle(prefix + ' ' + antList[ant_index])
        figAnt.text(0.45, 0.05, 'Frequency [GHz]')
        figAnt.text(0.03, 0.45, 'Tsys (solid) and Trec (dotted) [K]', rotation=90)
    #
    #-------- Plot BP
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index)
        for spw_index in range(spwNum):
            AntSpwIndex = ant_index* spwNum + spw_index
            chNum = len(freqList[spw_index]); chRange = range(int(0.05*chNum), int(0.95*chNum))
            for scan_index in range(scanNum):
                TsysPL = figAnt.add_subplot(scanNum, spwNum, spwNum* scan_index + spw_index + 1 )
                timeLabel = qa.time('%fs' % (atmTime[scan_index]), form='fits')[0]
                for pol_index in range(2):
                    plotTrx  = TrxList[AntSpwIndex][pol_index][chRange]
                    plotTsys = TskyList[AntSpwIndex][chRange, scan_index] + plotTrx
                    TsysPL.plot( freqList[spw_index][chRange], plotTsys, ls='steps-mid', label = 'Tsys Pol '+ PolList[pol_index])
                    TsysPL.plot( freqList[spw_index][chRange], plotTrx,  ls=':', label = 'Trec Pol ' + PolList[pol_index])
                #
                TsysPL.axis([np.min(freqList[spw_index]), np.max(freqList[spw_index]), 0.0, plotMax])
                if scan_index == 0: TsysPL.set_title('SPW ' + `spwList[spw_index]`)
                if scan_index < scanNum - 1: TsysPL.set_xticklabels([])
                if spw_index == 0: TsysPL.text(np.min(freqList[spw_index]), 0.8* plotMax, timeLabel, fontsize='8')
                else: TsysPL.set_yticklabels([])
                #
            #
        #
        figAnt.savefig(pp, format='pdf')
        #
    #
    plt.close('all')
    pp.close()
    return
#
