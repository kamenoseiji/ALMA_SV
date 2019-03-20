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
        TauPL.tick_params(axis='both', labelsize=6)
        TauPL.plot(freqList[spw_index][chRange], Tau0spec[spw_index][chRange], ls='steps-mid')
        text_sd = 'SPW = %d' % (spwList[spw_index])
        TauPL.text(np.min(freqList[spw_index]), 1.01* plotMax, text_sd, fontsize='8')
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
    #plotMax = np.max(TrxList) + np.max(TskyList)
    plotMax = 1.5* np.median(TrxList) + 1.5* np.median(TskyList)
    PolList = ['X', 'Y']
    #-------- Prepare Plots
    figAnt = plt.figure(figsize = (8, 11))
    figAnt.suptitle(prefix)
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'Tsys (solid) and Trec (dotted) [K]', rotation=90)
    #-------- Plot BP
    for ant_index in range(antNum):
        if ant_index > 0:
            for PL in TsysPL: figAnt.delaxes(PL)
        #
        TsysPL = []
        for spw_index in range(spwNum):
            #AntSpwIndex = ant_index* spwNum + spw_index
            chNum = len(freqList[spw_index]); chRange = range(int(0.05*chNum), int(0.95*chNum))
            for scan_index in range(scanNum):
                currentPL = figAnt.add_subplot(scanNum, spwNum, spwNum* scan_index + spw_index + 1 )
                TsysPL = TsysPL + [currentPL]
                timeLabel = qa.time('%fs' % (atmTime[scan_index]), form='fits')[0]
                for pol_index in range(2):
                    plotTrx  = TrxList[spw_index][pol_index, chRange, ant_index, scan_index]
                    plotTsys = TskyList[spw_index][chRange, ant_index, scan_index] + plotTrx
                    currentPL.plot( freqList[spw_index][chRange], plotTsys, ls='steps-mid', label = 'Tsys Pol '+ PolList[pol_index])
                    currentPL.plot( freqList[spw_index][chRange], plotTrx,  ls=':', label = 'Trec Pol ' + PolList[pol_index])
                #
                currentPL.axis([np.min(freqList[spw_index]), np.max(freqList[spw_index]), 0.0, plotMax])
                currentPL.tick_params(axis='both', labelsize=6)
                if scan_index == 0: currentPL.set_title('SPW ' + `spwList[spw_index]`)
                if scan_index == 0 and spw_index == 0: currentPL.text(1.2* np.min(freqList[0]) - 0.2* np.max(freqList[0]), (0.9 + 0.075*scanNum)*plotMax, antList[ant_index], fontsize='16')
                if scan_index < scanNum - 1: currentPL.set_xticklabels([])
                if spw_index == 0: currentPL.text(np.min(freqList[spw_index]), 0.8* plotMax, timeLabel, fontsize='8')
                else: currentPL.set_yticklabels([])
                #
            #
        #
        plt.show()
        figAnt.savefig(pp, format='pdf')
        #
        #
    #
    for PL in TsysPL: figAnt.delaxes(PL)
    plt.close('all')
    pp.close()
    del(TsysPL)
    del(figAnt)
    return
#
#-------- Plot Bandpass :w
def plotBP(pp, prefix, antList, spwList, BPscan, BPList):
    plotMax = 1.5
    msfile = prefix + '.ms'
    antNum, spwNum = len(antList), len(spwList)
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' Scan ' + `BPscan`)
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'Bandpass Amplitude and Phase', rotation=90)
    #-------- Plot BP
    for ant_index in range(antNum):
        if ant_index > 0:
            for PL in AmpList: figAnt.delaxes(PL)
            for PL in PhsList: figAnt.delaxes(PL)
        #
        AmpList, PhsList = [], []
        for spw_index in range(spwNum):
            chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* Freq  # GHz
            AmpPL = figAnt.add_subplot(2, spwNum, spw_index + 1 )
            PhsPL = figAnt.add_subplot(2, spwNum, spwNum + spw_index + 1 )
            AmpList = AmpList + [AmpPL]
            PhsList = PhsList + [PhsPL]
            for pol_index in range(ppolNum):
                plotBP = BPList[spw_index][ant_index,pol_index]
                AmpPL.plot(Freq, abs(plotBP), ls='steps-mid', label = 'Pol=' + PolList[pol_index])
                PhsPL.plot( Freq, np.angle(plotBP), '.', label = 'Pol=' + PolList[pol_index])
            #
            if spw_index == 0: AmpPL.set_title(antList[ant_index])
            AmpPL.axis([np.min(Freq), np.max(Freq), 0.0, 1.25*plotMax])
            AmpPL.tick_params(axis='both', labelsize=6)
            AmpPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            AmpPL.text( np.min(Freq), 1.1* plotMax, 'SPW=' + `spwList[spw_index]` + ' Amp')
            PhsPL.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
            PhsPL.tick_params(axis='both', labelsize=6)
            PhsPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            PhsPL.text( np.min(Freq), 2.5, 'SPW=' + `spwList[spw_index]` + ' Phase')
        #
        plt.show()
        figAnt.savefig(pp, format='pdf')
    #
    plt.close('all')
    pp.close()
    del(AmpList)
    del(PhsList)
    del(AmpPL)
    del(PhsPL)
    del(figAnt)
    return
#
