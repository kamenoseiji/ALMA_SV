import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
#
#-------- Set Color Map
lineCmap = plt.get_cmap('Set1')
#-------- Plot optical depth
def plotTau(prefix, spwList, freqList, Tau0spec):
    figTau = plt.figure(0, figsize = (11,8))
    figTau.suptitle(prefix + ' Zenith Opacity')
    figTau.text(0.45, 0.05, 'Frequency [GHz]')
    figTau.text(0.03, 0.45, 'Optical Depth', rotation=90)
    spwNum = len(spwList)
    plotMax = 0.01
    for spw_index in range(spwNum): plotMax = max(plotMax, np.max(Tau0spec[spw_index]))
    for spw_index in range(spwNum):
        chNum = len(freqList[spw_index]); chRange = range(int(0.05*chNum), int(0.95*chNum))
        TauPL = figTau.add_subplot(1, spwNum, spw_index + 1 )
        TauPL.axis([np.min(freqList[spw_index]), np.max(freqList[spw_index]), 0.0, 1.05* plotMax])
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
        TsysPL, TsysMax = [], []
        for spw_index in range(spwNum): TsysMax = TsysMax + [1.7* np.percentile(TrxList[spw_index], 75, axis=(0,1,3))[ant_index] + 1.5* np.median(TskyList[spw_index])]
        plotMax = max(TsysMax)
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
#-------- Plot autocorrelation power spectra
def plotAC(prefix, antList, spwList, freqList, AC):
    pp = PdfPages('AC_' + prefix + '.pdf')
    antNum, spwNum, polNum = len(antList), len(spwList), AC[0].shape[2]
    polName = ['X', 'Y']
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' Power Spectra')
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.5, 'Median amplitude and variation [dB]', rotation=90)
    #phaseCalIndex = range(55,60) + range(135,140) + range(160,165)
    #targetIndex = range(55,135) + range(140,160)
    #-------- Plot AC
    for ant_index in range(antNum):
        if ant_index > 0:
            for PL in ACList: figAnt.delaxes(PL)
            for PL in SDList: figAnt.delaxes(PL)
        #
        ACList, SDList = [], []
        ACMAX = 0.0
        for spw_index in range(spwNum): ACMAX = max(ACMAX, np.max(AC[spw_index][ant_index]))
        for spw_index in range(spwNum):
            Freq = freqList[spw_index]
            for pol_index in range(polNum):
                '''
                ACPL = figAnt.add_subplot(2, spwNum, spwNum* pol_index + spw_index + 1)
                ACList = ACList + [ACPL]
                plotAC = np.mean(AC[spw_index][ant_index, targetIndex, pol_index], axis=0) / np.mean(AC[spw_index][ant_index, phaseCalIndex, pol_index], axis=0)
                maxAC, minAC, maxFreq = np.max(plotAC), np.min(plotAC), Freq[np.argmax(plotAC)]
                text_sd = 'Peak = %.3f at %.2f GHz' % (maxAC, maxFreq)
                plotMax, plotMin = max(1.1, maxAC), min(0.99, minAC)
                if spw_index == 0 and pol_index == 0: ACPL.text(1.2* np.min(Freq) - 0.2* np.max(Freq), 1.1*plotMax-0.1*plotMin, antList[ant_index], fontsize='10')
                ACPL.get_xaxis().get_major_formatter().set_useOffset(False)
                #ACPL.xaxis.set_major_locator(ptick.MultipleLocator(0.01))
                #ACPL.xaxis.set_major_formatter(ptick.FormatStrFormatter('%0.2f'))
                ACPL.axis([np.min(Freq), np.max(Freq), plotMin, plotMax])
                ACPL.tick_params(axis='both', labelsize=6)
                #ACPL.set_xticklabels([])
                ACPL.set_xticks(np.arange( int(np.min(Freq)*100)/100.0, int(np.max(Freq)*100)/100.0, 0.01))
                ACPL.text( np.min(Freq), 0.92*plotMax + 0.08*plotMin, 'AC SPW=' + `spwList[spw_index]` + ' Pol-' + polName[pol_index], fontsize=7)
                ACPL.text( np.min(Freq), 0.85*plotMax + 0.15*plotMin, text_sd, fontsize=7)
                ACPL.plot(Freq, plotAC, ls='steps-mid')
                '''
                ACPL = figAnt.add_subplot(4, spwNum, spwNum* pol_index + spw_index + 1)
                SDPL = figAnt.add_subplot(4, spwNum, spwNum* (2+pol_index) + spw_index + 1)
                ACList = ACList + [ACPL]; SDList = SDList + [SDPL]
                plotAC = 10.0* np.log10(np.median(AC[spw_index][ant_index, :, pol_index], axis=0) / ACMAX)
                maxAC, minAC, maxFreq = np.max(plotAC), np.min(plotAC), Freq[np.argmax(plotAC)]
                text_sd = 'Peak = %.1f dB at %.2f GHz' % (maxAC, maxFreq)
                plotMax, plotMin = max(5.0, maxAC), min(-5.0, minAC)
                if spw_index == 0 and pol_index == 0: ACPL.text(1.2* np.min(Freq) - 0.2* np.max(Freq), 1.1*plotMax+0.5, antList[ant_index], fontsize='10')
                ACPL.axis([np.min(Freq), np.max(Freq), plotMin, plotMax])
                ACPL.tick_params(axis='both', labelsize=6)
                ACPL.set_xticklabels([])
                ACPL.text( np.min(Freq), 0.92*plotMax + 0.08*plotMin, 'AC SPW=' + `spwList[spw_index]` + ' Pol-' + polName[pol_index], fontsize=7)
                ACPL.text( np.min(Freq), 0.85*plotMax + 0.15*plotMin, text_sd, fontsize=7)
                ACPL.plot(Freq, plotAC, ls='steps-mid')
                #
                plotSD = 10.0* np.log10(np.std(AC[spw_index][ant_index, :, pol_index], axis=0) / np.median(AC[spw_index][ant_index, :, pol_index]))
                maxSD, minSD, maxFreq = np.max(plotSD), np.min(plotSD), Freq[np.argmax(plotSD)]
                text_sd = '%.1f at %.2f GHz' % (maxSD, maxFreq)
                bgcolor = 'green'
                if maxSD > -30.0: bgcolor = 'orange'
                if maxSD > -20.0: bgcolor = 'red'
                plotMax, plotMin = max(-20.0, maxSD), min(-40.0, minSD)
                SDPL.axis([np.min(Freq), np.max(Freq), 1.05*plotMin - 0.05*plotMax, 1.05*plotMax - 0.05*plotMin])
                SDPL.axhspan(ymin=-30.0, ymax=plotMax, color=bgcolor, alpha=0.1) 
                SDPL.tick_params(axis='both', labelsize=6)
                SDPL.get_xaxis().get_major_formatter().set_useOffset(False)
                if pol_index == 0: SDPL.set_xticklabels([])
                SDPL.text( np.min(Freq), 0.92*plotMax + 0.08*plotMin, 'SD SPW=' + `spwList[spw_index]` + ' Pol-' + polName[pol_index], fontsize=7)
                SDPL.text( np.min(Freq), 0.85* plotMax + 0.15*plotMin, text_sd, fontsize=7)
                SDPL.plot(Freq, plotSD, ls='steps-mid')
            #
        #
        plt.show()
        figAnt.savefig(pp, format='pdf')
    #
    plt.close('all')
    pp.close()
    del(ACList)
    del(SDList)
    del(ACPL)
    #del(SDPL)
    return
#
#-------- Plot Bandpass
def plotBP(pp, prefix, antList, spwList, BPscan, BPList, bunchNum=1):
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
            chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* bunchVec(Freq, bunchNum)  # GHz
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
#-------- Plot XY-phase spectra
def plotXYP(pp, prefix, spwList, XYspec, bunchNum=1):
    spwNum = len(spwList)
    figXYP = plt.figure(figsize = (11, 8))
    figXYP.suptitle(prefix + ' XY Phase')
    figXYP.text(0.45, 0.05, 'Frequency [GHz]')
    figXYP.text(0.03, 0.45, 'XY Phase [deg]', rotation=90)
    for spw_index in range(spwNum):
        spw = spwList[spw_index]
        chNum, chWid, Freq = GetChNum(prefix + '.ms', spw); Freq = 1.0e-9* bunchVec(Freq, bunchNum)  # GHz
        PhsPL = figXYP.add_subplot(1, spwNum, spw_index + 1)
        XYP  = XYspec[spw_index]
        PhsPL.plot( Freq, np.angle(XYP)*180.0/pi, '.', label = 'SPW ' + `spw`)
        PhsPL.axis([np.min(Freq), np.max(Freq), -180.0, 180.0])
        PhsPL.tick_params(axis='both', labelsize=6)
        PhsPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
    #
    plt.show()
    figXYP.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#
#-------- Plot D-term spectra
def plotDSpec(pp, prefix, antList, spwList, DxList, DyList):
    plotMax = 0.12
    antNum, spwNum = len(antList), len(spwList)
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' D-term spectra')
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'D-term Spectra (Real and Imaginary)', rotation=90)
    for ant_index in range(antNum):
        if ant_index > 0:
            for PL in DxPList: figAnt.delaxes(PL)
            for PL in DyPList: figAnt.delaxes(PL)
        #
        DxPList, DyPList = [], []
        for spw_index in range(spwNum):
            DxPL = figAnt.add_subplot( 2, spwNum, spw_index + 1 )
            DyPL = figAnt.add_subplot( 2, spwNum, spw_index + spwNum + 1 )
            DxPList = DxPList + [DxPL]
            DyPList = DyPList + [DyPL]
            #
            plotDx, plotDy = DxList[ant_index, spw_index], DyList[ant_index, spw_index]
            DxPL.plot( FreqList[spw_index], plotDx.real, ls='steps-mid', label = 'reDx')
            DxPL.plot( FreqList[spw_index], plotDx.imag, ls='steps-mid', label = 'imDx')
            DxPL.axis([np.min(FreqList[spw_index]), np.max(FreqList[spw_index]), -plotMax, plotMax])
            DyPL.plot( FreqList[spw_index], plotDy.real, ls='steps-mid', label = 'reDy')
            DyPL.plot( FreqList[spw_index], plotDy.imag, ls='steps-mid', label = 'imDy')
            DyPL.axis([np.min(FreqList[spw_index]), np.max(FreqList[spw_index]), -plotMax, plotMax])
            #
            if spw_index == 0: DxPL.set_title(antList[ant_index])
            DxPL.tick_params(axis='both', labelsize=6)
            DyPL.tick_params(axis='both', labelsize=6)
        #
        DxPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
        DyPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
        plt.show()
        figAnt.savefig(pp, format='pdf')
    #
    plt.close('all'); pp.close()
    del(DxPList); del(DyPList); del(DxPL); del(DyPL); del(figAnt)
    return
#
