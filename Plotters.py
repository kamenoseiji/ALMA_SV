import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import matplotlib.cm as cm
#
#-------- Plot optical depth
def plotTau(prefix, antList, spwList, secZ, plotTsky, tempAtm, Tau0, TrxFlag, plotMax, PLOTFMT='png'):
    PolList = ['X', 'Y']
    airmass = np.arange( 1.0, 1.25*np.max(secZ), 0.01)
    figTau = plt.figure(0, figsize = (11,8))
    figTau.suptitle(prefix + ' Optical Depth')
    figTau.text(0.45, 0.05, 'Airmass')
    figTau.text(0.03, 0.45, 'Sky Temperature [K]', rotation=90)
    antNum, spwNum  = len(antList), len(spwList)
    for spw_index in range(spwNum):
        for pol_index in range(2):
            TskyPL = figTau.add_subplot(2, spwNum, spwNum* pol_index + spw_index + 1 )
            TskyPL.axis([1.0, 1.25*np.max(secZ), 0.0, plotMax])
            TskyPL.plot( airmass, 2.713* np.exp(-Tau0[spw_index]* airmass) + tempAtm* (1.0 - np.exp(-Tau0[spw_index]* airmass)), '-')
            for ant_index in range(antNum):
                # plotTsky = chAvgTsky[ant_index, spw_index, pol_index] - TantN[ant_index, spw_index, pol_index]
                TskyPL.scatter( secZ[ant_index], plotTsky[ant_index, spw_index, pol_index], s=15* TrxFlag[ant_index, spw_index, pol_index] + 1, color=cm.gist_ncar( float(ant_index) / antNum ), alpha=0.25, label = antList[ant_index])
            #
            text_sd = 'Pol %s Tau(zenith)=%6.4f' % (PolList[pol_index], Tau0[spw_index])
            TskyPL.text(1.01, 0.95* plotMax, text_sd, fontsize='9')
            if pol_index == 0: TskyPL.set_title('SPW ' + `spwList[spw_index]`)
        #
    #
    TskyPL.legend(loc = 'lower right', prop={'size' :7}, numpoints = 1)
    if PLOTFMT == 'png': figTau.savefig('TAU_' + prefix + '.png')
    else : figTau.savefig('TAU_' + prefix + '.pdf')
    plt.close('all')
    return
#
#-------- Plot Tsys and Trx Spectrum
def plotTsys(prefix, antList, ambTime, spwList, TrxList, TskyList, PLOTFMT='png', timeThresh=30.0):
    #-------- Plots for Tsys spectra
    TimePlot = ambTime[np.where( diff(ambTime) < timeThresh)[0].tolist()]
    numTimePlot = len(TimePlot)
    antNum, spwNum  = len(antList), len(spwList)
    plotMax = 1.5 * np.median(TrxList[0] + TskyList[0])
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
            chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* Freq  # GHz
            for scan_index in range(numTimePlot):
                TsysPL = figAnt.add_subplot(numTimePlot, spwNum, spwNum* scan_index + spw_index + 1 )
                timeRange = np.where( abs( ambTime - TimePlot[scan_index]) < timeThresh )[0].tolist()
                timeLabel = qa.time('%fs' % (TimePlot[scan_index]), form='fits')[0]
                for pol_index in range(2):
                    plotTrx  = np.mean(TrxList[AntSpwIndex][pol_index][:, timeRange], axis=1)
                    plotTsys = np.mean(TskyList[AntSpwIndex][pol_index][:, timeRange], axis=1) + plotTrx
                    TsysPL.plot( Freq, plotTsys, ls='steps-mid', label = 'Tsys_Pol=' + PolList[pol_index])
                    TsysPL.plot( Freq, plotTrx,  ls=':', label = 'Trec_Pol=' + PolList[pol_index])
                #
                TsysPL.axis([np.min(Freq), np.max(Freq), 0.0, plotMax])
                if scan_index == 0: TsysPL.set_title('SPW ' + `spwList[spw_index]`)
                if scan_index < numTimePlot - 1: TsysPL.set_xticklabels([])
                if spw_index == 0: TsysPL.text(np.min(Freq), 0.8* plotMax, timeLabel, fontsize='8')
                else: TsysPL.set_yticklabels([])
                #
            #
        #
        if PLOTFMT == 'png': figAnt.savefig('TSYS_' + prefix + '_' + antList[ant_index] + '.png')
        else : figAnt.savefig('TSYS_' + prefix + '_' + antList[ant_index] + '.pdf')
        #
    #
    plt.close('all')
    return
#
