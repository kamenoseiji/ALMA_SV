#-------- Script to compare ACA power and BB power
execfile(SCR_DIR + 'interferometry.py')
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
#-------- Script to compare ACA power and BB power
def BB_filter(timeBB, dataBB):
	index = np.where( abs(dataBB[0]) > 0.0)[0]
	return timeBB[index], abs(dataBB[0,index])
#
#-------- Scan Time in ACA data
def timeMatch( timeBB, timeACA ):
	return np.where( (timeACA < np.max(timeBB)) & (timeACA > np.min(timeBB)) )[0] 
#

def timeRange(refTime, timeBB, timeWidth):
	return np.where( abs(timeBB - refTime) < 0.5*timeWidth)[0]
#
def TRX(Phot, Pamb, Thot, Tamb):
    return (Thot* Pamb - Phot* Tamb) / (Phot - Pamb)
#
def TSYS(Psky, Pamb, Tamb):
    return (Psky* Tamb) / (Pamb - Psky)
#
msfile = prefix + '.ms'
antList = GetAntName(msfile)
antNum  = len(antList)
polNum  = len(pol)
spwNum  = len(spw_ACA)

#-------- TsysScan 
print '                         | BB_detector     | Correlator CHAV | BB_detector     | Correlator CHAV'
print 'ANT BB PL  Tamb   Thot   | Phot/Pamb   Trx | Phot/Pamb   Trx | Phot/Psky  Tsys | Phot/Psky  Tsys'
print '-------------------------+-----------------+-----------------+-----------------+----------------'
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index, figsize = (11, 8))
    figAnt.suptitle('Tsys and Trec ' + prefix + ' ' + antList[ant_index] + ' Scan = ' + `TsysScan`)
    figAnt.text(0.5, 0.03, 'Frequency [GHz]')
    figAnt.text(0.03, 0.5, 'Tsys and Trec [K]', rotation=90)
    for spw_index in range(spwNum):
        tempAmb, tempHot = GetLoadTemp(msfile, ant_index, spw_BB[spw_index])
        chNum, chWid, freq = GetChNum(msfile, spw_ACA[spw_index]); freq = freq * 1.0e-9
        TsysPL = figAnt.add_subplot( int(ceil(sqrt(spwNum))), int(ceil(sqrt(spwNum))), spw_index + 1 )
        for pol_index in range(polNum):
            timeACA, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_ACA[spw_index], TsysScan)
            timeBB,  dataBB = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_BB[spw_index], TsysScan)
            timeBB,  dataBB = BB_filter( timeBB, dataBB )
            refSpec = dataXY[:,0].real
            calSpec = dataXY.real
            chNum = calSpec.shape[0]
            chRange = range( int(0.05*chNum), int(0.95*chNum))
            if chNum == 1:
                powerACA = abs(calSpec[0,:])
            else:
                powerACA = abs(np.mean(calSpec[chRange,:], axis=0))
            #
            #-------- States of Hot/Amb/Sky
            edge = np.where( diff(dataBB) > 0.1*min(dataBB) )[0]
            skyRange = range(0, edge[0]); ambRange = range(edge[0]+1, edge[1]); hotRange = range(edge[1]+1, len(timeBB))
            #-------- Time integtartion in each state
            PskyBB,  PambBB,  PhotBB  = np.mean(dataBB[skyRange]), np.mean(dataBB[ambRange]), np.mean(dataBB[hotRange])
            PskyACA, PambACA, PhotACA = np.mean( powerACA[timeMatch(timeBB[skyRange], timeACA)]), np.mean( powerACA[timeMatch(timeBB[ambRange], timeACA)]), np.mean( powerACA[timeMatch(timeBB[hotRange], timeACA)])
            SskyACA, SambACA, ShotACA = np.mean(calSpec[:,timeMatch(timeBB[skyRange], timeACA)], axis=1), np.mean(calSpec[:,timeMatch(timeBB[ambRange], timeACA)], axis=1), np.mean(calSpec[:,timeMatch(timeBB[hotRange], timeACA)], axis=1)
            #-------- Tsys and Trx calculation
            TrxBB, TrxACA, TsysBB, TsysACA = TRX(PhotBB, PambBB, tempHot, tempAmb), TRX(PhotACA, PambACA, tempHot, tempAmb), TSYS(PskyBB, PambBB, tempAmb), TSYS(PskyACA, PambACA, tempAmb)
            TrxACAspec, TsysACAspec = TRX(ShotACA, SambACA, tempHot, tempAmb), TSYS(SskyACA, SambACA, tempAmb)
            print '%s %d %s: %6.2f %6.2f | %7.4f %7.2f | %7.4f %7.2f | %7.4f %7.2f | %7.4f %7.2f' % (antList[ant_index], spw_index, pol[pol_index], tempAmb, tempHot, PhotBB/PambBB, TrxBB, PhotACA/PambACA, TrxACA, PhotBB/PskyBB, TsysBB, PhotACA/PskyACA, TsysACA)
            TsysPL.plot( freq, TsysACAspec, ls='steps-mid', label = 'Tsys Pol=' + pol[pol_index])
            TsysPL.plot( freq, TrxACAspec,  ls='steps-mid', label = 'Trec Pol=' + pol[pol_index])
            #TsysPL.ticklabel_format(axis='y',scilimits=(0,0))
        #
        TsysPL.axis([np.min(freq), np.max(freq), 0.0, 1.25* np.max(TsysACAspec)])
        TsysPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
        TsysPL.yaxis.offsetText.set_fontsize(10)
        TsysPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        TsysPL.text(np.min(freq), 1.1* np.max(TsysACAspec), 'SPW=' + `spw_ACA[spw_index]`)
    #
    figAnt.savefig('Tsys_' + prefix + '_' + antList[ant_index] + '_Scan' + `TsysScan` + '.png')
#
plt.close('all')
