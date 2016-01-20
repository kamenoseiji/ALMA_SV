#-------- Script to compare ACA power and BB power
execfile(SCR_DIR + 'interferometry.py')

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
    for spw_index in range(spwNum):
        tempAmb, tempHot = GetLoadTemp(msfile, ant_index, spw_BB[spw_index])
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
            edge = np.where( diff(dataBB) > 0.1*min(dataBB) )[0]
            skyRange = range(0, edge[0])
            ambRange = range(edge[0]+1, edge[1])
            hotRange = range(edge[1]+1, len(timeBB))
            PskyBB,  PambBB,  PhotBB  = np.mean(dataBB[skyRange]), np.mean(dataBB[ambRange]), np.mean(dataBB[hotRange])
            PskyACA, PambACA, PhotACA = np.mean( powerACA[timeMatch(timeBB[skyRange], timeACA)]), np.mean( powerACA[timeMatch(timeBB[ambRange], timeACA)]), np.mean( powerACA[timeMatch(timeBB[hotRange], timeACA)])
            TrxBB, TrxACA, TsysBB, TsysACA = TRX(PhotBB, PambBB, tempHot, tempAmb), TRX(PhotACA, PambACA, tempHot, tempAmb), TSYS(PskyBB, PambBB, tempAmb), TSYS(PskyACA, PambACA, tempAmb)
            print '%s %d %s: %6.2f %6.2f | %7.4f %7.2f | %7.4f %7.2f | %7.4f %7.2f | %7.4f %7.2f' % (antList[ant_index], spw_index, pol[pol_index], tempAmb, tempHot, PhotBB/PambBB, TrxBB, PhotACA/PambACA, TrxACA, PhotBB/PskyBB, TsysBB, PhotACA/PskyACA, TsysACA)
#
