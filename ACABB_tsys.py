#-------- Script to compare ACA power and BB power
execfile('/Users/kameno/Programs/ALMA_SV/interferometry.py')

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

msfile = prefix + '.ms'
antList = GetAntName(msfile)
antNum  = len(antList)
polNum  = len(pol)
spwNum  = len(spw_ACA)

#-------- TsysScan 
for spw_index in range(spwNum):
	for ant_index in range(antNum):
		for pol_index in range(polNum):
			timeACA, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_ACA[spw_index], TsysScan)
			timeBB,  dataBB = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_BB[spw_index], TsysScan)
			timeBB,  dataBB = BB_filter( timeBB, dataBB )
			refSpec = dataXY[:,0].real
			calSpec = np.zeros( dataXY.shape )
			for index in range(dataXY.shape[1]):
				#temp = VanvQ4(dataXY[:, index].real / refSpec)
				#calSpec[:, index] = Polynomial( temp, 4.57*temp, coeff)* refSpec
				#calSpec[:, index] = Polynomial( temp, temp*2.4, coeff)* refSpec
				calSpec[:, index] = dataXY[:, index].real
			#
			powerACA = abs(np.mean(calSpec, axis=0))
			edge = np.where( diff(dataBB) > 0.1*min(dataBB) )[0]
			skyRange = range(0, edge[0])
			ambRange = range(edge[0]+1, edge[1])
			hotRange = range(edge[1]+1, len(timeBB))
			PskyBB,  PambBB,  PhotBB  = np.mean(dataBB[skyRange]), np.mean(dataBB[ambRange]), np.mean(dataBB[hotRange])
			PskyACA, PambACA, PhotACA = np.mean( powerACA[timeMatch(timeBB[skyRange], timeACA)]), np.mean( powerACA[timeMatch(timeBB[ambRange], timeACA)]), np.mean( powerACA[timeMatch(timeBB[hotRange], timeACA)])
			print 'SPW index=%d  ant_index=%d pol_index=%d' % (spw_index, ant_index, pol_index)
			print 'AMB/SKY: BB=%e  ACA=%e' % (PambBB/PskyBB, PambACA/PskyACA)
			print 'HOT/SKY: BB=%e  ACA=%e' % (PhotBB/PskyBB, PhotACA/PskyACA)
#
