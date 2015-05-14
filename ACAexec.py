#-------- Script to compare ACA power and BB power
execfile('/users/skameno/Scripts/interferometry.py')
VanvQ4 = loadVanvQ4('/users/skameno/Scripts/VanvQ4.data')	# 4-bit Van Vleck
coeff  = loadAcorrCoeff('/users/skameno/Scripts/ACAVanvCoeff.data')

#-------- Script to compare ACA power and BB power
def BB_filter(timeBB, dataBB):
	index = np.where( abs(dataBB[0]) > 0.0)[0]
	return timeBB[index], abs(dataBB[0,index])
#
#-------- Scan Time in ACA data
def timeMatch( timeBB, timeACA ):
	return np.where( (timeACA < np.max(timeBB)) & (timeACA > np.min(timeBB)) )[0] 
#

msfile = prefix + '.ms'
antList = GetAntName(msfile)
antNum  = len(antList)
polNum  = len(pol)
scanNum  = len(scan)
spwNum  = len(spw_ACA)

#-------- MarsScan 
for scan_index in range(scanNum):
	for spw_index in range(spwNum):
		fig = plt.figure(figsize = (8,11))
		text_sd = '%s Scan=%d BB=%d SPW=%d' % (prefix, scan[scan_index], spw_BB[spw_index], spw_ACA[spw_index])
		for ant_index in range(antNum):
			for pol_index in range(polNum):
				timeBB, dataBB = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_BB[spw_index], scan[scan_index])
				timeBB, dataBB = BB_filter( timeBB, dataBB )
				timeACA, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_ACA[spw_index], scan[scan_index])
				refSpec = dataXY[:,0].real
				calSpec = np.zeros( dataXY.shape )
				for index in range(dataXY.shape[1]):
					calSpec[:, index] = dataXY[:, index].real / refSpec
				#	temp = VanvQ4(dataXY[:, index].real / refSpec)
				#	calSpec[:, index] = Polynomial( temp, 4.57*temp, coeff)* refSpec
				#
				plt.subplot(antNum, polNum, ant_index* polNum + pol_index+1)
				plt.plot( timeBB, abs(dataBB)/abs(dataBB[0]), '.', label="BB power")
				plt.plot(timeACA, np.mean(calSpec, 0)/np.mean(calSpec,0)[0], ls='steps-mid', label="ACA power (corrected)")
				plt.ylim(0.9, 1.1* np.max(abs(dataBB)/abs(dataBB[0])))
				#plt.legend(loc='upper left')
				plt.text(0.7*np.max(timeACA)+0.3*np.min(timeACA), 1.05*np.max(abs(dataBB)/abs(dataBB[0])), 'Ant=' + antList[ant_index] + ', Pol=' + pol[pol_index], size='x-small')
				#plt.xlabel('Time', fontsize=9)
				#plt.ylabel('Scaled Power [a.u.]', fontsize=9)
				#plt.suptitle('Ant=' + antList[ant_index] + ', Pol=' + pol[pol_index])
			#
		#
		plt.suptitle(prefix + 'Scan=' + `scan[scan_index]`)
		#plt.savefig(prefix + 'Scan' + `scan[scan_index]` + '.pdf', form='pdf')
		#plt.close()
	#
#
