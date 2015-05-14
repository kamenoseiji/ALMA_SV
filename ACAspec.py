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
#-------- Scan Pattern
def onSource( spec, lineRange ):
	power = np.mean(abs(spec[lineRange,:]), 0)
	return  np.where(power > 0.5*(np.median(power) + np.max(power)))[0]
#
def offSource( spec, lineRange ):
	power = np.mean(abs(spec[lineRange,:]), 0)
	return  np.where(power < 0.5*(np.median(power) + np.min(power)))[0]
#
def lineSkip( spec, lineRange):
	temp = spec
	relChRange = lineRange - np.min(lineRange)
	numRange   = float(len(relChRange))
	temp[lineRange] = spec[lineRange[0]-1]* (numRange - relChRange)/(numRange+1)+  spec[lineRange[1]]* (relChRange + 1.0)/(numRange+1)
	return temp
#

msfile = prefix + '.ms'
antList = GetAntName(msfile)
antNum  = len(antList)
polNum  = len(pol)
spwNum  = len(spw_TDM)
chNum, chWid, TDMfreq = GetChNum(msfile, spw_TDM[0])
TDMfreq = 1.0e-9* TDMfreq	# GHz
lineRange = range( min(np.where(TDMfreq > freqRange[0])[0]), max(np.where(TDMfreq < freqRange[1])[0]))
TrxXY  = np.zeros([antNum, spwNum, polNum, chNum])
TsysXY = np.zeros([antNum, spwNum, polNum, chNum])
TrxBB  = np.zeros([antNum, spwNum, polNum])
TsysBB = np.zeros([antNum, spwNum, polNum])
#-------- Tsys Scan
for spw_index in range(spwNum):
	for ant_index in range(antNum):
		tempAmb, tempHot = GetLoadTemp(msfile, ant_index, spw_BB[spw_index])
		for pol_index in range(polNum):
			timeXY, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_TDM[spw_index], TsysScan)
			edge = np.where( diff(np.mean(abs(dataXY), 0)) > 0.02* np.mean(abs(dataXY)) )[0]
			skyRange = range(0, edge[0]-1)
			ambRange = range(edge[0]+2, edge[1]-1)
			hotRange = range(edge[1]+2, len(timeXY))
			refSpec = np.mean(dataXY[:, ambRange].real, axis=1)
			calSpec = np.zeros(dataXY.shape)
			for index in range(dataXY.shape[1]):
				#temp = dataXY[:, index].real / refSpec
				#calSpec[:, index] = Polynomial( temp, np.mean(temp), coeff)* refSpec
				calSpec[:,index] = dataXY[:, index].real / refSpec
			#
			PskyXY = lineSkip(np.median( calSpec[:, skyRange], 1), range(min(lineRange)-2, max(lineRange)+5))
			PambXY = np.median( calSpec[:, ambRange], 1)
			PhotXY = np.median( calSpec[:, hotRange], 1)
			TrxXY[ant_index, spw_index, pol_index, :]  = (tempHot* PambXY - PhotXY* tempAmb) / (PhotXY - PambXY)
			TsysXY[ant_index, spw_index, pol_index, :] = (PskyXY* tempAmb) / (PambXY - PskyXY)
			#
			#timeBB,  dataBB = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_BB[spw_index], TsysScan)
			#timeBB,  dataBB = BB_filter( timeBB, dataBB )
			#edge = np.where( diff(dataBB) > 1000.0 )[0]
			#skyRange = range(0, edge[0]-1)
			#ambRange = range(edge[0]+2, edge[1]-1)
			#hotRange = range(edge[1]+2, len(timeBB))
			#PskyBB,  PambBB,  PhotBB  = np.median(dataBB[skyRange]), np.median(dataBB[ambRange]), np.median(dataBB[hotRange])
			#TrxBB[ant_index, spw_index, pol_index]  = (tempHot* PambBB - PhotBB* tempAmb) / (PhotBB - PambBB)
			#TsysBB[ant_index, spw_index, pol_index] = (PskyBB* tempAmb) / (PambBB - PskyBB)
#			print 'TRX_BB=%5.1f  TRX_XY=%5.1f  TSYS_BB=%5.1f TSYS_XY=%5.1f' % (TrxBB[ant_index, spw_index, pol_index], TrxXY[ant_index, spw_index, pol_index], TsysBB[ant_index, spw_index, pol_index], TsysXY[ant_index, spw_index, pol_index])
		#
	#
#

#-------- OnOff Scan
for spw_index in range(spwNum):
	chNum, chWid, freq = GetChNum(msfile, spw_ACA[spw_index])
	freq = freq* 1.0e-9
	lineRange = range( min(np.where(freq > freqRange[0])[0]), max(np.where(freq < freqRange[1])[0]))
	fig = plt.figure(figsize = (8,11))
	text_sd = '%s Scan=%d SPW=%d' % (prefix, scan, spw_ACA[spw_index])
	for ant_index in range(antNum):
		for pol_index in range(polNum):
			timeACA, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_ACA[spw_index], scan)
			onIndex = onSource( dataXY, lineRange)
			offIndex = offSource( dataXY, lineRange)
			refSpec = np.mean(dataXY[:,offIndex].real, 1)
			calSpec = np.zeros( dataXY.shape )
			for index in range(dataXY.shape[1]):
				#calSpec[:, index] = dataXY[:, index].real / refSpec
				#temp = VanvQ4(dataXY[:, index].real / refSpec)		# assuming VQ4 is done online
				temp = dataXY[:, index].real / refSpec
				calSpec[:, index] = Polynomial( temp, np.mean(temp), coeff)* refSpec
			#
			TsysFreq = interp1d(TDMfreq, TsysXY[ant_index, spw_index, pol_index, :], kind='cubic')
			lineSpec = np.mean(calSpec[:, onIndex], 1)/np.mean(calSpec[:, offIndex], 1) - 1.0
			lineSpec = TsysFreq(freq)* lineSpec
			integTA = np.sum(lineSpec[lineRange])* abs(chWid[0]) * 1.0e-6
			plt.subplot(antNum, polNum, ant_index* polNum + pol_index+1)
			plt.plot(freq, lineSpec, ls='steps-mid', label="ACA Spectrum (corrected)")
			plt.text(0.8*max(freq) + 0.2*min(freq), 0.9*max(lineSpec) + 0.1*min(lineSpec) , antList[ant_index] + ' ' + pol[pol_index], size='x-small')
			text_TA = 'Integ Ta = %5.1f K MHz' % integTA
			plt.text(0.8*max(freq) + 0.2*min(freq), 0.8*max(lineSpec) + 0.2*min(lineSpec) , text_TA, size='x-small')
			print antList[ant_index] + ' ' + pol[pol_index] + ' %6.2f' % integTA
			#plt.xlabel('Frequency [GHz]', fontsize=9)
			#plt.ylabel('Scaled Power [a.u.]', fontsize=9)
		#
	#
	plt.suptitle(prefix + ' SPW=' + `spw_ACA[spw_index]` + ' Scan=' + `scan` + ' RAW')
	plt.savefig(prefix + 'SPW' + `spw_ACA[spw_index]` + 'Scan' + `scan` + 'RAW.pdf', form='pdf')
	plt.close()
#
