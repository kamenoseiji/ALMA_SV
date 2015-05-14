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

def timeRange(refTime, timeBB, timeWidth):
	return np.where( abs(timeBB - refTime) < 0.5*timeWidth)[0]
#
def antIndex(prefix, antList):
	msfile = prefix + '.ms'
	antListInMS = GetAntName(msfile)
	antIndexInMS = []
	for ant_index in range(len(antList)):
		antIndexInMS.append( antListInMS.tolist().index(antList[ant_index]))
	#
	return antIndexInMS
#
def scanPattern(scanTime, scanGap):
	gap = np.where( diff(scanTime) > scanGap )[0]
	scanNum = len(gap) + 1
	ST_index = append(0, gap+1)
	ED_index = append(gap, len(scanTime)-1)
	return ST_index, ED_index
#
def TsysSpec(prefix, TsysScan, TsysSPW):
	msfile = prefix + '.ms'
	antList = GetAntName(msfile)
	antNum  = len(antList)
	polNum  = len(pol)
	chNum, chWid, freq = GetChNum(msfile, TsysSPW); Tsysfreq = freq* 1.0e-9 # GHz
	TrxSpec = np.zeros([antNum, polNum, chNum]); TsysSpec = np.zeros([antNum, polNum, chNum])
	#
	#-------- Get Physical Temperature of loads
	for ant_index in range(antNum):
		tempAmb, tempHot = GetLoadTemp(msfile, ant_index, TsysSPW)
		#
		for pol_index in range(polNum):
			timeXY, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, TsysSPW, TsysScan)
			#
			#-------- Time range of Sky/Amb/Hot
			edge = np.where( diff(timeXY) > 1.0 )[0]
			skyRange = range(0, edge[0])
			ambRange = range(edge[0]+1, edge[1])
			hotRange = range(edge[1]+1, len(timeXY))
			#
			#-------- Calc. Tsys Spectrum
			Psky, Pamb, Phot = np.mean(dataXY[:,skyRange].real, 1), np.mean(dataXY[:,ambRange].real, 1), np.mean(dataXY[:,hotRange].real, 1)
			TrxSpec[ant_index, pol_index]  = (tempHot* Pamb - Phot* tempAmb) / (Phot - Pamb)
			TsysSpec[ant_index, pol_index] = (Psky* tempAmb) / (Pamb - Psky)
			print '%s SPW=%d pol=%s: Trx=%5.1f Tsys=%5.1f' % (antList[ant_index], TsysSPW, pol[pol_index], np.median(TrxSpec[ant_index, pol_index]), np.median(TsysSpec[ant_index, pol_index]))
		#
	#
	return TrxSpec, TsysSpec
#
#
msfile = prefix + '.ms'
antList = GetAntName(msfile)
antNum  = len(antList)
polNum  = len(pol)
scanNum  = len(scan)
spwNum  = len(spw_ACA)
logfile = open(prefix + '_BBLOG.log', 'w')

#-------- Tsys spectrum for specified antennas
Trx, Tsys = TsysSpec( prefix, TsysScan, TsysSPW )
antListInACA = antIndex(prefix, antList)
TsysACA  = np.median(Tsys[antListInACA], axis=2)

#-------- SourceScan
for scan_index in range(scanNum):
	print 'SCAN=%d' % scan[scan_index]
	print 'ANT POL SPW SS  TaACA TaBB Ratio'
	for spw_index in range(spwNum):
		fig = plt.figure(figsize = (11,8))
		fig.text(0.05, 0.45, 'Total Power [scaled by Tsys]', rotation=90)
		fig.text(0.45, 0.05, 'Scan Relative Time [sec]')
		text_sd = '%s Scan=%d BB=%d SPW=%d' % (prefix, scan[scan_index], spw_BB[spw_index], spw_ACA[spw_index])
		for ant_index in range(antNum):
			for pol_index in range(polNum):
				timeBB, dataBB = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_BB[spw_index], scan[scan_index])
				timeBB, dataBB = BB_filter( timeBB, dataBB )
				timeACA, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_ACA[spw_index], scan[scan_index])
				refSpec = dataXY[:,0].real
				calSpec = np.zeros( dataXY.shape )
				for index in range(dataXY.shape[1]):
					calSpec[:, index] = dataXY[:, index].real
				#
				plotBB = abs(dataBB)/abs(dataBB[0])
				plotACA= np.mean(calSpec, 0)/np.mean(calSpec,0)[0]
				timeMatchedBB = np.zeros(len(timeACA))
				ACA_BB_ratio  = np.zeros(len(timeACA))
				for index in range(len(timeACA)):
					indexRange = timeRange(timeACA[index], timeBB, 0.14)
					if( len(indexRange) > 0):
						timeMatchedBB[index] = np.mean(plotBB[indexRange])
					#
				#
				indexAvail = np.where( timeMatchedBB > 0.0)[0]
				ST_index, ED_index = scanPattern(timeACA, 1.0)
				#-------- Loop for subscan
				SSnum = len(SS)
				for ss_index in range(SSnum):
					plt.subplot(antNum, polNum* SSnum, (ant_index* polNum + pol_index) * SSnum+ ss_index + 1)
					indexRange = list(set(indexAvail) & set(range(ST_index[SS[ss_index]], ED_index[SS[ss_index]])))
					relTimeACA = timeACA[indexRange] - np.mean(timeACA[indexRange])
					resultACA, errorACA = fitGauss( relTimeACA, plotACA[indexRange], 1.0e-4* np.ones([len(relTimeACA)]))
					resultBB, errorBB = fitGauss( relTimeACA, timeMatchedBB[indexRange], 1.0e-4* np.ones([len(relTimeACA)]))
					plt.plot(relTimeACA, plotACA[indexRange], '.', label="ACA-power")
					plt.plot(relTimeACA, timeMatchedBB[indexRange], 'o', label="BB-power")
					text_sd = '%s %s %d %d: %5.3f %5.3f %5.3f' % (antList[ant_index], pol[pol_index], spw_ACA[spw_index], SS[ss_index], TsysACA[ant_index, pol_index]*resultACA[0], TsysACA[ant_index, pol_index]*resultBB[0], resultACA[0]/resultBB[0])
					print text_sd
					plotx = arange( np.min(relTimeACA), np.max(relTimeACA), 0.01)
					ploty = resultACA[0]*np.exp(-0.5* ((plotx - resultACA[1])/ resultACA[2])**2 ) + resultACA[3] + resultACA[4]*plotx
					plt.plot( plotx, ploty )
					plt.text(np.min(plotx), np.max(ploty), text_sd, size='x-small')
					ploty = resultBB[0]*np.exp(-0.5* ((plotx - resultBB[1])/ resultBB[2])**2 ) + resultBB[3] + resultBB[4]*plotx
					plt.plot( plotx, ploty )
					logfile.write(text_sd + '\n')
				#
				#
				#for index in indexAvail:
				#	text_sd = '%s %s %d %6.3f %6.3f' % (antList[ant_index], pol[pol_index], index,  TsysACA[ant_index, pol_index]* (plotACA[index] - 1.0), TsysACA[ant_index, pol_index]*(timeMatchedBB[index] - 1.0))
				#	logfile.write(text_sd + '\n')
				#
				#plt.subplot(antNum, polNum, ant_index* polNum + pol_index+1)
				#plt.plot(timeACA[indexAvail], timeMatchedBB[indexAvail], '.', label="BB power")
				#plt.plot(timeACA[indexAvail], plotACA[indexAvail], ls='steps-mid', label="ACA power")
				#plt.plot(timeACA[indexAvail], ACA_BB_ratio, ls='steps-mid', label="ACA/BB ratio")
				#plt.text(0.6*np.max(timeACA)+0.4*np.min(timeACA), 0.9*np.max(plotBB) + 0.1*np.min(plotBB), 'Ant=' + antList[ant_index] + ', Pol=' + pol[pol_index], size='x-small')
				#print '%s, %s, %d, %5.3f, %5.3f' % (antList[ant_index], pol[pol_index], spw_ACA[spw_index], np.max(ACA_BB_ratio-1.0)*100.0, np.min(ACA_BB_ratio-1.0)*100.0)
				#plt.xlabel('Time', fontsize=9)
				#plt.ylabel('Scaled Power [a.u.]', fontsize=9)
				#plt.suptitle('Ant=' + antList[ant_index] + ', Pol=' + pol[pol_index])
			#
		#
		plt.suptitle(prefix + 'Scan=' + `scan[scan_index]` + 'Spw=' + `spw_ACA[spw_index]`)
		plt.savefig(prefix + '.Scan' + `scan[scan_index]` + '.SPW' + `spw_ACA[spw_index]` + '.raw.pdf', form='pdf')
		plt.close()
	#
#
logfile.close()

"""
timeMatchedBB = np.zeros(len(timeACA))
for index in range(len(timeACA)):
	indexRange = timeRange(timeACA[index], timeBB, 0.14)
	if( len(indexRange) > 0):
		timeMatchedBB[index] = np.mean(plotBB[indexRange])
	#
#
indexAvail = np.where( timeMatchedBB > 0.0)[0]

plt.plot(timeMatchedBB[indexAvail], '.')
plt.plot(plotACA[indexAvail], ls='steps-mid')
"""
