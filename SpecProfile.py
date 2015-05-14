#-------- Script to compare ACA power and BB power
execfile('/users/skameno/Scripts/interferometry.py')
#VanvQ4 = loadVanvQ4('/users/skameno/Scripts/VanvQ4.data')	# 4-bit Van Vleck
#coeff  = loadAcorrCoeff('/users/skameno/Scripts/ACAVanvCoeff.data')

#-------- Scan Time in ACA data
def timeMatch( timeBB, timeACA ):
	return np.where( (timeACA < np.max(timeBB)) & (timeACA > np.min(timeBB)) )[0] 
#

def timeRange(refTime, timeBB, timeWidth):
	return np.where( abs(timeBB - refTime) < 0.5*timeWidth)[0]
#
def specInterp(freq, spec):
	return interp1d(freq, spec, kind='cubic')
#
def chRange(freq, range):
	return np.where( (freq > range[0]) & (freq < range[1]))[0]
#
def lineBL(freq, spec, lineRange, BLrange):
	BLch = np.where(((freq > BLrange[0]) & (freq < BLrange[1])) | ((freq > BLrange[2]) & (freq < BLrange[3])))[0]
	SPch = np.where( (freq > lineRange[0]) & (freq < lineRange[1]) )[0]
	corrSpec = spec[SPch] - np.mean(spec[BLch])
	return np.max(corrSpec), np.sum(corrSpec)
#
def scanPattern(scanTime, scanGap):
	gap = np.where( diff(scanTime) > scanGap )[0]
	scanNum = len(gap) + 1
	ST_index = append(0, gap+1)
	ED_index = append(gap, len(scanTime)-1)
	return ST_index, ED_index
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
def TsysSpec(prefix, TsysScan, TsysSPW):
	msfile = prefix + '.ms'
	antList = GetAntName(msfile)
	antNum  = len(antList)
	polNum  = len(pol)
	chNum, chWid, freq = GetChNum(msfile, TsysSPW); Tsysfreq = freq* 1.0e-9	# GHz
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
			# print '%s SPW=%d pol=%s: Trx=%5.1f Tsys=%5.1f' % (antList[ant_index], TsysSPW, pol[pol_index], np.median(TrxSpec[ant_index, pol_index]), np.median(TsysSpec[ant_index, pol_index]))
		#
	#
	return TrxSpec, TsysSpec
#
def lineCalib(spec, Tsys, offRange, onRange):
	offSpec = np.mean(spec[:, offRange], axis=1)
	onSpec  = np.mean(spec[:, onRange], axis=1)
	return Tsys* (onSpec - offSpec)/offSpec
#

#-------- Tsys spectrum for specified antennas
TrxBLC, TsysBLC = TsysSpec( prefixBLC, TsysScan, TsysSPW )
TrxACA, TsysACA = TsysSpec( prefixACA, TsysScan, TsysSPW )
antListInBLC = antIndex(prefixBLC, antList)
antListInACA = antIndex(prefixACA, antList)
TrxBLC, TsysBLC  = TrxBLC[antListInBLC], TsysBLC[antListInBLC]
TrxACA, TsysACA  = TrxACA[antListInACA], TsysACA[antListInACA]

#-------- Plot Tsys spectrum
fig = plt.figure(figsize = (11, 8))
plt.suptitle(prefixBLC + '  Scan ' + `TsysScan`)
chNum, chWid, freq = GetChNum(prefixBLC+'.ms', TsysSPW); Tsysfreq = freq* 1.0e-9	# GHz
antNum, polNum = TrxBLC.shape[0], TrxBLC.shape[1]
for ant_index in range(antNum):
	for pol_index in range(polNum):
		#-------- Plot. Tsys Spectrum
		plt.subplot(polNum,  antNum, antNum* pol_index + ant_index + 1)
		xlim=[np.min(Tsysfreq), np.max(Tsysfreq)]
		ylim=[0.0, 1.2* np.max(TsysACA[ant_index, pol_index, 4:chNum])]
		plt.plot(Tsysfreq, TrxBLC[ant_index, pol_index], ls='steps-mid')
		plt.plot(Tsysfreq, TsysBLC[ant_index, pol_index], ls='steps-mid')
		plt.plot(Tsysfreq, TrxACA[ant_index, pol_index], ls='steps-mid', linestyle='--')
		plt.plot(Tsysfreq, TsysACA[ant_index, pol_index], ls='steps-mid', linestyle='--')
		text_sd = antList[ant_index] + ' SPW=' + `TsysSPW` + ' Pol=' + pol[pol_index]
		plt.text(0.9*xlim[0]+0.1*xlim[1], 0.1*ylim[0]+0.9*ylim[1], text_sd, size='x-small')
		text_sd = 'BLC: Trx= %5.1f K Tsys= %5.1f K' % (np.median(TrxBLC[ant_index, pol_index]), np.median(TsysBLC[ant_index, pol_index]))
		plt.text(0.9*xlim[0]+0.1*xlim[1], 0.15*ylim[0]+0.85*ylim[1], text_sd, size='x-small')
		text_sd = 'ACA: Trx= %5.1f K Tsys= %5.1f K' % (np.median(TrxACA[ant_index, pol_index]), np.median(TsysACA[ant_index, pol_index]))
		plt.text(0.9*xlim[0]+0.1*xlim[1], 0.20*ylim[0]+0.80*ylim[1], text_sd, size='x-small')
		plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]], fontsize=3)
	#
#
#plt.savefig( prefix + 'Scan' + `TsysScan` + '.pdf')
#plt.close('all')
#
#
logfile = open(prefixACA + '.' + `onSubScan[0]` + '.log', 'w')
fig = plt.figure(figsize = (11, 8))
plt.suptitle(prefixBLC + '  Ta Scan=' + `Scan`)
chNumBLC, chWidBLC, freq = GetChNum(prefixBLC+'.ms', spw_BLC); freqBLC = freq* 1.0e-9	# GHz
chNumACA, chWidACA, freq = GetChNum(prefixACA+'.ms', spw_ACA); freqACA = freq* 1.0e-9	# GHz
print 'ANT Scan SPW Time POL      ACA       BLC Ratio'
for ant_index in range(antNum):
	for pol_index in range(polNum):
		plt.subplot(polNum,  antNum, antNum* pol_index + ant_index + 1)
		smoothTsysACA = interp1d(Tsysfreq, TsysACA[ant_index, pol_index])
		ACAant = antListInACA[ant_index]
		timeACA, dataACA = GetVisibility(prefixACA+'.ms', ACAant, ACAant, pol_index, spw_ACA, Scan)
		#-------- subScan Pattern Analysis
		ST_index, ED_index = scanPattern(timeACA, 1.0)
		#for SS_index in range(len(onSubScan)):
		#for SS_index in [11, 32]:
		for SS_index in [SS]:
			offRange = range(ST_index[offSubScan[SS_index]], ED_index[offSubScan[SS_index]])
			onRange  = range(ST_index[onSubScan[SS_index]], ED_index[onSubScan[SS_index]])
			xlim, ylim = plotRange, [-0.5, 30.0]
			chRangeACA = np.where( (freqACA > plotRange[0]) & (freqACA < plotRange[1] ))[0]
			#-------- Scaling for raster scan
			p = np.ones([5, 3])
			x = np.zeros([5])
			y = np.zeros([5])
			weight = np.zeros([5])
			stackACA = np.zeros([len(chRangeACA)])
			for timeIndex in range(14, 19):
				relIndex = timeIndex - 14
				#
				calSpecACA = lineCalib( dataACA[chRangeACA,:].real, smoothTsysACA(freqACA[chRangeACA]), offRange, [onRange[timeIndex]] )
				peakTA_ACA, integTA_ACA = lineBL(freqACA[chRangeACA], calSpecACA, lineRangeACA, freeRangeACA)
				y[relIndex] = integTA_ACA
				x[relIndex] = relIndex - 2.0
				stackACA = stackACA + calSpecACA
			#
			p[:, 1] = x; p[:, 2] = x**2
			param = np.linalg.solve(np.dot(p.T,  p),  np.dot(p.T, y))
			ymax = param[0] - param[1]**2 / (4.0* param[2])
			weight = (param[0] + param[1]* x + param[2]* x**2) / ymax
			calSpecACA = stackACA / np.sum(weight)
			#integTA_BLC = integTA_BLC * np.median(chWidBLC) * 1.0e-6
			#integTA_ACA = integTA_ACA * np.median(chWidACA) * 1.0e-6
			#text_sd = '%s  %02d  %02d  %02d %s  %8.3f %8.3f %5.3f' % (antList[ant_index], onSubScan[SS_index], TsysSPW, timeIndex, pol[pol_index], integTA_ACA, integTA_BLC, integTA_ACA/integTA_BLC)
			#logfile.write(text_sd + '\n')
			#print text_sd
			#
			plt.plot(freqACA[chRangeACA], calSpecACA* scaleFact, ls='steps-mid') 
			plt.axis([345.85, 345.725, -1, 30], fontsize=6)
			plt.legend(('SubScan 11', 'SubScan 32'), 'upper left', prop={'size' :9})
			text_sd = antList[ant_index] + ' SPW=' + `spw_BLC` + ' Pol=' + pol[pol_index]
			plt.text(345.76, 29, text_sd, size='x-small')
			text_sd = 'Peak = %5.3f K' % (np.max(calSpecACA* scaleFact))
			plt.text(345.76, 27, text_sd, size='x-small')
		#
	#
#
logfile.close()
plt.savefig( prefixACA + '_' + `onSubScan[0]` + '_' + `SS_index` + '.ACASpec.pdf')
plt.close('all')
