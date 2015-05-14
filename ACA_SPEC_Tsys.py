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
#
def lineCalib(spec, Tsys, offRange, onRange):
	offSpec = np.mean(spec[:, offRange], axis=1)
	onSpec  = np.mean(spec[:, onRange], axis=1)
	return Tsys* (onSpec - offSpec)/offSpec
#

#-------- Tsys spectrum for specified antennas
logfile = open(prefixACA + '.SC' + `Scan` + '.SPW' + `spw_ACA` + '.log', 'w')
TDMmsfile = prefixBLC + '.ms'
ACAmsfile = prefixACA + '.ms'
TrxBLC, TsysBLC = TsysSpec( TDMmsfile, pol, TsysScan,  TsysSPWBLC, logfile, False )
TrxACA, TsysACA = TsysSpec( ACAmsfile, pol, TsysScan,  TsysSPWACA, logfile, False )
antListInBLC = antIndex(prefixBLC, antList)
antListInACA = antIndex(prefixACA, antList)
TrxBLC, TsysBLC  = TrxBLC[antListInBLC], TsysBLC[antListInBLC]
TrxACA, TsysACA  = TrxACA[antListInACA], TsysACA[antListInACA]
#
#-------- Plot Tsys spectrum
fig = plt.figure(figsize = (11, 8))
plt.suptitle(prefixBLC + '  Scan ' + `TsysScan`)
chNum, chWid, freq = GetChNum(prefixBLC+'.ms', TsysSPWBLC); Tsysfreq = freq* 1.0e-9	# GHz
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
		text_sd = antList[ant_index] + ' SPW=' + `TsysSPWBLC` + ' Pol=' + pol[pol_index]
		plt.text(0.9*xlim[0]+0.1*xlim[1], 0.1*ylim[0]+0.9*ylim[1], text_sd, size='x-small')
		text_sd = 'BLC: Trx= %5.1f K Tsys= %5.1f K' % (np.median(TrxBLC[ant_index, pol_index]), np.median(TsysBLC[ant_index, pol_index]))
		plt.text(0.9*xlim[0]+0.1*xlim[1], 0.15*ylim[0]+0.85*ylim[1], text_sd, size='x-small')
		text_sd = 'ACA: Trx= %5.1f K Tsys= %5.1f K' % (np.median(TrxACA[ant_index, pol_index]), np.median(TsysACA[ant_index, pol_index]))
		plt.text(0.9*xlim[0]+0.1*xlim[1], 0.20*ylim[0]+0.80*ylim[1], text_sd, size='x-small')
		plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]], fontsize=3)
	#
#
plt.savefig( prefixBLC + '_TsysScan' + `TsysScan` + '.pdf')
plt.close('all')
#
#
fig = plt.figure(figsize = (11, 8))
plt.suptitle(prefixBLC + '  Ta Scan=' + `Scan`)
chNumBLC, chWidBLC, freq = GetChNum(prefixBLC+'.ms', spw_BLC); freqBLC = freq* 1.0e-9	# GHz
chNumACA, chWidACA, freq = GetChNum(prefixACA+'.ms', spw_ACA); freqACA = freq* 1.0e-9	# GHz
print 'ANT Scan SPW Time POL      ACA       BLC Ratio'
for ant_index in range(antNum):
	for pol_index in range(polNum):
		plt.subplot(polNum,  antNum, antNum* pol_index + ant_index + 1)
		if chWid[0] < 0:		# LSB
			smoothTsysBLC = interp1d(Tsysfreq[::-1], TsysBLC[ant_index, pol_index, ::-1])
			smoothTsysACA = interp1d(Tsysfreq[::-1], TsysACA[ant_index, pol_index, ::-1])
		else:					# USB
			smoothTsysBLC = interp1d(Tsysfreq, TsysBLC[ant_index, pol_index])
			smoothTsysACA = interp1d(Tsysfreq, TsysACA[ant_index, pol_index])
		#
		BLCant, ACAant = antListInBLC[ant_index], antListInACA[ant_index]
		timeBLC, dataBLC = GetVisibility(prefixBLC+'.ms', BLCant, BLCant, pol_index, spw_BLC, Scan)
		timeACA, dataACA = GetVisibility(prefixACA+'.ms', ACAant, ACAant, pol_index, spw_ACA, Scan)
		#-------- subScan Pattern Analysis
		ST_index, ED_index = scanPattern(timeBLC, 1.0)
		for SS_index in range(len(onSubScan)):
			offRange = range(ST_index[offSubScan[SS_index]], ED_index[offSubScan[SS_index]])
			onRange  = range(ST_index[onSubScan[SS_index]], ED_index[onSubScan[SS_index]])
			xlim, ylim = plotRange, [-0.5, 30.0]
			chRangeBLC = np.where( (freqBLC > plotRange[0]) & (freqBLC < plotRange[1] ))[0]
			chRangeACA = np.where( (freqACA > plotRange[0]) & (freqACA < plotRange[1] ))[0]
			integTA_BLC = np.zeros(len(onRange))
			integTA_ACA = np.zeros(len(onRange))
			for timeIndex in range(len(onRange)):
				calSpecTDM = lineCalib( dataBLC[chRangeBLC,:].real, smoothTsysBLC(freqBLC[chRangeBLC]), offRange, [onRange[timeIndex]] )
				peakTA_BLC, integTA_BLC[timeIndex] = lineBL(freqBLC[chRangeBLC], calSpecTDM, lineRangeTDM, freeRangeTDM)
				#
				calSpecACA = lineCalib( dataACA[chRangeACA,:].real, smoothTsysACA(freqACA[chRangeACA]), offRange, [onRange[timeIndex]] )
				peakTA_ACA, integTA_ACA[timeIndex] = lineBL(freqACA[chRangeACA], calSpecACA, lineRangeACA, freeRangeACA)
				integTA_BLC[timeIndex] = integTA_BLC[timeIndex] * abs(np.median(chWidBLC)) * 1.0e-6
				integTA_ACA[timeIndex] = integTA_ACA[timeIndex] * abs(np.median(chWidACA)) * 1.0e-6
				text_sd = '%s  %02d  %02d  %02d %s  %8.3f %8.3f %5.3f' % (antList[ant_index], onSubScan[SS_index], TsysSPWBLC, timeIndex, pol[pol_index], integTA_ACA[timeIndex], integTA_BLC[timeIndex], integTA_ACA[timeIndex]/integTA_BLC[timeIndex])
				logfile.write(text_sd + '\n')
				print text_sd
			#
			if SS_index == 11:
				plt.plot(integTA_BLC, ls='steps-mid')
				plt.plot(integTA_ACA, ls='steps-mid')
				plt.legend(('TDM', 'ACA'), 'upper left', prop={'size' :9})
				text_sd = antList[ant_index] + ' SPW=' + `spw_BLC` + ' Pol=' + pol[pol_index]
				plt.text(22, 0.95*np.max(integTA_BLC), text_sd, size='x-small')
				text_sd = 'ACA/TDM = %5.3f ' % (np.sum(integTA_ACA)/np.sum(integTA_BLC))
				plt.text(22, 0.9*np.max(integTA_BLC), text_sd, size='x-small')
			#
		#
	#
#
logfile.close()
plt.savefig( prefixACA + '.' + `onSubScan[0]` + '.ACA.pdf')
plt.close('all')
