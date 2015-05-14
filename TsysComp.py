#-------- Script to compare ACA power and BB power
execfile('/users/skameno/Scripts/interferometry.py')
VanvQ3, Qeff = loadVanvQ3('/users/skameno/Scripts/ACAPowerCoeff.data')	# 3-bit Van Vleck
VanvQ4 = loadVanvQ4('/users/skameno/Scripts/VanvQ4.data')	# 4-bit Van Vleck
#coeff  = loadAcorrCoeff('/users/skameno/Scripts/ACAVanvCoeff.data')

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
"""
def TsysSpec(prefix, TsysScan, TsysSPW, logfile, vanvSW):
	msfile = prefix + '.ms'
	antList = GetAntName(msfile)
	antNum  = len(antList)
	polNum  = len(pol)
	chNum, chWid, freq = GetChNum(msfile, TsysSPW); Tsysfreq = freq* 1.0e-9 # GHz
	TrxSp = np.zeros([antNum, polNum, chNum]); TsysSp = np.zeros([antNum, polNum, chNum])
	if (chNum == 128):
		chRange = range(10, 118)
	else:
		chRange = range(8, 116)
	#
	#-------- Get Physical Temperature of loads
	text_sd =  '%s SPW=%d SCAN=%d' % (prefix, TsysSPW, TsysScan)
	print text_sd
	logfile.write(text_sd + '\n')
	for ant_index in range(antNum):
		tempAmb, tempHot = GetLoadTemp(msfile, ant_index, TsysSPW)
		#
		for pol_index in range(polNum):
			timeXY, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, TsysSPW, TsysScan)
			#
			#-------- Time range of Sky/Amb/Hot
			edge = np.where( diff(timeXY) > 1.0 )[0]
			skyRange = range(2, edge[0])
			ambRange = range(edge[0]+3, edge[1])
			hotRange = range(edge[1]+3, len(timeXY)-1)
			#
			if vanvSW:
				#-------- Van Vleck Correction for 3-bit in Lag Domain
				#refZeroLag = np.sum(dataXY[:, skyRange].real) / len(skyRange)
				#refZeroLag = np.sum(dataXY[:, ambRange].real) / len(ambRange)
				#refZeroLag = np.mean(dataXY[:, ambRange].real)
				#print '%5.3f ' % (refZeroLag *  2**(-36))
				#refZeroLag = 2**36
				#refZeroLag = np.mean(dataXY[:, skyRange].real)* 3.2
				refZeroLag = np.mean(dataXY[:, ambRange].real)
				for index in range(dataXY.shape[1]):
					#temp = fft(dataXY[:, index].real)
					#ZeroLag = temp.real[0]
					#VanvCorrect = Qeff( ZeroLag / refZeroLag)
					##print '%5.3f %5.3f %5.3f %5.3f ' % (ZeroLag, ZeroLag/refZeroLag, VanvQ3(ZeroLag/refZeroLag), VanvCorrect)
					##ReTemp, ImTemp = temp.real / VanvCorrect, temp.imag / VanvCorrect;  ReTemp[0] = ZeroLag
					##temp = ReTemp + (0+1j)*ImTemp
					##temp = temp* VanvQ3(ZeroLag/refZeroLag)
					##temp[0] = temp[0] / VanvCorrect
					#temp = temp / VanvCorrect
					#dataXY[:,index] = ifft(temp).real
					#temp = dataXY[:, index].real
					#ZeroLag = np.mean(temp)
					#VanvCorrect = Qeff( ZeroLag / refZeroLag)
					#dataXY[:,index] = temp / VanvCorrect
				#
				#
			#
			#-------- Calc. Tsys Spectrum
			Psky, Pamb, Phot = np.mean(dataXY[:,skyRange].real, 1), np.mean(dataXY[:,ambRange].real, 1), np.mean(dataXY[:,hotRange].real, 1)
			TrxSp[ant_index, pol_index]  = (tempHot* Pamb - Phot* tempAmb) / (Phot - Pamb)
			TsysSp[ant_index, pol_index] = (Psky* tempAmb) / (Pamb - Psky)
			#TrxMean  = (tempHot* np.mean(Pamb[chRange]) - np.mean(Phot[chRange])* tempAmb) / (np.mean(Phot[chRange]) - np.mean(Pamb[chRange]))
			#TsysMean = (np.mean(Psky[chRange])* tempAmb) / (np.mean(Pamb[chRange]) - np.mean(Psky[chRange]))
			#print '%s %s: %5.1f %5.1f' % (antList[ant_index], pol[pol_index], np.median(TrxSp[ant_index, pol_index]), np.median(TsysSp[ant_index, pol_index]))
			#text_sd = '%s %s: %5.1f %5.1f' % (antList[ant_index], pol[pol_index], TrxMean, TsysMean)
			text_sd = '%s %s: %5.1f %5.1f' % (antList[ant_index], pol[pol_index], np.median(TrxSp[ant_index, pol_index]), np.median(TsysSp[ant_index, pol_index]))
			print text_sd
			logfile.write(text_sd + '\n')
		#
	#
	return TrxSp, TsysSp
#
def PowerSpec(prefix, TsysScan, TsysSPW, vanvSW):
	msfile = prefix + '.ms'
	antList = GetAntName(msfile)
	antNum  = len(antList)
	polNum  = len(pol)
	chNum, chWid, freq = GetChNum(msfile, TsysSPW); Tsysfreq = freq* 1.0e-9 # GHz
	PskySpec = np.zeros([antNum, polNum, chNum]); PambSpec = np.zeros([antNum, polNum, chNum]); PhotSpec = np.zeros([antNum, polNum, chNum])
	#-------- Get Physical Temperature of loads
	for ant_index in range(antNum):
		tempAmb, tempHot = GetLoadTemp(msfile, ant_index, TsysSPW)
		#
		for pol_index in range(polNum):
			timeXY, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, TsysSPW, TsysScan)
			#
			#-------- Time range of Sky/Amb/Hot
			edge = np.where( diff(timeXY) > 1.0 )[0]
			skyRange = range(2, edge[0])
			ambRange = range(edge[0]+3, edge[1])
			hotRange = range(edge[1]+3, len(timeXY)-1)
			#
			if vanvSW:
				#-------- Van Vleck Correction for 3-bit in Lag Domain
				#refZeroLag = np.sum(dataXY[:, skyRange].real) / len(skyRange)
				refZeroLag = np.sum(dataXY[:, skyRange].real) / len(skyRange)* (2.9/2.4)**2
				for index in range(dataXY.shape[1]):
					temp = fft(dataXY[:, index].real)
					ZeroLag = temp.real[0]
					VanvCorrect = Qeff(ZeroLag / refZeroLag)
					ReTemp, ImTemp = temp.real / VanvCorrect, temp.imag / VanvCorrect;  ReTemp[0] = ZeroLag
					dataXY[:,index] = ifft(ReTemp + (0 + 1.0j)* ImTemp ).real
				#
			#
			#-------- Calc. Tsys Spectrum
			PambSpec[ant_index, pol_index] = np.mean(dataXY[:,ambRange].real, 1)
			scaleFact = tempAmb / np.mean(PambSpec[ant_index, pol_index])
			PambSpec[ant_index, pol_index] = scaleFact* PambSpec[ant_index, pol_index] 
			PskySpec[ant_index, pol_index] = scaleFact* np.mean(dataXY[:,skyRange].real, 1)
			PhotSpec[ant_index, pol_index] = scaleFact* np.mean(dataXY[:,hotRange].real, 1)
		#
	#
	return Tsysfreq, PskySpec, PambSpec, PhotSpec
#
"""
TDMmsfile = TDMprefix + '.ms'
ACAmsfile = ACAprefix + '.ms'
logfile = open(TDMprefix + '_' + `TsysScan` + '_' + `TsysSPW` + '_TsysLOG.log', 'w')
antList = GetAntName(TDMmsfile)
antNum  = len(antList)
polNum  = len(pol)
#-------- Tsys spectrum for specified antennas
TDMfreq, TDMPsky, TDMPamb, TDMPhot = PowerSpec( TDMmsfile, pol, TsysScan, TsysSPW, False )
ACAfreq, ACAPsky, ACAPamb, ACAPhot = PowerSpec( ACAmsfile, pol, TsysScan, TsysSPW, vanvSW )
TDMTrx, TDMTsys = TsysSpec( TDMmsfile, pol, TsysScan, TsysSPW, logfile, False )
ACATrx, ACATsys = TsysSpec( ACAmsfile, pol, TsysScan, TsysSPW, logfile, vanvSW  )
fig = plt.figure(figsize = (11, 8))
plt.suptitle(TDMprefix + ' / ' + ACAprefix)
fig.text(0.05, 0.45, 'Tsys, Trx [K]', rotation=90)
fig.text(0.45, 0.05, 'Frequency [GHz]')
xlim = [np.min(ACAfreq), np.max(ACAfreq)]; ylim= [0.0, np.max(TDMTsys)]
for ant_index in range(antNum):
	for pol_index in range(polNum):
		plt.subplot(polNum,  antNum, antNum* pol_index + ant_index + 1)
		plt.plot(TDMfreq, TDMTrx[ant_index, pol_index], ls='steps-mid')
		plt.plot(TDMfreq, TDMTsys[ant_index, pol_index], ls='steps-mid')
		plt.plot(ACAfreq, ACATrx[ant_index, pol_index], '.')
		plt.plot(ACAfreq, ACATsys[ant_index, pol_index], '.')
		plt.legend(('Trx', 'Tsys'), 'upper left', prop={'size' :9})
		plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]], fontsize=3) 
		text_sd = '%s %s' % (antList[ant_index], pol[pol_index]); plt.text(0.25*xlim[0]+0.75*xlim[1], 0.95*ylim[1], text_sd, size='x-small')
	#
#
#plt.savefig(TDMprefix + '_' + `TsysScan` + '_' + `TsysSPW` + '_TsysSpec.pdf', form='pdf'); plt.close()
#
fig = plt.figure(figsize = (11, 8))
plt.suptitle(TDMprefix + ' / ' + ACAprefix)
fig.text(0.05, 0.45, 'Power [K]', rotation=90)
fig.text(0.45, 0.05, 'Frequency [GHz]')
xlim = [np.min(ACAfreq), np.max(ACAfreq)]; ylim= [0.0, np.max(TDMPhot)]
for ant_index in range(antNum):
	for pol_index in range(polNum):
		plt.subplot(polNum,  antNum, antNum* pol_index + ant_index + 1)
		plt.plot(TDMfreq, TDMPsky[ant_index, pol_index], ls='steps-mid')
		plt.plot(TDMfreq, TDMPamb[ant_index, pol_index], ls='steps-mid')
		plt.plot(TDMfreq, TDMPhot[ant_index, pol_index], ls='steps-mid')
		plt.plot(ACAfreq, ACAPsky[ant_index, pol_index], '.')
		plt.plot(ACAfreq, ACAPamb[ant_index, pol_index], '.')
		plt.plot(ACAfreq, ACAPhot[ant_index, pol_index], '.')
		plt.legend(('Sky', 'Amb', 'Hot'), 'upper left', prop={'size' :9})
		plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]], fontsize=3) 
		text_sd = '%s %s' % (antList[ant_index], pol[pol_index]); plt.text(0.25*xlim[0]+0.75*xlim[1], 0.95*ylim[1], text_sd, size='x-small')
	#
#
#plt.savefig(TDMprefix + '_' + `TsysScan` + '_' + `TsysSPW` + '_PowerSpec.pdf', form='pdf'); plt.close()
#
guide = np.array([0.0, 550.0])
fig = plt.figure(figsize = (11, 8))
plt.suptitle(TDMprefix + ' / ' + ACAprefix)
fig.text(0.05, 0.45, 'ACA Power [K]', rotation=90)
fig.text(0.45, 0.05, 'TDM Power [K]')
if (ACAPsky.shape[2] == 124):
	TDMchRange = range(2,126)
else:
	TDMchRange = range(128)
#
for ant_index in range(antNum):
	for pol_index in range(polNum):
		plt.subplot(polNum,  antNum, antNum* pol_index + ant_index + 1)
		plt.plot(TDMPsky[ant_index, pol_index, TDMchRange], ACAPsky[ant_index, pol_index], '.')
		plt.plot(TDMPamb[ant_index, pol_index, TDMchRange], ACAPamb[ant_index, pol_index], '.')
		plt.plot(TDMPhot[ant_index, pol_index, TDMchRange], ACAPhot[ant_index, pol_index], '.')
		plt.legend(('Sky', 'Amb', 'Hot'), 'upper left', prop={'size' :9})
		text_sd = '%s %s' % (antList[ant_index], pol[pol_index]); plt.text(0.75*guide[1], 0.99*guide[1], text_sd, size='x-small')
		plt.plot(guide, guide); plt.plot(guide, 0.95*guide); plt.plot(guide, 1.05*guide)
	#
#
#plt.savefig(TDMprefix + '_' + `TsysScan` + '_' + `TsysSPW` + '_PowerCorr.pdf', form='pdf'); plt.close()
logfile.close()
