#-------- Script to compare ACA power and BB power
execfile('/users/skameno/Scripts/interferometry.py')
VanvQ3, Qeff = loadVanvQ3('/users/skameno/Scripts/ACAPowerCoeff.data')	# 3-bit Van Vleck
#VanvQ4 = loadVanvQ4('/users/skameno/Scripts/VanvQ4.data')	# 4-bit Van Vleck
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
def TsysSpec(prefix, TsysScan, TsysSPW, vanvSW=False):
	msfile = prefix + '.ms'
	antList = GetAntName(msfile)
	antNum  = len(antList)
	polNum  = len(pol)
	chNum, chWid, freq = GetChNum(msfile, TsysSPW); Tsysfreq = freq* 1.0e-9 # GHz
	TrxSp = np.zeros([antNum, polNum, chNum]); TsysSp = np.zeros([antNum, polNum, chNum])
	TrxMean = np.zeros([antNum, polNum]); TsysMean = np.zeros([antNum, polNum])
	if (chNum == 128):
		chRange = range(10, 118)
	else:
		chRange = range(8, 116)
	#
	#-------- Get Physical Temperature of loads
	#text_sd =  '%s SPW=%d SCAN=%d' % (prefix, TsysSPW, TsysScan)
	#print text_sd
	#logfile.write(text_sd + '\n')
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
			if vanvSW:
				#-------- Van Vleck Correction for 3-bit in Lag Domain
				refZeroLag = np.sum(dataXY[:, ambRange].real) / len(skyRange)
				for index in range(dataXY.shape[1]):
					temp = fft(dataXY[:, index].real)
					ZeroLag = temp.real[0]
					VanvCorrect = Qeff(ZeroLag / refZeroLag)
					ReTemp, ImTemp = temp.real / VanvCorrect, temp.imag / VanvCorrect;  ReTemp[0] = ZeroLag
					dataXY[:,index] = ifft(ReTemp + (0 + 1.0j)* ImTemp ).real
				#
				#
			#
			#-------- Calc. Tsys Spectrum
			Psky, Pamb, Phot = np.mean(dataXY[:,skyRange].real, 1), np.mean(dataXY[:,ambRange].real, 1), np.mean(dataXY[:,hotRange].real, 1)
			TrxSp[ant_index, pol_index]  = (tempHot* Pamb - Phot* tempAmb) / (Phot - Pamb)
			TsysSp[ant_index, pol_index] = (Psky* tempAmb) / (Pamb - Psky)
			TrxMean[ant_index, pol_index]  = (tempHot* np.mean(Pamb[chRange]) - np.mean(Phot[chRange])* tempAmb) / (np.mean(Phot[chRange]) - np.mean(Pamb[chRange]))
			TsysMean[ant_index, pol_index] = (np.mean(Psky[chRange])* tempAmb) / (np.mean(Pamb[chRange]) - np.mean(Psky[chRange]))
		#
	#
	return TrxSp, TsysSp, TrxMean, TsysMean
#
def TsysPM(prefix, TsysScan, TsysSPW):
	msfile = prefix + '.ms'
	antList = GetAntName(msfile)
	antNum  = len(antList)
	polNum  = len(pol)
	PskyBB = np.zeros([antNum, polNum])
	PambBB = np.zeros([antNum, polNum])
	PhotBB = np.zeros([antNum, polNum])
	TrxBB = np.zeros([antNum, polNum])
	TsysBB = np.zeros([antNum, polNum])
	for ant_index in range(antNum):
		tempAmb, tempHot = GetLoadTemp(msfile, ant_index, TsysSPW)
		for pol_index in range(polNum):
			timeBB, dataBB = GetVisibility(msfile, ant_index, ant_index, pol_index, TsysSPW, TsysScan)
			timeBB,  dataBB = BB_filter( timeBB, dataBB )
			edge = np.where( diff(timeBB) > 1.0 )[0]
			skyRange = range(0, edge[0])
			ambRange = range(edge[0]+1, edge[1])
			hotRange = range(edge[1]+1, len(timeBB))
			PskyBB[ant_index, pol_index],  PambBB[ant_index, pol_index],  PhotBB[ant_index, pol_index]  = np.mean(dataBB[skyRange]), np.mean(dataBB[ambRange]), np.mean(dataBB[hotRange])
			TrxBB[ant_index, pol_index]  = (tempHot* PambBB[ant_index, pol_index] - PhotBB[ant_index, pol_index]* tempAmb) / (PhotBB[ant_index, pol_index] - PambBB[ant_index, pol_index])
			TsysBB[ant_index, pol_index] = (PskyBB[ant_index, pol_index]* tempAmb) / (PambBB[ant_index, pol_index] - PskyBB[ant_index, pol_index])
		#
	#
	return PskyBB, PambBB, PhotBB, TrxBB, TsysBB
#
def PowerSpec(prefix, TsysScan, TsysSPW, vanvSW=False):
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
			skyRange = range(0, edge[0])
			ambRange = range(edge[0]+1, edge[1])
			hotRange = range(edge[1]+1, len(timeXY))
			#
			if vanvSW:
				#-------- Van Vleck Correction for 3-bit in Lag Domain
				refZeroLag = np.sum(dataXY[:, ambRange].real) / len(skyRange)
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
msfile = TDMprefix + '.ms'
logfile = open(TDMprefix + '_' + `TsysScan` + '_' + `TsysSPW` + '_TsysLOG.log', 'w')
antList = GetAntName(msfile)
antNum  = len(antList)
polNum  = len(pol)
#-------- Tsys spectrum for specified antennas
TDMfreq, TDMPsky, TDMPamb, TDMPhot = PowerSpec( TDMprefix, TsysScan, TsysSPW, False )
ACAfreq, ACAPsky, ACAPamb, ACAPhot = PowerSpec( ACAprefix, TsysScan, TsysSPW, True )
TDMTrx, TDMTsys, TDMTrxMean, TDMTsysMean = TsysSpec( TDMprefix, TsysScan, TsysSPW, False )
TDM_BB_Psky, TDM_BB_Pamb, TDM_BB_Phot, TDM_BB_Trx, TDM_BB_Tsys = TsysPM(TDMprefix, TsysScan, BBSPW)
ACATrx, ACATsys, ACATrxMean, ACATsysMean = TsysSpec( ACAprefix, TsysScan, TsysSPW, True  )
ACA_BB_Psky, ACA_BB_Pamb, ACA_BB_Phot, ACA_BB_Trx, ACA_BB_Tsys = TsysPM(ACAprefix, TsysScan, BBSPW)
text_sd =  '%s SPW=%d SCAN=%d' % (prefix, TsysSPW, TsysScan)
print text_sd
#logfile.write(text_sd + '\n')
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
		text_sd = '%s %d %s %s %d %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f' % (TDMprefix, TsysScan, antList[ant_index], pol[pol_index], TsysSPW, TDM_BB_Trx[ant_index, pol_index], TDM_BB_Tsys[ant_index, pol_index], TDMTrxMean[ant_index, pol_index], TDMTsysMean[ant_index, pol_index], ACA_BB_Trx[ant_index, pol_index], ACA_BB_Tsys[ant_index, pol_index], ACATrxMean[ant_index, pol_index], ACATsysMean[ant_index, pol_index])
		print text_sd
		#logfile.write(text_sd + '\n')
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
