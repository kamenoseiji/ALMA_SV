import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fftpack

#---- Date & Time
def md2doy(year, month, day):			# Calculate Day of Year
	doyMonth = [[0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334], [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]]
	leap_year = 0 if (year % 4) else 1	# Leap Year Flag
	return( doyMonth[leap_year][(month - 1)] + day )

def doy2mjd(year, doy, hh, mm, ss):		# Calculate MJD (in unit of second)
	MJD_1901 = 15384
	DAY_PER_4YEAR = 1461
	DAY_PER_YEAR  = 365
	SEC_PER_DAY   = 86400
	year -= 1901
	mjd = (year / 4)* DAY_PER_4YEAR + (year % 4)* DAY_PER_YEAR + doy + MJD_1901 
	sec = ((hh* 60) + mm)* 60 + ss
	return( mjd* SEC_PER_DAY + sec )

#---- Calculate Allan Variance 
def allanVar(x, lag):
	vecSize = len(x);   diffSize = vecSize - lag;   avSize = diffSize - lag
	temp = x[lag:vecSize] - x[0:diffSize]
	temp2= temp[lag:diffSize] - temp[0:avSize]
	return np.dot( temp2, temp2) / (2* avSize* lag* lag)

#---- Functions for data handling in Measurement Set
def GetAntName(msfile):
	tb.open(msfile+'/'+'ANTENNA')
	namelist = tb.getcol("NAME")
	tb.close()
	return namelist

#-- Read Visibilities and Store into an Array
def GetVisibility(msfile, ant1, ant2, pol, spwID, field):
	Out='ANTENNA1 == '+`ant1`+' && ANTENNA2 == '+`ant2` + ' && DATA_DESC_ID == '+`spwID` + '&& FIELD_ID == ' + `field`
	tb.open(msfile)
	antXantYspw = tb.query(Out)
	timeXY = antXantYspw.getcol('TIME')
	try:
		dataXY = antXantYspw.getcol('FLOAT_DATA')[pol]
	except:
		dataXY = antXantYspw.getcol('DATA')[pol]
	tb.close()
	return timeXY, dataXY

#-- Read Frequency Setups
def GetChNum(msfile, spwID):
	tb.open(msfile + '/' + 'SPECTRAL_WINDOW')
	chNum = tb.getcell("NUM_CHAN", spwID)
	chWid = tb.getcell("CHAN_WIDTH", spwID)
	freq  = tb.getcell("CHAN_FREQ", spwID)
	tb.close()
	return chNum, chWid, freq

#-- Correlation Function to remove standing wave
def spec2corr(spec):
	nspec = len(spec)
	tmpspec = np.append(spec, np.zeros(nspec))
	corr = fftpack.ifft(tmpspec)
	return corr

def hpf(nspec, width):
	filter = np.ones(nspec)
	for index in range(2*width):
		atten = np.exp( -0.5* index* index / (width* width))
		filter[index + 1] = 1.0 - atten
		filter[nspec - index - 1] = 1.0 - atten
		filter[0] = 0
	return filter

def corr2spec(corr):
	nspec = len(corr)/2
	spec = fftpack.fft(corr)[:nspec]
	return spec[:nspec].real

#-- Measure variable components of spur
def avSpec(timeXY, spec, chNum, chWid):
	timeNum = len(timeXY)
	time_incr = min(np.diff(timeXY))
	av_free = 3.0 / abs(2.0* chWid[0]* time_incr)	# Spurious-free AV, Hanning-window is assumed

	#-- Gain Calibration
	Gain = np.mean(spec[(0.1*chNum):(0.9*chNum)], 0)
	timeSpec = spec / Gain

	#-- High-pass filter to discriminate narrow-band spurious from bandpass instability 
	filter = hpf(2*chNum, int(chNum* 0.15))			# 5% bandwidth high-pass filter
	scaleFact = 1.0 / np.mean(filter* filter)
	for time_index in range(timeNum):
		timeSpec[:,time_index] = corr2spec(filter* spec2corr(timeSpec[:,time_index]))* scaleFact

	#-- Caclulate time-based Allan Variances for every spectral channel
	avSpec = np.zeros(chNum)
	for ch_index in range(chNum):
		avSpec[ch_index] = allanVar(timeSpec[ch_index,:], 1)

	return avSpec, av_free

def VarSpec(timeXY, spec, chNum, chWid):
	timeNum = len(timeXY)
	time_incr = min(np.diff(timeXY))
	var_free = 1.0 / abs(2.0* chWid[0]* time_incr)	# Spurious-free AV, Hanning-window is assumed

	#-- Gain Calibration
	Gain = np.mean(spec[(0.1*chNum):(0.9*chNum)], 0)
	timeSpec = spec / Gain

	#-- High-pass filter to discriminate narrow-band spurious from bandpass instability 
	filter = hpf(2*chNum, chNum* 0.05)			# 5% bandwidth high-pass filter
	scaleFact = 1.0 / np.mean(filter* filter)
	for time_index in range(timeNum):
		timeSpec[:,time_index] = corr2spec(filter* spec2corr(timeSpec[:,time_index]))* scaleFact

	#-- Caclulate time-based Allan Variances for every spectral channel
	varSpec = np.var(timeSpec, 1)
	return varSpec, var_free

def OnOffSpec(spec, integPP):
	numCh, numScan = spec.shape[0], spec.shape[1]
	weight = np.tile( np.array(np.repeat([1, -1], integPP)), int(numScan/(2* integPP)) )
	modSpec = spec[:, range(len(weight))]* weight
	integOnOffSpec = np.mean(modSpec, axis=1)* 2.0
	return integOnOffSpec

#------------ AV method to detect spurious signals
def plotSpurAV(prefix, spwList, polList, fieldID, sp_thresh):
	msfile = prefix + '.ms'
	antList = GetAntName(msfile)
	antNum, spwNum, polNum  = len(antList), len(spwList), len(polList)
	file_log = open(prefix + '-SP.log', 'w')		# Prepare Log File
	for ant_index in range(antNum):
		fig = plt.figure(figsize = (8, 11)) 
		fig.text(0.05, 0.45, 'Allan Variance', rotation=90)
		for pol_index in range(polNum):
			for spw_index in range(spwNum):
				plt.subplot( spwNum, polNum, spw_index* polNum + pol_index + 1 )
				chNum, chWid, freq = GetChNum(msfile, spwList[spw_index])		# Frequency Setting
				timeXY, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, spwList[spw_index], fieldID)
				av_sp, av_th = avSpec(timeXY, abs(dataXY), chNum, chWid)		# Allan Variance at sp channel and thermal noise
				smp_num = len(timeXY)

				freq /= 1.0e9									# Hz -> GHz
				snr3 = av_th* (1.0 + 1.0/sqrt(smp_num/9.0 - 0.5))		# Spur-to-Noisr Ratio
				av_thresh = 3.0* sp_thresh + av_th		# -56 dB after 16-hour integration 

				#-------- List NCR channels
				fill_color = 'g'			# Safe color
				ncr_index = np.where(av_sp[(chNum*0.05):(chNum*0.95)] > av_thresh)[0] + (chNum*0.05)
				if( len(ncr_index) > 0):
					for index in ncr_index:
						text_sp = '%s %s %s%d %s%04d %7.3f %s %s%5.1e%s%5.1e' % (antList[ant_index], polList[pol_index], 'SPW=', spwList[spw_index], 'CH=', index, freq[index], 'GHz', 'AV=', av_sp[index], '/', av_thresh)
						print text_sp
						file_log.write(text_sp + '\n')	# Logging
					fill_color = 'm'		# If detect spurious, warning color

				#-------- Plot AV
				xlim = [min(freq), max(freq)]; ylim = [av_th/3, av_th* 300]				# Plot Range
				plt.fill([xlim[0], xlim[0], xlim[1], xlim[1]], [snr3, av_thresh, av_thresh, snr3], 'b', alpha=0.2)
				plt.text(0.25*xlim[0] + 0.75*xlim[1], snr3, 'SNR = 3', color='b', size='x-small')
				plt.fill([xlim[0], xlim[0], xlim[1], xlim[1]], [ylim[1], av_thresh, av_thresh , ylim[1]], fill_color, alpha=0.2)
				plt.text(0.25*xlim[0] + 0.75*xlim[1], av_thresh, 'Requirement', color='r', size='x-small')
				plt.semilogy(freq, av_sp, ls='steps-mid')
				plt.xticks(fontsize=6); plt.yticks(fontsize=6)
				plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]], fontsize=6)
				plt.text(0.95*xlim[0] + 0.05*xlim[1], ylim[1]*0.3, 'Pol=' + polList[pol_index] + ' SPW=' + `spwList[spw_index]`, size='x-small')

			plt.xlabel('Frequency [GHz]', size='small')
		plt.suptitle('Spurious Scan Ant=' + antList[ant_index])
		plt.savefig(prefix + '_Ant-'+ antList[ant_index] + 'SPW' + `min(spw)` + '-' + `max(spw)` + '-SP.pdf', form='pdf')
		plt.close()
	file_log.close()	# Save Log File
	return timeXY, dataXY

#------------ On - off bandpass stability
def plotOnOffSpec(prefix, timeRange, bwf, spwList, polList, fieldID, integ, thresh):
	msfile = prefix + '.ms'
	antList = GetAntName(msfile)
	antNum, spwNum, polNum  = len(antList), len(spwList), len(polList)
	for ant_index in range(antNum):
		fig = plt.figure(figsize = (8, 11)) 
		fig.text(0.05, 0.45, 'On - Off Spectrum [scalsed by Tsys]', rotation=90)
		for pol_index in range(polNum):
			for spw_index in range(spwNum):
				plt.subplot( spwNum,  polNum, spw_index* polNum + pol_index + 1 )
				chNum, chWid, freq = GetChNum(msfile, spwList[spw_index])		# Frequency Setting
				chRange = range( int(chNum* (1.0 - bwf)), int(chNum* bwf))
				timeXY, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, spwList[spw_index], fieldID)
				integPP = int(max(1.0, round(integ / np.median(np.diff(timeXY)))))
				spec = abs(dataXY[:, timeRange])		# Power Spectrum
				BP = np.mean(spec, axis=1)				# Bandpass Shape
				spec = (spec.T / BP).T					# Bandpass Calibration
				OnOff = OnOffSpec( spec, integPP )		# On - Off scan
				sd_spec, pe_spec = np.std( OnOff[chRange] ), max( abs(OnOff[chRange] - np.mean(OnOff[chRange])) )
				#-------- Color of filled area
				fill_color = 'g'			# Safe color
				if(pe_spec > thresh[spw_index]):
					fill_color = 'm'		# If pe > thresh, warning color

				#-------- Plot On-Off spectrum
				freq /= 1.0e9									# Hz -> GHz
				xlim, ylim = [min(freq), max(freq)], [-2.5*thresh[spw_index], 2.5*thresh[spw_index]]						# Plot Range
				xrange, yrange = [min(freq[chRange]), max(freq[chRange])], [-thresh[spw_index], thresh[spw_index]]	# Fill Range
				yrange += np.mean(OnOff[chRange])
				plt.fill([xrange[0], xrange[0], xrange[1], xrange[1]], [yrange[1], yrange[0], yrange[0] , yrange[1]], fill_color, alpha=0.1)
				plt.plot(freq, OnOff, ls='steps-mid')
				plt.xticks(fontsize=6); plt.yticks(fontsize=6)
				plt.axis([xlim[0], xlim[1], min(ylim[0], min(OnOff[chRange])), max(ylim[1], max(OnOff[chRange]))], fontsize=6)
				text_sd = '%s%d %s%d  %s %5.1e %s %5.1e' % ('Pol=', pol_index, 'SPW=', spwList[spw_index], 'SD=', sd_spec, 'PE=', pe_spec)
				plt.text(xrange[0], yrange[1]*0.8 + yrange[0]*0.2, text_sd, size='x-small')
			plt.xlabel('Frequency [GHz]', size='small')
		plt.suptitle('On-Off Spectrum Ant=' + antList[ant_index])
		plt.savefig(prefix + '_Ant-'+ antList[ant_index] + 'SPW' + `min(spw)` + '-' + `max(spw)` + '-BP.pdf', form='pdf')
		plt.close()
	return spec

#-------- Procedures
startMJD = doy2mjd(scanStart[0], md2doy(scanStart[0], scanStart[1], scanStart[2]), scanStart[3], scanStart[4], scanStart[5])
endMJD   = doy2mjd(scanEnd[0], md2doy(scanEnd[0], scanEnd[1], scanEnd[2]), scanEnd[3], scanEnd[4], scanEnd[5])
for file_index in range(len(prefix)):
	print 'Processing ' + prefix[file_index]
	timeXY, dataXY = plotSpurAV(prefix[file_index], spw, pol, fieldID, sp_thresh)
	#listobs(prefix[file_index]+'.ms', listfile=prefix[file_index]+'.listobs')
	#timeRange = range( min(np.where(timeXY > startMJD)[0]), max(np.where(timeXY < endMJD)[0]) )
	#spec = plotOnOffSpec(prefix[file_index], timeRange, bwf, spw, pol, fieldID, integPP, np.tile([thresh], len(spw)))

