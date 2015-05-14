import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fftpack

#-------- Definitions for Functions Used in the Allan Varianve Analysis

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
def GetVisibity(msfile, ant1, ant2, pol, spwID, field):
	Out='ANTENNA1 == '+`ant1`+' && ANTENNA2 == '+`ant2` + ' && DATA_DESC_ID == '+`spwID` + '&& FIELD_ID == ' + `field`
	tb.open(msfile)
	antXantYspw = tb.query(Out)
	timeXY = antXantYspw.getcol('TIME')
	try:
		dataXY = antXantYspw.getcol('DATA')[pol]
	except:
		dataXY = antXantYspw.getcol('FLOAT_DATA')[pol]
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

#	return np.append(corr[nspec:(2*nspec)], corr[:nspec])


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
	return abs(spec[:nspec])

#-- Measure variable components of spur
def avSpec(msfile, ant_index, pol, spw, field):
	#-- Read data from msfile
	timeXY, dataXY = GetVisibity(msfile, ant_index, ant_index, pol, spw, field)
	chNum, chWid, freq = GetChNum(msfile, spw)
	timeNum = len(timeXY)
	time_incr = min(np.diff(timeXY))
	av_free = 1.0 / abs(chWid[0]* time_incr)	# Spurious-free AV

	#-- Gain Calibration
	Gain = np.mean(abs(dataXY[(0.1*chNum):(0.9*chNum)]), 0)
	timeSpec = dataXY / Gain
	filter = hpf(2*chNum, 96)
	for time_index in range(timeNum):
		timeSpec[:,time_index] = corr2spec(filter* spec2corr(timeSpec[:,time_index]))

	#-- Caclulate time-based Allan Variances for every spectral channel
	avSpec = np.zeros(chNum)
	for ch_index in range(chNum):
		avSpec[ch_index] = allanVar(timeSpec[ch_index,:], 1)

	var_sp = (avSpec - av_free) / (3.0* timeNum)
	sigma_sp = np.sqrt(abs(var_sp))
	return sigma_sp, av_free/(3.0* timeNum)

#-- Find SPWs at highest spectral resolutions
def spurSPW(msfile):
	spwList = []
	tb.open(msfile+'/'+'SOURCE')
	spwID=tb.getcol("SPECTRAL_WINDOW_ID")
	for spw_index in range(len(spwID)):
		chNum, chWid, freq = GetChNum(msfile, spwID[spw_index])
		if(chNum > 3000):
			spwList.append(spwID[spw_index]) 
	tb.close()
	return(spwList)

#------------ Procedures
def plotSpurAV(prefix):
	msfile = prefix + '.ms'
	vis = msfile
	pol = ['XX', 'YY']
	antList = GetAntName(msfile)
	antNum  = len(antList)
	spwList = spurSPW(msfile)
	for spw in spwList:
		for pol_index in range(len(pol)):
			fig = plt.figure(figsize = (8, 11)) 
			print('Processing SPW=' + `spw` + ' Pol=' + pol[pol_index])
			for ant_index in range(antNum):
				plt.subplot( antNum, 1, ant_index+1 )
				sigma_sp, var_sp = avSpec(msfile, ant_index, pol_index, spw, 0)	# sqrt(Allan variance excess)
				chNum, chWid, freq = GetChNum(msfile, spw)		# Frequency Setting
				freq /= 1.0e9									# Hz -> GHz
				snr3, snr7 = np.sqrt(var_sp/3.0)* 3.0, np.sqrt(var_sp/3.0)* 7.0			# Spur-to-Noisr Ratio
				xlim = [min(freq), max(freq)]; ylim = [3.0e-6, 3.0e-4]			# Plot Range
				plt.fill([xlim[0], xlim[0], xlim[1], xlim[1]], [snr3, snr7, snr7, snr3], 'y', alpha=0.2)
				plt.fill([xlim[0], xlim[0], xlim[1], xlim[1]], [ylim[1], snr7, snr7 , ylim[1]], 'm', alpha=0.2)
				plt.semilogy(freq, sigma_sp, ls='steps-mid')
				plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]])
				plt.text(0.95*xlim[0] + 0.05*xlim[1], ylim[0]*0.4 + ylim[1]*0.4, antList[ant_index] + ' Pol=' + pol[pol_index])
				if(ant_index == antNum - 1):
					plt.xlabel('Frequency [GHz]')
			plt.suptitle('Spurious Scan for SPW=' + `spw`)
			plt.savefig(prefix + '_Pol-'+ pol[pol_index] + '_SPW-' + `spw` + '.pdf', form='pdf')
			plt.close()
	return
#
#-- Setup
for file_index in range(len(prefix)):
	plotSpurAV(prefix[file_index])

