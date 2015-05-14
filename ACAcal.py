import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import scipy.fftpack as fftpack
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
#---- Functions for data handling in Measurement Set
def GetAntName(msfile):
	tb.open(msfile+'/'+'ANTENNA')
	namelist = tb.getcol("NAME")
	tb.close()
	return namelist

def Ant2Bl(ant1, ant2):
	antenna1 = max(ant1, ant2); antenna2 = min(ant1, ant2)
	return antenna1* (antenna1 - 1)/2 + antenna2

def Bl2Ant(baseline):
	ant1 = int(sqrt(2.0* baseline))
	while( (ant1* (ant1 + 1)/2 ) > baseline):
		ant1 -= 1
	return [(ant1+1), int(baseline - ant1*(ant1 + 1)/2)]

def GetTimerecord(msfile, ant1, ant2, spwID, fieldID):
	Out='ANTENNA1 == '+`ant1`+' && ANTENNA2 == '+`ant2` + ' && DATA_DESC_ID == '+`spwID` + ' && FIELD_ID == ' + `fieldID`
	tb.open(msfile)
	antXantYspw = tb.query(Out)
	interval=tb.getcol('INTERVAL')
	timeXY = antXantYspw.getcol('TIME')
	tb.close()
	return interval, timeXY

def GetCalTable(msfile, ant, spwID):
	Out='ANTENNA_ID == '+`ant` + ' && SPECTRAL_WINDOW_ID == ' + `spwID`
	tb.open(msfile + '/' + 'CALDEVICE')
	antXantYspw = tb.query(Out)
	calTemp = antXantYspw.getcol("TEMPERATURE_LOAD")
	tb.close()
	return calTemp[0,0], calTemp[1,0]

#-- Read Visibilities and Store into an Array
def GetVisibility(msfile, ant1, ant2, pol, spwID, fieldID):
	Out='ANTENNA1 == '+`ant1`+' && ANTENNA2 == '+`ant2` + ' && DATA_DESC_ID == '+`spwID` + ' && FIELD_ID == ' + `fieldID`
	tb.open(msfile)
	antXantYspw = tb.query(Out)
	timeXY = antXantYspw.getcol('TIME')
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

def corr2spec(corr):
	nspec = len(corr)/2
	spec = fftpack.fft(corr)[:nspec]
	return abs(spec[:nspec])

def corr2Xspec(corr):
	nspec = len(corr)/2
	spec = fftpack.fft(corr)[:nspec]
	return spec[:nspec]

#-- Measure variable components of spur

#-- Find SPWs at highest spectral resolutions
def spwList(msfile):
	spwList = []
	chNumList = []
	tb.open(msfile+'/'+'SOURCE')
	spwID=tb.getcol("SPECTRAL_WINDOW_ID")
	for spw_index in range(len(spwID)):
		chNum, chWid, freq = GetChNum(msfile, spwID[spw_index])
		if(chNum > 30):
			spwList.append(spwID[spw_index]) 
			chNumList.append(chNum) 
	tb.close()
	return spwList, chNumList

def rawCrossSpec(prefix, pol):
	polNum = len(pol)
	msfile = prefix + '.ms'
	antList = GetAntName(msfile)
	antNum = len(antList)
	blNum  = antNum* (antNum - 1)/2
	spwID, chNumList = spwList(msfile)
	spwNum = len(spwID)
	chNum  = max(chNumList)
	rawSpec = np.zeros([blNum, polNum, spwNum, chNum], dtype=complex)
	freq    = np.zeros([spwNum, chNum])
	for spw_index in range(spwNum):
		chNum, chWid, freq[spw_index] = GetChNum(msfile, spwID[spw_index])

	for bl_index in range(blNum):
		ants = Bl2Ant(bl_index)
		for pol_index in range(polNum):
			for spw_index in range(spwNum):
				timeXY, dataXY = GetVisibility(msfile, ants[1], ants[0], pol_index, spwID[spw_index], 0)
				rawSpec[bl_index, pol_index, spw_index] = np.mean(dataXY, axis=1)
	return freq*1.0e-9, abs(chWid[0])*1.0e-9, rawSpec

def rawPowerSpec(prefix, pol, spwID, fieldID):
	polNum = len(pol)
	spwNum = len(spwID)
	msfile = prefix + '.ms'
	spw_List, chNumList = spwList(msfile)
	antList = GetAntName(msfile)
	antNum = len(antList)
	spwNum = len(spwID)
	chNum  = chNumList[spw_List.index(spwID[0])]
	rawSpec = np.zeros([antNum, polNum, spwNum, chNum])
	freq    = np.zeros([spwNum, chNum])
	for spw_index in range(spwNum):
		chNum, chWid, freq_tmp = GetChNum(msfile, spwID[spw_index])
		freq[spw_index, 0:chNum] = freq_tmp

	for ant_index in range(antNum):
		for pol_index in range(polNum):
			for spw_index in range(spwNum):
				timeXY, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, spwID[spw_index], fieldID)
				rawSpec[ant_index, pol_index, spw_index, 0:chNum] = abs(np.mean(dataXY, axis=1))
	#return freq*1.0e-9, abs(chWid[0])*1.0e-9, rawSpec
	return freq*1.0e-9, abs(chWid[0])*1.0e-9, dataXY

def loadAcorrCoeff( calFile ):
	ACA_Caldata = loadtxt(calFile)
	#bitDist = np.transpose(ACA_Caldata[2:10])
	#analogPower  = ACA_Caldata[0]
	digitalPower = ACA_Caldata[1]
	fitCoeff = ACA_Caldata[10:12]
	return interp1d(digitalPower, fitCoeff[0], kind='cubic'), interp1d(digitalPower, fitCoeff[1], kind='cubic')

def ACA_Acorr_Correct( spec, coeff0, coeff1, coeff2):
	outSpec = np.zeros(len(spec))
	for index in range(len(spec)):
		outSpec[index] = coeff0 + spec[index]* (coeff1 + coeff2* spec[index])
	return(outSpec)
	
def VanvCorrect(  spec ):
	corrFn = 2.0* spec2corr(spec)
	vanvCorr = np.zeros(len(corrFn), dtype=complex)
	vanvSpec = np.zeros(len(spec))
	for index in range(len(corrFn)):
		vanvCorr[index] = abs(corrFn[index].real)
		if(corrFn[index].real < 0.0):
			vanvCorr[index] = -vanvCorr[index]
#		vanvCorr[index] = vanvCorr[index] + (0+1j)* interP(thresh, abs(corrFn[index].imag))[0]
		vanvCorr[index] = vanvCorr[index] + (0+1j)* abs(corrFn[index].imag)
		if(corrFn[index].imag < 0.0):
			vanvCorr[index] = np.conjugate(vanvCorr[index])
	vanvSpec = corr2spec( vanvCorr )
	return vanvSpec

