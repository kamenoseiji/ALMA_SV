import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import scipy
from scipy import constants
from scipy.fftpack import fft, ifft
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import LSQUnivariateSpline
from scipy.interpolate import griddata
import scipy.optimize
import time
import datetime
#-------- Muller Matrix
def MullerMatrix(Dx0, Dy0, Dx1, Dy1):
    return np.array([
        [1.0, Dx1.conjugate(), Dx0, Dx0* Dx1.conjugate()],
        [Dy1.conjugate(), 1.0, Dx0* Dy1.conjugate(), Dx0],
        [Dy0, Dy0* Dx1.conjugate(), 1.0, Dx1.conjugate()],
        [Dy0* Dy1.conjugate(), Dy0, Dy1.conjugate(), 1.0]])
#
def InvMullerMatrix(Dx0, Dy0, Dx1, Dy1):
    return np.array([
        [1.0, -Dx1.conjugate(), -Dx0, Dx0* Dx1.conjugate()],
        [-Dy1.conjugate(), 1.0, Dx0* Dy1.conjugate(), -Dx0],
        [-Dy0, Dy0* Dx1.conjugate(), 1.0, -Dx1.conjugate()],
        [Dy0* Dy1.conjugate(), -Dy0, -Dy1.conjugate(), 1.0]]) / ((1.0 - Dx0* Dy0)*(1.0 - Dx1.conjugate()* Dy1.conjugate()))
#
def PAMatrix(PA):
    cs = math.cos(2.0* PA)
    sn = math.sin(2.0* PA)
    return np.array([
        [1.0,  cs,  sn,  0.0],
        [0.0, -sn,  cs,  1.0j],
        [0.0, -sn,  cs, -1.0j],
        [1.0, -cs, -sn,  0.0]])
#
def InvPAMatrix(PA):
    cs = math.cos(2.0* PA)
    sn = math.sin(2.0* PA)
    return 0.5*np.array([
        [1.0, 0.0, 0.0, 1.0],
        [ cs, -sn, -sn, -cs],
        [ sn,  cs,  cs, -sn],
        [0.0,-1.0j,1.0j, 0.0]])
#
def AzEl2PA(az, el, lat):        # Azimuth, Elevation, Latitude in [rad]
    cos_lat = np.cos(lat)
    return np.arctan( -cos_lat* np.sin(az) / (np.sin(lat)* np.cos(el) - cos_lat* np.sin(el)* np.cos(az)) )
#
#-------- Greenwidge Mean Sidereal Time
def mjd2gmst( mjd, ut1utc ):        # mjd in [day], ut1utc in [sec]
    FACT = [24110.54841, 8640184.812866, 0.093104, 0.0000062]
    MJD_EPOCH = 51544.5             # MJD at 2000 1/1 12:00:00 UT
    TU_UNIT   = 36525.0
    SEC_PER_DAY = 86400.0
    tu = (mjd - MJD_EPOCH) / TU_UNIT
    ut1 = modf(mjd)[0]* SEC_PER_DAY + ut1utc
    gmst = (ut1 + FACT[0] + ((FACT[3]* tu + FACT[2])* tu + FACT[1])* tu) / SEC_PER_DAY
    return 2.0* pi* modf(gmst)[0]
#
def gst2lst( gst, longitude ):      # gst, lambda in [rad]
    return( gst + longitude)
#
def azel2radec( az, el, lst, latitude):
    dec = np.arcsin( np.sin(el)* np.sin(latitude) + np.cos(el)* np.cos(latitude)* np.cos(az) )
    ha  = np.arctan2( -np.sin(az)* np.cos(el)/np.cos(dec), (np.sin(el) - np.sin(dec)* np.sin(latitude))/(np.cos(dec)* np.cos(latitude)) )
    ra  = lst - ha
    return ra, dec
#
def gst2ha( gst, longitude, ra ):      # gst, lambda, ra in [rad]
    lst = gst + longitude
    ha  = lst - ra
    return 2.0* pi* modf((ha + pi)/ (2.0* pi))[0] - pi
#
def radec2ecliptic( ra, dec, mjd ):          # J2000 -> ecliptic, mjd in [day]
    MJD_EPOCH = 51544.5             # MJD at 2000 1/1 12:00:00 UT
    TU_UNIT   = 36525.0             # Julian Century
    tu = (mjd - MJD_EPOCH) / TU_UNIT    # Julian Century from J2000.0
    inclination = 0.4090926006005829 + ((((-2.104091376015386e-13* tu - 2.792526803190927e-12)* tu + 9.712757287348442e-09)* tu - 8.876938501115603e-10)* tu - 1.9833368961184175e-06)* tu
    cs, sn = cos(inclination), sin(inclination)
    Xa, Ya, Za = np.cos(dec)* np.cos(ra), np.cos(dec)* np.sin(ra), np.sin(dec)
    Xb, Yb, Zb = Xa, cs* Ya + sn* Za, -sn* Ya + cs* Za
    return np.arctan2(Yb, Xb), np.arcsin(Zb)
#
def ecliptic2radec( longitude, latitude, mjd ):          # ecliptic -> J2000, mjd in [day]
    MJD_EPOCH = 51544.5             # MJD at 2000 1/1 12:00:00 UT
    TU_UNIT   = 36525.0             # Julian Century
    tu = (mjd - MJD_EPOCH) / TU_UNIT    # Julian Century from J2000.0
    inclination = 0.4090926006005829 + ((((-2.104091376015386e-13* tu - 2.792526803190927e-12)* tu + 9.712757287348442e-09)* tu - 8.876938501115603e-10)* tu - 1.9833368961184175e-06)* tu
    cs, sn = cos(inclination), sin(inclination)
    Xa, Ya, Za = np.cos(latitude)* np.cos(longitude), np.cos(latitude)* np.sin(longitude), np.sin(latitude)
    Xb, Yb, Zb = Xa, cs* Ya - sn* Za, sn* Ya + cs* Za
    return np.arctan2(Yb, Xb), np.arcsin(Zb)
#
#-------- Time-based matching between time tags in visibilities and in scan pattern 
def indexList( refArray, motherArray ):     # Compare two arrays and return matched index
    IL = []
    for currentItem in refArray:
        IL = IL + np.where( motherArray == currentItem )[0].tolist()
    #
    return IL
#
def AzElMatch( refTime, scanTime, thresh, Az, El ):
    index = np.where( abs(scanTime - refTime) < thresh)[0]
    return np.median(Az[index]), np.median(El[index])
#
def GetAntD(antName):
    antD = 12.0
    if antName.find('C') > -1:
        antD = 7.0
    #
    return(antD)
#
def GetFWHM(msfile, spw, antD ):    # antD is the antenna diameter [m]
    Num, chWid, Freq = GetChNum(msfile, spw)
    wavelength = constants.c / np.median(Freq)
    return 1.13* 180.0* 3600.0* wavelength / (antD* pi) # Gaussian beam, in unit of arcsec
#
def GetSourceList(msfile):              # Source List
    sourceList, posList = [], []
    tb.open( msfile + '/SOURCE')
    SourceID   = tb.getcol('SOURCE_ID')
    SourceName = tb.getcol('NAME')
    SourcePos  = tb.getcol('DIRECTION')
    tb.close()
    sourceNum = len(np.unique(SourceID))
    for source_index in range(sourceNum):
        IDindex = np.where( SourceID == source_index)[0][0]
        sourceList.append(SourceName[IDindex])
        posList.append(SourcePos[:,IDindex])
    #
    return sourceList, posList
#
def GetAzEl(msfile):
	Out = msfile + '/' + 'POINTING'
	tb.open(Out)
	Direction=tb.getcol("DIRECTION")
	Time     = tb.getcol("TIME")
	AntID    = tb.getcol("ANTENNA_ID")
	tb.close()
	return Time, AntID, Direction[0,0], Direction[1,0]
#
def GetAzOffset(msfile):
	Out = msfile + '/' + 'POINTING'
	tb.open(Out)
	Offset   = tb.getcol("POINTING_OFFSET")
	Time     = tb.getcol("TIME")
	AntID    = tb.getcol("ANTENNA_ID")
	tb.close()
	return Time, AntID, Offset[:,0]*180*3600/pi
#
def TimeExtract(AntID, scanTime, keyTime):		# Find index in scanTime with the keyTime
	time_index = range(len(keyTime))
	for scan_index in time_index:
		index = np.where( AntID == 0 )[0]
		time_index[scan_index] = index[np.argmin( abs(scanTime[index] - keyTime[scan_index]))]
	#
	return time_index

def AzElExtract(antNum, AntID, timeXY, scanTime, Offset):
	timeIndex = range(len(timeXY))
	for scanIndex in range(len(timeXY)):
		index = np.where( AntID == 0 )[0]
		timeIndex[scanIndex] = np.argmin( abs(scanTime[index] - timeXY[scanIndex]))
	return Offset[0, timeIndex], Offset[1, timeIndex]

def isoDateTime( integerTime ):
	TU = str(integerTime).split('.')
	return( qa.time(TU[0]+'s', form="fits")[0] + '.' + TU[1])

def AllanVar(x, lag):
	vecSize = len(x);	diffSize = vecSize - lag;	avSize = diffSize - lag
	temp = x[lag:vecSize] - x[0:diffSize]
	temp2= temp[lag:diffSize] - temp[0:avSize]
	return np.dot( temp2, temp2) / (2* avSize* lag* lag)

def GetLoadTemp(msfile, AntID, spw):
	Out = msfile + '/' + 'CALDEVICE'
	tb.open(Out)
	Condition = 'ANTENNA_ID == ' + `AntID` + ' && SPECTRAL_WINDOW_ID == ' + `spw`
	temp = tb.query(Condition).getcol('TEMPERATURE_LOAD')
	tb.close()
	return np.median(temp[0]), np.median(temp[1])
#

def GetAntName(msfile):
	tb.open(msfile+'/'+'ANTENNA')
	namelist = tb.getcol("NAME")
	tb.close()
	return namelist

def GetBasePair(AntNum):
	BaseNum=int((float(AntNum))*(float(AntNum)-1)/2)
	BasePair=[]
	for i in range(AntNum):
		for j in range(i+1,AntNum):
			Test=[]
			Test.append(i)
			Test.append(j)
			BasePair.append(Test)
	return BasePair

def GetTimerecord(msfile, ant1, ant2, pol, spwID, PScan):
	Out='ANTENNA1 == '+`ant1`+' && ANTENNA2 == '+`ant2` + ' && DATA_DESC_ID == '+`spwID`  + ' && SCAN_NUMBER == ' + `PScan`
	tb.open(msfile)
	antXantYspw = tb.query(Out)
	interval=antXantYspw.getcol('INTERVAL')
	timeXY = antXantYspw.getcol('TIME')
	tb.close()
	return interval, timeXY
#
def GetVisAllBL(msfile, spwID, scanID):
	antNum = len(GetAntName(msfile))
	corrNum= antNum* (antNum + 1)/2		# Number of correlations (with autocorr)
	blNum  = corrNum - antNum
	#
	Out='DATA_DESC_ID == '+`spwID` + ' && SCAN_NUMBER == ' + `scanID`
	tb.open(msfile)
	antXantYspw = tb.query(Out, sortlist='noduplicates ANTENNA1, ANTENNA2, TIME')
	timeXY = antXantYspw.getcol('TIME')
	timeNum = len(timeXY) / corrNum
	try:
		dataXY = antXantYspw.getcol('DATA')		# dataXY in array[pol, ch, baselinextime]
	except:
		dataXY = antXantYspw.getcol('FLOAT_DATA')		# dataXY in array[pol, ch, baselinextime]
	tb.close()
	polNum, chNum = dataXY.shape[0], dataXY.shape[1]
	timeStamp = timeXY.reshape(corrNum, timeNum)[0]
	acorr_index = range(antNum)
	xcorr_index = range(blNum)
	for ant_index in range(antNum):
		acorr_index[ant_index] = Ant2Bla_RevLex(ant_index, ant_index, antNum)
	for bl_index in range(blNum):
		ant1, ant0 = Bl2Ant(bl_index)
		xcorr_index[bl_index] = Ant2Bla_RevLex(ant0, ant1, antNum)
	Pspec = dataXY.reshape(polNum, chNum, corrNum, timeNum)[:,:,acorr_index,:]
	Xspec = dataXY.reshape(polNum, chNum, corrNum, timeNum)[:,:,xcorr_index,:]
	return timeStamp, Pspec, Xspec
#
def GetUVW(msfile, spwID, scanID):
    antNum = len(GetAntName(msfile))
    corrNum= antNum* (antNum + 1)/2		# Number of correlations (with autocorr)
    blNum  = corrNum - antNum
    #
    Out='DATA_DESC_ID == '+`spwID` + ' && SCAN_NUMBER == ' + `scanID`
    tb.open(msfile)
    antXantYspw = tb.query(Out, sortlist='noduplicates ANTENNA1, ANTENNA2, TIME')
    timeXY = antXantYspw.getcol('TIME')
    timeNum = len(timeXY) / corrNum
    uvw = antXantYspw.getcol('UVW')
    tb.close()
    timeStamp = timeXY.reshape(corrNum, timeNum)[0]
    xcorr_index = range(blNum)
    for bl_index in range(blNum):
        ant1, ant0 = Bl2Ant(bl_index)
        xcorr_index[bl_index] = Ant2Bla_RevLex(ant0, ant1, antNum)
    #
    UVW = uvw.reshape(3, corrNum, timeNum)[:,xcorr_index,:]
    return timeStamp, UVW
#
def GetVisibility(msfile, ant1, ant2, pol, spwID, scanID):
	nameList = GetAntName(msfile)
	NumAnt   = len(nameList)
	BasePair = GetBasePair(NumAnt)
	NumBase  = len(BasePair)
	Out='ANTENNA1 == '+`ant1`+' && ANTENNA2 == '+`ant2` + ' && DATA_DESC_ID == '+`spwID`  + ' && SCAN_NUMBER == ' + `scanID`
	tb.open(msfile)
	antXantYspw = tb.query(Out)
	timeXY = antXantYspw.getcol('TIME')
	try:
		dataXY = antXantYspw.getcol('DATA')[pol]
	except:
		dataXY = antXantYspw.getcol('FLOAT_DATA')[pol]
	tb.close()
	return timeXY, dataXY
#
def GetPSpec(msfile, ant, spwID):
    Out='ANTENNA1 == '+`ant`+' && ANTENNA2 == '+`ant` + ' && DATA_DESC_ID == '+`spwID`
    tb.open(msfile)
    antXantYspw = tb.query(Out)
    timeXY = antXantYspw.getcol('TIME')
    try:
        dataXY = antXantYspw.getcol('DATA')
    except:
        dataXY = antXantYspw.getcol('FLOAT_DATA')
    tb.close()
    return timeXY, dataXY.real
#
#-------- Mapping antList in refList
def antIndex(refList, antList): 
    antMap = []
    for ant_index in range(len(antList)):
        if antList[ant_index] in refList:
            antMap.append(refList.index(antList[ant_index]))
        else:
            antMap.append( -1 )
        #
    #
    return antMap
#
def GetChNum(msfile, spwID):
	tb.open(msfile + '/' + 'SPECTRAL_WINDOW')
	chNum = tb.getcell("NUM_CHAN", spwID)
	chWid = tb.getcell("CHAN_WIDTH", spwID)
	freq  = tb.getcell("CHAN_FREQ", spwID)
	tb.close()
	return chNum, chWid, freq

def Ant2Bl(ant1, ant2):		# Antenna -> baseline index (without autocorr)
	antenna1 = max(ant1, ant2); antenna2 = min(ant1, ant2)
	return antenna1* (antenna1 - 1)/2 + antenna2
#
def Ant2BlD(ant1, ant2):		# Antenna -> baseline index (without autocorr)
	antenna1 = max(ant1, ant2); antenna2 = min(ant1, ant2)
	return antenna1* (antenna1 - 1)/2 + antenna2, (ant1 < ant2)
#
def Bl2Ant(baseline):
	ant1 = int(sqrt(2.0* baseline))
	while( (ant1* (ant1 + 1)/2 ) > baseline):
		ant1 -= 1
	return [(ant1+1), int(baseline - ant1*(ant1 + 1)/2)]
#
def Ant2Bla_RevLex(ant0, ant1, antNum):		# Reverse Lexical, with autcorr
	antenna0 = min(ant0, ant1); antenna1 = max(ant0, ant1)
	kernel = antNum* antenna0 - antenna0* (antenna0 - 1)/2
	return kernel + antenna1 - antenna0
#
def Ant2Bl_RevLex(ant1, ant2, antnum):
	antenna1 = max(ant1, ant2); antenna2 = min(ant1, ant2)
	return int(antnum* antenna2 - (antenna2 + 1)* (antenna2 + 2) / 2  + antenna1)
#
#-------- BL-ANT indexing
ANT0 = []; ANT1 = []
for bl_index in range(2016):    # Maximum number of baseline
    ants = Bl2Ant(bl_index)
    ANT0.append(ants[0])        # bl -> ant0 (baseline-end antenna) mapping
    ANT1.append(ants[1])        # bl -> ant1 (baseline-begin antenna) mapping
#
#ANT0 = np.array(ANT0); ANT1 = np.array(ANT1)
#
def BlAmpMatrix(antNum):
    blNum = antNum* (antNum - 1) / 2
    blamp_matrix = np.zeros([blNum, antNum])
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        blamp_matrix[bl_index, ants[0]] = 1
        blamp_matrix[bl_index, ants[1]] = 1
    return blamp_matrix
#
def BlPhaseMatrix(antNum):
    blNum = antNum* (antNum - 1) / 2
    blphs_matrix = np.zeros([blNum, (antNum - 1)])
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        blphs_matrix[bl_index, (ants[0] - 1)] = 1
        if(ants[1] > 0):
            blphs_matrix[bl_index, (ants[1] - 1)] = -1
    return blphs_matrix
#
def DxMatrix(antNum):
    blNum = antNum* (antNum - 1) /2
    Dx_matrix = np.zeros([blNum, antNum])
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        Dx_matrix[bl_index, ants[1]] = 1
    #
    return Dx_matrix
#
def DyMatrix(antNum):
    blNum = antNum* (antNum - 1) /2
    Dy_matrix = np.zeros([blNum, antNum])
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        Dy_matrix[bl_index, ants[0]] = 1
    #
    return Dy_matrix
#
def DxImMatrix(antNum):
    blNum = antNum* (antNum - 1) / 2
    blphs_matrix = np.zeros([blNum, (antNum - 1)])
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        if(ants[1] > 0):
            blphs_matrix[bl_index, (ants[1] - 1)] = 1
    return blphs_matrix
#
def DyImMatrix(antNum):
    blNum = antNum* (antNum - 1) / 2
    blphs_matrix = np.zeros([blNum, (antNum - 1)])
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        blphs_matrix[bl_index, (ants[0] - 1)] = 1
    return blphs_matrix
#
def ant2blphs( ant_phase, antphs_error ):
	antnum = len(ant_phase)					# Number of antennas
	blnum  = antnum* (antnum - 1) / 2		# Number of baselines
	bl_phase = np.zeros(blnum); bl_phase_err = np.zeros(blnum)
	for bl_index in range(blnum):
		ants =  Bl2Ant(bl_index)		# antennas used in the baseline
		bl_phase[bl_index] = ant_phase[ants[0]] - ant_phase[ants[1]]
		bl_phase_err[bl_index] = sqrt( antphs_error[ants[0]]* antphs_error[ants[0]] + antphs_error[ants[0]] * antphs_error[ants[0]])

	return arctan2(sin(bl_phase), cos(bl_phase)), bl_phase_err
#
def ant2blamp(ant_amp, antamp_error):
	antnum = len(ant_amp)				# Number of antennas
	blnum  = antnum* (antnum - 1) / 2		# Number of baselines
	bl_amp = np.zeros(antnum); bl_amp_err = np.zeros(antnum)			# Prepare output vectors
	for bl_index in range(blnum):
		ants = Bl2Ant(bl_index)
		bl_amp[bl_index] = sqrt(ant_amp[ants[0]] * ant_amp[ants[1]])
		bl_amp_err[bl_index] = sqrt( antamp_error[ants[0]]* antamp_error[ants[0]] + antamp_error[ants[1]]* antamp_error[ants[1]])
	return  bl_amp, bl_amp_err


def clamp_solve(bl_amp, bl_error):
	blnum  =  len(bl_amp)
	antnum =  Bl2Ant(blnum)[0]		# Number of baselines and antennas
	weight =  np.divide(1.0, np.multiply(bl_error, bl_error))
	log_bl =  np.log(bl_amp)	# weights from standard error
	p_matrix = BlAmpMatrix(antnum)									# Elements of the matrix
	ptwp     = np.dot(p_matrix.T, np.dot(np.diag(weight), p_matrix))
	ptwp_inv = scipy.linalg.inv(ptwp)
	solution = np.exp(np.dot(ptwp_inv,  np.dot(p_matrix.T,  np.dot(np.diag(weight), log_bl))))
	return solution, solution* np.sqrt(np.diag(ptwp_inv))

def cldelay_solve(bl_delay, bl_error):
	blnum  =  len(bl_delay)
	antnum =  Bl2Ant(blnum)[0]
	weight =  np.divide(1.0, np.multiply(bl_error, bl_error))
	solution = np.zeros(antnum-1)
	
	#---- Partial Matrix
	p_matrix   = BlPhaseMatrix(antnum)
	ptwp       = np.dot(p_matrix.T, np.dot(np.diag(weight), p_matrix))
	ptwp_inv   = scipy.linalg.inv(ptwp)
	solution   = np.dot(ptwp_inv,  np.dot(p_matrix.T,  np.dot(np.diag(weight), bl_delay)))
	return np.append(0, solution), np.append(0, np.sqrt(np.diag(ptwp_inv)))

def clcomplex_solve(bl_vis, bl_error):
	blnum  =  len(bl_vis)
	antnum =  Bl2Ant(blnum)[0]
	weight =  np.divide(1.0, np.multiply(bl_error, bl_error))
	weight = np.append(weight, weight)
	#
	resid  =  np.zeros(2* blnum)
	niter  = 0
	correction = np.ones(2* antnum - 1)
	solution   = np.zeros(2* antnum - 1)
	#
	#---- Initial solution
	solution[0] = sqrt(abs(bl_vis[0]))		# Refant has only real part
	#solution[0] = sqrt(np.median(abs(bl_vis)))		# Refant has only real part
	for ant_index in range(1, antnum) :
		solution[ant_index]			= bl_vis[Ant2Bl(0, ant_index )].real / solution[0]
		solution[antnum + ant_index - 1]= bl_vis[Ant2Bl(0, ant_index )].imag / solution[0]
	#
	while (np.dot(correction, correction) > 1e-12) and (niter < 10) :
		complex_matrix = np.zeros((2*blnum, 2*antnum - 1))
		#-------- Residual Vector
		for bl_index in range(blnum):
			ants = Bl2Ant(bl_index)
			if ants[1] != 0:
				resid[bl_index]			= bl_vis[bl_index].real - (solution[ants[0]]* solution[ants[1]] + solution[antnum + ants[0] - 1]* solution[antnum + ants[1] - 1])	# Real part
				resid[blnum + bl_index] = bl_vis[bl_index].imag - (solution[ants[1]]* solution[antnum + ants[0] - 1] - solution[ants[0]]* solution[antnum + ants[1] - 1])	# Imag part
			else:
				resid[bl_index]			= bl_vis[bl_index].real - (solution[ants[0]]* solution[0])	# Real part
				resid[blnum + bl_index] = bl_vis[bl_index].imag - (solution[antnum + ants[0] - 1]* solution[0])	# Imag part
		#print 'Iteration ' + `niter` + '  Resid = ' + `np.dot(resid, resid)`
		#
		#---- Partial Matrix
		for bl_index in range(blnum):
			ants = Bl2Ant(bl_index)
			complex_matrix[bl_index, ants[0]] = solution[ants[1]]
			complex_matrix[bl_index, ants[1]] = solution[ants[0]]
			if ants[1] != 0:
				complex_matrix[bl_index, antnum + ants[0] - 1] = solution[antnum + ants[1] - 1]
				complex_matrix[bl_index, antnum + ants[1] - 1] = solution[antnum + ants[0] - 1]
				complex_matrix[blnum + bl_index, ants[1]] =  solution[antnum + ants[0] - 1]
				complex_matrix[blnum + bl_index, ants[0]] = -solution[antnum + ants[1] - 1]
				complex_matrix[blnum + bl_index, antnum + ants[1] - 1] = -solution[ants[0]]
				complex_matrix[blnum + bl_index, antnum + ants[0] - 1] =  solution[ants[1]]
			else:		# No ants[1]
				complex_matrix[blnum + bl_index, 0]		= solution[antnum + ants[0] - 1]
				complex_matrix[blnum + bl_index, antnum + ants[0] - 1]= solution[0]
		#
		ptwp = np.dot( complex_matrix.T, np.dot(np.diag(weight), complex_matrix))
		ptwp_inv   = scipy.linalg.inv(ptwp)
		correction = np.dot(ptwp_inv,  np.dot(complex_matrix.T, np.dot(np.diag(weight), resid)))
		solution   = np.add(solution, correction)
		niter      =  niter + 1
	#
	return solution[range(antnum)] + 1j* np.append(0, solution[range(antnum, 2*antnum-1)])
#
def clphase_solve(bl_phase, bl_error):
	blnum  =  len(bl_phase)
	antnum =  Bl2Ant(blnum)[0]
	weight =  np.divide(1.0, np.multiply(bl_error, bl_error))

	resid  =  np.zeros(blnum)
	niter  = 0
	correction = np.ones(antnum - 1)
	solution   = np.zeros(antnum - 1)
	
	#---- Initial solution
	for ant_index in range(1,antnum) :
		solution[ant_index - 1] = bl_phase[Ant2Bl(0, ant_index )]
	
	while  (np.dot(correction, correction) > 1e-12) and (niter < 10):	# Loop for iterations
		#---- Residual Vector
		for bl_index in range(blnum):
			ants = Bl2Ant(bl_index)
			phs_diff  =  bl_phase[bl_index] - solution[ants[0] - 1]
			if  ants[1] != 0 :
				phs_diff = phs_diff + solution[ants[1]-1]  	# Baselines without refant
			resid[bl_index] = atan2( sin(phs_diff), cos(phs_diff) )* weight[bl_index]	# 2pi ambibuity removal
	
		#---- Partial Matrix
		p_matrix   = BlPhaseMatrix(antnum)
		ptwp       = np.dot(p_matrix.T, np.dot(np.diag(weight), p_matrix))
		ptwp_inv   = scipy.linalg.inv(ptwp)
		correction = np.dot(ptwp_inv,  np.dot(p_matrix.T, resid))
		solution   = np.add(solution, correction)
		niter      =  niter + 1
    #
	return np.append(0, solution), np.append(0, np.sqrt(np.diag(ptwp_inv)))
#
def Vis2solveDD(Vis, PS):
    blNum  = len(Vis) / 4                   # (I, Q, U, V)
    antNum = Bl2Ant(blNum)[0]               # Number of tracking antennas
    Dx = np.zeros(antNum, dtype=complex)    # Dx, Dy solutions for scanning antenna
    Dy = np.zeros(antNum, dtype=complex)    # Dx, Dy solutions for scanning antenna
    W = np.diag( np.r_[0.1*np.ones(blNum), np.ones(blNum), np.ones(blNum), 0.1*np.ones(blNum), 0.1*np.ones(blNum), np.ones(blNum), np.ones(blNum), 0.1*np.ones(blNum)] )
    ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
    for loop_index in range(2):
        residVis = np.zeros(8* blNum)    # Residual (real and imaginary) vector (Obs - Model)
        P = np.zeros([8*blNum, 4* antNum]) # [ReXX, ReXY, ReYX, ReYY, ImXX, ImXY, ImYX, ImYY] x [ReDx, ImDx, ReDy, ImDy], No Im part in refant
        #-------- Determine P matrix
        for bl_index in range(blNum):
            # ants = Bl2Ant(bl_index)
            stokesReIndex = range(bl_index, 4* blNum, blNum)
            stokesImIndex = range(bl_index + 4*blNum, 8*blNum, blNum)
            ModelVis = np.dot(MullerMatrix(Dx[ant1[bl_index]], Dy[ant1[bl_index]], Dx[ant0[bl_index]], Dy[ant0[bl_index]]), PS)
            residVis[stokesReIndex] = Vis[stokesReIndex].real  - ModelVis.real
            residVis[stokesImIndex] = Vis[stokesReIndex].imag  - ModelVis.imag
            #-------- Derivative by  ReDx0
            DeltaP = 100.0*(np.dot(MullerMatrix(Dx[ant1[bl_index]] + 0.01, Dy[ant1[bl_index]], Dx[ant0[bl_index]], Dy[ant0[bl_index]]), PS) - ModelVis)
            P[stokesReIndex, ant1[bl_index]] += DeltaP.real
            P[stokesImIndex, ant1[bl_index]] += DeltaP.imag
            #-------- Derivative by ImDx0
            DeltaP = 100.0*(np.dot(MullerMatrix(Dx[ant1[bl_index]] + 0.01j, Dy[ant1[bl_index]], Dx[ant0[bl_index]], Dy[ant0[bl_index]]), PS) - ModelVis)
            P[stokesReIndex, antNum + ant1[bl_index]] += DeltaP.real
            P[stokesImIndex, antNum + ant1[bl_index]] += DeltaP.imag
            #-------- Derivative by  ReDx1
            DeltaP = 100.0*(np.dot(MullerMatrix(Dx[ant1[bl_index]], Dy[ant1[bl_index]], Dx[ant0[bl_index]] + 0.01, Dy[ant0[bl_index]]), PS) - ModelVis)
            P[stokesReIndex, ant0[bl_index]] += DeltaP.real
            P[stokesImIndex, ant0[bl_index]] += DeltaP.imag
            #-------- Derivative by ImDx1
            DeltaP = 100.0*(np.dot(MullerMatrix(Dx[ant1[bl_index]], Dy[ant1[bl_index]], Dx[ant0[bl_index]] + 0.01j, Dy[ant0[bl_index]]), PS) - ModelVis)
            P[stokesReIndex,antNum + ant0[bl_index]] += DeltaP.real
            P[stokesImIndex,antNum + ant0[bl_index]] += DeltaP.imag
            #-------- Derivative by  ReDy0
            DeltaP = 100.0*(np.dot(MullerMatrix(Dx[ant1[bl_index]], Dy[ant1[bl_index]] + 0.01, Dx[ant0[bl_index]], Dy[ant0[bl_index]]), PS) - ModelVis)
            P[stokesReIndex, 2*antNum + ant1[bl_index]] += DeltaP.real
            P[stokesImIndex, 2*antNum + ant1[bl_index]] += DeltaP.imag
            #-------- Derivative by ImDy0
            DeltaP = 100.0*(np.dot(MullerMatrix(Dx[ant1[bl_index]], Dy[ant1[bl_index]] + 0.01j, Dx[ant0[bl_index]], Dy[ant0[bl_index]]), PS) - ModelVis)
            P[stokesReIndex, 3*antNum + ant1[bl_index]] += DeltaP.real
            P[stokesImIndex, 3*antNum + ant1[bl_index]] += DeltaP.imag
            #-------- Derivative by  ReDy1
            DeltaP = 100.0*(np.dot(MullerMatrix(Dx[ant1[bl_index]], Dy[ant1[bl_index]], Dx[ant0[bl_index]], Dy[ant0[bl_index]] + 0.01), PS) - ModelVis)
            P[stokesReIndex, 2*antNum + ant0[bl_index]] += DeltaP.real
            P[stokesImIndex, 2*antNum + ant0[bl_index]] += DeltaP.imag
            #-------- Derivative by ImDy1
            DeltaP = 100.0*(np.dot(MullerMatrix(Dx[ant1[bl_index]], Dy[ant1[bl_index]], Dx[ant0[bl_index]], Dy[ant0[bl_index]] + 0.01j), PS) - ModelVis)
            P[stokesReIndex, 3*antNum + ant0[bl_index]] += DeltaP.real
            P[stokesImIndex, 3*antNum + ant0[bl_index]] += DeltaP.imag
        #
        #P = P[:, range(3*antNum) + range(3*antNum + 1, 4*antNum)]
        P = P[:, range(antNum) + range(antNum + 1, 4*antNum)]
        PTWP = np.dot(P.T, np.dot(W, P))
        D_vec = np.dot( scipy.linalg.inv( PTWP ), np.dot(P.T, np.dot(W, residVis)))
        # print 'Iter%d : %e %e' % (loop_index, np.dot(residVis, residVis), np.dot(D_vec, D_vec))
        Dx += D_vec[0:antNum] + 1.0j* np.r_[0, D_vec[(antNum):(2*antNum-1)]]
        Dy += D_vec[(2*antNum-1):(3*antNum-1)] + 1.0j* D_vec[(3*antNum-1):(4*antNum-1)]
        #Dy += D_vec[(2*antNum):(3*antNum)] + 1.0j* np.r_[0, D_vec[(3*antNum):(4*antNum-1)]]
    #
    return Dx, Dy
#
def Vis2solveD(Vis, DtX, DtY, PS ):
    trkAntNum = len(DtX)  # Number of tracking antennas
    W = np.diag( np.r_[0.1*np.ones(trkAntNum), np.ones(trkAntNum), np.ones(trkAntNum), 0.1*np.ones(trkAntNum), 0.1*np.ones(trkAntNum), np.ones(trkAntNum), np.ones(trkAntNum), 0.1*np.ones(trkAntNum)] )
    Dx = 0.0 + 0.0j
    Dy = 0.0 + 0.0j
    for loop_index in range(2):
        residVis = np.zeros(4* trkAntNum, dtype=complex)    # Residual (real and imaginary) vector (Obs - Model)
        P = np.zeros([8*trkAntNum, 4]) # [ReXX, ReXY, ReYX, ReYY, ImXX, ImXY, ImYX, ImYY] x [ReDx, ImDx, ReDy, ImDy]
        #-------- Determine P matrix
        for index in range(trkAntNum):
            stokesReIndex = range(index, 4* trkAntNum, trkAntNum)
            stokesImIndex = range(index + 4* trkAntNum, 8*trkAntNum, trkAntNum)
            ModelVis = np.dot(MullerMatrix(DtX[index], DtY[index], Dx, Dy), PS)
            residVis[stokesReIndex] = Vis[stokesReIndex] - ModelVis
            #-------- Derivative by  ReDx
            DeltaP = 100.0*(np.dot(MullerMatrix(DtX[index], DtY[index], Dx + 0.01,  Dy), PS) - ModelVis)
            P[stokesReIndex, 0] += DeltaP.real
            P[stokesImIndex, 0] += DeltaP.imag
            #-------- Derivative by ImDx
            DeltaP = 100.0*(np.dot(MullerMatrix(DtX[index], DtY[index], Dx + 0.01j, Dy), PS) - ModelVis)
            P[stokesReIndex, 1] += DeltaP.real
            P[stokesImIndex, 1] += DeltaP.imag
            #-------- Derivative by  ReDy
            DeltaP = 100.0*(np.dot(MullerMatrix(DtX[index], DtY[index], Dx, Dy + 0.01),  PS) - ModelVis)
            P[stokesReIndex, 2] += DeltaP.real
            P[stokesImIndex, 2] += DeltaP.imag
            #-------- Derivative by ImDx
            DeltaP = 100.0*(np.dot(MullerMatrix(DtX[index], DtY[index], Dx, Dy + 0.01j), PS) - ModelVis)
            P[stokesReIndex, 3] += DeltaP.real
            P[stokesImIndex, 3] += DeltaP.imag
        #
        PTWP = np.dot(P.T, np.dot(W, P))
        D_vec = np.dot( np.linalg.inv( PTWP ), np.dot(P.T, np.dot(W, np.r_[residVis.real, residVis.imag])))
        Dx += (D_vec[0] + 1.0j*D_vec[1])
        Dy += (D_vec[2] + 1.0j*D_vec[3])
        # print 'Iter%d : Residual = %e' % (loop_index, np.dot(np.r_[residVis.real, residVis.imag], np.r_[residVis.real, residVis.imag]))
    #
    return Dx, Dy
#
def beamF(disk2FWHMratio):     # diskR / FWHM ratio
    disk2sigma = disk2FWHMratio * 2.3548200450309493   # FWHM / (2.0* sqrt(2.0* log(2.0)))
    return( 2.0* (1.0 - exp(-0.5* (disk2sigma)**2)) / (disk2sigma**2) )
#
def Tb2Flux(Tb, Freq, diskR):   # Tb [K], Freq [GHz], diskR [arcsec]
    c = 299792458               # m/s
    hPlanck = 6.62606957e-7     # scaled by 1e27
    h_over_kb = 0.04799243      # scaled by 1e9
    solidAngle = 7.384135e15* diskR**2  # scaled by 1e26
    intensity = 2.0* hPlanck* Freq**3 / ( c**2 * (exp(h_over_kb* Freq / Tb) - 1.0))
    return(intensity* solidAngle)   # Jy
#
def corr2spec( corr ):
	nspec = len(corr)/2
	spec  = fft(corr)[1:nspec]
	return spec[:nspec]

def spec2Acorr(spec):
	nspec = len(spec)
	tmpspec = np.append(spec, np.zeros(1, dtype=complex))
	tmpspec = np.append(tmpspec, spec[255:1:-1])
	corr = ifft(tmpspec).real
	return(corr)
#

def spec2corr(spec):
	nspec = len(spec)
	tmpspec = np.append(spec, np.zeros(nspec, dtype=complex))
	corr = ifft(tmpspec)
	return np.append(corr[nspec:(2*nspec)], corr[:nspec])


def delay_cal( spec, delay ):
	# spec : input spectrum (complex)
	# delay : delay[1] = initial phase, delay[2] = delay
	# delay_cal() returns delay-calibrated spectrum
	#
	nspec = len( spec )
	twiddle = np.exp( pi* delay* np.multiply(range(-nspec/2, nspec/2), 1j) / nspec )
	return np.multiply(spec, twiddle)

def delay_search( spec ):
	nspec = len( spec )
	#-------- Search for delay
	corrAbs = abs(spec2corr(spec))
	if( max(corrAbs) == 0.0 ):	# Unavailable baseline
		return 0.0, 1.0e-20	
	#
	delay = np.where(corrAbs == max(corrAbs))[0][0] - nspec # Coarse Delay
	trial_delay = np.multiply(range(-2,3), 0.5)
	trial_amp   = np.zeros(5) 
	for i in range(5):
		# trial_amp[i] = abs(sum(delay_cal(spec, delay + trial_delay[i])))
		trial_amp[i] = abs(np.mean(delay_cal(spec, delay + trial_delay[i])))
	fit = np.polyfit(trial_delay, trial_amp, 2)
	return delay - fit[1]/(2.0*fit[0]), fit[2] - fit[1]**2/(4*fit[0])	# return delay and amp

def blGain( blSpec ):				# Average in spectral channels
	return np.mean(blSpec, 0)

def blBp( blSpec ):					# Average in time
	return np.mean(blSpec, 1)

def bunchVec( spec, bunch, bunchNum ):
	chNum = len(spec)/bunchNum
	totalNum = chNum* bunchNum
	return(np.mean( spec[0:totalNum].reshape(chNum, bunchNum), axis=1))

def specBunch( blSpec, chBunch, timeBunch ):
	chNum, timeNum = blSpec.shape[0], blSpec.shape[1]
	totalNum = chNum* timeNum
	tmp = np.mean(blSpec.reshape(chBunch, totalNum/chBunch), 0).reshape((chNum/chBunch), timeNum).T.reshape(timeBunch, totalNum/chBunch/timeBunch)
	return np.mean(tmp, 0).reshape((timeNum/timeBunch), (chNum/chBunch)).T 

def gainAnt(vis, viserr):		# vis[blNum]
	blNum = len(vis);  antNum = Bl2Ant(blNum)[0]
	GainAmp_ant = clamp_solve( abs(vis), viserr)[0]
	GainPhs_ant = clphase_solve( np.angle(vis), viserr/abs(vis))[0]
	return GainAmp_ant* np.exp(1j * GainPhs_ant)

def bpPhsAnt(spec):			# Determine phase-only antenna-based BP
	blnum, chNum, timeNum = len(spec), len(spec[0]), len(spec[0,0])
	antnum = Bl2Ant(blnum)[0]
	BP_bl  = np.zeros([blnum, chNum], dtype=complex)
	BP_ant = np.zeros([antnum, chNum], dtype=complex)
	for bl_index in range(blnum):
		BP_bl[bl_index,:] = np.exp(1.0j* np.angle(blBp(spec[bl_index])))
	for ch_index in range(chNum):
		BP_ant[:,ch_index] = np.exp(1.0j* clphase_solve(np.angle(BP_bl[:,ch_index]) , [np.std(np.angle(BP_bl[:,ch_index]))]*blnum)[0])
	return BP_ant
#
def delayCalSpec( Xspec, chRange ):     # chRange = [startCH:stopCH] specifies channels to determine delay 
    blNum, chNum = Xspec.shape[0], Xspec.shape[1]
    delay_bl = np.zeros(blNum)
    delayCalXspec = np.zeros([blNum, chNum], dtype=complex)
    delayResults = np.apply_along_axis( delay_search, 1, Xspec[:, chRange] )    # BL delays in delayResults[:,0], Amps in delayResults[:,1]
    #
    delay_ant = cldelay_solve(delayResults[:,0], 1.0/delayResults[:,1])[0]
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        delay_bl[bl_index] = delay_ant[ants[0]] - delay_ant[ants[1]]
    #
    delayCalXspec = np.apply_along_axis( delay_cal, 0, Xspec, delay_bl)
    return delay_ant, delayCalXspec
#
def delayCalSpec2( Xspec, chRange, sigma ):  # chRange = [startCH:stopCH] specifies channels to determine delay 
	blNum, chNum = Xspec.shape[0], Xspec.shape[1]
	delay_bl, amp_bl = np.zeros(blNum), np.zeros(blNum)
	delayCalXspec = np.zeros([blNum, chNum], dtype=complex)
	for bl_index in range(blNum):
		try:
			delay_bl[bl_index], amp_bl[bl_index] = delay_search(Xspec[bl_index, chRange])
		except:
			delay_bl[bl_index] = 0.0
			amp_bl[bl_index] = sigma
		#
	#
	delay_ant, delay_err = cldelay_solve(delay_bl, sigma/amp_bl)
	for bl_index in range(blNum):
		ants = Bl2Ant(bl_index); ant1, ant2 = ants[1], ants[0]
		delayCalXspec[bl_index] = delay_cal(Xspec[bl_index], delay_ant[ant2] - delay_ant[ant1])
	#   
	return delay_ant, delay_err, delayCalXspec
#
def BPtable(msfile, spw, BPScan):   # 
    pPol, cPol = [0,1], []  # parallel and cross pol
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, BPScan)    # Xspec[pol, ch, bl, time]
    polNum, chNum = Pspec.shape[0], Pspec.shape[1]
    if polNum == 4:
        pPol, cPol = [0,3], [1,2]  # parallel and cross pol
    #
    chRange = range(int(0.05*chNum), int(0.95*chNum))                   # Trim band edge
    BP_ant  = np.ones([antNum, ppolNum, chNum], dtype=complex)          # BP_ant[ant, pol, ch]
    XPspec  = np.zeros([polNum, chNum, blNum], dtype=complex)           # XPspec[pol, ch, bl]
    BPXY = np.ones([chNum, blNum], dtype=complex)
    BPYX = np.ones([chNum, blNum], dtype=complex)
    #---- Phase inversion for inverted baselines
    Xspec = Xspec[:,:,blMap]                                            # Canonical ordering 
    Ximag = Xspec.transpose(0,1,3,2).imag* (-2.0* np.array(blInv) + 1.0)
    Xreal = Xspec.transpose(0,1,3,2).real
    #---- Phase inversion for cross-pol
    if cpolNum > 0: # Full polarization pairs
        Xspec[0].imag = Ximag[0].transpose(0,2,1)   # XX
        Xspec[3].imag = Ximag[3].transpose(0,2,1)   # YY
        Xspec[1].real = (Xreal[1]*(1.0 - np.array(blInv)) + Xreal[2]* np.array(blInv)).transpose(0,2,1) # ReXY
        Xspec[1].imag = (Ximag[1]*(1.0 - np.array(blInv)) + Ximag[2]* np.array(blInv)).transpose(0,2,1) # ImXY
        Xspec[2].real = (Xreal[2]*(1.0 - np.array(blInv)) + Xreal[1]* np.array(blInv)).transpose(0,2,1) # ReYX
        Xspec[2].imag = (Ximag[2]*(1.0 - np.array(blInv)) + Ximag[1]* np.array(blInv)).transpose(0,2,1) # ImYX
        chAvgXX = np.mean(Xspec[0,chRange], axis=0 )
        chAvgYY = np.mean(Xspec[3,chRange], axis=0 )
    else:   # parallel polarization only
        Xspec[0].imag = Ximag[0].transpose(0,2,1)
        Xspec[1].imag = Ximag[1].transpose(0,2,1)
        chAvgXX = np.mean(Xspec[0,chRange], axis=0 )
        chAvgYY = np.mean(Xspec[1,chRange], axis=0 )
    #-------- Antenna-based Gain Cal
    GainX = np.apply_along_axis( gainComplex, 0, chAvgXX)
    GainY = np.apply_along_axis( gainComplex, 0, chAvgYY)
    if cpolNum > 0: # Full polarization pairs
        for ch_index in range(chNum):
            Xspec[0, ch_index] = gainCalVis( Xspec[0,ch_index], GainX, GainX)
            Xspec[1, ch_index] = gainCalVis( Xspec[1,ch_index], GainX, GainY)
            Xspec[2, ch_index] = gainCalVis( Xspec[2,ch_index], GainY, GainX)
            Xspec[3, ch_index] = gainCalVis( Xspec[3,ch_index], GainY, GainY)
        #
        XCspec = np.mean(Xspec, axis=3)[cPol]                         # Time Average and Select Pol
    else:
        for ch_index in range(chNum):
            Xspec[0, ch_index] = gainCalVis(Xspec[0,ch_index], GainX, GainX)
            Xspec[1, ch_index] = gainCalVis(Xspec[1,ch_index], GainY, GainY)
        #
    #
    #-------- Time Average
    XPspec[pPol[0]] = np.mean(Xspec, axis=3)[pPol[0]]  # Time Average and Select Pol
    XPspec[pPol[1]] = np.mean(Xspec, axis=3)[pPol[1]]  # Time Average and Select Pol
    if cpolNum > 0: # Full polarization pairs
        XPspec[cPol[0]] = np.mean(Xspec, axis=3)[cPol[0]]   # Time Average and Select Pol
        XPspec[cPol[1]] = np.mean(Xspec, axis=3)[cPol[1]]   # Time Average and Select Pol
    #
    #-------- Antenna-based bandpass spectra
    for pol_index in range(ppolNum):
        #-------- Solution (BL -> Ant)
        BP_ant[:,pol_index] = np.apply_along_axis(gainComplex, 0, XPspec[pPol[pol_index]].T)
    #
    #-------- Bandpass Correction for Cross-pol
    XYdelay = 0.0;
    if cpolNum > 0: # Full polarization pairs
        BPXY = (BP_ant[ant0, 0]* BP_ant[ant1, 1].conjugate()).T
        BPYX = (BP_ant[ant0, 1]* BP_ant[ant1, 0].conjugate()).T
        XC = np.mean( (XCspec[0] / BPXY), axis=1 ) + np.mean( (XCspec[1] / BPYX), axis=1 ).conjugate()
        XYdelay, amp = delay_search( XC[chRange] )
        XYdelay = (float(chNum) / float(len(chRange)))* XYdelay
    #
    return BP_ant, XYdelay
#
def bpCal(spec, BP0, BP1):      # spec[blNum, chNum, timeNum]
    blnum, chNum, timeNum = len(spec), len(spec[0]), len(spec[0,0])
    ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
    return (spec.transpose(2,0,1) / (BP1[ant0]* BP0[ant1].conjugate())).transpose(1, 2, 0)
	#bpCalXX = np.zeros([blnum, chNum, timeNum], dtype=complex)
	#BP_bl = np.zeros([blnum, chNum], dtype=complex)
	#for bl_index in range(blnum):
	#	ants = Bl2Ant(bl_index)
	#	BP_bl[bl_index] = BP1[ants[0]]* BP0[ants[1]].conjugate()
	#	bpCalXX[bl_index] = (spec[bl_index].T / BP_bl[bl_index]).T
	#return bpCalXX
#
def phaseCal(spec, Gain):   # spec[blNum, chNum, timeNum], Gain[antNum, chNum, timeNum]
    blnum, chNum, timeNum = spec.shape[0], spec.shape[1], spec.shape[2]
    ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
    return spec * abs(Gain[ant0]* Gain[ant1].conjugate()) / (Gain[ant0]* Gain[ant1].conjugate())
	#gainCalXX = np.zeros([blnum, chNum, timeNum], dtype=complex)
	#for bl_index in range(blnum):
	#	ants = Bl2Ant(bl_index); ant1, ant2 = ants[1], ants[0]
	#	Gain_bl = Gain[ant2].conjugate() *  Gain[ant1]
	#	Gain_bl = Gain_bl / abs(Gain_bl)
	#	gainCalXX[bl_index] = spec[bl_index] * Gain_bl
	#return gainCalXX
#
def gainCal(spec, Gain):   # spec[blNum, chNum, timeNum], Gain[antNum, chNum, timeNum]
    blnum, chNum, timeNum = spec.shape[0], spec.shape[1], spec.shape[2]
    ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
    return spec / (Gain0[ant0]* Gain1[ant1].conjugate())
	#gainCalXX = np.zeros([blnum, chNum, timeNum], dtype=complex)
	#for bl_index in range(blnum):
	#	ants = Bl2Ant(bl_index); ant1, ant2 = ants[1], ants[0]
	#	Gain_bl = Gain[ant2] *  Gain[ant1].conjugate()
	#	gainCalXX[bl_index] = spec[bl_index] / Gain_bl
	#return gainCalXX
#
def gainCalVis(vis, Gain1, Gain0):      # vis[blNum, timeNum], Gain[antNum, timeNum]
    blNum, timeNum = vis.shape[0], vis.shape[1]
    ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
    return vis / (Gain0[ant0]* Gain1[ant1].conjugate())
#
def P2P(vector):
	return np.max(vector) - np.min(vector)

def PeakExcess(vector):
	meanVec = np.mean(vector)
	return max((np.max(vector) - meanVec), (meanVec - np.min(vector)))

def GetBunchedVis(msfile, ant1, ant2, pol, spw, field, chBunch, timeBunch):
	timeXY, dataXY = GetVisibity(msfile, ant1, ant2, pol, spw, field)
	return specBunch(dataXY, chBunch, timeBunch)


def bandpassCorrection(msfile, antnum, pol, spw, field, chBunch, timeBunch):
	blnum = antnum* (antnum - 1) / 2
#-------- First baseline to get time index
	timeXY = GetTimerecord(msfile, 0, 1, pol, spw, field)
	chNum, chWid  = GetChNum(msfile, spw)

	scanEnd   = np.where(np.diff(timeXY) > 5.0)[0]
	scanStart = np.append(0, scanEnd + 1)
	scanEnd   = np.append(scanEnd, 2* max(scanEnd) - scanEnd[len(scanEnd)-2])
	timeNum = scanEnd.max() + 1

	XX = np.zeros([blnum, chNum/chBunch, timeNum/timeBunch], dtype=complex)
	print ' -- Reading Visibility Data for SPW=' + `spw`
	for bl_index in range(blnum):
		ants = Bl2Ant(bl_index)
		if ants[1] == 0:
			print '    Baseline ' + `ants[0]`,
		if ants[1] == ants[0] - 1 :
			print '-' + `ants[1]`
		else:
			print '-' + `ants[1]`,
		bl_index = Ant2Bl(ants[1], ants[0])
		XX[bl_index] = GetBunchedVis(msfile, ants[1], ants[0], pol, spw, field, chBunch, timeBunch)

	timeXY  = np.mean(timeXY[range(timeNum)].reshape(timeBunch, timeNum/timeBunch), 0)
	timeNum /= timeBunch
	chNum   /= chBunch
#-------- Delay Determination for Antenna
	print ' -- Delay Calibration'
	delayCalXX = np.zeros([blnum, chNum, timeNum], dtype=complex)
	delay_bl, delay_ant = np.zeros(blnum), np.zeros([antnum, timeNum])
	startCH, stopCH = int(0.1* chNum), int(0.9* chNum) 
	for time_index in range(timeNum):
	 	for bl_index in range(blnum):
 			delay_bl[bl_index] = delay_search(XX[bl_index, startCH:stopCH, time_index])
		delay_ant[:,time_index] = cldelay_solve(delay_bl, [np.std(delay_bl)]*blnum)[0]
		#
		#---- Apply Delay Correction for visibilities
 		for bl_index in range(blnum):
			ants = Bl2Ant(bl_index); ant1, ant2 = ants[1], ants[0]
			delayCalXX[bl_index, :, time_index] = delay_cal(XX[bl_index, :, time_index], delay_ant[ant2, time_index] - delay_ant[ant1, time_index])

#-------- Overall Bandpass Calibration
	print ' -- Bandpass Calibration'
	BP_ant = bpPhsAnt(delayCalXX)
	bpCalXX = bpCal(delayCalXX, BP_ant)

#-------- Gain Calibration
	print ' -- Gain Calibration'
	Gain_ant = gainAnt(bpCalXX[:, startCH:stopCH,:])
	gainCalXX = gainCal(bpCalXX, Gain_ant)

#-------- Integrated Spectra
	print ' -- Time Integration'
	BP_bl = np.zeros([blnum, chNum], dtype=complex)
	for bl_index in range(blnum):
		BP_bl[bl_index] = blBp(gainCalXX[bl_index])
	BPAmp_ant = np.zeros([antnum, chNum])
	BPPhs_ant = np.zeros([antnum, chNum])
	BP_ant = np.zeros([antnum, chNum], dtype=complex)
	for ch_index in range(chNum):
		BPAmp_ant[:,ch_index] = clamp_solve(abs(BP_bl[:,ch_index]), [np.std(abs(BP_bl[:,ch_index]))]*blnum)[0]
		BPPhs_ant[:,ch_index] = clphase_solve(np.angle(BP_bl[:,ch_index]), [np.std(np.angle(BP_bl[:,ch_index]))]*blnum)[0]
	BP_ant = BPAmp_ant* np.exp(1j* BPPhs_ant)
	bpCalXX = bpCal(gainCalXX, BP_ant)
	return delay_ant, Gain_ant, BP_ant, bpCalXX
#
def bandpassStability(bpCalXX, segNum):
#-------- Segment Spectra
	blnum   = len(bpCalXX)
	antnum  = Bl2Ant(blnum)[0]
	chNum   = len(bpCalXX[0])
	timeNum = len(bpCalXX[0,0])
	startCH, stopCH = int(0.1* chNum), int(0.9* chNum) 
	seg_len = int(timeNum / segNum)					# Number of time in a segment
	seg_blXX = np.zeros([blnum, chNum], dtype=complex)			# Baseline-based integrated spectrum
	segAmp_ant = np.zeros([antnum, chNum])						# Amplitude of antenna-based spectrum
	segPhs_ant = np.zeros([antnum, chNum])						# Phase of antenna-based spectrum
	segXX = np.zeros([segNum, antnum, chNum], dtype=complex)	# Antenna-based integrated spectrum
	sd_spec = np.zeros([antnum, segNum]); pp_spec = np.zeros([antnum, segNum]); pe_spec = np.zeros([antnum, segNum])	# For statistics
	sd_phas = np.zeros([antnum, segNum]); pp_phas = np.zeros([antnum, segNum]); pe_phas = np.zeros([antnum, segNum])	# For statistics
	#
	for seg_index in range(segNum):
		for bl_index in range(blnum):
			seg_blXX[bl_index] = blBp(bpCalXX[bl_index,:,(seg_index*seg_len):((seg_index + 1)*seg_len)])
		for ch_index in range(chNum):
			segAmp_ant[:,ch_index] = clamp_solve(abs(seg_blXX[:,ch_index]), [np.std(abs(seg_blXX[:,ch_index]))]*blnum)[0]
			segPhs_ant[:,ch_index] = clphase_solve(np.angle(seg_blXX[:,ch_index]), [np.std(np.angle(seg_blXX[:,ch_index]))]*blnum)[0]
		segXX[seg_index] = segAmp_ant* np.exp(1j* segPhs_ant)
		#
		for ant_index in range(antnum):
			sd_spec[ant_index, seg_index] = np.std(abs(segXX[seg_index,ant_index,startCH:stopCH]))
			pp_spec[ant_index, seg_index] = P2P(abs(segXX[seg_index,ant_index,startCH:stopCH]))
			pe_spec[ant_index, seg_index] = PeakExcess(abs(segXX[seg_index,ant_index,startCH:stopCH]))
			sd_phas[ant_index, seg_index] = np.std(np.angle(segXX[seg_index,ant_index,startCH:stopCH]))
			pp_phas[ant_index, seg_index] = P2P(np.angle(segXX[seg_index,ant_index,startCH:stopCH]))
			pe_phas[ant_index, seg_index] = PeakExcess(np.angle(segXX[seg_index,ant_index,startCH:stopCH]))
		#
	return sd_spec, pp_spec, pe_spec, sd_phas, pp_phas, pe_phas
#
#-------- Smoothing complex vector
def splineComplex( axis, samplePoints, vector ):
    Samp = np.std( abs( vector ))             # Smoothing factor for amplitude
    Sphs = Samp / mean(abs(vector))   # Smoothing factor for phase
    vec_abs = abs( vector )
    vec_phs = np.angle( vector )
    vec_real = np.cos(vec_phs)
    vec_imag = np.sin(vec_phs)
    SPL_amp = UnivariateSpline(samplePoints, abs(vector), s=0.01*Samp)
    SPL_real = UnivariateSpline(samplePoints, vector.real, s=0.01*Sphs)
    SPL_imag = UnivariateSpline(samplePoints, vector.imag, s=0.01*Sphs)
    smth_phs = np.arctan2( SPL_imag(axis), SPL_real(axis) )
    return( SPL_amp(axis)* exp(1.0j* smth_phs) )
#
#-------- Van Vleck Correction
def loadVanvQ4( File ):
	VanvData = loadtxt(File)
	analogPower = VanvData[0]
	Q4Power     = VanvData[1]
	return interp1d(Q4Power, analogPower) 
#
#
def loadVanvQ3( File ):
	VanvData = loadtxt(File)
	Q3Power     = VanvData[0]
	Ratio 		= VanvData[1]
	Qeff 		= VanvData[2]
	return interp1d(Q3Power, Ratio), interp1d(Q3Power, Qeff) 
#
#
def loadAcorrCoeff( calFile ):
	ACA_Caldata = loadtxt(calFile)
	analogPower = ACA_Caldata[0]
	scaledPower = ACA_Caldata[1]
	Q4VanvPower = ACA_Caldata[2]
	scaleCoeff  = ACA_Caldata[3:8]
	vanvCoeff   = ACA_Caldata[8:13]
	scaleFact = Q4VanvPower[10]
	Q4VanvPower = Q4VanvPower / scaleFact
	return( [interp1d(Q4VanvPower, vanvCoeff[0], kind='cubic'), interp1d(Q4VanvPower, vanvCoeff[1], kind='cubic'), interp1d(Q4VanvPower, vanvCoeff[2], kind='cubic'), interp1d(Q4VanvPower, vanvCoeff[3], kind='cubic'), interp1d(Q4VanvPower, vanvCoeff[4], kind='cubic')] )
#
def Polynomial( Qpower, Ppower, coeff ):
	order = len(coeff) - 1
	result = coeff[order]( Ppower )
	for index in range( (order - 1), -1, -1):
		result = coeff[index](Ppower) + Qpower* result
	#
	return result
#
def residSimpleGauss(param, x, y):
    return y - param[0]* np.exp( -0.5* ((x - param[1])/param[2])**2 )
#
def simpleGaussFit( x, y ):
    param = [np.max(y), 0.0, np.std(x)]
    result = scipy.optimize.leastsq( residSimpleGauss, param, args=(x, y))
    return result[0]
#
#-------- Residual from Gaussian, used in fitGauss
def residGauss(param, x, y, yerr):
	# param[3]: amplitude, mean, sd, bias, rate of Gaussian function
	return (y - param[0]*np.exp( -0.5* ((x - param[1])/ param[2])**2 ) - param[3] - param[4]*x) / yerr

#-------- Gauss fit
def fitGauss( x, y, yerr ):
	param = [np.max(y), 0.0, 0.2*np.max(x), np.median(y), (y[len(y)-1] - y[0])/(x[len(x)-1] - x[0]) ]
	result = scipy.optimize.leastsq( residGauss, param, args=(x, y, yerr), full_output=True)
	return result[0], np.sqrt([result[1][0,0], result[1][1,1], result[1][2,2], result[1][3,3], result[1][4,4]])
#
#-------- 2-D Gauss fit
def resid2DGauss(param, z, x, y):
    argGauss = ((x - param[1])/param[3])**2 + ((y - param[2])/param[4])**2
    return z - param[0]* np.exp( -0.5* argGauss )
#
def simple2DGaussFit( z, x, y ):
    param = [np.max(y), 0.0, 0.0, np.std(x), np.std(y)]
    result = scipy.optimize.leastsq(resid2DGauss, param, args=(z, x, y))
    return result[0]
#
#-------- ACD edge pattern
def ACDedge( timeXY ):
    timeSKY = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE")
    timeAMB = msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT")
    timeHOT = msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT")
    #
    skyRange = range(2, edge[0]-1)
    ambRange = range(edge[0]+3, edge[1]-1)
    hotRange = range(edge[1]+3, len(timeXY)-1)
    return skyRange, ambRange, hotRange
#
#-------- 3-bit VanVleck Correction
def Vanv3bitCorr( dataXY, refRange, Qeff ):
	refZeroLag = np.mean(dataXY[:, refRange].real)
	for index in range(dataXY.shape[1]):
		temp = dataXY[:, index].real
		ZeroLag = np.mean(temp)
		VanvCorrect = Qeff( ZeroLag / refZeroLag)
		dataXY[:,index] = temp / VanvCorrect
	#
	return dataXY
#
#
#-------- Tsys from ACD
def TsysSpec(msfile, pol, TsysScan, spw, vanvSW):
    #if vanvSW:
    #	VanvQ3, Qeff = loadVanvQ3('/users/skameno/Scripts/ACAPowerCoeff.data')  # 3-bit Van Vleck
    #
	antList = GetAntName(msfile)
	antNum  = len(antList)
	polNum  = len(pol)
	chNum, chWid, freq = GetChNum(msfile, spw)
	TrxSp = np.zeros([antNum, polNum, chNum]); TsysSp = np.zeros([antNum, polNum, chNum])
	chRange = range(int(0.1*chNum), int(0.9*chNum))
	#
	text_sd =  '%s SPW=%d SCAN=%d' % (msfile, spw, TsysScan)
	print text_sd
	#-------- Loop for Antenna
	for ant_index in range(antNum):
		#-------- Get Physical Temperature of loads
		tempAmb, tempHot = GetLoadTemp(msfile, ant_index, spw)
		#
		for pol_index in range(polNum):
			timeXY, dataXY = GetVisibility(msfile, ant_index, ant_index, pol[pol_index], spw, TsysScan)
			#-------- Time range of Sky/Amb/Hot
			skyRange, ambRange, hotRange = ACDedge(timeXY)
			#-------- Van Vleck Correction
			if vanvSW:
				dataXY = Vanv3bitCorr(dataXY, ambRange, Qeff)
			#
			#-------- Calc. Tsys Spectrum
			Psky, Pamb, Phot = np.mean(dataXY[:,skyRange].real, 1), np.mean(dataXY[:,ambRange].real, 1), np.mean(dataXY[:,hotRange].real, 1)
			TrxSp[ant_index, pol_index]  = (tempHot* Pamb - Phot* tempAmb) / (Phot - Pamb)
			TsysSp[ant_index, pol_index] = (Psky* tempAmb) / (Pamb - Psky)
			text_sd = '%s PL%s %5.1f %5.1f' % (antList[ant_index], pol[pol_index], np.median(TrxSp[ant_index, pol_index]), np.median(TsysSp[ant_index, pol_index]))
			print text_sd
			# logfile.write(text_sd + '\n')
		#
	#
	return TrxSp, TsysSp
#
#-------- Power Spectra
def PowerSpec(msfile, pol, TsysScan, TsysSPW, vanvSW):
	VanvQ3, Qeff = loadVanvQ3('/users/skameno/Scripts/ACAPowerCoeff.data')  # 3-bit Van Vleck
	antList = GetAntName(msfile)
	antNum  = len(antList)
	polNum  = len(pol)
	chNum, chWid, freq = GetChNum(msfile, TsysSPW); Tsysfreq = freq* 1.0e-9 # GHz
	PskySpec = np.zeros([antNum, polNum, chNum])
	PambSpec = np.zeros([antNum, polNum, chNum])
	PhotSpec = np.zeros([antNum, polNum, chNum])
	#-------- Loop for antennas
	for ant_index in range(antNum):
		#-------- Get Physical Temperature of loads
		tempAmb, tempHot = GetLoadTemp(msfile, ant_index, TsysSPW)
		#
		#-------- Loop for polarizations
		for pol_index in range(polNum):
			timeXY, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, TsysSPW, TsysScan)
			#
			#-------- Time range of Sky/Amb/Hot
			skyRange, ambRange, hotRange = ACDedge(timeXY)
			#-------- Van Vleck Correction
			if vanvSW:
				dataXY = Vanv3bitCorr(dataXY, ambRange, Qeff)
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
#-------- Bandpass Plog
def PlotBP(msfile, antList, spwList, BPList):
    spwNum, antNum = len(spwList), len(antList)
    #-------- Prepare Plots
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index, figsize = (11, 8))
        figAnt.suptitle(prefix + ' ' + antList[ant_index])
        figAnt.text(0.45, 0.05, 'Frequency [GHz]')
        figAnt.text(0.03, 0.45, 'Bandpass Amplitude and Phase', rotation=90)
    #
    #-------- Plot BP
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index)
        for spw_index in range(spwNum):
            chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* Freq  # GHz
            BPampPL = figAnt.add_subplot( 2, spwNum, spw_index + 1 )
            BPphsPL = figAnt.add_subplot( 2, spwNum, spw_index + spwNum + 1 )
            for pol_index in range(ppolNum):
                plotBP = BPList[spw_index][ant_index, pol_index]
                BPampPL.plot( Freq, abs(plotBP), ls='steps-mid', label = 'Pol=' + PolList[pol_index])
                BPampPL.axis([np.min(Freq), np.max(Freq), 0.0, 1.25])
                BPampPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
                BPampPL.yaxis.offsetText.set_fontsize(10)
                BPampPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
                BPphsPL.plot( Freq, np.angle(plotBP), '.', label = 'Pol=' + PolList[pol_index])
                BPphsPL.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
            #
            BPampPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            BPphsPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
            BPampPL.text( np.min(Freq), 1.1, 'SPW=' + `spwList[spw_index]` + ' Amp')
            BPphsPL.text( np.min(Freq), 2.5, 'SPW=' + `spwList[spw_index]` + ' Phase')
        #
    #
    return figAnt
#
#
#-------- Amp and Phase Plot
def plotAmphi(fig, freq, spec):
	amp = abs(spec);	phase = np.angle(spec)
	#ampAxis = fig.add_axes((0, 0.2, 1, 0.8))
	#Panel = fig.add_subplot(111)
	ampAxis = fig.add_axes((0.1, 0.3, 0.89, 0.6))	# left, bottom, width, height
	phsAxis = fig.add_axes((0.1, 0.1, 0.89, 0.2), sharex = ampAxis)
	ampAxis.tick_params(labelbottom = 'off')
	ampAxis.plot(freq, amp, ls='steps-mid')
	ampAxis.axis( [min(freq), max(freq), 0, 1.1* max(amp)], size='x-small' )
	phsAxis.plot(freq, phase, '.')
	phsAxis.axis( [min(freq), max(freq), -pi, pi], size='x-small' )
	return
#
def gainComplex( bl_vis ):
    blnum  =  len(bl_vis)
    antnum =  Bl2Ant(blnum)[0]
    #
    resid  =  np.zeros(2* blnum)
    correction = np.ones(2* antnum - 1)
    solution   = np.zeros(2* antnum - 1)
    #
    #---- Initial solution
    solution[0] = sqrt(abs(bl_vis[0]))		# Refant has only real part
    for ant_index in range(1, antnum) :
        solution[ant_index]			= bl_vis[Ant2Bl(0, ant_index )].real / solution[0]
        solution[antnum + ant_index - 1]= bl_vis[Ant2Bl(0, ant_index )].imag / solution[0]
    #
    for iter_index in range(2):
        complex_matrix = np.zeros((2*blnum, 2*antnum - 1))
        #-------- Residual Vector
        for bl_index in range(blnum):
            ants = Bl2Ant(bl_index)
            if ants[1] != 0:
                resid[bl_index]			= bl_vis[bl_index].real - (solution[ants[0]]* solution[ants[1]] + solution[antnum + ants[0] - 1]* solution[antnum + ants[1] - 1])	# Real part
                resid[blnum + bl_index] = bl_vis[bl_index].imag - (solution[ants[1]]* solution[antnum + ants[0] - 1] - solution[ants[0]]* solution[antnum + ants[1] - 1])	# Imag part
            else:
                resid[bl_index]			= bl_vis[bl_index].real - (solution[ants[0]]* solution[0])	# Real part
                resid[blnum + bl_index] = bl_vis[bl_index].imag - (solution[antnum + ants[0] - 1]* solution[0])	# Imag part
        #
        #---- Partial Matrix
        for bl_index in range(blnum):
            ants = Bl2Ant(bl_index)
            complex_matrix[bl_index, ants[0]] = solution[ants[1]]
            complex_matrix[bl_index, ants[1]] = solution[ants[0]]
            if ants[1] != 0:
                complex_matrix[bl_index, antnum + ants[0] - 1] = solution[antnum + ants[1] - 1]
                complex_matrix[bl_index, antnum + ants[1] - 1] = solution[antnum + ants[0] - 1]
                complex_matrix[blnum + bl_index, ants[1]] =  solution[antnum + ants[0] - 1]
                complex_matrix[blnum + bl_index, ants[0]] = -solution[antnum + ants[1] - 1]
                complex_matrix[blnum + bl_index, antnum + ants[1] - 1] = -solution[ants[0]]
                complex_matrix[blnum + bl_index, antnum + ants[0] - 1] =  solution[ants[1]]
            else:		# No ants[1]
                complex_matrix[blnum + bl_index, 0]		= solution[antnum + ants[0] - 1]
                complex_matrix[blnum + bl_index, antnum + ants[0] - 1]= solution[0]
        #
        ptp = np.dot(complex_matrix.T, complex_matrix)
        ptp_inv   = scipy.linalg.inv(ptp)
        correction = np.dot(ptp_inv,  np.dot(complex_matrix.T, resid))
        solution   = np.add(solution, correction)
    #
    return solution[range(antnum)] + 1j* np.append(0, solution[range(antnum, 2*antnum-1)])
#
#-------- Function to calculate visibilities
def polariVis( Xspec ):     # Xspec[polNum, blNum, chNum, timeNum]
    blNum, chNum, timeNum   = Xspec.shape[1], Xspec.shape[2], Xspec.shape[3]
    chRange = range( int(chNum*0.06), int(chNum* 0.96))
    #-------- Visibilities
    XX = Xspec[0]     # XX[BL, CH, TIME]
    XY = Xspec[1]     # XY[BL, CH, TIME]
    YX = Xspec[2]     # YX[BL, CH, TIME]
    YY = Xspec[3]     # YY[BL, CH, TIME]
    #-------- Bandpass table
    print '--- Making Antenna-based Bandpass Table'
    BPX = bpPhsAnt(XX)          # BPX[ANT, CH] (phase only)
    BPY = bpPhsAnt(YY)          # BPY[ANT, CH] (phase only)
    #-------- Bandpass phase correction
    print '--- Applying Bandpass Calibration'
    XXbpcal = bpCal(XX, BPX, BPX)    # XXbpcal[BL, CH, TIME]
    YYbpcal = bpCal(YY, BPY, BPY)    # YYbpcal[BL, CH, TIME]
    XYbpcal = bpCal(XY, BPX, BPY)    # XYbpcal[BL, CH, TIME]
    YXbpcal = bpCal(YX, BPY, BPX)    # YXbpcal[BL, CH, TIME]
    #-------- channel average
    print '--- Channel-averaging visibilities'
    if len(chRange) == 0:
        XXchav = XXbpcal[:,0]
        XYchav = XYbpcal[:,0]
        YXchav = YXbpcal[:,0]
        YYchav = YYbpcal[:,0]
    else:
        XXchav = np.mean( XXbpcal[:,chRange], axis=1)  # XXchav[BL, TIME]
        XYchav = np.mean( XYbpcal[:,chRange], axis=1)  # XXchav[BL, TIME]
        YXchav = np.mean( YXbpcal[:,chRange], axis=1)  # XXchav[BL, TIME]
        YYchav = np.mean( YYbpcal[:,chRange], axis=1)  # XXchav[BL, TIME]
    #
    #-------- Gain Calibration
    print '--- Solution for antenna-based gain'
    GainX = np.apply_along_axis( gainComplex, 0, XXchav )
    GainY = np.apply_along_axis( gainComplex, 0, YYchav )
    VisXX = np.mean(gainCalVis( XXchav, GainX, GainX ), axis = 0)
    VisYY = np.mean(gainCalVis( YYchav, GainY, GainY ), axis = 0)
    VisXY = np.mean(gainCalVis( XYchav, GainX, GainY ), axis = 0)
    VisYX = np.mean(gainCalVis( YXchav, GainY, GainX ), axis = 0)
    #
    return GainX, GainY, VisXX, VisXY, VisYX, VisYY
#
#-------- Determine antenna-based gain with polarized source
def polariGain( XX, YY, PA, StokesQ, StokesU):
    blNum, timeNum = XX.shape[0], XX.shape[1]
    csPA = np.cos(2.0* PA)
    snPA = np.sin(2.0* PA)
    Xscale = 1.0 / (1.0 + StokesQ* csPA + StokesU* snPA)
    Yscale = 1.0 / (1.0 - StokesQ* csPA - StokesU* snPA)
    #
    ScaleXX = np.dot(XX, np.diag(Xscale))
    ScaleYY = np.dot(YY, np.diag(Yscale))
    #
    GainX = np.apply_along_axis( gainComplex, 0, ScaleXX)
    GainY = np.apply_along_axis( gainComplex, 0, ScaleYY)
    return GainX, GainY
#
def XY2Stokes(PA, VisXY, VisYX):
    #-------- Least-Square fit for polarizatino parameters (Q, U, XYphase, Dx, Dy)
    timeNum = len(VisXY)
    sinPA2 = np.sin(2.0*PA)
    cosPA2 = np.cos(2.0*PA)
    P = np.zeros([7, 4* timeNum])       # 7 parameters to solve, 4 (ReXY, ImXY, ReYX, ImYX) * timeNum measurements
    solution = np.array([0.1, 0.1, np.angle( np.mean(VisXY[1])), 0.0, 0.0, 0.0, 0.0])   # Initial parameters : StokesQ, StokesU, XYphase, Re(Dx+Dy*), Im(Dx+Dy*)
    W = np.diag(4.0* np.ones([4* timeNum])/ np.var(VisXY.imag + VisYX.imag))
    #-------- Iteration loop
    for index in range(10):
        sinPhi = np.sin(solution[2])
        cosPhi = np.cos(solution[2])
        UC_QS = cosPA2* solution[1] - sinPA2* solution[0]   # U cosPA2 - Q sinPA2
        modelVis = np.r_[
             cosPhi* UC_QS + solution[3],       # Re XY*
             sinPhi* UC_QS + solution[4],       # Im XY*
             cosPhi* UC_QS + solution[5],       # Re YX*
            -sinPhi* UC_QS + solution[6] ]      # Im YX*
        #-------- Partial matrix
        P[0] = np.r_[-sinPA2* cosPhi, -sinPA2* sinPhi, -sinPA2* cosPhi,  sinPA2* sinPhi]
        P[1] = np.r_[ cosPA2* cosPhi,  cosPA2* sinPhi,  cosPA2* cosPhi, -cosPA2* sinPhi]
        P[2] = np.r_[-sinPhi* UC_QS, cosPhi* UC_QS, -sinPhi* UC_QS, -cosPhi* UC_QS]
        P[3] = np.r_[np.ones([timeNum]),  np.zeros([timeNum]), np.zeros([timeNum]), np.zeros([timeNum])]
        P[4] = np.r_[np.zeros([timeNum]), np.ones([timeNum]),  np.zeros([timeNum]), np.zeros([timeNum])]
        P[5] = np.r_[np.zeros([timeNum]), np.zeros([timeNum]), np.ones([timeNum]),  np.zeros([timeNum])]
        P[6] = np.r_[np.zeros([timeNum]), np.zeros([timeNum]), np.zeros([timeNum]), np.ones([timeNum])]
        PTWP_inv = scipy.linalg.inv(np.dot(P, np.dot(W, P.T)))
        vecVis = np.r_[ VisXY.real, VisXY.imag, VisYX.real, VisYX.imag ]
        residual = vecVis - modelVis
        correction = np.dot( PTWP_inv, np.dot (P, np.dot(W, residual)))
        # print 'Iteration %d : correction = %e' % (index, np.dot(correction,correction))
        solution   = solution + correction
        if np.dot(correction,correction) < 1.0e-15:
            break
        #
    #
    solution[2] = np.arctan2( np.sin(solution[2]), np.cos(solution[2]) )    # Remove 2pi ambiguity
    return(solution, np.sqrt(np.diag(PTWP_inv)))
#
#-------- GridPoint
def GridPoint( value, samp_x, samp_y, point_x, point_y, kernel ):
    #---- Check NaN and replace with 0
    nan_index = np.where( value != value )[0]
    value[nan_index] = 0.0
    #---- Distance from gridding points
    dist_sq = (samp_x - point_x)**2 + (samp_y - point_y)**2
    dist_thresh = 9.0 * kernel**2
    index = np.where( dist_sq < dist_thresh)[0]
    wt = exp( -0.5* dist_sq[index] / kernel**2 )
    nan_index = np.where( value[index] != value[index] )[0]
    wt[nan_index] = 0.0
    sumWt = np.sum(wt)
    if sumWt < 1.0e-3:
        return 0.0
    #
    return np.sum(value[index]* wt) / sumWt
#
#-------- GridData
def GridData( value, samp_x, samp_y, grid_x, grid_y, kernel ):
    gridNum = len(grid_x)
    results = np.zeros(gridNum)
    for index in range(gridNum):
     results[index] = GridPoint( value, samp_x, samp_y, grid_x[index], grid_y[index], kernel)
    #
    return results
#
#-------- Disk Visibility
def diskVis(diskRadius, u):
    # diskRadius : radius of planet disk [rad]
    # u          : spatial frequency (= baseline / wavelength)
    argument = 2.0* pi* u* diskRadius
    return 2.0* scipy.special.jn(1, argument) / argument
#
#-------- Disk Visibility with primary beam correction, u must be smaller than 0.3/diskRadius
def diskVisBeam(diskRadius, u, primaryBeam):
    # diskRadius : radius of planet disk [rad]
    # u          : spatial frequency (= baseline / wavelength)
    # primaryBeam: FWHM of primary beam [rad]
    disp_u = (0.30585 / diskRadius)**2 + 2.0* log(2.0)/(pi* primaryBeam)**2
    return beamF(diskRadius/primaryBeam)* np.exp(-0.5* (u**2 / disp_u))
#
