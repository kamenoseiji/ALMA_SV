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
from scipy.sparse import lil_matrix
import scipy.optimize
import time
import datetime
#======== Baseline and Antenna Indexing
KERNEL_BL = arange(64)*arange(1,65)/2
def indexList( refArray, motherArray ):     # Compare two arrays and return matched index
    IL = []
    for currentItem in refArray: IL = IL + np.where( motherArray == currentItem )[0].tolist()
    return IL
#
def timeMatch( refTime, scanTime, thresh): # Time-based matching
    match = np.where( abs(scanTime - refTime) < thresh)[0].tolist()
    return len(match)
#
def Ant2Bl(ant1, ant2):	    # Antenna -> baseline index (without autocorr)
    antenna1 = max(ant1, ant2); antenna2 = min(ant1, ant2)
    return antenna1* (antenna1 - 1)/2 + antenna2
#
def Ant2BlD(ant1, ant2):    # Antenna -> baseline index and direction (True if inverted)
    antenna1 = max(ant1, ant2); antenna2 = min(ant1, ant2)
    return antenna1* (antenna1 - 1)/2 + antenna2, (ant1 < ant2)
#
def Bl2Ant(bl_index):     # Baseline -> antenna indexing (canonical ordering)
    ant1 = max(np.where(KERNEL_BL<= bl_index)[0]) + 1
    return ant1, bl_index - KERNEL_BL[ant1 - 1]
#
def revList(inList):
    listLen = len(inList)
    outList = []
    for index in range(listLen):
        outList.append( inList.index(index) )
    #
    return outList
#
ANT0 = []; ANT1 = []     # List the BL -> antenna indexing
for bl_index in range(2016):    # Maximum number of baseline
    ants = Bl2Ant(bl_index)
    ANT0.append(ants[0])        # bl -> ant0 (baseline-end antenna) mapping
    ANT1.append(ants[1])        # bl -> ant1 (baseline-begin antenna) mapping
#
def Ant2Bla_RevLex(ant0, ant1, antNum):    # Reverse Lexical, with autcorr
    antenna0 = min(ant0, ant1); antenna1 = max(ant0, ant1)
    kernel = antNum* antenna0 - antenna0* (antenna0 - 1)/2
    return kernel + antenna1 - antenna0
#
def Ant2Bl_RevLex(ant1, ant2, antnum):    # Reverse Lexical, without autcorr
    antenna1 = max(ant1, ant2); antenna2 = min(ant1, ant2)
    return int(antnum* antenna2 - (antenna2 + 1)* (antenna2 + 2) / 2  + antenna1)
#
def subArrayIndex(Flag):          #-------- SubArray Indexing
    blNum = len(Flag); antNum = Bl2Ant(blNum)[0]; kernelBL = KERNEL_BL[range(antNum-1)].tolist()
    ant0, ant1 = ANT0[0:blNum], ANT1[0:blNum]
    flagIndex = np.where(Flag == 1)[0].tolist()       # Baselines: uvDist < UVlimit
    flagRefIndex = list( set(flagIndex) & set(kernelBL))            # Baselines including refant and uvDist < UVlimit 
    SAantennas = [0] + list(np.array(ant0)[flagRefIndex]); SAantennas.sort()
    SAantNum = len(SAantennas); SAblNum = SAantNum* (SAantNum - 1)/2
    SAblMap = []
    for bl_index in range(SAblNum):
        SAblMap = SAblMap + [Ant2Bl(SAantennas[ant0[bl_index]], SAantennas[ant1[bl_index]])]
    #
    SAblFlag = np.zeros([SAblNum]); SAblFlag[indexList(np.array(flagIndex), np.array(SAblMap))] = 1.0
    SAant0, SAant1 = np.array(ant0)[SAblMap].tolist(), np.array(ant1)[SAblMap].tolist()
    return SAantennas, SAblMap, SAblFlag, SAant0, SAant1
#
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
def AzEl2PA(az, el, lat=-23.029/180.0*pi): # Azimuth, Elevation, Latitude (default=ALMA) in [rad]
    cos_lat = np.cos(lat)
    #return np.arctan( -cos_lat* np.sin(az) / (np.sin(lat)* np.cos(el) - cos_lat* np.sin(el)* np.cos(az)) )
    return np.arctan2( -cos_lat* np.sin(az), (np.sin(lat)* np.cos(el) - cos_lat* np.sin(el)* np.cos(az)) )
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
#======== MS data interface
def AzElMatch( refTime, scanTime, AntID, targetAnt, thresh, Az, El ):
    antTimeIndex = np.where(AntID == targetAnt)[0].tolist()
    time_pointer = np.array(antTimeIndex)[np.where( abs(scanTime[antTimeIndex] - refTime) < thresh)[0].tolist()].tolist()
    return np.median(Az[time_pointer]), np.median(El[time_pointer])
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
#
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
def GetVisCross(msfile, spwID, scanID):
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
	for bl_index in range(blNum):
		ant1, ant0 = Bl2Ant(bl_index)
		xcorr_index[bl_index] = Ant2Bla_RevLex(ant0, ant1, antNum)
	Xspec = dataXY.reshape(polNum, chNum, corrNum, timeNum)[:,:,xcorr_index,:]
	return timeStamp, Xspec
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
def antRefScan( msfile, timeRange ):    # Check scanning and tracking antennas
    antList = GetAntName(msfile)
    antNum = len(antList)
    scanRange = np.zeros(antNum)
    Time, AntID, Offset = GetAzOffset(msfile)
    for ant_index in range(antNum):
        time_index = np.where( (AntID == ant_index) )[0]
        time_index = time_index[np.where( (Time[time_index] >= timeRange[0]))[0]]
        time_index = time_index[np.where( (Time[time_index] <= timeRange[1]))[0]]
        scanRange[ant_index] = max( Offset[0, time_index] ) - min( Offset[0, time_index] )
    #
    trkAntIndex  = np.where( scanRange == 0.0 )[0]
    scanAntIndex = np.where( scanRange >  0.0 )[0]
    return trkAntIndex.tolist(), scanAntIndex.tolist(), Time, Offset
#
def GetChNum(msfile, spwID):
	tb.open(msfile + '/' + 'SPECTRAL_WINDOW')
	chNum = tb.getcell("NUM_CHAN", spwID)
	chWid = tb.getcell("CHAN_WIDTH", spwID)
	freq  = tb.getcell("CHAN_FREQ", spwID)
	tb.close()
	return chNum, chWid, freq

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
def BlPhaseMatrix(antPhs):
    antNum = len(antPhs)
    blNum = antNum* (antNum - 1) / 2
    antGain = np.exp( 1.0j * antPhs)
    blphs_matrix = np.zeros([blNum, antNum], dtype=complex)
    for bl_index in range(blNum):
        blphs_matrix[bl_index, ANT1[bl_index]] = 1.0j* antGain[ANT1[bl_index]]* antGain[ANT0[bl_index]].conjugate()
        blphs_matrix[bl_index, ANT0[bl_index]] =-1.0j* antGain[ANT1[bl_index]]* antGain[ANT0[bl_index]].conjugate()
    #
    return blphs_matrix[:,1:antNum]
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
#
def logamp_solve(bl_amp):
    blnum  =  len(bl_amp); log_bl =  np.log(bl_amp)
    antNum =  Bl2Ant(blnum)[0]
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    PTP_inv = ((2.0* antNum - 2.0)* np.diag(np.ones(antNum)) - 1.0) / (2.0* (antNum - 1.0)* (antNum - 2.0))
    PTV = np.zeros(antNum)
    for ant_index in range(antNum):
        index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
        PTV[ant_index] += np.sum(log_bl[index0]) + np.sum(log_bl[index1])
    #
    return np.exp(PTP_inv.dot(PTV))
#
def PTPmatrix(Gain):  # Gain is a vector of antenna-based gain amplitude (real)
    antNum = len(Gain); normG = Gain.dot(Gain)
    PTP = np.zeros([antNum, antNum]) + Gain
    for ant_index in range(antNum): 
        PTP[ant_index,:] *= Gain[ant_index]
        PTP[ant_index, ant_index] = normG - Gain[ant_index]**2
    #
    return PTP
#
def clamp_solve(bl_amp, niter=2):
    blnum  =  len(bl_amp)
    antGain = logamp_solve(bl_amp); antNum = len(antGain)
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    for iter_index in range(niter):
        resid = bl_amp - antGain[ant0]* antGain[ant1]
        y = np.zeros(antNum)
        for ant_index in range(antNum):
            index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
            y[ant_index] += antGain[ant1[index0]].dot(resid[index0])
            y[ant_index] += antGain[ant0[index1]].dot(resid[index1])
        #
        L = np.linalg.cholesky(PTPmatrix(antGain))
        t = np.linalg.solve(L, y)
        correction = np.linalg.solve(L.T, t)
        antGain += correction; antGain = abs(antGain)
    #
    return antGain
#

def cldelay_solve(bl_delay):    # see http://qiita.com/kamenoseiji/items/782031a0ce8bbc1dc99c
    blNum = len(bl_delay); antNum = Bl2Ant(blNum)[0]
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    PTP_inv = (np.diag(np.ones(antNum - 1)) + 1.0) / antNum
    PTY = np.zeros(antNum - 1)
    for ant_index in range(1, antNum):
        index0 = np.where(ant0 == ant_index)[0].tolist()
        index1 = np.where(ant1 == ant_index)[0].tolist()
        PTY[ant_index - 1] += np.sum(bl_delay[index0]) - np.sum(bl_delay[index1])
    #
    return np.array( [0.0] + (PTP_inv.dot(PTY)).tolist() )
#
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
def clphase_solve(Vis, iterNum = 2):
    Vis = Vis / abs(Vis)    # Normalization
    blNum = len(Vis); antNum = Bl2Ant(blNum)[0]
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    PTP_inv = (np.diag(np.ones(antNum - 1)) + 1.0) / antNum
    antGain = np.append(1.0 + 0.0j, Vis[KERNEL_BL[0:antNum-1]])
    #
    for iter_index in range(iterNum):
        resid = Vis - (antGain[ant0] * antGain[ant1].conjugate())
        print 'Iter %d: resid = %f' % (iter_index, np.sum(abs(resid)**2))
        PTY = np.zeros(antNum - 1, dtype=complex)
        for ant_index in range(1,antNum):
            index0 = np.where(ant0 == ant_index)[0].tolist()
            index1 = np.where(ant1 == ant_index)[0].tolist()
            Y = np.zeros(blNum, dtype=complex)
            Y[index0] = -1.0j* antGain[ant_index].conjugate()* antGain[ant1[index0]]
            Y[index1] =  1.0j* antGain[ant_index]* antGain[ant0[index1]].conjugate()
            PTY[ant_index - 1] += Y.dot(resid)
        #
        antPhs  = np.append(0, np.angle(antGain[1:antNum]) + PTP_inv.dot(PTY.real))
        antGain = np.cos(antPhs) + 1.0j* np.sin(antPhs)
    #
    return antPhs
#
def MullerVector(Dx0, Dy0, Dx1, Dy1, Unity):
    P = np.array([[Unity,               Dx1.conjugate(),      Dx0,                  Dx0* Dx1.conjugate()],
                  [Dy1.conjugate(),     Unity,                Dx0* Dy1.conjugate(), Dx0                 ],
                  [Dy0,                 Dx1.conjugate()* Dy0, Unity,                Dx1.conjugate()     ],
                  [Dy0*Dy1.conjugate(), Dy0,                  Dy1.conjugate(),      Unity               ]]) #.transpose(2,0,1)
    return P
#
def dMdDVec(Dx1, Dy1, Unity):
    return np.array([
        [0.0*Unity, 0.0*Unity, Unity,  Dx1.conjugate()],
        [0.0*Unity, 0.0*Unity, Dy1.conjugate(), Unity],
        [Unity, Dx1.conjugate(), 0.0*Unity, 0.0*Unity],
        [Dy1.conjugate(), Unity, 0.0*Unity, 0.0*Unity]])
#
def KMvec(Dx, Dy, Unity):
    return np.array([[ Unity, Dx ], [Dy, Unity]] )
#
def Vis2solveDD(Vis, PS):
    blNum  = Vis.shape[1]; antNum = Bl2Ant(blNum)[0]                   # (I, Q, U, V)
    ant0, ant1 = ANT0[0:blNum], ANT1[0:blNum]
    def Dresid( D ):
        antNum = len(D)/4; blNum = antNum* (antNum - 1) /2
        Dx = D[0:antNum] + (0.0+1.0j)* D[antNum:(2*antNum)]
        Dy = D[(2*antNum):(3*antNum)] + (0.0+1.0j)*D[(3*antNum):(4*antNum)]
        resid   = (Vis - np.dot(MullerVector( Dx[ant0], Dy[ant0], Dx[ant1], Dy[ant1], np.ones(blNum, dtype=complex) ).transpose(2,0,1), PS).T).reshape(4* blNum)
        resReal = np.r_[resid.real, resid.imag]
        return np.dot(resReal, resReal)
    #
    fit = scipy.optimize.minimize(Dresid, x0=np.zeros(4*antNum), method="tnc")['x']
    return fit[0:antNum] + (0.0+1.0j)*fit[antNum:(2*antNum)], fit[(2*antNum):(3*antNum)] + (0.0+1.0j)*fit[(3*antNum):(4*antNum)]
#
def VisPA_solveDM(Vis, PA, Stokes):
    blNum  = Vis.shape[1]; antNum = Bl2Ant(blNum)[0]; PAnum = len(PA)
    P = np.zeros([4*antNum, 2*blNum*PAnum])
    def Dresid(D):
        antNum = len(D)/4; blNum = antNum* (antNum - 1) /2
        Dx = D[0:antNum] + (0.0+1.0j)* D[antNum:(2*antNum)]
        Dy = D[(2*antNum):(3*antNum)] + (0.0+1.0j)*D[(3*antNum):(4*antNum)]
        resid = np.array([], dtype=complex)
        for PA_index in range(PAnum):
            PS = np.dot(PAMatrix(PA[PA_index]), Stokes)
            resid = np.append(resid, (Vis[:,:,PA_index] - np.dot(MullerVector( Dx[ant0], Dy[ant0], Dx[ant1], Dy[ant1], np.ones(blNum, dtype=complex) ).transpose(2,0,1), PS).T).reshape(4* blNum))
        #
        resReal = np.r_[resid.real, resid.imag]
        return resReal
    #
    P[0] = 1.0e3* (Dresid() - Dresid(np.zeros(4*antNum)))
#
def VisPA_solveD(Vis, PA, Stokes):
    PAnum, blNum = len(PA), Vis.shape[1]; antNum = Bl2Ant(blNum)[0]; PABLnum = PAnum* blNum
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    CS, SN = np.cos(2.0* PA), np.sin(2.0*PA)
    QCpUS = Stokes[1]*CS + Stokes[2]*SN
    UCmQS = Stokes[2]*CS - Stokes[1]*SN
    Unity = np.ones(PAnum); Zeroty = np.zeros(PAnum)
    resid = np.zeros(4* PABLnum, dtype=complex)
    resid[0:PABLnum]          = (Vis[0] - Unity - QCpUS).reshape(PABLnum)
    resid[PABLnum:2*PABLnum]  = (Vis[1] - UCmQS).reshape(PABLnum)
    resid[2*PABLnum:3*PABLnum]= (Vis[2] - UCmQS).reshape(PABLnum)
    resid[3*PABLnum:4*PABLnum]= (Vis[3] - Unity + QCpUS).reshape(PABLnum)
    #
    P = np.zeros([4, antNum, 4, blNum, PAnum], dtype=complex)   # [(Dx0,Dy0,Dx1,Dy1), ant, (XX,XY,YX,YY), bl, PA]
    for ant_index in range(antNum):
        index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
        P[0, ant_index, 0, index0] = UCmQS      # XX / Dx0
        P[2, ant_index, 0, index1] = UCmQS      # XX / Dx1

        P[1, ant_index, 1, index0] = Unity + QCpUS      # XY / Dy0
        P[2, ant_index, 1, index1] = Unity - QCpUS      # XY / Dx1

        P[0, ant_index, 2, index0] = Unity - QCpUS      # YX / Dx0
        P[3, ant_index, 2, index1] = Unity + QCpUS      # YX / Dy1

        P[2, ant_index, 3, index0] = UCmQS      # YX / Dx0
        P[3, ant_index, 3, index1] = UCmQS      # YX / Dy1
    #
    PM = P.reshape(4*antNum, 4*blNum*PAnum)
    PTP = PM.dot(PM.T)
    #def Dresid(D):
    #    antNum = len(D)/4; blNum = antNum* (antNum - 1) /2
    #    Dx = D[0:antNum] + (0.0+1.0j)* D[antNum:(2*antNum)]
    #    Dy = D[(2*antNum):(3*antNum)] + (0.0+1.0j)*D[(3*antNum):(4*antNum)]
    #    resid = np.array([], dtype=complex)
    #    for PA_index in range(PAnum):
    #        PS = np.dot(PAMatrix(PA[PA_index]), Stokes)
    #        resid = np.append(resid, (Vis[:,:,PA_index] - np.dot(MullerVector( Dx[ant0], Dy[ant0], Dx[ant1], Dy[ant1], np.ones(blNum, dtype=complex) ).transpose(2,0,1), PS).T).reshape(4* blNum))
    #    #
    #    resReal = np.r_[resid.real, resid.imag]
    #    return np.dot(resReal, resReal)
    ##
    #fit = scipy.optimize.minimize(Dresid, x0=np.zeros(4*antNum), method="tnc")['x']
    return fit[0:antNum] + (0.0+1.0j)*fit[antNum:(2*antNum)], fit[(2*antNum):(3*antNum)] + (0.0+1.0j)*fit[(3*antNum):(4*antNum)]
#
def Vis2solveDDD(Vis, PS):
    blNum  = Vis.shape[1]; antNum = Bl2Ant(blNum)[0]                   # (I, Q, U, V)
    ant0, ant1 = ANT0[0:blNum], ANT1[0:blNum]
    dSolMap = range(2* antNum) + range(2*antNum+1, 4*antNum)
    Dx, Dy, Unity = np.zeros(antNum, dtype=complex), np.zeros(antNum, dtype=complex), np.ones(blNum, dtype=complex)     # Dx, Dy solutions for scanning antenna
    #weight = np.r_[0.01*np.ones(blNum), np.ones(blNum), np.ones(blNum), 0.01*np.ones(blNum), 0.01*np.ones(blNum), np.ones(blNum), np.ones(blNum), 0.1*np.ones(blNum)]
    weight = np.ones(8*blNum)
    for loop_index in range(2):
        resid   = (Vis - np.dot(MullerVector( Dx[ant0], Dy[ant0], Dx[ant1], Dy[ant1], Unity ).transpose(2,0,1), PS).T).reshape(4* blNum)
        resReal = np.r_[resid.real, resid.imag]
        #
        P = np.zeros( [4*antNum, 8*blNum] )
        K0, K1 = KMvec(Dx[ant0], Dy[ant0], Unity).transpose(2,0,1), KMvec(Dx[ant1], Dy[ant1], Unity).transpose(2,0,1)
        dVisX0 = np.dot(K1, PS[2:4])  # dM/dX0
        dVisY0 = np.dot(K1, PS[0:2])  # dM/dY0
        dVisX1 = np.dot(K0, PS[[1,3]])# dM/dX1
        dVisY1 = np.dot(K0, PS[[0,2]])# dM/dX1
        #
        for bl_index in range(blNum):
            a0, a1 = ant0[bl_index], ant1[bl_index]
            #-------- Derivative by real
            P[0*antNum + a0, 0*blNum + bl_index] += dVisX0[bl_index,0].real # /reX0
            P[0*antNum + a0, 1*blNum + bl_index] += dVisX0[bl_index,1].real # /reX0
            P[1*antNum + a0, 2*blNum + bl_index] += dVisY0[bl_index,0].real # /reY0
            P[1*antNum + a0, 3*blNum + bl_index] += dVisY0[bl_index,1].real # /reY0
            P[0*antNum + a1, 0*blNum + bl_index] += dVisX1[bl_index,0].real # /reX1
            P[0*antNum + a1, 2*blNum + bl_index] += dVisX1[bl_index,1].real # /reX1
            P[1*antNum + a1, 1*blNum + bl_index] += dVisY1[bl_index,0].real # /reY1
            P[1*antNum + a1, 3*blNum + bl_index] += dVisY1[bl_index,1].real # /reY1
            P[0*antNum + a0, 4*blNum + bl_index] += dVisX0[bl_index,0].imag # /reX0
            P[0*antNum + a0, 5*blNum + bl_index] += dVisX0[bl_index,1].imag # /reX0
            P[1*antNum + a0, 6*blNum + bl_index] += dVisY0[bl_index,0].imag # /reY0
            P[1*antNum + a0, 7*blNum + bl_index] += dVisY0[bl_index,1].imag # /reY0
            P[0*antNum + a1, 4*blNum + bl_index] += dVisX1[bl_index,0].imag # /reX1
            P[0*antNum + a1, 6*blNum + bl_index] += dVisX1[bl_index,1].imag # /reX1
            P[1*antNum + a1, 5*blNum + bl_index] += dVisY1[bl_index,0].imag # /reY1
            P[1*antNum + a1, 7*blNum + bl_index] += dVisY1[bl_index,1].imag # /reY1
            P[2*antNum + a0, 0*blNum + bl_index] -= dVisX0[bl_index,0].imag # /imX0
            P[2*antNum + a0, 1*blNum + bl_index] -= dVisX0[bl_index,1].imag # /imX0
            P[3*antNum + a0, 2*blNum + bl_index] -= dVisY0[bl_index,0].imag # /imY0
            P[3*antNum + a0, 3*blNum + bl_index] -= dVisY0[bl_index,1].imag # /imY0
            P[2*antNum + a1, 0*blNum + bl_index] += dVisX1[bl_index,0].imag # /imX1
            P[2*antNum + a1, 2*blNum + bl_index] += dVisX1[bl_index,1].imag # /imX1
            P[3*antNum + a1, 1*blNum + bl_index] += dVisY1[bl_index,0].imag # /imY1
            P[3*antNum + a1, 3*blNum + bl_index] += dVisY1[bl_index,1].imag # /imY1
            P[2*antNum + a0, 4*blNum + bl_index] += dVisX0[bl_index,0].real # /imX0
            P[2*antNum + a0, 5*blNum + bl_index] += dVisX0[bl_index,1].real # /imX0
            P[3*antNum + a0, 6*blNum + bl_index] += dVisY0[bl_index,0].real # /imY0
            P[3*antNum + a0, 7*blNum + bl_index] += dVisY0[bl_index,1].real # /imY0
            P[2*antNum + a1, 4*blNum + bl_index] -= dVisX1[bl_index,0].real # /imX1
            P[2*antNum + a1, 6*blNum + bl_index] -= dVisX1[bl_index,1].real # /imX1
            P[3*antNum + a1, 5*blNum + bl_index] -= dVisY1[bl_index,0].real # /imY1
            P[3*antNum + a1, 7*blNum + bl_index] -= dVisY1[bl_index,1].real # /imY1
        #
        #P = P[dSolMap]
        PWtP = np.dot(P, np.dot(np.diag(weight), P.T))
        # correction = np.dot(scipy.linalg.inv(PWtP), np.dot(P, weight* resReal))
        correction = scipy.linalg.solve(PWtP, np.dot(P, weight* resReal))
        #Dx += correction[range(antNum)]           + 1.0j* np.append(0, correction[range(2*antNum, 3*antNum-1)])
        #Dy += correction[range(antNum, 2*antNum)] + 1.0j* correction[range(3*antNum-1, 4*antNum-1)]
        Dx += correction[range(antNum)]           + 1.0j* correction[range(2*antNum, 3*antNum)]
        Dy += correction[range(antNum, 2*antNum)] + 1.0j* correction[range(3*antNum, 4*antNum)]
    #
    return Dx, Dy
#
def Vis2solveDS(Vis, DtX, DtY, PS ):
    def DTresid( D ):
        Dx, Dy = D[0] + (0.0 + 1.0j)*D[1], D[2] + (0.0 + 1.0j)*D[3]
        ModelVis = np.array([
           PS[0] + PS[1]*DtX.conjugate() + PS[2]*Dx + PS[3]*Dx*DtX.conjugate(),
           PS[0]*DtY.conjugate() + PS[1] + PS[2]*Dx*DtY.conjugate() + PS[3]*Dx,
           PS[0]*Dy + PS[1]*Dy*DtX.conjugate() + PS[2] + PS[3]*DtX.conjugate(),
           PS[0]*Dy*DtY.conjugate() + PS[1]*Dy + PS[2]*DtY.conjugate() + PS[3]])
        residVis = (Vis - ModelVis).reshape(4*trkAntNum)
        return np.dot( residVis, residVis.conjugate() ).real
    #
    fit = scipy.optimize.minimize(DTresid, x0=np.zeros(4), method="tnc")['x']
    return fit[0]+(0.0+1.0j)*fit[1], fit[2]+(0.0+1.0j)*fit[3]
#
def Vis2solveD(Vis, DtX, DtY, PS ):
    trkAntNum = len(DtX)  # Number of tracking antennas
    #weight = np.r_[0.1*np.ones(trkAntNum), np.ones(trkAntNum), np.ones(trkAntNum), 0.1*np.ones(trkAntNum), 0.1*np.ones(trkAntNum), np.ones(trkAntNum), np.ones(trkAntNum), 0.1*np.ones(trkAntNum)]
    weight = np.ones(8* trkAntNum)
    Dx, Dy = 0.0 + 0.0j, 0.0 + 0.0j
    Unity = np.ones(trkAntNum)
    Zeros = 0.0* Unity
    for loop_index in range(2):
        #-------- Determine P matrix
        ModelVis = np.array([
            PS[0] + PS[1]*DtX.conjugate() + PS[2]*Dx + PS[3]*Dx*DtX.conjugate(),
            PS[0]*DtY.conjugate() + PS[1] + PS[2]*Dx*DtY.conjugate() + PS[3]*Dx,
            PS[0]*Dy + PS[1]*Dy*DtX.conjugate() + PS[2] + PS[3]*DtX.conjugate(),
            PS[0]*Dy*DtY.conjugate() + PS[1]*Dy + PS[2]*DtY.conjugate() + PS[3]])
        residVis = (Vis - ModelVis).reshape(4*trkAntNum)
        M = dMdDVec(DtX, DtY, Unity).transpose(2,0,1)
        dVis = np.dot( M, PS )
        #
        P = np.c_[
            np.r_[ dVis[:, 0].real,  dVis[:, 1].real, Zeros, Zeros, dVis[:, 0].imag, dVis[:, 1].imag, Zeros, Zeros],
            np.r_[-dVis[:, 0].imag, -dVis[:, 1].imag, Zeros, Zeros, dVis[:, 0].real, dVis[:, 1].real, Zeros, Zeros],
            np.r_[Zeros, Zeros,  dVis[:, 2].real,  dVis[:, 3].real, Zeros, Zeros, dVis[:, 2].imag, dVis[:, 3].imag],
            np.r_[Zeros, Zeros, -dVis[:, 2].imag, -dVis[:, 3].imag, Zeros, Zeros, dVis[:, 2].real, dVis[:, 3].real]]
        #
        tPWP = np.dot(P.T, np.dot(np.diag(weight), P))
        correction = scipy.linalg.solve( tPWP, np.dot(P.T, weight* np.r_[residVis.real, residVis.imag]))
        Dx += (correction[0] + 1.0j*correction[1])
        Dy += (correction[2] + 1.0j*correction[3])
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
def BPtable(msfile, spw, BPScan, blMap, blInv):   # 
    blNum = len(blMap); antNum = Bl2Ant(blNum)[0]
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, BPScan)    # Xspec[pol, ch, bl, time]
    ant0, ant1, polNum, chNum, timeNum = ANT0[0:blNum], ANT1[0:blNum], Pspec.shape[0], Pspec.shape[1], Pspec.shape[3]
    chRange = range(int(0.05*chNum), int(0.95*chNum))                   # Trim band edge
    BP_ant  = np.ones([antNum, 2, chNum], dtype=complex)          # BP_ant[ant, pol, ch]
    kernel_index = KERNEL_BL[0:(antNum-1)]
    if polNum == 4:
        polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
        Xspec  = CrossPolBL(Xspec[:,:,blMap], blInv)
        #---- Delay Cal
        timeAvgSpecX, timeAvgSpecY = np.mean(Xspec[0,chRange][:,kernel_index], axis=2), np.mean(Xspec[3,chRange][:,kernel_index], axis=2)
        antDelayX = np.append(np.array([0.0]), len(chRange)* np.apply_along_axis(delay_search, 0, timeAvgSpecX)[0]/chNum)
        antDelayY = np.append(np.array([0.0]), len(chRange)* np.apply_along_axis(delay_search, 0, timeAvgSpecY)[0]/chNum)
        delayCalTable = np.ones([2, antNum, chNum], dtype=complex)
        for ant_index in range(antNum):
	        delayCalTable[0,ant_index] = np.exp(pi* antDelayX[ant_index]* np.multiply(range(-chNum/2, chNum/2), 1j) / chNum )
	        delayCalTable[1,ant_index] = np.exp(pi* antDelayY[ant_index]* np.multiply(range(-chNum/2, chNum/2), 1j) / chNum )
        #
		delayCaledXspec = (Xspec.transpose(3,0,2,1) * delayCalTable[polYindex][:,ant0] / delayCalTable[polXindex][:,ant1]).transpose(1, 3, 2, 0)
        #---- Gain Cal
        Gain = np.array([gainComplexVec(np.mean(delayCaledXspec[0,chRange], axis=0)), gainComplexVec(np.mean(delayCaledXspec[3,chRange], axis=0))])
        CaledXspec = (Xspec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())).transpose(1,0,2,3)
        #---- Coherent time-averaging
        XPspec = np.mean(CaledXspec, axis=3)  # Time Average
        #---- Antenna-based bandpass table
        BP_ant[:,0], BP_ant[:,1] = gainComplexVec(XPspec[0].T), gainComplexVec(XPspec[3].T)
        #---- XY delay
        BPCaledXspec = XPspec.transpose(2, 0, 1) /(BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate())
        BPCaledXYSpec = np.mean(BPCaledXspec[:,1], axis=0) +  np.mean(BPCaledXspec[:,2], axis=0).conjugate()
        XYdelay, amp = delay_search( BPCaledXYSpec[chRange] )
        XYdelay = (float(chNum) / float(len(chRange)))* XYdelay
        BPCaledXYSpec = BPCaledXYSpec / abs(BPCaledXYSpec)
    else:
        Xspec  = ParaPolBL(Xspec[:,:,blMap], blInv)
        #---- Delay Cal
        timeAvgSpecX, timeAvgSpecY = np.mean(Xspec[0,chRange][:,kernel_index], axis=2), np.mean(Xspec[1,chRange][:,kernel_index], axis=2)
        antDelayX = np.append(np.array([0.0]), len(chRange)* np.apply_along_axis(delay_search, 0, timeAvgSpecX)[0]/chNum)
        antDelayY = np.append(np.array([0.0]), len(chRange)* np.apply_along_axis(delay_search, 0, timeAvgSpecY)[0]/chNum)
        delayCalTable = np.ones([2, antNum, chNum], dtype=complex)
        for ant_index in range(antNum):
	        delayCalTable[0,ant_index] = np.exp(pi* antDelayX[ant_index]* np.multiply(range(-chNum/2, chNum/2), 1j) / chNum )
	        delayCalTable[1,ant_index] = np.exp(pi* antDelayY[ant_index]* np.multiply(range(-chNum/2, chNum/2), 1j) / chNum )
        #
		delayCaledXspec = (Xspec.transpose(3,0,2,1) * delayCalTable[:,ant0] / delayCalTable[:,ant1]).transpose(1, 3, 2, 0)
        #---- Gain Cal
        Gain = np.array([gainComplexVec(np.mean(delayCaledXspec[0,chRange], axis=0)), gainComplexVec(np.mean(delayCaledXspec[1,chRange], axis=0))])
        CaledXspec = (Xspec.transpose(1,0,2,3) / (Gain[:,ant0]* Gain[:,ant1].conjugate())).transpose(1,0,2,3)
        #---- Coherent time-averaging
        XPspec = np.mean(CaledXspec, axis=3)  # Time Average
        #---- Antenna-based bandpass table
        BP_ant[:,0], BP_ant[:,1] = gainComplexVec(XPspec[0].T), gainComplexVec(XPspec[1].T)
        BPCaledXYSpec = np.ones(chNum, dtype=complex)
        XYdelay = 0.0   # No XY correlations
    #
    return BP_ant, BPCaledXYSpec, XYdelay
#
def bpCal(spec, BP0, BP1):      # spec[blNum, chNum, timeNum]
    blnum, chNum, timeNum = len(spec), len(spec[0]), len(spec[0,0])
    ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
    return (spec.transpose(2,0,1) / (BP1[ant0]* BP0[ant1].conjugate())).transpose(1, 2, 0)
#
def phaseCal(spec, Gain):   # spec[blNum, chNum, timeNum], Gain[antNum, chNum, timeNum]
    blnum, chNum, timeNum = spec.shape[0], spec.shape[1], spec.shape[2]
    ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
    return spec * abs(Gain[ant0]* Gain[ant1].conjugate()) / (Gain[ant0]* Gain[ant1].conjugate())
#
def gainCal(spec, Gain):   # spec[blNum, chNum, timeNum], Gain[antNum, chNum, timeNum]
    blnum, chNum, timeNum = spec.shape[0], spec.shape[1], spec.shape[2]
    ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
    return spec / (Gain0[ant0]* Gain1[ant1].conjugate())
#
#-------- SmoothGain
def smoothGain( timeValue, complexValue ):
    amp    = np.abs(complexValue); meanAmp = np.mean(amp)
    meanPhs = np.angle(np.mean(complexValue)); meanTwiddle = np.cos(meanPhs) + (0.0 - 1.0j)*np.sin(meanPhs) 
    biasedPhs =  np.angle(complexValue * meanTwiddle)
    ampSP  = UnivariateSpline( timeValue, amp, w=amp/meanAmp, s=0.1)
    phsSP  = UnivariateSpline( timeValue, biasedPhs + meanPhs, w=amp/meanAmp, s=0.1)
    return ampSP, phsSP
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


#-------- First baseline to get time index
def bandpassCorrection(msfile, antnum, pol, spw, field, chBunch, timeBunch):
	blnum = antnum* (antnum - 1) / 2
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
def PMatrix(CompSol):
    antNum = len(CompSol)
    PM = np.zeros([2*antNum-1, 2*antNum-1])
    SumSqr = CompSol.dot(CompSol.conjugate()).real
    for ant_index in range(antNum):
        PM[ant_index, ant_index] = SumSqr - abs(CompSol[ant_index])**2      # Upper left diagnoal
        for ant_index2 in range(ant_index+1, antNum):                       # Upper left non-diagonal
            PM[ant_index, ant_index2] = CompSol[ant_index].real * CompSol[ant_index2].real - CompSol[ant_index].imag * CompSol[ant_index2].imag
            PM[ant_index2, ant_index] = PM[ant_index, ant_index2]
        #
        #
    #
    for ant_index in range(1, antNum):
        PM[antNum + ant_index - 1, antNum + ant_index - 1] = PM[ant_index, ant_index]   # Lower right diagnal
        for ant_index2 in range(ant_index+1, antNum):                                   # Lower right non-diagonal
            PM[antNum + ant_index - 1, antNum + ant_index2 - 1] = CompSol[ant_index].imag * CompSol[ant_index2].imag - CompSol[ant_index].real * CompSol[ant_index2].real
            PM[antNum + ant_index2 - 1, antNum + ant_index - 1] = PM[antNum + ant_index - 1, antNum + ant_index2 - 1]
        #
        #   Lower left diagnoal = 0
        for ant_index2 in range(ant_index+1, antNum):                       # Lower left non-diagonal
            PM[antNum + ant_index - 1, ant_index2] = CompSol[ant_index].imag* CompSol[ant_index2].real + CompSol[ant_index].real* CompSol[ant_index2].imag
        for ant_index2 in range(0, ant_index):                              # Upper right non-diagonal
            PM[antNum + ant_index - 1, ant_index2] = CompSol[ant_index].imag* CompSol[ant_index2].real + CompSol[ant_index].real* CompSol[ant_index2].imag
        #
    #
    for ant_index in range(antNum):  # Upper right
        PM[ant_index, range(antNum, 2*antNum-1)] = PM[range(antNum, 2*antNum-1), ant_index]
    #
    return PM
#
def PTdotR(CompSol, Cresid):
    antNum = len(CompSol)
    blNum = antNum* (antNum-1) / 2
    ant0, ant1= np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    PTR = np.zeros(2*antNum)
    for ant_index in range(antNum):
        index0 = np.where(ant0 == ant_index)[0].tolist()
        index1 = np.where(ant1 == ant_index)[0].tolist()
        PTR[range(ant_index)]          += (CompSol[ant_index].real* Cresid[index0].real + CompSol[ant_index].imag* Cresid[index0].imag)
        PTR[range(ant_index+1,antNum)] += (CompSol[ant_index].real* Cresid[index1].real - CompSol[ant_index].imag* Cresid[index1].imag)
        PTR[range(antNum, antNum+ant_index)]    += (CompSol[ant_index].imag* Cresid[index0].real - CompSol[ant_index].real* Cresid[index0].imag)
        PTR[range(antNum+ant_index+1,2*antNum)] += (CompSol[ant_index].imag* Cresid[index1].real + CompSol[ant_index].real* Cresid[index1].imag)
    #
    return PTR[range(antNum) + range(antNum+1, 2*antNum)]
#
def gainComplexVec( bl_vis, niter=2 ):       # bl_vis[baseline, channel]
    ChavVis = np.median(bl_vis.real, axis=1) + (0.0+1.0j)*np.median(bl_vis.imag, axis=1)
    blNum, chNum  =  bl_vis.shape[0], bl_vis.shape[1]
    antNum =  Bl2Ant(blNum)[0]
    ant0, ant1, kernelBL = ANT0[0:blNum], ANT1[0:blNum], KERNEL_BL[range(antNum-1)].tolist()
    Solution, CompSol = np.zeros([antNum, chNum], dtype=complex), np.zeros(antNum, dtype=complex)
    #---- Initial solution
    CompSol[0] = sqrt(abs(ChavVis[0])) + 0j
    CompSol[1:antNum] = ChavVis[kernelBL] / CompSol[0]
    #---- Global iteration
    PTP = PMatrix(CompSol)
    L = np.linalg.cholesky(PTP)
    for ch_index in range(chNum):
        Cresid = bl_vis[:,ch_index] - CompSol[ant0]* CompSol[ant1].conjugate()
        t = np.linalg.solve(L, PTdotR(CompSol, Cresid))
        correction = np.linalg.solve(L.T, t)
        Solution[:,ch_index] = CompSol + correction[range(antNum)] + 1.0j* np.append(0, correction[range(antNum, 2*antNum-1)])
    #
    #---- Local iteration
    for iter_index in range(niter):
        for ch_index in range(chNum):
            PTP = PMatrix(Solution[:,ch_index])
            Cresid = bl_vis[:,ch_index] - Solution[ant0, ch_index]* Solution[ant1, ch_index].conjugate()
            L = np.linalg.cholesky(PTP)
            t = np.linalg.solve(L, PTdotR(Solution[:,ch_index], Cresid))
            correction = np.linalg.solve(L.T, t)
            Solution[:,ch_index] = Solution[:,ch_index] + correction[range(antNum)] + 1.0j* np.append(0, correction[range(antNum, 2*antNum-1)])
        #
    #
    return Solution
#
def gainComplex( bl_vis, niter=2 ):
    blNum  =  len(bl_vis)
    antNum =  Bl2Ant(blNum)[0]
    ant0, ant1, kernelBL = ANT0[0:blNum], ANT1[0:blNum], KERNEL_BL[range(antNum-1)].tolist()
    CompSol = np.zeros(antNum, dtype=complex)
    #---- Initial solution
    CompSol[0] = sqrt(abs(bl_vis[0])) + 0j
    CompSol[1:antNum] = bl_vis[kernelBL] / CompSol[0]
    #----  Iteration
    for iter_index in range(niter):
        PTP        = PMatrix(CompSol)
        L          = np.linalg.cholesky(PTP)         # Cholesky decomposition
        Cresid     = bl_vis - CompSol[ant0]* CompSol[ant1].conjugate()
        t          = np.linalg.solve(L, PTdotR(CompSol, Cresid))
        correction = np.linalg.solve(L.T, t)
        CompSol    = CompSol + correction[range(antNum)] + 1.0j* np.append(0, correction[range(antNum, 2*antNum-1)])
    #
    return CompSol
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
    csPA, snPA = np.cos(2.0* PA), np.sin(2.0* PA)
    Xscale = 1.0 / (1.0 + StokesQ* csPA + StokesU* snPA)
    Yscale = 1.0 / (1.0 - StokesQ* csPA - StokesU* snPA)
    #
    ScaledXX, ScaledYY = XX * Xscale, YY* Yscale
    return gainComplexVec(ScaledXX), gainComplexVec(ScaledYY)
    #return np.apply_along_axis( gainComplex, 0, ScaledXX), np.apply_along_axis( gainComplex, 0, ScaledYY)
#
def XXYY2QU(PA, Vis):       # <XX*>, <YY*> to determine Q and U
    timeNum, sinPA2, cosPA2 = len(PA),np.sin(2.0*PA), np.cos(2.0*PA)
    W = np.ones(timeNum) / (np.var(Vis[0].imag) + np.var(Vis[1].imag))   # weight
    XX_YY = Vis[0].real - Vis[1].real
    P = np.array(np.c_[cosPA2, sinPA2]).T
    return scipy.linalg.solve(np.dot(P, np.dot(np.diag(W), P.T)), np.dot(P, W* XX_YY))
#
def XY2Phase(PA, Q, U, Vis):       # XY*, YX* to determine XYphase
    UC_QS = U* np.cos(2.0* PA) - Q* np.sin(2.0* PA)
    correlation = np.dot(Vis[0], UC_QS) + np.dot(Vis[1].conjugate(), UC_QS)
    return np.angle(correlation)
#
def XY2Stokes(PA, Vis, solution):       # XY*, YX* to determine Q, U, XYphase, Dx, Dy
    #-------- Least-Square fit for polarizatino parameters (Q, U, XYphase, Dx, Dy)
    timeNum = len(Vis[0])
    Unity, Zeroty = np.ones([timeNum]), np.zeros([timeNum])
    sinPA2, cosPA2 = np.sin(2.0*PA), np.cos(2.0*PA)
    P = np.zeros([7, 4* timeNum])       # 7 parameters to solve, 4 (ReXY, ImXY, ReYX, ImYX) * timeNum measurements
    W = np.ones(4* timeNum)/ (np.var(Vis[0].imag + Vis[1].imag))
    #-------- Iteration loop
    for index in range(10):
        sinPhi, cosPhi = np.sin(solution[2]), np.cos(solution[2])
        UC_QS = cosPA2* solution[1] - sinPA2* solution[0]   # U cosPA2 - Q sinPA2
        modelVis = np.r_[
             cosPhi* UC_QS + solution[3],       # Re XY*
             sinPhi* UC_QS + solution[4],       # Im XY*
             cosPhi* UC_QS + solution[5],       # Re YX*
            -sinPhi* UC_QS + solution[6] ]      # Im YX*
        vecVis = np.r_[ Vis[0].real, Vis[0].imag, Vis[1].real, Vis[1].imag ]
        residual = vecVis - modelVis
        #-------- Partial matrix
        #                      ReXY             ImXY             ReYX             ImYX
        P[0] = np.r_[-sinPA2* cosPhi, -sinPA2* sinPhi, -sinPA2* cosPhi,  sinPA2* sinPhi]    # dQ
        P[1] = np.r_[ cosPA2* cosPhi,  cosPA2* sinPhi,  cosPA2* cosPhi, -cosPA2* sinPhi]    # du
        P[2] = np.r_[-sinPhi* UC_QS, cosPhi* UC_QS, -sinPhi* UC_QS, -cosPhi* UC_QS]         # dPhi
        P[3] = np.r_[ Unity, Zeroty,  Zeroty,  Zeroty]  # Leakage
        P[4] = np.r_[Zeroty,  Unity,  Zeroty,  Zeroty]  # Leakage
        P[5] = np.r_[Zeroty, Zeroty,   Unity,  Zeroty]  # Leakage
        P[6] = np.r_[Zeroty, Zeroty,  Zeroty,   Unity]  # Leakage
        PTWP_inv = scipy.linalg.inv(np.dot(P, np.dot(np.diag(W), P.T)))
        correction = np.dot( PTWP_inv, np.dot (P, W* residual))
        # print 'Iteration %d : correction = %e' % (index, np.dot(correction,correction))
        solution   = solution + correction
        if np.dot(correction,correction) < 1.0e-15: break
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
def diskVisBeam(diskShape, u, v, primaryBeam):
    # diskShape  : planet disk diameter [MajorAxis, MinorAxis, PA] (rad)
    # u,v        : spatial frequency (= baseline / wavelength) 
    # primaryBeam: FWHM of primary beam [rad]
    cs, sn = np.cos(diskShape[2]), np.sin(diskShape[2])
    diskRadius = 0.5* np.sqrt(diskShape[0]* diskShape[1])
    DSmaj = 1.0 / np.sqrt( (0.30585 / diskShape[0])**2 + 2.0* log(2.0)/(pi* primaryBeam)**2 )    # Primary-beam correction
    DSmin = 1.0 / np.sqrt( (0.30585 / diskShape[1])**2 + 2.0* log(2.0)/(pi* primaryBeam)**2 )    # Primary-beam correction
    uvDisp = (DSmin*(u* cs - v* sn))**2 + (DSmaj*(u* sn + v* cs))**2 
    return beamF(diskRadius/primaryBeam)* np.exp(-0.5* uvDisp)
#
#-------- ArrayCenterAntenna
def bestRefant(uvDist):
    blNum = len(uvDist)
    antNum, ant0, ant1 = Bl2Ant(blNum)[0], ANT0[0:blNum], ANT1[0:blNum]
    blCounter = np.zeros([antNum])
    distOrder = np.argsort(uvDist)
    for bl_index in distOrder:
        blCounter[ant0[bl_index]] += 1
        blCounter[ant1[bl_index]] += 1
        if np.max(blCounter) > 3: break
    #
    return np.argmax(blCounter)
#
#-------- ParallelPol Visibility
def ParaPolBL(Xspec, blInv):
    Neg, Pos = (0.0 + np.array(blInv)), (1.0 - np.array(blInv))
    Tspec = Xspec.copy()
    Tspec[0]   = (Xspec[0].transpose(0,2,1)* Pos + Xspec[0].conjugate().transpose(0,2,1)* Neg).transpose(0,2,1) # XX
    Tspec[1]   = (Xspec[1].transpose(0,2,1)* Pos + Xspec[1].conjugate().transpose(0,2,1)* Neg).transpose(0,2,1) # XY
    return Tspec
#
#-------- CrossPol Visibility
def CrossPolBL(Xspec, blInv):
    Neg, Pos = (0.0 + np.array(blInv)), (1.0 - np.array(blInv))
    Tspec = Xspec.copy()
    Tspec[0]   = (Xspec[0].transpose(0,2,1)* Pos + Xspec[0].conjugate().transpose(0,2,1)* Neg).transpose(0,2,1) # XX
    Tspec[1]   = (Xspec[1].transpose(0,2,1)* Pos + Xspec[2].conjugate().transpose(0,2,1)* Neg).transpose(0,2,1) # XY
    Tspec[2]   = (Xspec[2].transpose(0,2,1)* Pos + Xspec[1].conjugate().transpose(0,2,1)* Neg).transpose(0,2,1) # YX
    Tspec[3]   = (Xspec[3].transpose(0,2,1)* Pos + Xspec[3].conjugate().transpose(0,2,1)* Neg).transpose(0,2,1) # YY
    return Tspec
#
#-------- Tool 
def get_progressbar_str(progress):
    MAX_LEN = 48
    BAR_LEN = int(MAX_LEN * progress)
    return ('[' + '=' * BAR_LEN + ('>' if BAR_LEN < MAX_LEN else '') + ' ' * (MAX_LEN - BAR_LEN) + '] %.1f%%' % (progress * 100.))
#
