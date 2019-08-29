import re
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
BANDPA = [0.0, 45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0]   # X-pol orientation for Band-1, 2, 3, 4, 5, 6, 7, 8, 9, and 10
BANDFQ = [0.0, 43.2, 75.0, 97.5, 132.0, 183.0, 233.0, 343.5, 460.0, 650.0, 870.0]   # Standard frequency [GHz]
Tcmb = 2.725    # CMB temperature
kb        = 1.38064852e3 # Boltzman constant (* 1e26 for Jy scaling)
RADDEG = 180.0 / math.pi
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
def subArrayIndex(Flag, refant):          #-------- SubArray Indexing
    blNum = len(Flag); antNum = Bl2Ant(blNum)[0]
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    kernelBL = np.where(ant0 == refant)[0].tolist() + np.where(ant1 == refant)[0].tolist()
    flagIndex = np.where(Flag[kernelBL] > 0.5)[0].tolist()
    useKernelBL = np.array(kernelBL)[flagIndex].tolist()
    SAantennas = unique(np.append(ant0[useKernelBL], ant1[useKernelBL]))
    SAantMap = [refant] + sort(np.array(list(set(SAantennas) - set([refant])))).tolist()
    SAantNum = len(SAantennas); SAblNum = SAantNum* (SAantNum - 1)/2
    SAblMap, SAblInv = range(SAblNum), range(SAblNum)
    for bl_index in range(SAblNum): SAblMap[bl_index], SAblInv[bl_index] = Ant2BlD(SAantMap[ant0[bl_index]], SAantMap[ant1[bl_index]])
    return SAantMap, SAblMap, SAblInv
#
def linearRegression( x, y, err):
    weight = 1.0 / err**2
    Sw, Swxx, Swx, Swy, Swxy = np.sum(weight), np.sum( weight* x**2), weight.dot(x), weight.dot(y), np.sum(weight* x* y)
    det    = Sw* Swxx - Swx**2
    sol    = np.array([Swxx* Swy - Swx* Swxy, Sw* Swxy - Swx* Swy]) / det
    solerr = np.sqrt(np.array([Swxx, Sw]) / det)
    return sol, solerr
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
def MullerVector(Dx0, Dy0, Dx1, Dy1, Unity):
    P = np.array([[Unity,               Dx1.conjugate(),      Dx0,                  Dx0* Dx1.conjugate()],
                  [Dy1.conjugate(),     Unity,                Dx0* Dy1.conjugate(), Dx0                 ],
                  [Dy0,                 Dx1.conjugate()* Dy0, Unity,                Dx1.conjugate()     ],
                  [Dy0*Dy1.conjugate(), Dy0,                  Dy1.conjugate(),      Unity               ]]) #.transpose(2,0,1)
    return P
#
def InvMullerVector(Dx0, Dy0, Dx1, Dy1, Unity):
    return np.array([
        [Unity, -Dx1.conjugate(), -Dx0, Dx0* Dx1.conjugate()],
        [-Dy1.conjugate(), Unity, Dx0* Dy1.conjugate(), -Dx0],
        [-Dy0, Dy0* Dx1.conjugate(), Unity, -Dx1.conjugate()],
        [Dy0* Dy1.conjugate(), -Dy0, -Dy1.conjugate(), Unity]]) / ((Unity - Dx0* Dy0)*(Unity - Dx1.conjugate()* Dy1.conjugate()))
#
def PSvector(PA, Stokes):
    """
    PSvector returns a vector of Stokes parameters observced at specified parallactic angle
    """
    cs, sn = np.cos(2.0* PA), np.sin(2.0* PA)
    QCpUS, UCmQS = Stokes[1]* cs + Stokes[2]* sn, Stokes[2]* cs - Stokes[1]* sn
    return np.array([ Stokes[0] + QCpUS, UCmQS + (1.0j* Stokes[3]), UCmQS - (1.0j* Stokes[3]), Stokes[0] - QCpUS])
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
def PAVector(PA, Unity):
    cs, sn = np.cos(2.0* PA), np.sin(2.0* PA)
    Zeroty = 0.0* Unity
    return np.array([
        [Unity, cs, sn, Zeroty],
        [Zeroty, -sn, cs,  1.0j*Unity],
        [Zeroty, -sn, cs, -1.0j*Unity],
        [Unity, -cs, -sn, Zeroty]])
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
def InvPAVector(PA, Unity):
    cs, sn = np.cos(2.0* PA), np.sin(2.0* PA)
    Zeroty = 0.0* Unity
    return 0.5*np.array([
        [Unity, Zeroty, Zeroty, Unity],
        [ cs, -sn, -sn, -cs],
        [ sn,  cs,  cs, -sn],
        [Zeroty,-1.0j*Unity,1.0j*Unity, Zeroty]])
#
def AzEl2PA(az, el, lat=-23.029/180.0*pi): # Azimuth, Elevation, Latitude (default=ALMA) in [rad]
    cos_lat, sin_lat = np.cos(lat), np.sin(lat)
    #return np.arctan( -cos_lat* np.sin(az) / (np.sin(lat)* np.cos(el) - cos_lat* np.sin(el)* np.cos(az)) )
    return np.arctan2( -cos_lat* np.sin(az), (sin_lat* np.cos(el) - cos_lat* np.sin(el)* np.cos(az)) )
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
def AzElMatch( refTime, scanTime, AntID, targetAnt, Az, El ):
    antTimeIndex = np.where(AntID == targetAnt)[0].tolist()
    timeNum = len(refTime)
    az, el = np.zeros(timeNum), np.zeros(timeNum)
    for time_index in range(timeNum):
        time_ptr = np.argmin( abs(scanTime[antTimeIndex] - refTime[time_index]) )
        az[time_index], el[time_index] = np.median(Az[antTimeIndex[time_ptr]]), np.median(El[antTimeIndex[time_ptr]])
    return az, el
#
#-------- Get atmCal SPWs
def GetAtmSPWs(msfile):
    msmd.open(msfile)
    atmSPWs = list( (set(msmd.tdmspws()) | set(msmd.fdmspws())) & set(msmd.spwsforintent("CALIBRATE_ATMOSPHERE*"))); atmSPWs.sort()
    if len(atmSPWs) == 0:
        atmSPWList = msmd.spwsforintent("CALIBRATE_ATMOSPHERE*").tolist()
        tb.open(msfile + '/' + 'SPECTRAL_WINDOW')
        for spwID in atmSPWList:
            if tb.getcell("NUM_CHAN", spwID) > 60: atmSPWs = atmSPWs + [spwID]
        #
        tb.close()
    #
    msmd.close()
    return atmSPWs
#
#-------- Get Bandpass SPWs
def GetBPcalSPWs(msfile):
    msmd.open(msfile)
    bpSPWs  = msmd.spwsforintent("CALIBRATE_PHASE*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_FLUX*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_BANDPASS*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_DELAY*").tolist(); bpSPWs.sort()
    msmd.close()
    return bpSPWs
#
def GetBandNames(msfile, atmSPWs=[]):
    if len(atmSPWs) < 1: atmSPWs = GetAtmSPWs(msfile)
    msmd.open(msfile)
    atmspwNames = msmd.namesforspws(atmSPWs)
    atmBandNames, atmPattern = [], r'RB_..'
    for spwName in atmspwNames : atmBandNames = atmBandNames + re.findall(atmPattern, spwName)
    msmd.close(); msmd.done()
    return atmBandNames
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
    #tb.open( msfile + '/SOURCE')
    tb.open( msfile + '/FIELD')
    SourceID = tb.getcol('SOURCE_ID')
    SourceName = tb.getcol('NAME')
    #SourcePos  = tb.getcol('DIRECTION')
    SourcePos  = tb.getcol('PHASE_DIR')
    tb.close()
    IDNum = np.max(SourceID) + 1
    sourceList, posList = ['']*IDNum, np.zeros([IDNum, 2])
    IDList = np.unique(SourceID).tolist()
    for ID in IDList:
        IDindex = SourceID.tolist().index(ID) 
        sourceList[ID] = SourceName[IDindex]
        #posList[ID] = SourcePos[:,IDindex]
        posList[ID] = SourcePos[:,0,IDindex]
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
def GetPolQuery(sourceName, mjdSec, Freq, SCR_DIR, R_DIR = '' ):
    # Freq (input) : frequency in [GHz]
    os.system('rm -rf CalQU.data')
    text_sd = R_DIR + 'Rscript %spolQuery.R -D%s -F%f %s' % (SCR_DIR, qa.time('%fs' % mjdSec, form='ymd')[0], Freq, sourceName)
    os.system(text_sd)
    fp = open('CalQU.data')
    lines = fp.readlines()
    fp.close()
    catI, catQ, catU = {}, {}, {}
    for eachLine in lines:
        catI[eachLine.split()[0]] = float(eachLine.split()[1])
        catQ[eachLine.split()[0]] = float(eachLine.split()[2])
        catU[eachLine.split()[0]] = float(eachLine.split()[3])
    #
    return catI, catQ, catU
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
def GetTemp(msfile):
    Out = msfile + '/WEATHER'
    tb.open(Out)
    temp = tb.getcol('TEMPERATURE')
    tb.close()
    return np.median(temp)
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

def GetTimerecord(msfile, ant1, ant2, spwID, PScan):
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
def GetVisAllBL(msfile, spwID, scanID, fieldID=-1, AC=True):
    antNum = len(GetAntName(msfile))
    corrNum= antNum* (antNum + 1)/2		# Number of correlations (with autocorr)
    blNum  = corrNum - antNum
    Out = 'DATA_DESC_ID == ' + `spwID` + ' && SCAN_NUMBER == ' + `scanID`
    if fieldID >= 0: Out = Out + ' && FIELD_ID == ' + `fieldID`
    tb.open(msfile)
    antXantYspw = tb.query(Out, sortlist='noduplicates ANTENNA1, ANTENNA2, TIME')
    timeXY = antXantYspw.getcol('TIME')
    if AC:
        timeNum = len(timeXY) / corrNum
    else:
        timeNum = len(timeXY) / blNum
    try:
        dataXY = antXantYspw.getcol('DATA')		# dataXY in array[pol, ch, baselinextime]
    except:
        dataXY = antXantYspw.getcol('FLOAT_DATA')		# dataXY in array[pol, ch, baselinextime]
    tb.close()
    polNum, chNum = dataXY.shape[0], dataXY.shape[1]
    if AC:
        timeStamp = timeXY.reshape(corrNum, timeNum)[0]
        acorr_index = range(antNum)
        xcorr_index = range(blNum)
        for ant_index in range(antNum):
            acorr_index[ant_index] = Ant2Bla_RevLex(ant_index, ant_index, antNum)
        for bl_index in range(blNum):
            ant1, ant0 = Bl2Ant(bl_index)
            xcorr_index[bl_index] = Ant2Bla_RevLex(ant0, ant1, antNum)
        #
        Pspec = dataXY.reshape(polNum, chNum, corrNum, timeNum)[:,:,acorr_index,:]
        Xspec = dataXY.reshape(polNum, chNum, corrNum, timeNum)[:,:,xcorr_index,:]
    else: 
        timeStamp = timeXY.reshape(blNum, timeNum)[0]
        xcorr_index = range(blNum)
        for bl_index in range(blNum):
            ant1, ant0 = Bl2Ant(bl_index)
            xcorr_index[bl_index] = Ant2Bl_RevLex(ant0, ant1, antNum)
        #
        Xspec = dataXY.reshape(polNum, chNum, blNum, timeNum)[:,:,xcorr_index,:]
        Pspec = np.ones([polNum, chNum, antNum, timeNum])
    #
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
def GetPSpecScan(msfile, ant, spwID, scanID):
    Out='ANTENNA1 == '+`ant`+' && ANTENNA2 == '+`ant` + ' && DATA_DESC_ID == '+`spwID` + ' && SCAN_NUMBER == ' + `scanID`
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
def antRefScan( msfile, timeRange, antFlag=[] ):    # Check scanning and tracking antennas
    antList = GetAntName(msfile)
    antNum = len(antList)
    UseAntID = set(range(antNum))
    for flagAnt in antFlag:
        UseAntID = UseAntID - set(np.where(antList == flagAnt)[0])
    #
    UseAntID = list(UseAntID)
    scanRange = np.ones(antNum)
    Time, AntID, Offset = GetAzOffset(msfile)
    for ant_index in UseAntID:
        time_index = np.where( (AntID == ant_index) )[0]
        if len(time_index) == 0:
            scanRange[ant_index] = 0.0
        else:
            time_index = time_index[np.where( (Time[time_index] - timeRange[0]) > -24e-3)[0]]
            time_index = time_index[np.where( (Time[time_index] - timeRange[1]) <  24e-3)[0]]
            if len(time_index) == 0:
                scanRange[ant_index] = 0.0
            else:
                scanRange[ant_index] = max( Offset[0, time_index] ) - min( Offset[0, time_index] )
            #
        #
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
    blNum  =  len(bl_amp); log_bl =  np.log(bl_amp)
    antNum =  Bl2Ant(blNum)[0]
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    PTP_inv = ((2.0* antNum - 2.0)* np.diag(np.ones(antNum)) - 1.0) / (2.0* (antNum - 1.0)* (antNum - 2.0))
    PTV = np.zeros(antNum)
    for ant_index in range(antNum):
        index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
        PTV[ant_index] += np.sum(log_bl[index0]) + np.sum(log_bl[index1])
    #
    return np.exp(PTP_inv.dot(PTV))
#
def ATAmatrix(Gain):  # Gain is a vector of antenna-based gain amplitude (real)
    antNum = len(Gain); normG = Gain.dot(Gain)
    PTP = np.zeros([antNum, antNum]) + Gain
    for ant_index in range(antNum): 
        PTP[ant_index,:] *= Gain[ant_index]
        PTP[ant_index, ant_index] = normG - Gain[ant_index]**2
    #
    return PTP
#
def CTCmatrix(Gain):  # Gain is a vector of antenna-based gain amplitude (imag)
    antNum = len(Gain); normG = Gain.dot(Gain)
    PTP = np.zeros([antNum, antNum]) + Gain
    for ant_index in range(antNum): 
        PTP[ant_index,:] *= (-Gain[ant_index])
        PTP[ant_index, ant_index] = normG - Gain[ant_index]**2
    #
    return PTP
#
def ATBmatrix(Gain):  # Gain is a vector of antenna-based complex gain
    antNum = len(Gain); normG = Gain.real.dot(Gain.imag)
    PTP = np.zeros([antNum, antNum]) + Gain.imag
    for ant_index in range(antNum): 
        PTP[ant_index,:] = Gain.real[ant_index]* Gain.imag + Gain.imag[ant_index]* Gain.real
        PTP[ant_index, ant_index] = 0.0
    #
    return PTP[1:antNum]
#
def clamp_solve(bl_amp, niter=2):
    blNum  =  len(bl_amp)
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
        L = np.linalg.cholesky(ATAmatrix(antGain))
        t = np.linalg.solve(L, y)
        correction = np.linalg.solve(L.T, t)
        antGain += correction; antGain = abs(antGain)
    #
    return antGain
#
def cldelay_solve(bl_delay):    # see http://qiita.com/kamenoseiji/items/782031a0ce8bbc1dc99c
    blNum = len(bl_delay); antNum = Bl2Ant(blNum)[0]
    PTY = np.zeros(antNum)
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    for ant_index in range(antNum):
        index0 = np.where(ant0 == ant_index)[0].tolist()
        index1 = np.where(ant1 == ant_index)[0].tolist()
        PTY[ant_index] += (np.sum(bl_delay[index0]) - np.sum(bl_delay[index1]))
    #
    return PTY / antNum
#
def GFS_delay(Vis, Freq, niter=2):             # Vis[ch, bl]
    blNum = Vis.shape[1]
    antNum = Bl2Ant(blNum)[0]
    antDelay = np.zeros(antNum)
    BBFreq = Freq - np.mean(Freq)
    PM = GFSmatrix(antDelay, BBFreq)
    antDelay = antDelay + np.sum(PM * Vis.T, axis=(1,2))
'''
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
'''
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
    WeightBL = abs(Vis)
    Vis = Vis / WeightBL    # Normalization
    WeightAnt = clamp_solve(WeightBL)
    blNum, antNum = len(Vis), len(WeightAnt)
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    antGain = np.append(1.0 + 0.0j, Vis[KERNEL_BL[0:antNum-1]])
    for iter_index in range(iterNum):
        # resid = WeightBL* (Vis - (antGain[ant0] * antGain[ant1].conjugate()))
        resid = Vis - (antGain[ant0] * antGain[ant1].conjugate())
        # print 'Iter %d : resid = %e : sum(antGain) = %e' % (iter_index, resid.dot(resid.conjugate()).real, sum(antGain).real)
        sumRealAntGain = np.sum(antGain).real
        alpha = 1.0 - sumRealAntGain**2
        PTY = np.ones(antNum, dtype=complex)* sumRealAntGain
        PTP_inv = (np.diag(((3.0 - antNum)* alpha - antNum)*np.ones(antNum)) + alpha) / ((2.0* antNum - 4.0)*alpha**2 + (4.0* antNum - antNum**2)* alpha - antNum**2)
        for ant_index in range(antNum):
            index0 = np.where(ant0 == ant_index)[0].tolist()
            index1 = np.where(ant1 == ant_index)[0].tolist()
            Y = np.zeros(blNum, dtype=complex)
            Y[index0] = -1.0j* antGain[ant_index].conjugate()* antGain[ant1[index0]]
            Y[index1] =  1.0j* antGain[ant_index]* antGain[ant0[index1]].conjugate()
            PTY[ant_index] += Y.dot(resid)
        #
        antPhs  = np.angle(antGain) + PTP_inv.dot(PTY.real)
        antGain = np.cos(antPhs) + 1.0j* np.sin(antPhs)
    #
    #return antGain * np.exp(-1.0j* np.angle(np.sum(antGain)))
    return antGain
#
'''
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
    return antGain
#
'''
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
#-------- D-term determination using full-polarization visibilities to a calibrator whose Stokes parameters are known
def VisPA_solveD(Vis, PA, Stokes, Dx=[], Dy=[]):
    # Vis (input) : full-polarization visibilities [pol, bl, PA]
    # PA  (input) : parallactic angle + BandPA
    # Stokes (input): Stokes parameters of the source
    PAnum, blNum = len(PA), Vis.shape[1]; antNum = Bl2Ant(blNum)[0]; PABLnum = PAnum* blNum
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    CS, SN = np.cos(2.0* PA), np.sin(2.0*PA)
    QCpUS = Stokes[1]*CS + Stokes[2]*SN
    UCmQS = Stokes[2]*CS - Stokes[1]*SN
    ssqQCpUS = QCpUS.dot(QCpUS)     # sum( QCpUS^2 )
    sumQCpUS = np.sum(QCpUS)        # sum( QCpUS )
    if len(Dx) == 0 :                    # Start over the initial D-term value
        #-------- <XX*> to determine Dx (initial value)
        PTP_inv = np.zeros([2*antNum-1, 2*antNum-1])
        PTP_inv[0:antNum][:,0:antNum] = ((2.0* antNum - 2.0)* np.diag(np.ones(antNum)) - 1.0) / (2.0* (antNum - 1.0)* (antNum - 2.0))
        PTP_inv[antNum:(2*antNum-1)][:,antNum:(2*antNum-1)] = (np.diag(np.ones(antNum-1)) + 1.0) / antNum
        PTP_inv /= ssqQCpUS
        PTY = np.zeros([2*antNum-1])
        resid = Vis[0] - (1.0 + QCpUS)
        for ant_index in range(1,antNum):
            index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
            PTY[ant_index] = np.sum(resid[index0].real.dot(QCpUS)) + np.sum(resid[index1].real.dot(QCpUS))
            PTY[antNum + ant_index - 1] = np.sum(resid[index0].imag.dot(QCpUS)) - np.sum(resid[index1].imag.dot(QCpUS))
        #
        index1 = np.where(ant1 == 0)[0].tolist(); PTY[0] = np.sum(resid[index1].real.dot(QCpUS))
        Solution = PTP_inv.dot(PTY); Dx = Solution[0:antNum] + (1.0j)* np.append(0.0, Solution[antNum:2*antNum-1])
        #
        #-------- <YY*> to determine Dy (initial value)
        resid = Vis[3] - (1.0 - QCpUS)
        for ant_index in range(1,antNum):
            index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
            PTY[ant_index] = np.sum(resid[index0].real.dot(QCpUS)) + np.sum(resid[index1].real.dot(QCpUS))
            PTY[antNum + ant_index - 1] = np.sum(resid[index0].imag.dot(QCpUS)) - np.sum(resid[index1].imag.dot(QCpUS))
        #
        index1 = np.where(ant1 == 0)[0].tolist(); PTY[0] = np.sum(resid[index1].real.dot(QCpUS))
        Solution = PTP_inv.dot(PTY); Dy = Solution[0:antNum] + (1.0j)* np.append(0.0, Solution[antNum:2*antNum-1])
    #
    #-------- <XY*> and <YX*> to determine Dx and Dy
    PTP = np.zeros([4*antNum, 4*antNum])            # (Dx.real, Dx.imag, Dy.real, Dy.imag)^2
    PTP[0:2*antNum][:,0:2*antNum] = (ssqQCpUS - 2.0*sumQCpUS + PAnum)* (antNum - 1.0)* np.identity(2*antNum)
    PTP[2*antNum:4*antNum][:,2*antNum:4*antNum] = (ssqQCpUS + 2.0*sumQCpUS + PAnum)* (antNum - 1.0)* np.identity(2*antNum)
    PTP[2*antNum:3*antNum][:,0:antNum] = (PAnum - ssqQCpUS) * (1.0 - np.identity(antNum))
    PTP[0:antNum][:,2*antNum:3*antNum] = PTP[2*antNum:3*antNum][:,0:antNum]
    PTP[3*antNum:4*antNum][:,antNum:2*antNum] = -PTP[2*antNum:3*antNum][:,0:antNum]
    PTP[antNum:2*antNum][:,3*antNum:4*antNum] = -PTP[2*antNum:3*antNum][:,0:antNum]
    #-------- Residual Vector
    residXY, residYX = Vis[1] - UCmQS, Vis[2] - UCmQS
    PTRX, PTRY = np.zeros(antNum, dtype=complex), np.zeros(antNum, dtype=complex)
    for ant_index in range(antNum):
        index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
        residXY[index0] -= Dx[ant_index]* (Stokes[0] - QCpUS)
        residXY[index1] -= Dy[ant_index].conjugate()* (Stokes[0] + QCpUS)
        residYX[index0] -= Dy[ant_index]* (Stokes[0] + QCpUS)
        residYX[index1] -= Dx[ant_index].conjugate()* (Stokes[0] - QCpUS)
    #
    #-------- PTR vector
    for ant_index in range(antNum):
        index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
        PTRX[ant_index] += (Stokes[0] - QCpUS).dot(np.sum(residXY[index0], axis=0))
        PTRX[ant_index] += (Stokes[0] - QCpUS).dot(np.sum(residYX[index1], axis=0).conjugate())
        PTRY[ant_index] += (Stokes[0] + QCpUS).dot(np.sum(residYX[index0], axis=0))
        PTRY[ant_index] += (Stokes[0] + QCpUS).dot(np.sum(residXY[index1], axis=0).conjugate())
    #
    resid = np.array([PTRX.real, PTRX.imag, PTRY.real, PTRY.imag]).reshape(4*antNum)
    #-------- Solution
    L = np.linalg.cholesky(PTP)
    t = np.linalg.solve(L, resid)
    Solution = np.linalg.solve(L.T, t)
    #-------- Correction
    Dx += Solution[0:antNum] + (1.0j)* Solution[antNum:2*antNum]
    Dy += Solution[2*antNum:3*antNum] + (1.0j)* Solution[3*antNum:4*antNum]
    return Dx, Dy
#
#-------- D-term determination using full-polarization visibilities to a calibrator whose Stokes parameters are known
def VisMuiti_solveD(Vis, QCpUS, UCmQS, Dx=[], Dy=[], I=1.0):
    # Vis (input) : full-polarization visibilities [pol, bl, PA]
    # QCpUS  (input) : Q cos(2PA) + U sin(2PA)
    # UCmQS  (input) : U cos(2PA) - Q sin(2PA)
    PAnum, blNum = len(QCpUS), Vis.shape[1]; antNum = Bl2Ant(blNum)[0]; PABLnum = PAnum* blNum
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    ssqUCmQS = UCmQS.dot(UCmQS)     # sum( UCmQS^2 )
    ImmQCpUS = (I - QCpUS).dot(I - QCpUS)  # sum (I - QCpUS)^2)
    IppQCpUS = (I + QCpUS).dot(I + QCpUS)  # sum (I + QCpUS)^2)
    IpmQCpUS = (I + QCpUS).dot(I - QCpUS)  # sum (I + QCpUS)(I - QCpUS)
    if len(Dx) == 0 :                    # Start over the initial D-term value
        #-------- <XX*> to determine Dx (initial value)
        PTP_inv = np.zeros([2*antNum-1, 2*antNum-1])
        PTP_inv[0:antNum][:,0:antNum] = ((2.0* antNum - 2.0)*np.identity(antNum) - 1.0) / (2.0* (antNum - 2.0)* (antNum - 3.0))
        PTP_inv[antNum:(2*antNum-1)][:,antNum:(2*antNum-1)] = (np.identity(antNum - 1) + 1.0) / antNum
        PTP_inv *= (blNum/ ssqUCmQS)
        PTY = np.zeros([2*antNum-1])
        resid = Vis[0] - (I + QCpUS)
        for ant_index in range(1, antNum):
            index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
            PTY[ant_index] = np.sum((resid[index0].real).dot(QCpUS)) + np.sum((resid[index1].real).dot(QCpUS))
            PTY[antNum + ant_index - 1] = np.sum((resid[index0].imag).dot(QCpUS)) - np.sum((resid[index1].imag).dot(QCpUS))
        #
        index1 = np.where(ant1 == 0)[0].tolist(); PTY[0] = np.sum((resid[index1].real).dot(QCpUS))
        Solution = PTP_inv.dot(PTY); newDx = Solution[0:antNum] + (1.0j)* np.append(0.0, Solution[antNum:(2*antNum-1)])
        #-------- <YY*> to determine Dy (initial value)
        resid = Vis[3] - (I - QCpUS)
        for ant_index in range(antNum):
            index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
            PTY[ant_index] = np.sum((resid[index0].real).dot(QCpUS)) + np.sum((resid[index1].real).dot(QCpUS))
            PTY[antNum + ant_index - 1] = np.sum((resid[index0].imag).dot(QCpUS)) - np.sum((resid[index1].imag).dot(QCpUS))
        #
        index1 = np.where(ant1 == 0)[0].tolist(); PTY[0] = np.sum((resid[index1].real).dot(QCpUS))
        Solution = PTP_inv.dot(PTY); newDy = Solution[0:antNum] + (1.0j)* np.append(0.0, Solution[antNum:(2*antNum-1)])
    else:
        newDx, newDy = Dx.copy(), Dy.copy()
    #
    #-------- <XY*> and <YX*> to determine Dx and Dy
    PTP = np.zeros([4*antNum, 4*antNum])            # (Dx.real, Dx.imag, Dy.real, Dy.imag)
    PTP[0:2*antNum][:,0:2*antNum] = ImmQCpUS* (antNum - 1.0)* np.identity(2*antNum)    # P00, P01, P10, and P11
    PTP[2*antNum:4*antNum][:,2*antNum:4*antNum] = IppQCpUS* (antNum - 1.0)* np.identity(2*antNum)  # P22, P23, P32, and P33
    PTP[2*antNum:3*antNum][:,0:antNum] = IpmQCpUS* (1.0 - np.identity(antNum))         # P02
    PTP[3*antNum:4*antNum][:,antNum:2*antNum] = -PTP[2*antNum:3*antNum][:,0:antNum]    # P13
    PTP[0:antNum][:,2*antNum:3*antNum] = PTP[2*antNum:3*antNum][:,0:antNum]            # P20
    PTP[antNum:2*antNum][:,3*antNum:4*antNum] = -PTP[2*antNum:3*antNum][:,0:antNum]    # P31
    PTY = np.zeros(4* antNum)
    #-------- Residual Vector
    residXY, residYX = Vis[1] - UCmQS, Vis[2] - UCmQS
    for ant_index in range(antNum):
        index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
        residXY[index0] -= newDx[ant_index]* (I - QCpUS)
        residXY[index1] -= newDy[ant_index].conjugate()* (I + QCpUS)
        residYX[index0] -= newDy[ant_index]* (I + QCpUS)
        residYX[index1] -= newDx[ant_index].conjugate()* (I - QCpUS)
    #
    #-------- PTR vector
    for ant_index in range(antNum):
        index0, index1 = np.where(ant0 == ant_index)[0].tolist(), np.where(ant1 == ant_index)[0].tolist()
        PTY[ant_index] += (I - QCpUS).dot(np.sum(residXY[index0].real, axis=0))
        PTY[ant_index] += (I - QCpUS).dot(np.sum(residYX[index1].real, axis=0))
        PTY[antNum + ant_index] += (I - QCpUS).dot(np.sum(residXY[index0].imag, axis=0))
        PTY[antNum + ant_index] -= (I - QCpUS).dot(np.sum(residYX[index1].imag, axis=0))
        PTY[2*antNum + ant_index] += (I + QCpUS).dot(np.sum(residXY[index1].real, axis=0))
        PTY[2*antNum + ant_index] += (I + QCpUS).dot(np.sum(residYX[index0].real, axis=0))
        PTY[3*antNum + ant_index] -= (I + QCpUS).dot(np.sum(residXY[index1].imag, axis=0))
        PTY[3*antNum + ant_index] += (I + QCpUS).dot(np.sum(residYX[index0].imag, axis=0))
    #
    #-------- Solution
    L = np.linalg.cholesky(PTP)
    t = np.linalg.solve(L, PTY)
    Solution = np.linalg.solve(L.T, t)
    #-------- Correction
    newDx += Solution[0:antNum] + (1.0j)* Solution[antNum:2*antNum]
    newDy += Solution[2*antNum:3*antNum] + (1.0j)* Solution[3*antNum:4*antNum]
    return newDx, newDy
#
def TransferD(Vis, DtX, DtY, PS):
    refAntNum, PAnum = len(DtX), PS.shape[1]
    #
    A0 = np.repeat(PS[2], refAntNum)  +  np.outer(PS[3], DtX.real).reshape(PAnum* refAntNum)
    A1 = np.repeat(PS[3], refAntNum)  +  np.outer(PS[2], DtY.real).reshape(PAnum* refAntNum)
    A2 = -np.outer(PS[3], DtX.imag).reshape(PAnum* refAntNum)
    A3 = -np.outer(PS[2], DtY.imag).reshape(PAnum* refAntNum)
    #
    B0 = np.repeat(PS[0], refAntNum)  +  np.outer(PS[1], DtX.real).reshape(PAnum* refAntNum)
    B1 = np.repeat(PS[1], refAntNum)  +  np.outer(PS[0], DtY.real).reshape(PAnum* refAntNum)
    B2 = -np.outer(PS[1], DtX.imag).reshape(PAnum* refAntNum)
    B3 = -np.outer(PS[0], DtY.imag).reshape(PAnum* refAntNum)
    #
    resid = Vis.transpose(0,2,1).reshape(4, PAnum* refAntNum)
    resid[0] -= (np.repeat(PS[0], refAntNum) + np.outer(PS[1], DtX.conjugate()).reshape(PAnum* refAntNum))
    resid[1] -= (np.repeat(PS[1], refAntNum) + np.outer(PS[0], DtY.conjugate()).reshape(PAnum* refAntNum))
    resid[2] -= (np.repeat(PS[2], refAntNum) + np.outer(PS[3], DtX.conjugate()).reshape(PAnum* refAntNum))
    resid[3] -= (np.repeat(PS[3], refAntNum) + np.outer(PS[2], DtY.conjugate()).reshape(PAnum* refAntNum))
    #
    PTP_diag = np.array([
        A0.dot(A0) + A1.dot(A1) + A2.dot(A2) + A3.dot(A3), 
        A0.dot(A0) + A1.dot(A1) + A2.dot(A2) + A3.dot(A3), 
        B0.dot(B0) + B1.dot(B1) + B2.dot(B2) + B3.dot(B3),
        B0.dot(B0) + B1.dot(B1) + B2.dot(B2) + B3.dot(B3)])
    #
    PTdotR = np.array([
        A0.dot(resid[0].real) + A1.dot(resid[1].real) + A2.dot(resid[0].imag) + A3.dot(resid[1].imag),
       -A2.dot(resid[0].real) - A3.dot(resid[1].real) + A0.dot(resid[0].imag) + A1.dot(resid[1].imag),
        B0.dot(resid[2].real) + B1.dot(resid[3].real) + B2.dot(resid[2].imag) + B3.dot(resid[3].imag),
       -B2.dot(resid[2].real) - B3.dot(resid[3].real) + B0.dot(resid[2].imag) + B1.dot(resid[3].imag)])
    #
    Solution = PTdotR / PTP_diag
    return Solution[0] + 1.0j* Solution[1], Solution[2] + 1.0j* Solution[3]
#
def Vis2solveD(Vis, DtX, DtY, PS):
    refAntNum = len(DtX)
    A = np.array([ PS[3]* DtX.real + PS[2], PS[2]* DtY.real + PS[3], -PS[3]* DtX.imag, -PS[2]* DtY.imag])
    B = np.array([ PS[1]* DtX.real + PS[0], PS[0]* DtY.real + PS[1], -PS[1]* DtX.imag, -PS[0]* DtY.imag])
    #
    resid = Vis - np.array([ PS[0] + PS[1]* DtX.conjugate(), PS[1] + PS[0]* DtY.conjugate(), PS[2] + PS[3]* DtX.conjugate(), PS[3] + PS[2]* DtY.conjugate()])
    #
    PTP_diag = np.array([np.sum(np.diag(A.dot(A.T))), np.sum(np.diag(B.dot(B.T)))]).repeat(2)
    PTdotR = np.array([
        A.reshape(4*refAntNum).dot( np.r_[resid[0].real, resid[1].real, resid[0].imag, resid[1].imag]),
        A.reshape(4*refAntNum).dot( np.r_[resid[0].imag, resid[1].imag,-resid[0].real,-resid[1].real]),
        B.reshape(4*refAntNum).dot( np.r_[resid[2].real, resid[3].real, resid[2].imag, resid[3].imag]),
        B.reshape(4*refAntNum).dot( np.r_[resid[2].imag, resid[3].imag,-resid[2].real,-resid[2].real])])
    #
    Solution = PTdotR / PTP_diag
    return Solution[0] + 1.0j* Solution[1], Solution[2] + 1.0j* Solution[3]
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
	delay = argmax(corrAbs) - nspec # Coarse Delay
	trial_delay = np.multiply(range(-2,3), 0.5)
	trial_amp   = np.zeros(5) 
	for i in range(5):
		trial_amp[i] = abs(np.mean(delay_cal(spec, delay + trial_delay[i])))
	fit = np.polyfit(trial_delay, trial_amp, 2)
	return delay - fit[1]/(2.0*fit[0]), fit[2] - fit[1]**2/(4*fit[0])	# return delay and amp
#
def blGain( blSpec ):				# Average in spectral channels
	return np.mean(blSpec, 0)

def blBp( blSpec ):					# Average in time
	return np.mean(blSpec, 1)

def bunchVec( spec, bunchNum ):
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
    ant0, ant1 = np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    delayResults = np.apply_along_axis( delay_search, 1, Xspec[:, chRange] )    # BL delays in delayResults[:,0], Amps in delayResults[:,1]
    delay_ant = cldelay_solve(delayResults[:,0])* ((chNum + 0.0)/len(chRange)); antNum = len(delay_ant) # Antenna-based delay solutions
    twiddleAnt = np.exp( pi* np.outer( delay_ant, np.multiply(range(-chNum/2, chNum/2), 1j) / chNum ) ) # Antenna-baased Twiddle factor
    return Xspec * twiddleAnt[ant0] / twiddleAnt[ant1]
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
def BPtable(msfile, spw, BPScan, blMap, blInv, bunchNum=1, FG=np.array([]), TS=np.array([])): 
    blNum = len(blMap); antNum = Bl2Ant(blNum)[0]
    #timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, BPScan, -1, False)    # Xspec[pol, ch, bl, time]
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, BPScan, -1, True)    # Xspec[pol, ch, bl, time]
    if bunchNum > 1:
        def bunchN(Spec): return bunchVec(Spec, bunchNum)
        Xspec = np.apply_along_axis(bunchN, 1, Xspec)
    #
    if len(FG) > 0: flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
    else: flagIndex = range(len(timeStamp))
    #
    ant0, ant1, polNum, chNum, timeNum = ANT0[0:blNum], ANT1[0:blNum], Pspec.shape[0], Xspec.shape[1], Xspec.shape[3]
    chRange = range(int(0.05*chNum), int(0.95*chNum))                   # Trim band edge
    kernel_index = KERNEL_BL[0:(antNum-1)]
    if polNum == 4:
        BP_ant  = np.ones([antNum, 2, chNum], dtype=complex)          # BP_ant[ant, pol, ch]
        polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
        #print '  -- Mapping baselines...'
        Xspec  = CrossPolBL(Xspec[:,:,blMap], blInv)[:,:,:,flagIndex]
        #---- Delay Cal
        #print '  -- Delay Cal...'
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
        del delayCaledXspec
        CaledXspec = (Xspec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())).transpose(1,0,2,3)
        del Xspec
        #---- Coherent time-averaging
        #print '  -- Time Averaging...'
        XPspec = np.mean(CaledXspec, axis=3)  # Time Average
        del CaledXspec
        #---- Antenna-based bandpass table
        #print '  -- Antenna-based Solution...'
        BP_ant[:,0], BP_ant[:,1] = gainComplexVec(XPspec[0].T), gainComplexVec(XPspec[3].T)
        #---- XY delay
        #print '  -- XY delay...'
        BPCaledXspec = XPspec.transpose(2, 0, 1) /(BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate())
        del XPspec
        BPCaledXYSpec = np.mean(BPCaledXspec[:,1], axis=0) +  np.mean(BPCaledXspec[:,2], axis=0).conjugate()
        XYdelay, amp = delay_search( BPCaledXYSpec[chRange] )
        XYdelay = (float(chNum) / float(len(chRange)))* XYdelay
        BPCaledXYSpec = BPCaledXYSpec / abs(BPCaledXYSpec)
    elif polNum == 2:
        BP_ant  = np.ones([antNum, 2, chNum], dtype=complex)          # BP_ant[ant, pol, ch]
        Xspec  = ParaPolBL(Xspec[:,:,blMap], blInv)[:,:,:,flagIndex]
        #---- Delay Cal
        timeAvgSpecX, timeAvgSpecY = np.mean(Xspec[0,chRange][:,kernel_index], axis=2), np.mean(Xspec[1,chRange][:,kernel_index], axis=2)
        antDelayX = np.append(np.array([0.0]), len(chRange)* np.apply_along_axis(delay_search, 0, timeAvgSpecX)[0]/chNum)
        antDelayY = np.append(np.array([0.0]), len(chRange)* np.apply_along_axis(delay_search, 0, timeAvgSpecY)[0]/chNum)
        delayCalTable = np.ones([2, antNum, chNum], dtype=complex)
        for ant_index in range(antNum):
            print `ant_index` + ' ',
            delayCalTable[0,ant_index] = np.exp(pi* antDelayX[ant_index]* np.multiply(range(-chNum/2, chNum/2), 1j) / chNum )
            delayCalTable[1,ant_index] = np.exp(pi* antDelayY[ant_index]* np.multiply(range(-chNum/2, chNum/2), 1j) / chNum )
        #
        delayCaledXspec = (Xspec.transpose(3,0,2,1) * delayCalTable[:,ant0] / delayCalTable[:,ant1]).transpose(1, 3, 2, 0)
        #---- Gain Cal
        Gain = np.array([gainComplexVec(np.mean(delayCaledXspec[0,chRange], axis=0)), gainComplexVec(np.mean(delayCaledXspec[1,chRange], axis=0))])
        del delayCaledXspec
        CaledXspec = (Xspec.transpose(1,0,2,3) / (Gain[:,ant0]* Gain[:,ant1].conjugate())).transpose(1,0,2,3)
        del Xspec
        #---- Coherent time-averaging
        XPspec = np.mean(CaledXspec, axis=3)  # Time Average
        del CaledXspec
        #---- Antenna-based bandpass table
        BP_ant[:,0], BP_ant[:,1] = gainComplexVec(XPspec[0].T), gainComplexVec(XPspec[1].T)
        del XPspec
        BPCaledXYSpec = np.ones(chNum, dtype=complex)
        XYdelay = 0.0   # No XY correlations
    #
    else :
        BP_ant  = np.ones([antNum, 1, chNum], dtype=complex)          # BP_ant[ant, pol, ch]
        Xspec  = SinglePolBL(Xspec[:,:,blMap], blInv)[:,:,:,flagIndex]
        #---- Delay Cal
        timeAvgSpecX = np.mean(Xspec[0,chRange][:,kernel_index], axis=2)
        antDelayX = np.append(np.array([0.0]), len(chRange)* np.apply_along_axis(delay_search, 0, timeAvgSpecX)[0]/chNum)
        delayCalTable = np.ones([1, antNum, chNum], dtype=complex)
        for ant_index in range(antNum):
	        delayCalTable[0,ant_index] = np.exp(pi* antDelayX[ant_index]* np.multiply(range(-chNum/2, chNum/2), 1j) / chNum )
        #
		delayCaledXspec = (Xspec.transpose(3,0,2,1) * delayCalTable[:,ant0] / delayCalTable[:,ant1]).transpose(1, 3, 2, 0)
        #---- Gain Cal
        Gain = np.array([gainComplexVec(np.mean(delayCaledXspec[0,chRange], axis=0))])
        CaledXspec = (Xspec.transpose(1,0,2,3) / (Gain[0,ant0]* Gain[0,ant1].conjugate())).transpose(1,0,2,3)
        #---- Coherent time-averaging
        XPspec = np.mean(CaledXspec, axis=3)  # Time Average
        #---- Antenna-based bandpass table
        BP_ant[:,0] = gainComplexVec(XPspec[0].T)
        BPCaledXYSpec = np.ones(chNum, dtype=complex)
        XYdelay = 0.0   # No XY correlations
    #
    return BP_ant, BPCaledXYSpec, XYdelay, Gain
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
	timeXY = GetTimerecord(msfile, 0, 1, spw, field)
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
def splineComplex( samplePoints, vector ):
    SPL_real = UnivariateSpline(samplePoints, vector.real)
    SPL_imag = UnivariateSpline(samplePoints, vector.imag)
    return( SPL_real(samplePoints) + (0.0 + 1.0j)* SPL_imag(samplePoints) )
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
#======== Amplitude calibrations
#-------- Residuals for Tsky - secz regresion (optical depth + intercept)
def residTskyTransfer( param, Tamb, secz, Tsky, weight ):
    exp_Tau = np.exp( -param[1]* secz )
    return weight* (Tsky - (param[0] + Tcmb* exp_Tau  + Tamb* (1.0 - exp_Tau)))
#
#-------- Residuals for Tsky - secz regresion (without intercept)
def residTskyTransfer0( param, Tamb, secz, Tsky, weight ):
    exp_Tau = np.exp( -param[0]* secz )
    return weight* (Tsky - (Tcmb* exp_Tau  + Tamb* (1.0 - exp_Tau)))
#
#-------- Residuals for Tsky - secz regresion (fixed zenith optical depth)
def residTskyTransfer2( param, Tamb, Tau0, secz, Tsky, weight ):
    exp_Tau = np.exp( -Tau0* secz )
    return weight* (Tsky - (param[0] + Tcmb* exp_Tau  + Tamb* (1.0 - exp_Tau)))
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
def GFSmatrix(antDelay, Frequency):
    chNum, antNum = len(Frequency), len(antDelay)
    blNum = antNum* (antNum - 1) / 2
    ant0, ant1= np.array(ANT0[0:blNum]), np.array(ANT1[0:blNum])
    PM = np.zeros([antNum, blNum, chNum], dtype=complex)
    for ant_index in range(antNum):
        index0 = np.where(ant0 == ant_index)[0].tolist()
        index1 = np.where(ant1 == ant_index)[0].tolist()
        PM[ant_index, index1] =-np.exp(2.0j* pi* np.outer((antDelay[ant0[index1]] - antDelay[ant_index]), Frequency))
        PM[ant_index, index0] = np.exp(2.0j* pi* np.outer((antDelay[ant_index] - antDelay[ant1[index0]]), Frequency))
    #
    return -4.0j* pi* Frequency* PM
#
def PMatrix(CompSol):
    antNum = len(CompSol); matSize = 2*antNum-1
    PM = np.zeros([matSize, matSize])
    PM[0:antNum][:,0:antNum]   = ATAmatrix(CompSol.real) + CTCmatrix(CompSol.imag)
    PM[antNum:matSize][:,antNum:matSize] = ATAmatrix(CompSol.imag)[1:antNum][:,1:antNum] + CTCmatrix(CompSol.real)[1:antNum][:,1:antNum]
    PM[antNum:matSize][:,0:antNum]   = ATBmatrix(CompSol)
    PM[0:antNum][:,antNum:matSize]   = PM[antNum:matSize][:,0:antNum].T
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
    #---- Initial solution
    CompSol = np.append(sqrt(abs(ChavVis[0])) + 0j, ChavVis[kernelBL])
    CompSol[1:antNum] /=  CompSol[0]
    #---- Global iteration
    PTP = PMatrix(CompSol)
    L = np.linalg.cholesky(PTP)
    def GlobalGainSolve(vis):
        Cresid = vis - CompSol[ant0]* CompSol[ant1].conjugate()
        t = np.linalg.solve(L, PTdotR(CompSol, Cresid))
        correction = np.linalg.solve(L.T, t)
        return CompSol + correction[range(antNum)] + 1.0j* np.append(0, correction[range(antNum, 2*antNum-1)])
    #
    Solution = np.apply_along_axis(GlobalGainSolve, 0, bl_vis)
    #---- Local iteration
    def LocalGainSolve(visSol):
        vis, sol = visSol[0:blNum], visSol[blNum:blNum+antNum]
        PTP = PMatrix(sol)
        L = np.linalg.cholesky(PTP)
        Cresid = vis - sol[ant0]* sol[ant1].conjugate()
        t = np.linalg.solve(L, PTdotR(sol, Cresid))
        correction = np.linalg.solve(L.T, t)
        return sol + correction[range(antNum)] + 1.0j* np.append(0, correction[range(antNum, 2*antNum-1)])
    #
    for iter_index in range(niter): Solution = np.apply_along_axis(LocalGainSolve, 0, np.concatenate([bl_vis, Solution]))
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
def gainComplexErr( bl_vis, niter=2 ):
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
    return CompSol, logamp_solve(abs(bl_vis - CompSol[ant0]* CompSol[ant1].conjugate()))
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
def QUscale(PA,  StokesQ, StokesU):
    csPA, snPA = np.cos(2.0* PA), np.sin(2.0* PA)
    Xscale = 1.0 / (1.0 + StokesQ* csPA + StokesU* snPA)
    Yscale = 1.0 / (1.0 - StokesQ* csPA - StokesU* snPA)
#
def polariGain( XX, YY, QCpUS):
    blNum, timeNum = XX.shape[0], XX.shape[1]
    Xscale = 1.0 / (1.0 + QCpUS)
    Yscale = 1.0 / (1.0 - QCpUS)
    #
    ScaledXX, ScaledYY = XX * Xscale, YY* Yscale
    return gainComplexVec(ScaledXX), gainComplexVec(ScaledYY)
def XXYY2QU(PA, XX_YY):       # <XX*>, <YY*> to determine Q and U
    PAnum, SN, CS = len(PA),np.sin(2.0*PA), np.cos(2.0*PA)
    W = np.ones(PAnum) / (np.var(XX_YY[0].imag) + np.var(XX_YY[1].imag))   # weight
    XXmYY = (XX_YY[0].real - XX_YY[1].real) / (XX_YY[0].real + XX_YY[1].real)
    P = np.array(np.c_[np.ones(PAnum), CS, SN]).T
    return scipy.linalg.solve(np.dot(P, np.dot(np.diag(W), P.T)), np.dot(P, W* XXmYY))[[1,2]]
#
def XY2Phase(UC_QS, Vis):       # XY*, YX* to determine XYphase
    correlation = np.dot(Vis[0], UC_QS) + np.dot(Vis[1].conjugate(), UC_QS)
    return np.angle(correlation)
#
def XY2PhaseVec(TS, UC_QS, Vis):    # XY*, YX* to measuere XYphase variation
    product = (Vis[0] + Vis[1].conjugate())* UC_QS
    SP_phas = UnivariateSpline(TS, np.arctan(product.imag/product.real), w=abs(product)**2, s=0.01* np.median(abs(product)**2))
    return SP_phas(TS), product
#
def XY2Stokes(PA, Vis):            # XY*, YX* to determine Q, U
    UCmQS = 0.5* (Vis[0] + Vis[1]).real
    sinPA2, cosPA2 = np.sin(2.0*PA), np.cos(2.0*PA)
    SS, CC, SC, SY, CY = sinPA2.dot(sinPA2), cosPA2.dot(cosPA2), sinPA2.dot(cosPA2), sinPA2.dot(UCmQS), cosPA2.dot(UCmQS)
    PTP_inv = np.array([[CC,SC],[SC,SS]]) / (SS*CC - SC**2)
    QU = PTP_inv.dot( np.array([-SY, CY]) )
    return QU[0], QU[1]
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
def bestRefant(uvDist, useantList=[]):
    blNum = len(uvDist)
    antNum, ant0, ant1 = Bl2Ant(blNum)[0], ANT0[0:blNum], ANT1[0:blNum]
    if len(useantList) == 0: useantList = range(antNum)
    blCounter = np.zeros([antNum]); blCounter[useantList] += blNum
    distOrder = np.argsort(uvDist)
    for bl_index in distOrder:
        blCounter[ant0[bl_index]] += (blNum - bl_index)
        blCounter[ant1[bl_index]] += (blNum - bl_index)
        if np.max(blCounter) > 4* blNum: break
    #
    return np.argmax(blCounter)
#
#-------- SinglePol Visibility
def SinglePolBL(Xspec, blInv):
    Neg, Pos = (0.0 + np.array(blInv)), (1.0 - np.array(blInv))
    Tspec = Xspec.copy()
    Tspec[0]   = (Xspec[0].transpose(0,2,1)* Pos + Xspec[0].conjugate().transpose(0,2,1)* Neg).transpose(0,2,1) # XX
    return Tspec
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
