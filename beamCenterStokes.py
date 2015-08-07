execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import griddata
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#-------- Read (dAz, dEl) offsets from MS file and return them
def antRefScan( msfile ):
    antList = GetAntName(msfile)
    antNum = len(antList)
    scanRange = np.zeros(antNum)
    Time, AntID, Offset = GetAzOffset(msfile)
    for ant_index in range(antNum):
        time_index = np.where( AntID == ant_index )[0]
        scanRange[ant_index] = max( Offset[0, time_index] ) - min( Offset[0, time_index] )
    #
    refAntIndex  = np.where( scanRange == 0.0 )[0]
    scanAntIndex = np.where( scanRange >  0.0 )[0]
    return refAntIndex.tolist(), scanAntIndex.tolist()
#
#-------- Scanning Offset < threshold
def scanThresh( msfile, ant_index, thresh ):
    Time, AntID, Offset = GetAzOffset(msfile)
    time_index = np.where( AntID == ant_index )[0]
    onAxisIndex = np.where( Offset[0, time_index]**2 + Offset[1, time_index]**2 < thresh**2 )[0]
    return time_index[onAxisIndex].tolist()
    #return onAxisIndex.tolist()
#
#-------- Time-based matching between time tags in visibilities and in scan pattern 
def AzElMatch( refTime, scanTime, thresh, Az, El ):
    index = np.where( abs(scanTime - refTime) < thresh)[0]
    return np.median(Az[index]), np.median(El[index])
#
def timeMatch( refTime, scanTime, thresh):
    match = np.where( abs(scanTime - refTime) < thresh)[0].tolist()
    return len(match)
#
#----------------------------------------- Procedures
fileNum = len(prefix)
mjdSec = []
Az     = []
El     = []
PA     = []
visXX  = []
visXY  = []
visYX  = []
visYY  = []
#
for file_index in range(fileNum):
    msfile = wd + prefix[file_index] + '.ms'
    #-------- Antenna List
    antList = GetAntName(msfile)
    if refantName not in antList:
        print refantName + ' does not exist in this MS.'
        sys.exit()
    #
    refant_index = np.where( antList == refantName )[0][0]
    #-------- FWHM of the refant
    antD = 12.0
    if refantName.find('C') > -1: 
        antD = 7.0
    #
    FWHM = GetFWHM(msfile, spw[0], antD)
    print('Checking the Array ....')
    #-------- Reference and Scan Antennas
    refAnt, scanAnt = antRefScan(msfile)
    print('-------- Reference Antennas ----')
    for ant_index in refAnt:
        text_sd = 'Ref[%d]  / %d: %s ' % (ant_index, len(refAnt), antList[ant_index])
        print text_sd
    #
    if len(scanAnt) > 0 :
        print('-------- Scanning Antennas ----')
        for ant_index in scanAnt:
            text_sd = 'Scan[%d] / %d: %s ' % (ant_index, len(scanAnt), antList[ant_index])
            print text_sd
        #
    #
    scanTime, AntID, az, el = GetAzEl(msfile)
    if refant_index in refAnt:
        print( refantName + ' is a referencing antenna')
        index = np.where( AntID == refAnt[0]); scanTime = scanTime[index]; az = az[index]; el = el[index]
        AntIndex = refAnt
        antNum = len(refAnt)
    if refant_index in scanAnt:
        print( refantName + ' is a scanning antenna')
        index = scanThresh( msfile, 0, FWHM/20); scanTime = scanTime[index]; az = az[index]; el = el[index]
        antNum = len(antList)
        AntIndex = range(antNum)
    #
    blNum  = antNum* (antNum - 1) / 2
    antList = antList[AntIndex]
    antWeight = np.ones(antNum)
    #-------- Visibility sampling points
    interval, timeStamp = GetTimerecord(msfile, 0, refAnt[0], refAnt[0], spw[0], scan[0])
    timeNum = len(timeStamp)
    ScanAz = np.zeros([timeNum]); ScanEl = np.zeros([timeNum]); ScanPA = np.zeros([timeNum])
    #-------- Time index at on axis
    matchNum = np.zeros([timeNum])
    for time_index in range(timeNum):
        matchNum[time_index] = timeMatch( timeStamp[time_index], scanTime, np.median(interval))
    #
    onAxisIndex = np.where( matchNum > 0 )[0].tolist()
    #for time_index in range(timeNum):
    for time_index in onAxisIndex:
        ScanAz[time_index], ScanEl[time_index] = AzElMatch( timeStamp[time_index], scanTime, np.median(interval), az, el)
        ScanPA[time_index] = AzEl2PA(ScanAz[time_index], ScanEl[time_index], ALMA_lat)
    #
    #-- baseline-based weights
    blMap = range(blNum)
    blInv = [False]* blNum      # True -> inverted baseline
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        blMap[bl_index], blInv[bl_index] = Ant2BlD( AntIndex[ants[0]], AntIndex[ants[1]])
    #
    #-------- StokesParams
    for spw_index in range(len(spw)):
        chNum, chWid, Freq = GetChNum(msfile, spw[spw_index])
        vis_sigma = np.ones(blNum) / sqrt(2.0e9* np.median(diff(timeStamp)))
        print '--- Loading visibilities from MS'
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], scan[0])  # Xspec[POL, CH, BL, TIME]
        VisXX, VisXY, VisYX, VisYY = polariVis( Xspec[:,:,blMap].swapaxes(1,2) )# VisXX[TIME] --- avaraged in CH and BL
        #
        #for time_index in range(timeNum):
        for time_index in onAxisIndex:
            print '%d %f %f %f %f %f %f %f %f %f %f %f' % (timeStamp[time_index], ScanAz[time_index], ScanEl[time_index], ScanPA[time_index], VisXX[time_index].real, VisXX[time_index].imag, VisXY[time_index].real, VisXY[time_index].imag, VisYX[time_index].real, VisYX[time_index].imag, VisYY[time_index].real, VisYY[time_index].imag)
        #
    #
    mjdSec = np.append(mjdSec, timeStamp[onAxisIndex])
    Az     = np.append(Az, ScanAz[onAxisIndex])
    El     = np.append(El, ScanEl[onAxisIndex])
    PA     = np.append(PA, ScanPA[onAxisIndex])
    visXX  = np.append(visXX, VisXX[onAxisIndex])
    visXY  = np.append(visXY, VisXY[onAxisIndex])
    visYX  = np.append(visYX, VisYX[onAxisIndex])
    visYY  = np.append(visYY, VisYY[onAxisIndex])
#
np.save(prefix[0] + '.Ant.npy', antList)
np.save(prefix[0] + '.Azel.npy', np.array([mjdSec, Az, El, PA]))
np.save(prefix[0] + '.Vis.npy', np.array([visXX, visXY, visYX, visYY]))
