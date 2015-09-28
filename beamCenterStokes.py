execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import griddata
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#-------- Read (dAz, dEl) offsets from MS file and return them
def antRefScan( msfile, timeRange ):
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
    return trkAntIndex.tolist(), scanAntIndex.tolist()
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
#
for file_index in range(fileNum):
    msfile = wd + prefix[file_index] + '.ms'
    #-------- Time Range
    interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spw[file_index], scan[file_index])
    #-------- Antenna List
    antList = GetAntName(msfile)
    if refantName not in antList:
        print refantName + ' does not exist in this MS.'
        sys.exit()
    #
    refant_index = np.where( antList == refantName )[0][0]
    #-------- FWHM of the trkant
    antD = 12.0
    if refantName.find('C') > -1: 
        antD = 7.0
    #
    FWHM = GetFWHM(msfile, spw[file_index], antD)
    print('Checking the Array ....')
    #-------- Tracking and Scanning Antennas
    if len(trkAnt) == 0:
        trkAnt, scanAnt = antRefScan(msfile, [min(timeStamp), max(timeStamp)])
    #
    print('-------- Tracking Antennas ----')
    for ant_index in trkAnt:
        text_sd = 'Ref[%d]  / %d: %s ' % (ant_index, len(trkAnt), antList[ant_index])
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
    if refant_index in trkAnt:
        print( refantName + ' is a referencing antenna')
        index = np.where( AntID == trkAnt[0]); scanTime = scanTime[index]; az = az[index]; el = el[index]
        AntIndex = trkAnt
        antNum = len(trkAnt)
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
    interval, timeStamp = GetTimerecord(msfile, 0, trkAnt[0], trkAnt[0], spw[file_index], scan[file_index])
    timeNum = len(timeStamp)
    ScanAz = np.zeros([timeNum]); ScanEl = np.zeros([timeNum]); ScanPA = np.zeros([timeNum])
    #-------- Time index at on axis
    matchNum = np.zeros([timeNum])
    for time_index in range(timeNum):
        matchNum[time_index] = timeMatch( timeStamp[time_index], scanTime, np.median(interval))
    #
    onAxisIndex = np.where( matchNum > 0 )[0].tolist()
    for time_index in onAxisIndex:
        ScanAz[time_index], ScanEl[time_index] = AzElMatch( timeStamp[time_index], scanTime, np.median(interval), az, el)
        ScanPA[time_index] = AzEl2PA(ScanAz[time_index], ScanEl[time_index], ALMA_lat) - BANDPA
    #
    #-- baseline-based weights
    blMap = range(blNum)
    blInv = [False]* blNum      # True -> inverted baseline
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        blMap[bl_index], blInv[bl_index] = Ant2BlD( AntIndex[ants[0]], AntIndex[ants[1]])
    #
    #-------- Cross products (XX, YY, XY, YX)
    chNum, chWid, Freq = GetChNum(msfile, spw[file_index])
    if chNum > 1:
        chRange = range( int(0.05*chNum), int(0.95*chNum))
    #
    print '--- Loading visibilities from MS'
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[file_index], scan[file_index])  # Xspec[POL, CH, BL, TIME]
    if chNum == 1:
        temp = Xspec[:,0]
    else:
        temp = np.mean(Xspec[:,chRange], axis=1)
    #
    if file_index == 0:
        XX = temp[0, blMap]
        XY = temp[1, blMap]
        YX = temp[2, blMap]
        YY = temp[3, blMap]
    else:
        XX = hstack([XX, temp[0, blMap]])
        XY = hstack([XY, temp[1, blMap]])
        YX = hstack([YX, temp[2, blMap]])
        YY = hstack([YY, temp[3, blMap]])
    #
    mjdSec = np.append(mjdSec, timeStamp[onAxisIndex])
    Az     = np.append(Az, ScanAz[onAxisIndex])
    El     = np.append(El, ScanEl[onAxisIndex])
    PA     = np.append(PA, ScanPA[onAxisIndex])
#
solution = np.zeros([7])
for iter_index in range(3):
    print '---- Iteration ' + `iter_index` + ' for Stokes (Q, U) and Gain ----'
    GainX, GainY = polariGain(XX, YY, PA, solution[0], solution[1])
    VisXX = np.mean(gainCalVis( XX, GainX, GainX ), axis = 0)
    VisYY = np.mean(gainCalVis( YY, GainY, GainY ), axis = 0)
    VisXY = np.mean(gainCalVis( XY, GainX, GainY ), axis = 0)
    VisYX = np.mean(gainCalVis( YX, GainY, GainX ), axis = 0)
    solution, solerr = XY2Stokes(PA, VisXY, VisYX)
    text_sd = 'Q/I= %6.3f+-%6.4f  U/I= %6.3f+-%6.4f  XYphase= %6.3f+-%6.4f rad EVPA = %6.2f deg' % (solution[0], solerr[0], solution[1], solerr[1], solution[2], solerr[2], np.arctan(solution[1]/solution[0])*90.0/pi)
    print text_sd
#
plt.plot(PA, VisXY.real, '.', label = 'ReXY', color='cyan')
plt.plot(PA, VisXY.imag, '.', label = 'ImXY', color='darkblue')
plt.plot(PA, VisYX.real, '.', label = 'ReYX', color='magenta')
plt.plot(PA, VisYX.imag, '.', label = 'ImYX', color='darkred')
PArange = np.arange(min(PA), max(PA), 0.01)
plt.plot(PArange,  np.cos(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[3], '-', color='cyan')
plt.plot(PArange,  np.sin(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[4], '-', color='darkblue')
plt.plot(PArange,  np.cos(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[5], '-', color='magenta')
plt.plot(PArange, -np.sin(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[6], '-', color='darkred')
plt.xlabel('PA [rad]'); plt.ylabel('XY, YX (real and imaginary)')
plt.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
text_sd = 'Q/I=%6.3f+-%6.3f U/I=%6.3f+-%6.3f XYphase=%6.3f+-%6.3f rad (RefAnt:%s)' % (solution[0], solerr[0], solution[1], solerr[1], solution[2], solerr[2], antList[0]); plt.text(min(PA), min(VisXY.real), text_sd, size='x-small')
plt.savefig(prefix[0] + '-SPW' + `spw[0]` + '-' + refantName + 'QUXY.pdf', form='pdf')
#-------- Save Results
np.save(prefix[0] + '-SPW' + `spw[0]` + '-' + refantName + '.Ant.npy', antList)
np.save(prefix[0] + '-SPW' + `spw[0]` + '-' + refantName + '.Azel.npy', np.array([mjdSec, Az, El, PA]))
np.save(prefix[0] + '-SPW' + `spw[0]` + '-' + refantName + '.QUXY.npy', solution )
plt.close('all')
