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
#-------- Time-based matching between time tags in visibilities and in scan pattern 
def AzElMatch( refTime, scanTime, thresh, Az, El ):
    index = np.where( abs(scanTime - refTime) < thresh)[0]
    return np.median(Az[index]), np.median(El[index])
#
#----------------------------------------- Procedures
fileNum = len(prefix)
for file_index in range(fileNum):
    msfile = wd + prefix[file_index] + '.ms'
    #-------- Antenna List
    antList = GetAntName(msfile)
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
    scanTime, AntID, Az, El = GetAzEl(msfile)
    index = np.where( AntID == refAnt[0]); scanTime = scanTime[index]; Az = Az[index]; El = El[index]
    #
    AntIndex = refAnt
    antNum = len(refAnt)
    blNum  = antNum* (antNum - 1) / 2
    antList = antList[AntIndex]
    antWeight = np.ones(antNum)
    #-------- Visibility sampling points
    interval, timeStamp = GetTimerecord(msfile, 0, refAnt[0], refAnt[0], spw[0], scan[0])
    timeNum = len(timeStamp)
    ScanAz = np.zeros([timeNum]); ScanEl = np.zeros([timeNum]); ScanPA = np.zeros([timeNum])
    for time_index in range(timeNum):
        ScanAz[time_index], ScanEl[time_index] = AzElMatch( timeStamp[time_index], scanTime, np.median(interval), Az, El)
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
        VisXX, VisXY, VisYX, VisYY = polariVis( Xspec[:,:,blMap].swapaxes(1,2) )
        #
        for time_index in range(timeNum):
            print '%d %f %f %f %f %f %f %f %f %f %f %f' % (timeStamp[time_index], ScanAz[time_index], ScanEl[time_index], ScanPA[time_index], VisXX[time_index].real, VisXX[time_index].imag, VisXY[time_index].real, VisXY[time_index].imag, VisYX[time_index].real, VisYX[time_index].imag, VisYY[time_index].real, VisYY[time_index].imag)
        #
        plt.plot( ScanPA, abs( VisXX ), 'b.',)
        plt.plot( ScanPA, abs( VisXY ), 'r.',)
        plt.plot( ScanPA, abs( VisYX ), 'ro',)
        plt.plot( ScanPA, abs( VisYY ), 'bo',)
    #
#
