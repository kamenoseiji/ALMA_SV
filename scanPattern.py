execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
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
    time_index = np.where( AntID == scanAntIndex[0] )[0].tolist()
    return refAntIndex.tolist(), scanAntIndex.tolist(), Time[time_index], Offset[:, time_index]
#
#-------- Time-based matching between time tags in visibilities and in scan pattern 
def AzElMatch( refTime, scanTime, thresh, Offset ):
    index = np.where( abs(scanTime - refTime) < thresh)[0]
    return np.median(Offset[0, index]), np.median(Offset[1, index])
#
#-------- Divide scan into subscans, based on time continuity
def subScan(scanTime, gapThresh):
    gapIndex = np.where( diff(scanTime) > gapThresh )[0]
    subScanNum = len(gapIndex) + 1
    subScanStartIndex = np.append( array([0]),  gapIndex + 1).tolist()
    subScanEndIndex   = np.append( gapIndex, array([len(scanTime)-1])).tolist()
    return subScanStartIndex, subScanEndIndex
#
#-------- For drawing a circle
def circlePoints( x, y, radius ):
    angle = np.arange(-pi, (130/128)*pi, pi/128)
    return x + radius* np.cos(angle), y + radius* np.sin(angle)
#----------------------------------------- Procedures
msfile = wd + prefix + '.ms'
#-------- Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
#-------- Reference and Scan Antennas
refAnt, scanAnt, scanTime, Offset = antRefScan(msfile)
#-------- Visibility sampling points
interval, timeStamp = GetTimerecord(msfile, 0, refAnt[0], scanAnt[0], spw[0], scan[0])
timeNum = len(timeStamp)
ScanAz = np.zeros([timeNum]); ScanEl = np.zeros([timeNum])
for time_index in range(timeNum):
    ScanAz[time_index], ScanEl[time_index] = AzElMatch( timeStamp[time_index], scanTime, np.median(interval), Offset)
#
#-------- Frequency and Wavelength
chNum, chWid, Freq = GetChNum(msfile, spw[0])
wavelength = constants.c / Freq 
FWHM = 1.13* 180.0* 3600.0* wavelength / (12.0* pi) # Gaussian beam, in unit of arcsec
#-------- SubScan
subScanStartIndex, subScanEndIndex = subScan(timeStamp, 1.5*np.median(interval))
subScanNum = len(subScanStartIndex)
#
#-------- Prepare Figure
fig = plt.figure( figsize = (8,8))
for subScanIndex in range(subScanNum):
    time_index = range( subScanStartIndex[subScanIndex], subScanEndIndex[subScanIndex]+1 )
    plt.plot( ScanAz[time_index], ScanEl[time_index], '.', color="cyan" )
    plt.text( ScanAz[time_index[0]], ScanEl[time_index[0]], `subScanIndex`)
#
plt.title( prefix + ' ' + antList[scanAnt[0]] + ' SPW=' + `spw[0]` + ' SCAN=' + `scan[0]`)
plt.xlabel('Azimuth Offset [arcsec]')
plt.ylabel('Elevation Offset [arcsec]')
#-------- Draw FWHM beam
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM);   plt.plot( circle_x, circle_y )
plt.savefig(prefix + '_Scan.pdf', form='pdf')
#plt.close()
