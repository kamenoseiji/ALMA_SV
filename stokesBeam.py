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
    time_index = np.where( AntID == scanAntIndex[0] )[0].tolist()
    return refAntIndex.tolist(), scanAntIndex.tolist(), Time[time_index], Offset[:, time_index]
#
#-------- Read Grid Points
def readGrid( file ):
    file_ptr = open( file, 'r')
    #labels = file_ptr
    index = 0
    Az = []
    El = []
    for line in open( file ):
        items = line.split()
        if index == 0:
            label = items
        else:
            Az.append( float(items[1]))
            El.append( float(items[2]))
        #
        index = index + 1
    #
    return Az, El
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
#-------- Scanning Offset < threshold
def scanThresh( msfile, ant_index, thresh ):
    Time, AntID, Offset = GetAzOffset(msfile)
    time_index = np.where( AntID == ant_index )[0]
    onAxisIndex = np.where( Offset[0, time_index]**2 + Offset[1, time_index]**2 < thresh**2 )[0]
    return time_index[onAxisIndex].tolist()
#
#-------- Time-based matching
def timeMatch( refTime, scanTime, thresh):
    match = np.where( abs(scanTime - refTime) < thresh)[0].tolist()
    return len(match)
#
#-------- SmoothGain
def smoothGain( timeValue, complexValue ):
    realSP = UnivariateSpline( timeValue, complexValue.real )
    imagSP = UnivariateSpline( timeValue, complexValue.imag )
    return realSP, imagSP
#
#-------- For drawing a circle
def circlePoints( x, y, radius ):
    angle = np.arange(-pi, (130/128)*pi, pi/128)
    return x + radius* np.cos(angle), y + radius* np.sin(angle)
#
#-------- Antenna-based time-averated delay
def delayAnt( spec_BL_CH_TIME, chRange, timeRange ):
    chNum = spec_BL_CH_TIME.shape[1]
    blNum = spec_BL_CH_TIME.shape[0]
    antNum = Bl2Ant(blNum)[0]
    delay_bl= np.zeros([blNum])
    visamp  = np.zeros([blNum])
    delay_ant = np.zeros([antNum, len(timeRange)])
    delay_err = np.zeros([antNum, len(timeRange)])
    for time_index in range(len(timeRange)):
        for bl_index in range(blNum):
            delay_bl[bl_index], visamp[bl_index] = delay_search(XX[bl_index, chRange, timeRange[time_index]])
        #
        delay_ant[:, time_index], delay_err[:, time_index] = cldelay_solve( delay_bl* (1.0* chNum / len(chRange)), visamp )
    #
    return np.median(delay_ant, axis=1)
#
#-------- XY delay
def XYdel( spec_BL_CH_TIME, chRange, timeRange, delay_ant_XX, delay_ant_YY ):
    chNum = spec_BL_CH_TIME.shape[1]
    blNum = spec_BL_CH_TIME.shape[0]
    antNum = Bl2Ant(blNum)[0]
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
#----------------------------------------- Procedures
msfile = wd + prefix + '.ms'
solution = np.load(wd + QUXY + '.QUXY.npy')
AzEl = np.load(wd + QUXY + '.Azel.npy')
PA = UnivariateSpline(AzEl[0], AzEl[3], s=1.0e-5)
CalQ, CalU, XYphs = solution[0], solution[1], solution[2]
Dxy = solution[3] + (1.0j)* solution[4]
Dyx = solution[5] + (1.0j)* solution[6]
#-------- Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
blNum  = antNum* (antNum - 1) / 2
AntD = np.zeros([antNum])
for ant_index in range(antNum):
    AntD[ant_index] = GetAntD(antList[ant_index])
#
print('Checking the Array ....')
#-------- Reference and Scan Antennas
refAnt, scanAnt, scanTime, Offset = antRefScan(msfile)
print('-------- Reference Antennas ----')
for ant_index in refAnt:
    text_sd = 'Ref[%d]  / %d: %s ' % (ant_index, len(refAnt), antList[ant_index])
    print text_sd
#
print('-------- Scanning Antennas ----')
for ant_index in scanAnt:
    text_sd = 'Scan[%d] / %d: %s ' % (ant_index, len(scanAnt), antList[ant_index])
    print text_sd
#
AntIndex = refAnt + scanAnt
antWeight = np.ones(antNum)
antWeight[scanAnt] = 0.5
#-------- Visibility sampling points
interval, timeStamp = GetTimerecord(msfile, 0, refAnt[0], scanAnt[0], spw[0], scan[0])
vis_sigma = np.ones(blNum) / sqrt(2.0e9* np.median(diff(timeStamp)))
timeNum = len(timeStamp)
ScanAz = np.zeros([timeNum]); ScanEl = np.zeros([timeNum])
for time_index in range(timeNum):
    ScanAz[time_index], ScanEl[time_index] = AzElMatch( timeStamp[time_index], scanTime, np.median(interval), Offset)
#
#-------- Frequency and Wavelength
chNum, chWid, Freq = GetChNum(msfile, spw[0])
FWHM = GetFWHM(msfile, spw[0], AntD)
if chNum > 1:
    chRange = range( int(0.05*chNum), int(0.95*chNum))
#
#-------- Center position of scanning antennas
scanTime, AntID, az, el = GetAzEl(msfile)
index = scanThresh(msfile, scanAnt[0], FWHM[scanAnt[0]]/20)
scanTime = scanTime[index]
matchNum = np.zeros([timeNum])
for time_index in range(timeNum):
    matchNum[time_index] = timeMatch( timeStamp[time_index], scanTime, np.median(interval))
#
onAxisIndex = np.where( matchNum > 0 )[0].tolist()
#-------- Load Visibilities
BlMap  = range(blNum)
BlInv = [False]* blNum      # True -> inverted baseline
blWeight = np.ones([blNum])
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    BlMap[bl_index], BlInv[bl_index] = Ant2BlD(ants[0], ants[1])
    blWeight[bl_index] = antWeight[ants[0]]* antWeight[ants[1]]
    #BlMap[bl_index], BlInv[bl_index] = Ant2BlD( AntIndex[ants[0]], AntIndex[ants[1]])
    #blWeight[bl_index] = antWeight[AntIndex[ants[0]]]* antWeight[AntIndex[ants[1]]]
#
refBlIndex  = np.where(blWeight == 1.0)[0].tolist()
scanBlIndex = np.where(blWeight == 0.5)[0].tolist()
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], scan[0])
if chNum == 1:
    temp = Xspec[:,0]
else:
    temp = np.mean(Xspec[:,chRange], axis=1)
#
XX = temp[0,BlMap]
XY = temp[1,BlMap]
YX = temp[2,BlMap]
YY = temp[3,BlMap]
for bl_index in range(blNum):
    if BlInv[bl_index]:
        XX[bl_index] = XX[bl_index].conjugate()
        XY[bl_index] = XY[bl_index].conjugate()
        YX[bl_index] = YX[bl_index].conjugate()
        YY[bl_index] = YY[bl_index].conjugate()
    #
#
GainX, GainY = polariGain(XX, YY, PA(timeStamp), solution[0], solution[1])
Ucos_Qsin = solution[1]* np.cos(2.0*PA(timeStamp)) - solution[0]* np.sin(2.0*PA(timeStamp))
Qcos_Usin = solution[0]* np.cos(2.0*PA(timeStamp)) + solution[1]* np.sin(2.0*PA(timeStamp))
VisXX = gainCalVis( XX, GainX, GainX )
VisYY = gainCalVis( YY, GainY, GainY )
VisXY = gainCalVis( XY, GainX, GainY )* np.exp(0.0 - solution[2]*1.0j)
VisYX = gainCalVis( YX, GainY, GainX )* np.exp(0.0 + solution[2]*1.0j)
DxDx = VisXX - 1.0 - Qcos_Usin
DxDy = VisXY - Ucos_Qsin
DyDx = VisYX - Ucos_Qsin
DyDy = VisYY - 1.0 + Qcos_Usin

print('-------- Determining Antenna-based D-terms ----')
Pre = np.zeros([4*blNum, 2*antNum]) # [ReXX, ReXY, ReYX, ReYY] x [ReDx, ReDy]
Pim = np.zeros([4*blNum, 2*antNum]) # [ReXX, ReXY, ReYX, ReYY] x [ReDx, ReDy]
Dx = np.zeros([antNum, timeNum], dtype=complex)
Dy = np.zeros([antNum, timeNum], dtype=complex)
#
blGain = np.ones([blNum])
GainScaleX = 0.5 / np.max(abs(GainX))
#GainScaleY = 0.2 / np.max(abs(GainY))
for time_index in range(timeNum):
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        #blGain[bl_index] = GainScaleX* abs(GainX[ants[1], time_index])* abs(GainX[ants[0], time_index]) + GainScaleY* abs(GainY[ants[1], time_index])* abs(GainY[ants[0], time_index])
        blGain[bl_index] = GainScaleX* abs(GainX[ants[1], time_index])* abs(GainX[ants[0], time_index])
    #
    W  = np.diag(np.r_[blGain, blGain, blGain, blGain])
    #---- Real Part
    Pre[        0:   blNum ,      0:   antNum ] = Ucos_Qsin[time_index]* BlAmpMatrix(antNum)
    Pre[   blNum :(2*blNum),      0:   antNum ] = (1.0 - Qcos_Usin[time_index])* DxMatrix(antNum)
    Pre[   blNum :(2*blNum), antNum:(2*antNum)] = (1.0 + Qcos_Usin[time_index])* DyMatrix(antNum)
    Pre[(2*blNum):(3*blNum),      0:   antNum ] = (1.0 - Qcos_Usin[time_index])* DyMatrix(antNum)
    Pre[(2*blNum):(3*blNum), antNum:(2*antNum)] = (1.0 + Qcos_Usin[time_index])* DxMatrix(antNum)
    Pre[(3*blNum):(4*blNum), antNum:(2*antNum)] = Ucos_Qsin[time_index]* BlAmpMatrix(antNum)
    #---- Imag Part
    Pim[        0:   blNum ,      0:   antNum ] = Ucos_Qsin[time_index]* (DxMatrix(antNum) - DyMatrix(antNum))
    Pim[   blNum :(2*blNum),      0:   antNum ] = (1.0 - Qcos_Usin[time_index])* DxMatrix(antNum)
    Pim[   blNum :(2*blNum), antNum:(2*antNum)] = (1.0 - Qcos_Usin[time_index])* DyMatrix(antNum)
    Pim[(2*blNum):(3*blNum),      0:   antNum ] = (1.0 + Qcos_Usin[time_index])* DyMatrix(antNum)
    Pim[(2*blNum):(3*blNum), antNum:(2*antNum)] = (1.0 + Qcos_Usin[time_index])* DxMatrix(antNum)
    Pim[(3*blNum):(4*blNum), antNum:(2*antNum)] = Ucos_Qsin[time_index]* (DyMatrix(antNum) - DxMatrix(antNum))
    #
    PreWPre = np.dot( Pre.T, np.dot(W, Pre))
    PimWPim = np.dot( Pim.T, np.dot(W, Pim))
    DetPre = np.linalg.det( PreWPre )
    DetPim = np.linalg.det( PimWPim )
    ReXXYY = np.dot(W, np.r_[ DxDx[:,time_index].real, DxDy[:,time_index].real, DyDx[:,time_index].real, DyDy[:,time_index].real ])
    ImXXYY = np.dot(W, np.r_[ DxDx[:,time_index].imag, DxDy[:,time_index].imag, DyDx[:,time_index].imag, DyDy[:,time_index].imag ])
    print( `time_index` + '  det(PTP_real) = ' + `DetPre` + '  det(PTP_imag) = ' + `DetPim` )
    ReD = np.dot( np.linalg.inv( PreWPre ), np.dot( Pre.T, ReXXYY) )
    if DetPim != 0.0:
        ImD = np.dot( np.linalg.inv( PimWPim ), np.dot( Pim.T, ImXXYY) )
    else:
        ImD = np.dot( np.diag(1.0 / np.diag(PimWPim)), np.dot( Pim.T, ImXXYY) )
    #
    Dx[:, time_index] = ReD[0:antNum] + 1.0j* ImD[0:antNum]
    Dy[:, time_index] = ReD[antNum:(2*antNum)] + 1.0j* ImD[antNum:(2*antNum)]
#
"""
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    p[          bl_index,            ants[1]] =  1.0      # ReXY / ReDx
    p[  blNum + bl_index,   antNum + ants[1]] =  1.0      # ImXY / ImDx
    p[2*blNum + bl_index, 2*antNum + ants[1]] =  1.0      # ReYX / ReDy
    p[3*blNum + bl_index, 3*antNum + ants[1]] =  1.0      # ImYX / ImDy
    #
    p[          bl_index, 2*antNum + ants[0]] =  1.0      # ReXY / ReDy
    p[  blNum + bl_index, 3*antNum + ants[0]] = -1.0      # ImXY / ImDy
    p[2*blNum + bl_index,            ants[0]] =  1.0      # ReYX / ReDx
    p[3*blNum + bl_index,   antNum + ants[0]] = -1.0      # ImYX / ImDx
#
"""
"""
Dx = np.zeros([antNum, timeNum], dtype=complex)
Dy = np.zeros([antNum, timeNum], dtype=complex)
for time_index in range(timeNum):
    temp = np.r_[ DxDy[:,time_index].real, DxDy[:,time_index].imag, DyDx[:,time_index].real, DyDx[:,time_index].imag ]
    D = np.dot( np.linalg.inv(np.dot( p.T, p ) - np.diag(0.1*np.ones([4*antNum]))), np.dot( p.T, temp) )
    Dx[:,time_index] = D[0:antNum] + 1.0j* D[antNum:(2*antNum)]
    Dy[:,time_index] = D[(2*antNum):(3*antNum)] + 1.0j* D[(3*antNum):(4*antNum)]
#
"""


"""

#-------- PA and Phase
csPA = np.cos(2.0* PA(timeStamp))
snPA = np.sin(2.0* PA(timeStamp))
XX_YY = CalQ* csPA + CalU* snPA     # (XX - YY)/2
#--------
XX = Xspec[0,0,BlMap]
XY = Xspec[1,0,BlMap]
YX = Xspec[2,0,BlMap]
YY = Xspec[3,0,BlMap]
for bl_index in range(blNum):
    if BlInv[bl_index]:
        XX[bl_index] = XX[bl_index].conjugate()
        XY[bl_index] = XY[bl_index].conjugate()
        YX[bl_index] = YX[bl_index].conjugate()
        YY[bl_index] = YY[bl_index].conjugate()
    #
#
#def gainComplex( vis ):
#    #return clcomplex_solve( vis, 1.0e-8/abs(vis) )
#    return clcomplex_solve( vis, vis_sigma/blWeight )
#
GainX = np.apply_along_axis( gainComplex, 0, XX )/ sqrt(1.0 + XX_YY)
GainY = np.apply_along_axis( gainComplex, 0, YY )* exp( 1.0j* XYphs) / sqrt(1.0 - XX_YY)

VisXX = gainCalVis( XX, GainX, GainX)
VisXY = gainCalVis( XY, GainX, GainY)
VisYX = gainCalVis( YX, GainY, GainX)
VisYY = gainCalVis( YY, GainY, GainY)

refVisXX  = np.mean( VisXX[refBlIndex], axis=0)
#refVisXY  = (np.mean( VisXY[refBlIndex], axis=0) - Dxy)* exp(-1.0j* XYphs)
#refVisYX  = (np.mean( VisYX[refBlIndex], axis=0) + Dxy)* exp( 1.0j* XYphs)
refVisXY  = (np.mean( VisXY[refBlIndex], axis=0) - Dxy)
refVisYX  = (np.mean( VisYX[refBlIndex], axis=0) + Dxy)
refVisYY  = np.mean( VisYY[refBlIndex], axis=0)

ScanVisXX = np.mean( VisXX[scanBlIndex], axis=0)
#ScanVisXY = (np.mean( VisXY[scanBlIndex], axis=0) - Dxy)* exp(-1.0j* XYphs)
#ScanVisYX = (np.mean( VisYX[scanBlIndex], axis=0) + Dxy)* exp( 1.0j* XYphs)
ScanVisXY = (np.mean( VisXY[scanBlIndex], axis=0) - Dxy)
ScanVisYX = (np.mean( VisYX[scanBlIndex], axis=0) + Dxy)
ScanVisYY = np.mean( VisYY[scanBlIndex], axis=0)

refI = 0.5*(refVisXX + refVisYY).real
refQ = 0.5*(  csPA* (refVisXX - refVisYY).real - snPA* (refVisXY + refVisYX).real )
refU = 0.5*(  snPA* (refVisXX - refVisYY).real + csPA* (refVisXY + refVisYX).real )
refV = 0.5*( refVisXY - refVisYX).imag 
ScanI = 0.5*(ScanVisXX + ScanVisYY).real
ScanQ = 0.5*(  csPA* (ScanVisXX - ScanVisYY).real - snPA* (ScanVisXY + ScanVisYX).real )
ScanU = 0.5*(  snPA* (ScanVisXX - ScanVisYY).real + csPA* (ScanVisXY + ScanVisYX).real )
ScanV = 0.5*( ScanVisXY - ScanVisYX).imag 
logfile = open(prefix + '_Stokes.log', 'w')
text_sd = 'mjdSec       dAz    dEl   PA      I       Q       U       V'
logfile.write(text_sd + '\n')
for time_index in range(timeNum):
    text_sd = '%d %6.1f %6.1f %6.4f %7.4f %7.4f %7.4f %7.4f' % (timeStamp[time_index], ScanAz[time_index], ScanEl[time_index], PA(timeStamp[time_index]), ScanI[time_index], ScanQ[time_index], ScanU[time_index], ScanV[time_index])
    logfile.write(text_sd + '\n')
#
logfile.close()
#-------- SubScan
subScanStartIndex, subScanEndIndex = subScan(timeStamp, 1.5*np.median(interval))
subScanNum = len(subScanStartIndex)
beamCenterIndex = np.where( ScanAz**2 + ScanEl**2 < (FWHM/10.0)**2)[0].tolist()
#
#-------- Gain Cal at the beam center
for subScan_index in range(subScanNum):
    timeSpan = subScanEndIndex[subScan_index] - subScanStartIndex[subScan_index] + 1
    timeRange = range(subScanStartIndex[subScan_index] + int(timeSpan/5), subScanEndIndex[subScan_index] - int(timeSpan/5))
    subScanAngle = atan2( ScanEl[subScanEndIndex[subScan_index]], ScanAz[subScanEndIndex[subScan_index]])
    cs, sn = cos(subScanAngle), sin(subScanAngle)
    relativeScanPos = cs* ScanAz[timeRange] + sn* ScanEl[timeRange]
#
#-------- Plot
def circlePoints( x, y, radius ):
    angle = np.arange(-pi, (130/128)*pi, pi/128)
    return x + radius* np.cos(angle), y + radius* np.sin(angle)
#
fig = plt.figure( figsize = (10,8))
xi, yi = np.mgrid[ -floor(2.0*FWHM):floor(2.0*FWHM):128j, -floor(2.0*FWHM):floor(2.0*FWHM):128j]
#Imap = griddata( (ScanAz, ScanEl), ScanI, (xi, yi), method='cubic')
#Qmap = griddata( (ScanAz, ScanEl), ScanQ, (xi, yi), method='linear')
#Umap = griddata( (ScanAz, ScanEl), ScanU, (xi, yi), method='linear')
#Vmap = griddata( (ScanAz, ScanEl), ScanV, (xi, yi), method='linear')
Imap = GridData( ScanI, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Qmap = GridData( ScanQ, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Umap = GridData( ScanU, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Vmap = GridData( ScanV, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
plt.subplot( 2, 2, 1); plt.contourf(xi, yi, Imap, np.linspace(0.9, 1.1, 21), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesI')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 2); plt.contourf(xi, yi, Qmap, np.linspace(-0.12, 0.12, 25), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesQ')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 3); plt.contourf(xi, yi, Umap, np.linspace(-0.12, 0.12, 25), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesU')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 4); plt.contourf(xi, yi, Vmap, np.linspace(-0.05, 0.05, 21), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesV')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.plot( ScanAz, ScanEl, '.', color='k', alpha=0.1)
plt.axis([-2.0*FWHM, 2.0*FWHM, -2.0*FWHM, 2.0*FWHM])
plt.savefig( prefix + '-StokesMap.pdf', form='pdf'); plt.close()
plt.close()
#
#--------- Gridding for 11x11 sampling points
az, el = readGrid(gridFile)
logfile = open(prefix + '_StokesGrid.log', 'w')
text_sd = 'No. dAz      dEl        I        Q        U         V'
logfile.write(text_sd + '\n')
xi, yi = np.mgrid[ min(az):max(az):11j, max(el):min(el):11j]
Imap = GridData( ScanI, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Qmap = GridData( ScanQ, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Umap = GridData( ScanU, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Vmap = GridData( ScanV, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
for x_index in range(11):
    for y_index in range(11):
        text_sd = '%d %f %f %f %f %f %f' % (x_index*11 + y_index + 1, xi[x_index, y_index], yi[x_index, y_index], Imap[x_index, y_index], Qmap[x_index, y_index], Umap[x_index, y_index], Vmap[x_index, y_index])
        logfile.write(text_sd + '\n')
    #
#
logfile.close()
#-------- Plot within 2xFWHM
timeIndex = np.where( ScanAz**2 + ScanEl**2 < FWHM**2 )[0]
xi, yi = np.mgrid[ -floor(FWHM):floor(FWHM):128j, -floor(FWHM):floor(FWHM):128j]
fig = plt.figure( figsize = (10,8))
Qerr = GridData( ScanQ[timeIndex]-CalQ, ScanAz[timeIndex], ScanEl[timeIndex], xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Uerr = GridData( ScanU[timeIndex]-CalU, ScanAz[timeIndex], ScanEl[timeIndex], xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Perr = GridData( sqrt(ScanQ[timeIndex]**2 + ScanU[timeIndex]**2) - sqrt(CalQ**2 + CalU**2), ScanAz[timeIndex], ScanEl[timeIndex], xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
EVPAerr = GridData( np.arccos( (ScanQ[timeIndex]* CalQ + ScanU[timeIndex]* CalU) / sqrt( (ScanQ[timeIndex]**2 + ScanU[timeIndex]**2)* (CalQ**2 + CalU**2)))*90.0/pi, ScanAz[timeIndex], ScanEl[timeIndex], xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
#
plt.subplot( 2, 2, 1); plt.contourf(xi, yi, Qerr, np.linspace(-0.05, 0.05, 21, endpoint=True), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + ' Q_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 2); plt.contourf(xi, yi, Uerr, np.linspace(-0.05, 0.05, 21), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + ' U_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 3); plt.contourf(xi, yi, Perr, np.linspace(-0.05, 0.05, 21), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + ' P_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 4); plt.contourf(xi, yi, EVPAerr, np.linspace(-10, 10, 21), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + ' EVPA_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.savefig( prefix + '-StokesErr.pdf', form='pdf'); plt.close()
plt.close()
"""
