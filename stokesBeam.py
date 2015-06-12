execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import griddata
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
solution = np.load(wd + QUXY + 'QUXY.npy')
AzEl = np.load(wd + QUXY + '.Azel.npy')
PA = UnivariateSpline(AzEl[0], AzEl[3], s=1.0e-5)
CalQ, CalU, XYphs = solution[0], solution[1], solution[2]
Dxy = solution[3] + (1.0j)* solution[4]
#-------- Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
blNum  = antNum* (antNum - 1) / 2
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
wavelength = constants.c / Freq 
FWHM = 1.13* 180.0* 3600.0* wavelength / (12.0* pi) # Gaussian beam, in unit of arcsec
#-------- Gain Calibration for Refants
BlMap  = range(blNum)
BlInv = [False]* blNum      # True -> inverted baseline
blWeight = np.ones([blNum])
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    BlMap[bl_index], BlInv[bl_index] = Ant2BlD( AntIndex[ants[0]], AntIndex[ants[1]])
    blWeight[bl_index] = antWeight[AntIndex[ants[0]]]* antWeight[AntIndex[ants[1]]]
#
refBlIndex  = np.where(blWeight == 1.0)[0].tolist()
scanBlIndex = np.where(blWeight == 0.5)[0].tolist()
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], scan[0])
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
def gainComplex( vis ):
    #return clcomplex_solve( vis, 1.0e-8/abs(vis) )
    return clcomplex_solve( vis, vis_sigma/blWeight )
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
plt.subplot( 2, 2, 1); plt.contourf(xi, yi, Imap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesI')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 2); plt.contourf(xi, yi, Qmap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesQ')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 3); plt.contourf(xi, yi, Umap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesU')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 4); plt.contourf(xi, yi, Vmap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesV')
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
EVPAerr = GridData( np.arccos( (ScanQ[timeIndex]* CalQ + ScanU[timeIndex]* CalU) / (sqrt(ScanQ[timeIndex]**2 + ScanU[timeIndex]**2)* sqrt(CalQ**2 + CalU**2)))*90.0/pi, ScanAz[timeIndex], ScanEl[timeIndex], xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
plt.subplot( 2, 2, 1); plt.contourf(xi, yi, Qerr, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + ' Q_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 2); plt.contourf(xi, yi, Uerr, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + ' U_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 3); plt.contourf(xi, yi, Perr, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + ' P_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.subplot( 2, 2, 4); plt.contourf(xi, yi, EVPAerr, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + ' EVPA_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.savefig( prefix + '-StokesErr.pdf', form='pdf'); plt.close()
plt.close()
