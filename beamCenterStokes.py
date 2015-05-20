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
#-------- Time-based matching between time tags in visibilities and in scan pattern 
def AzElMatch( refTime, scanTime, thresh, Az, El ):
    index = np.where( abs(scanTime - refTime) < thresh)[0]
    return np.median(Az[index]), np.median(El[index])
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
def delayAnt( spec_BL_CH_TIME, chRange ):
    blNum = spec_BL_CH_TIME.shape[0]
    chNum = spec_BL_CH_TIME.shape[1]
    timeNum = spec_BL_CH_TIME.shape[2]
    antNum = Bl2Ant(blNum)[0]
    delay_bl= np.zeros([blNum, timeNum])
    visamp  = np.zeros([blNum, timeNum])
    delay_ant = np.zeros([antNum, timeNum])
    delay_err = np.zeros([antNum, timeNum])
    for time_index in range(timeNum)):
        for bl_index in range(blNum):
            delay_bl[bl_index], visamp[bl_index] = delay_search(XX[bl_index, chRange, time_index])
        #
        delay_ant[:, time_index], delay_err[:, time_index] = cldelay_solve( delay_bl* (1.0* chNum / len(chRange)), 1.0/visamp )
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
#-------- Antenna List
antList = GetAntName(msfile)
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
scanTime, AntID, Az, El = GetAzEl(msfile)
index = np.where( AntID == refAnt[0]); scanTime = scanTime[index]; Az = Az[index]; El = El[index]
#
AntIndex = refAnt
antNum = len(refAnt)
blNum  = antNum* (antNum - 1) / 2
antList = antList[AntIndex]
antWeight = np.ones(antNum)
#-------- Visibility sampling points
interval, timeStamp = GetTimerecord(msfile, 0, refAnt[0], scanAnt[0], spw[0], scan[0])
timeNum = len(timeStamp)
ScanAz = np.zeros([timeNum]); ScanEl = np.zeros([timeNum])
for time_index in range(timeNum):
    ScanAz[time_index], ScanEl[time_index] = AzElMatch( timeStamp[time_index], scanTime, np.median(interval), Az, El)
#
#-------- Frequency and Wavelength
chNum, chWid, Freq = GetChNum(msfile, spw[0])
wavelength = constants.c / np.mean(Freq)
FWHM = 1.13* 180.0* 3600.0* wavelength / (12.0* pi) # Gaussian beam, in unit of arcsec
vis_sigma = np.ones(blNum) / sqrt(2.0e9* np.median(diff(timeStamp)))
#-- baseline-based weights
blMap = range(blNum)
blInv = [False]* blNum      # True -> inverted baseline
blWeight = np.ones(blNum)
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    blMap[bl_index], blInv[bl_index] = Ant2BlD( AntIndex[ants[0]], AntIndex[ants[1]])
#
blWeight = np.ones([blNum])
#-------- Need Delay Cal?
if chNum > 31:
    chRange = range( int(chNum*0.05), int(chNum* 0.98))
    #-------- Visibilities
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], scan[0])
    timeNum = len(timeStamp); timeRange = range(timeNum)
    XX = Xspec[0, :, blMap]
    XY = Xspec[1, :, blMap]
    YX = Xspec[2, :, blMap]
    YY = Xspec[3, :, blMap]
    #
    #-------- Delay cal at beam center
    delay_XX = np.apply_along_axis( delayAnt, 2, XX[:,chRange])
    delay_YY = np.apply_along_axis( delayAnt, 2, YY[:,chRange])
    #delay_YY = delayAnt(YY, chRange, timeRange)
    #
    #delayXY = np.zeros([blNum, timeRange])
    #delayYX = np.zeros([blNum, timeRange])
    #delay_XY = np.zeros([blNum])
    #delay_YX = np.zeros([blNum])
    #for bl_index in range(blNum):
    #    for time_index in range( len(beamCenterIndex)):
    #        delayXY[bl_index, time_index], visamp = delay_search( XY[bl_index, chRange])
    #        delayYX[bl_index, time_index], visamp = delay_search( YX[bl_index, chRange])
    #    #
    #
    #fact = 1.0e9 / (2.0* np.median(abs(chWid))* len(chRange))
    #print '%s  SPW[%d] : XY_delay=%5.3f YX_delay=%5.3f [ns]' % (prefix, spw[0], np.mean(delayXY)* fact, np.mean(delayYX)* fact)
#
"""
#-------- Visibilities
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], scan[0])
XX = Xspec[0,0,blMap]
XY = Xspec[1,0,blMap]
YX = Xspec[2,0,blMap]
YY = Xspec[3,0,blMap]

#-------- Solve for antenna-based complex gains
def gainComplex( vis ):
    return clcomplex_solve( vis, vis_sigma/blWeight )
#
print '--- Solving Antenna-based Gain for XX'; GainXX = np.apply_along_axis( gainComplex, 0, XX )
print '--- Solving Antenna-based Gain for YY'; GainYY = np.apply_along_axis( gainComplex, 0, YY )
print '--- Solving Antenna-based Gain for XY'; GainXY = np.apply_along_axis( gainComplex, 0, XY )
print '--- Solving Antenna-based Gain for YX'; GainYX = np.apply_along_axis( gainComplex, 0, YX )
#
#-------- Gain Cal at the beam center
for subScan_index in range(subScanNum):
    timeSpan = subScanEndIndex[subScan_index] - subScanStartIndex[subScan_index] + 1
    timeRange = range(subScanStartIndex[subScan_index] + int(timeSpan/5), subScanEndIndex[subScan_index] - int(timeSpan/5))
    subScanAngle = atan2( ScanEl[subScanEndIndex[subScan_index]], ScanAz[subScanEndIndex[subScan_index]])
    cs, sn = cos(subScanAngle), sin(subScanAngle)
    relativeScanPos = cs* ScanAz[timeRange] + sn* ScanEl[timeRange]
    for ant_index in range(len(refAnt)):         #-------- Loop in ref ants
        GainXX[ant_index, timeRange] /= np.mean(GainXX[ant_index, timeRange])
        GainYY[ant_index, timeRange] /= np.mean(GainYY[ant_index, timeRange])
        GainXY[ant_index, timeRange] /= np.mean(GainXY[ant_index, timeRange])
        GainYX[ant_index, timeRange] /= np.mean(GainYX[ant_index, timeRange])
    #
    for ant_index in range(len(refAnt), antNum): #-------- Loop in scanning ants
        gaussParamX = simpleGaussFit(relativeScanPos, abs(GainXX[ant_index, timeRange]))
        gaussParamY = simpleGaussFit(relativeScanPos, abs(GainYY[ant_index, timeRange]))

        phaseX = np.angle(GainXX[ant_index, beamCenterIndex[subScan_index]])
        phaseY = np.angle(GainYY[ant_index, beamCenterIndex[subScan_index]])
        GainXX[ant_index, timeRange] /= (gaussParamX[0]* exp( (0.0+1.0j)* phaseX))
        GainYY[ant_index, timeRange] /= (gaussParamY[0]* exp( (0.0+1.0j)* phaseY))
        GainXY[ant_index, timeRange] /= ( sqrt(gaussParamX[0]* gaussParamY[0])* exp( (0.0+1.0j)* (phaseX - phaseY)) )
        GainYX[ant_index, timeRange] /= ( sqrt(gaussParamX[0]* gaussParamY[0])* exp( (0.0+1.0j)* (phaseY - phaseX)) )
        text_sd = '%d SubScan%d : Peak=%f  Offset=%f  Width=%f' % (ant_index, subScan_index, gaussParamX[0], gaussParamX[1], gaussParamX[2])
        print text_sd
    #
#
#-------- Plot
StokesI = 0.5*(abs(GainXX)**2 + abs(GainYY)**2)
StokesQ = 0.5*(abs(GainXX)**2 - abs(GainYY)**2)
StokesU = 0.5*(GainXY.real + GainYX.real)
StokesV = 0.5*(GainYX.imag - GainXY.imag)

for index in range(len(refAnt), antNum):
    fig = plt.figure( figsize = (10,8))
    xi, yi = np.mgrid[ -floor(2.0*FWHM):floor(2.0*FWHM):128j, -floor(2.0*FWHM):floor(2.0*FWHM):128j]
    #Imap = griddata( (ScanAz, ScanEl), StokesI[17], (xi, yi), method='linear')
    Imap = GridData( StokesI[index], ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5.0 ).reshape(len(xi), len(xi))
    Qmap = griddata( (ScanAz, ScanEl), StokesQ[index], (xi, yi), method='linear')
    Umap = griddata( (ScanAz, ScanEl), StokesU[index], (xi, yi), method='linear')
    Vmap = griddata( (ScanAz, ScanEl), StokesV[index], (xi, yi), method='linear')
    plt.subplot( 2, 2, 1); plt.contourf(xi, yi, Imap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-' + antList[AntIndex[index]] + '-StokesI')
    plt.subplot( 2, 2, 2); plt.contourf(xi, yi, Qmap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-' + antList[AntIndex[index]] + '-StokesQ')
    plt.subplot( 2, 2, 3); plt.contourf(xi, yi, Umap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-' + antList[AntIndex[index]] + '-StokesU')
    plt.subplot( 2, 2, 4); plt.contourf(xi, yi, Vmap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-' + antList[AntIndex[index]] + '-StokesV')
    plt.plot( ScanAz, ScanEl, '.', color='k', alpha=0.1)
    plt.axis([-2.0*FWHM, 2.0*FWHM, -2.0*FWHM, 2.0*FWHM])
    plt.savefig( prefix + '-' + antList[AntIndex[index]] + '.pdf', form='pdf'); plt.close()
#
"""
