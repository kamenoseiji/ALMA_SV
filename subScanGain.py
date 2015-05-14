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
print('Checking the Array ....')
#-------- Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
blNum  = antNum* (antNum - 1) / 2
#-------- Reference and Scan Antennas
refAnt, scanAnt, scanTime, Offset = antRefScan(msfile)
print('-------- Reference Antennas ----')
for ant_index in refAnt:
    text_sd = 'Ref[%d]  / %d: %s ' % (ant_index, len(refAnt), antList[ant_index])
    print text_sd
#
for ant_index in scanAnt:
    text_sd = 'Scan[%d] / %d: %s ' % (ant_index, len(scanAnt), antList[ant_index])
    print text_sd
#
AntIndex = refAnt + scanAnt
antWeight = np.ones(antNum)
antWeight[scanAnt] = 0.1
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
#-------- Pick central position
beamCenterIndex = []
for subScanIndex in range(subScanNum):
    time_index = range( subScanStartIndex[subScanIndex], subScanEndIndex[subScanIndex]+1 )
    beamCenterIndex.append( time_index[argmin( ScanAz[time_index]**2 + ScanEl[time_index]**2 )] )
#
vis_sigma = np.ones(blNum) / sqrt(2.0e9* np.median(diff(timeStamp)))
#-- baseline-based weights
blMap = range(blNum)
blInv = [False]* blNum      # True -> inverted baseline
blWeight = np.ones(blNum)
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    blMap[bl_index], blInv[bl_index] = Ant2BlD( AntIndex[ants[0]], AntIndex[ants[1]])
    blWeight[bl_index] = antWeight[AntIndex[ants[0]]]* antWeight[AntIndex[ants[1]]]
#
#-------- Need Delay Cal?
if len(delaySPW) != 0:
    #-------- Visibilities
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], scan[0])
    XX = Xspec[0, :, blMap]
    XY = Xspec[1, :, blMap]
    YX = Xspec[2, :, blMap]
    YY = Xspec[3, :, blMap]
    #
    #-------- Delay cal at beam center
    delay_XX = delayAnt(XX, chRange, beamCenterIndex)
    delay_YY = delayAnt(YY, chRange, beamCenterIndex)
    #
    delayXY = np.zeros([blNum, len(beamCenterIndex)])
    delayYX = np.zeros([blNum, len(beamCenterIndex)])
    delay_XY = np.zeros([blNum])
    delay_YX = np.zeros([blNum])
    for bl_index in range(blNum):
        for time_index in range( len(beamCenterIndex)):
            delayXY[bl_index, time_index], visamp = delay_search( XY[bl_index, chRange, beamCenterIndex[time_index]])
            delayYX[bl_index, time_index], visamp = delay_search( YX[bl_index, chRange, beamCenterIndex[time_index]])
        #
    #
    fact = 1.0e9 / (2.0* np.median(abs(chWid))* len(chRange))
    print '%s  SPW[%d] : XY_delay=%5.3f YX_delay=%5.3f [ns]' % (prefix, spw[0], np.mean(delayXY)* fact, np.mean(delayYX)* fact)
#
#-------- Visibilities
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], scan[0])
XX = Xspec[0,0,blMap]
XY = Xspec[1,0,blMap]
YX = Xspec[2,0,blMap]
YY = Xspec[3,0,blMap]
GainXX = np.zeros([antNum, timeNum], dtype=complex)
GainXY = np.zeros([antNum, timeNum], dtype=complex)
GainYX = np.zeros([antNum, timeNum], dtype=complex)
GainYY = np.zeros([antNum, timeNum], dtype=complex)
for time_index in range(timeNum):
    #try:
    GainXX[:, time_index] =  gainAnt( XX[:, time_index], vis_sigma/blWeight)
    GainYY[:, time_index] =  gainAnt( YY[:, time_index], vis_sigma/blWeight)
    GainXY[:, time_index] =  gainAnt( XY[:, time_index], vis_sigma/blWeight)
    GainXY[:, time_index] =  gainAnt( YX[:, time_index], vis_sigma/blWeight)
    #except:
    #    GainXX[:, time_index] = 0.0
    #
#
#-------- Gain Cal at the beam center
for subScan_index in range(subScanNum):
    timeSpan = subScanEndIndex[subScan_index] - subScanStartIndex[subScan_index] + 1
    timeRange = range(subScanStartIndex[subScan_index] + int(timeSpan/5), subScanEndIndex[subScan_index] - int(timeSpan/5))
    subScanAngle = atan2( ScanEl[subScanEndIndex[subScan_index]], ScanAz[subScanEndIndex[subScan_index]])
    cs, sn = cos(subScanAngle), sin(subScanAngle)
    relativeScanPos = cs* ScanAz[timeRange] + sn* ScanEl[timeRange]
    for ant_index in range(len(refAnt), antNum):
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
StokesI = abs(GainXX)**2 + abs(GainYY)**2
StokesQ = abs(GainXX)**2 - abs(GainYY)**2
StokesU = GainXY.real + GainYX.real
StokesV = GainXY.imag - GainYX.imag


for scanAntIndex in scanAnt:
    fig = plt.figure( figsize = (10,8))
    xi, yi = np.mgrid[ -floor(2.0*FWHM):floor(2.0*FWHM):128j, -floor(2.0*FWHM):floor(2.0*FWHM):128j]
    #Imap = griddata( (ScanAz, ScanEl), StokesI[17], (xi, yi), method='linear')
    Imap = GridData( StokesI[scanAntIndex], ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5.0 ).reshape(len(xi), len(xi))
    Qmap = griddata( (ScanAz, ScanEl), StokesQ[scanAntIndex], (xi, yi), method='linear')
    Umap = griddata( (ScanAz, ScanEl), StokesU[scanAntIndex], (xi, yi), method='linear')
    Vmap = griddata( (ScanAz, ScanEl), StokesV[scanAntIndex], (xi, yi), method='linear')
    plt.subplot( 2, 2, 1); plt.contourf(xi, yi, Imap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-' + antList[scanAntIndex] + '-StokesI')
    plt.subplot( 2, 2, 2); plt.contourf(xi, yi, Qmap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-' + antList[scanAntIndex] + '-StokesQ')
    plt.subplot( 2, 2, 3); plt.contourf(xi, yi, Umap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-' + antList[scanAntIndex] + '-StokesU')
    plt.subplot( 2, 2, 4); plt.contourf(xi, yi, Vmap, 25, cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-' + antList[scanAntIndex] + '-StokesV')
    plt.plot( ScanAz, ScanEl, '.', color='k', alpha=0.1)
    plt.axis([-2.0*FWHM, 2.0*FWHM, -2.0*FWHM, 2.0*FWHM])
    plt.savefig( prefix + '-' + antList[scanAntIndex] + '.pdf', form='pdf'); plt.close()
#

"""
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    delay_XY[bl_index] = np.median(delayXY, axis=1)[bl_index] # - delay_XX[ ants[0]] + delay_YY[ ants[1]]
    delay_YX[bl_index] = np.median(delayYX, axis=1)[bl_index] # + delay_XX[ ants[1]] - delay_YY[ ants[0]]
#

BLdelayXX = np.zeros([blNum, len(beamCenterIndex)])
BLdelayYY = np.zeros([blNum, len(beamCenterIndex)])
for bl_index in range(blNum):
    for time_index in range( len(beamCenterIndex)):
        BLdelayXX[bl_index, time_index], visamp = delay_search( XX[bl_index, chRange, beamCenterIndex[time_index]])
        BLdelayYY[bl_index, time_index], visamp = delay_search( YY[bl_index, chRange, beamCenterIndex[time_index]])
    #
#
"""



"""
#-------- Antenna-based Gain at the beam center
subScanPhaseXX = np.zeros
phs = np.zeros([antNum, timeNum]); phserr = np.zeros([antNum, timeNum])
for time_index in beamCenterIndex:
    phs[:, time_index], phserr[:, time_index] = clphase_solve( np.angle(XX[:,time_index]), abs(XX[:,time_index]) )
#
GreSP, GimSP = smoothGain( timeStamp[beamCenterIndex], np.exp(phs[1,beamCenterIndex]* 1.0j) )
GreSP.set_smoothing_factor(0.1); GimSP.set_smoothing_factor(0.1)
plt.plot( timeStamp, np.arctan2( GimSP(timeStamp), GreSP(timeStamp)), '-')
#-------- In-beam index
inBeamIndex = np.where( ScanAz**2 + ScanEl**2 < 4*FWHM**2 )[0]
inBeamNum   = len(inBeamIndex)
#-------- Prepare Figure
fig = plt.figure( figsize = (8,8))
plt.plot( ScanAz[inBeamIndex], ScanEl[inBeamIndex], '.', color="cyan" )
plt.axis([-1.5*FWHM, 1.5*FWHM, -1.5*FWHM, 1.5*FWHM])
plt.title( prefix + ' ' + antList[scanAnt[0]] + ' SPW=' + `spw[0]` + ' SCAN=' + `scan[0]`)
plt.xlabel('Azimuth Offset [arcsec]')
plt.ylabel('Elevation Offset [arcsec]')
#-------- Draw FWHM beam
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM);   plt.plot( circle_x, circle_y )
plt.savefig(prefix + '_ScanBeam.pdf', form='pdf')
plt.close()

#--------
GainXX = np.zeros([antNum, inBeamNum], dtype=complex)
GainXY = np.zeros([antNum, inBeamNum], dtype=complex)
GainYX = np.zeros([antNum, inBeamNum], dtype=complex)
GainYY = np.zeros([antNum, inBeamNum], dtype=complex)
for time_index in range(inBeamNum):
    print '%d / %d ' % (time_index, inBeamNum)
    #GainXX[:, time_index] = gainAnt( XX[:, inBeamIndex[time_index]], vis_sigma/blWeight)
    #GainYY[:, time_index] = gainAnt( YY[:, inBeamIndex[time_index]], vis_sigma/blWeight)
    GainXX[0, time_index] = XX[0, inBeamIndex[time_index]]
    GainXY[0, time_index] = XY[0, inBeamIndex[time_index]]
    GainYX[0, time_index] = YX[0, inBeamIndex[time_index]]
    GainYY[0, time_index] = YY[0, inBeamIndex[time_index]]
#
fig = plt.figure( figsize = (10,8))
xi, yi = np.mgrid[ -floor(FWHM):floor(FWHM):128j, -floor(FWHM):floor(FWHM):128j]
XXamp = griddata( (ScanAz[inBeamIndex], ScanEl[inBeamIndex]), abs( GainXX[0]) / max(abs( GainXX[0])), (xi, yi), method='cubic')
XYamp = griddata( (ScanAz[inBeamIndex], ScanEl[inBeamIndex]), abs( GainXY[0]) / max(abs( GainXY[0])), (xi, yi), method='cubic')
YXamp = griddata( (ScanAz[inBeamIndex], ScanEl[inBeamIndex]), abs( GainYX[0]) / max(abs( GainYX[0])), (xi, yi), method='cubic')
YYamp = griddata( (ScanAz[inBeamIndex], ScanEl[inBeamIndex]), abs( GainYY[0]) / max(abs( GainYY[0])), (xi, yi), method='cubic')
plt.subplot( 2, 2, 1) # XX amp
plt.contourf(xi, yi, XXamp, 25, cmap=plt.cm.jet); plt.colorbar()
plt.plot( ScanAz[inBeamIndex], ScanEl[inBeamIndex], '.', color='k', alpha=0.1)
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.title(prefix + '-XX amp')
plt.subplot( 2, 2, 2) # XX amp
plt.contourf(xi, yi, XYamp, 25, cmap=plt.cm.jet); plt.colorbar()
plt.plot( ScanAz[inBeamIndex], ScanEl[inBeamIndex], '.', color='k', alpha=0.1)
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.title(prefix + '-XY amp')
plt.subplot( 2, 2, 3) # XX amp
plt.contourf(xi, yi, YXamp, 25, cmap=plt.cm.jet); plt.colorbar()
plt.subplot( 2, 2, 4) # XX amp
plt.plot( ScanAz[inBeamIndex], ScanEl[inBeamIndex], '.', color='k', alpha=0.1)
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.title(prefix + '-YX amp')
plt.contourf(xi, yi, YYamp, 25, cmap=plt.cm.jet); plt.colorbar()
plt.plot( ScanAz[inBeamIndex], ScanEl[inBeamIndex], '.', color='k', alpha=0.1)
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.title(prefix + '-YY amp')
plt.savefig(prefix + '_BeamAmp.pdf', form='pdf')
plt.close()

fig = plt.figure( figsize = (10,8))
XXphs = griddata( (ScanAz[inBeamIndex], ScanEl[inBeamIndex]), np.angle( GainXX[0]), (xi, yi), method='cubic')
XYphs = griddata( (ScanAz[inBeamIndex], ScanEl[inBeamIndex]), np.angle( GainXY[0]), (xi, yi), method='cubic')
YXphs = griddata( (ScanAz[inBeamIndex], ScanEl[inBeamIndex]), np.angle( GainYX[0]), (xi, yi), method='cubic')
YYphs = griddata( (ScanAz[inBeamIndex], ScanEl[inBeamIndex]), np.angle( GainYY[0]), (xi, yi), method='cubic')
plt.subplot( 2, 2, 1) # XX amp
plt.contourf(xi, yi, XXphs, 25, cmap=plt.cm.jet); plt.colorbar()
plt.plot( ScanAz[inBeamIndex], ScanEl[inBeamIndex], '.', color='k', alpha=0.1)
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.title(prefix + '-XX phs')
plt.subplot( 2, 2, 2) # XX amp
plt.contourf(xi, yi, XYphs, 25, cmap=plt.cm.jet); plt.colorbar()
plt.plot( ScanAz[inBeamIndex], ScanEl[inBeamIndex], '.', color='k', alpha=0.1)
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.title(prefix + '-XY phs')
plt.subplot( 2, 2, 3) # XX amp
plt.contourf(xi, yi, YXphs, 25, cmap=plt.cm.jet); plt.colorbar()
plt.subplot( 2, 2, 4) # XX phs
plt.plot( ScanAz[inBeamIndex], ScanEl[inBeamIndex], '.', color='k', alpha=0.1)
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.title(prefix + '-YX phs')
plt.contourf(xi, yi, YYphs, 25, cmap=plt.cm.jet); plt.colorbar()
plt.plot( ScanAz[inBeamIndex], ScanEl[inBeamIndex], '.', color='k', alpha=0.1)
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.title(prefix + '-YY phs')
plt.savefig(prefix + '_BeamPhs.pdf', form='pdf')
plt.close()

"""
