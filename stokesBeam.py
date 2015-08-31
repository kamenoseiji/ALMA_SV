execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import griddata
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#
#-------- For Plot
def circlePoints( x, y, radius ):
    angle = np.arange(-pi, (130/128)*pi, pi/128)
    return x + radius* np.cos(angle), y + radius* np.sin(angle)
#
#-------- BL 2 Ant mapping
ANT0 = []
ANT1 = []
for bl_index in range(2016):
    ants = Bl2Ant(bl_index)
    ANT0.append(ants[0])
    ANT1.append(ants[1])
#
ANT0 = np.array(ANT0)
ANT1 = np.array(ANT1)
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
CalQ, CalU, GYphs = solution[0], solution[1], solution[2]
#-------- Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
blNum  = antNum* (antNum - 1) / 2
AntD = np.zeros([antNum])
for ant_index in range(antNum):
    AntD[ant_index] = GetAntD(antList[ant_index])
#
print('Checking the Array ....')
#-------- Tracking and Scanning Antennas
trkAnt, scnAnt, scanTime, Offset = antRefScan(msfile)
trkAntNum = len(trkAnt)
scnAntNum = len(scnAnt)
if refantName not in antList[trkAnt]:
   print refantName + ' does not exist in the tracking antenna list this MS.'
#
#-------- Antenna and BL Grouping
refAntIndex = np.where( antList == refantName )[0].tolist()       # Ref ant in all AntList
trkAnt.pop( trkAnt.index(refAntIndex[0])); trkAnt = refAntIndex + trkAnt
AntIndex = trkAnt + scnAnt     # Antenna List, refants are prior
antWeight = np.ones(antNum)
antWeight[scnAnt] = 0.5
#-- BL mapping for all baselines
BlMap  = range(blNum)
BlInv = [False]* blNum      # True -> inverted baseline
blWeight = np.ones([blNum])
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    BlMap[bl_index], BlInv[bl_index] = Ant2BlD(AntIndex[ants[0]], AntIndex[ants[1]])
    blWeight[bl_index] = antWeight[AntIndex[ants[0]]]* antWeight[AntIndex[ants[1]]]
#
trkBlIndex  = np.where(blWeight == 1.0)[0].tolist(); trkBlNum  = len(trkBlIndex)        # Ref-Ref baselines
ScTrBlIndex = np.where(blWeight == 0.5)[0].tolist(); ScTrBlNum = len(ScTrBlIndex)       # Ref-Scan baselines
ScScBlIndex = np.where(blWeight == 0.25)[0].tolist(); ScScBlNum = len(ScScBlIndex)      # Scan-Scan baselines
#-- Load all-baseline visibilities
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], scan[0])
#-------- Bunch/Sel
chNum = Xspec.shape[1]
if chNum == 1:
    temp = Xspec[:,0, BlMap]
else:
    # BPcal, Bunch/Sel process comes here!
    temp = np.mean(Xspec[:,chRange], axis=1)
#
Ximag = temp.transpose(0,2,1).imag * (-2.0* np.array(BlInv) + 1.0) 
temp.imag = Ximag.transpose(0,2,1)
XX = temp[0]; XY = temp[1]; YX = temp[2]; YY = temp[3]
#-- BL within antenna groups
TrkXX  = XX[trkBlIndex]; TrkXY  = XY[trkBlIndex]; TrkYX  = YX[trkBlIndex]; TrkYY  = YY[trkBlIndex]  # Ref-Ref 
ScTrXX = XX[ScTrBlIndex]; ScTrXY = XY[ScTrBlIndex]; ScTrYX = YX[ScTrBlIndex]; ScTrYY = YY[ScTrBlIndex] # Ref-Sc 
ScScXX = XX[ScScBlIndex]; ScScXY = XY[ScScBlIndex]; ScScYX = YX[ScScBlIndex]; ScScYY = YY[ScScBlIndex] # Ref-Sc 
#-------- Frequency and Wavelength
chNum, chWid, Freq = GetChNum(msfile, spw[0])
FWHM = GetFWHM(msfile, spw[0], AntD)
if chNum > 1:
    chRange = range( int(0.05*chNum), int(0.95*chNum))
#
#-------- Center position of scanning antennas
timeNum = len(timeStamp)
AZ = np.zeros([timeNum]); EL = np.zeros([timeNum]); PA = np.zeros([timeNum])
ScanAz = np.zeros([timeNum]); ScanEl = np.zeros([timeNum])
scanTime, AntID, az, el = GetAzEl(msfile)
index = np.where( AntID == refAntIndex[0])[0].tolist()
for time_index in range(timeNum):
    ScanAz[time_index], ScanEl[time_index] = AzElMatch( timeStamp[time_index], scanTime[index], np.median(diff(timeStamp)), Offset)
#
azel = np.r_[az[index], el[index]].reshape(2, len(index))
for time_index in range(timeNum):
    AZ[time_index], EL[time_index] = AzElMatch( timeStamp[time_index], scanTime[index], np.median(diff(timeStamp)), azel)
    PA[time_index] = AzEl2PA(AZ[time_index], EL[time_index], ALMA_lat)    #
#
Ucos_Qsin = solution[1]* np.cos(2.0*PA) - solution[0]* np.sin(2.0*PA)   # U cos - Q sin
Qcos_Usin = solution[0]* np.cos(2.0*PA) + solution[1]* np.sin(2.0*PA)   # Q cos + U sin
#-- Gain solutions using parallel-hand correlations
GainX, GainY = polariGain(XX, YY, PA, solution[0], solution[1])
GainY = GainY* exp(solution[2]* 1.0j)                                   # XY-phase correction
#
#-- Gain Correction
VisXX = gainCalVis( XX, GainX, GainX )
VisXY = gainCalVis( XY, GainX, GainY )
VisYX = gainCalVis( YX, GainY, GainX )
VisYY = gainCalVis( YY, GainY, GainY )
#-------- RefAnt D-term
DxDx = VisXX[trkBlIndex] - 1.0 - Qcos_Usin
DxDy = VisXY[trkBlIndex] - Ucos_Qsin
DyDx = VisYX[trkBlIndex] - Ucos_Qsin
DyDy = VisYY[trkBlIndex] - 1.0 + Qcos_Usin
print('-------- Determining Antenna-based D-terms (refants) ----')
phsNum = trkAntNum - 1     # Number of phase solutions (excluding refant)
Pre = np.zeros([4*trkBlNum, 2*trkAntNum])  # [ReXX, ReXY, ReYX, ReYY] x [ReDx, ReDy]
Pim = np.zeros([4*trkBlNum, 2*(phsNum)])   # [ReXX, ReXY, ReYX, ReYY] x [ReDx, ReDy], no refant component
Dx = np.zeros([antNum, timeNum], dtype=complex)
Dy = np.zeros([antNum, timeNum], dtype=complex)
blGain = np.ones([trkBlNum])
for time_index in range(timeNum):
    for bl_index in range(trkBlNum):
        ants = Bl2Ant(bl_index)
        blGain[bl_index] = abs(GainX[trkAnt[ants[1]], time_index])* abs(GainX[trkAnt[ants[0]], time_index]) + abs(GainY[trkAnt[ants[1]], time_index])* abs(GainY[trkAnt[ants[0]], time_index])
    #
    W  = np.diag(np.r_[blGain, blGain, blGain, blGain])
    #---- Real Part
    Pre[           0:   trkBlNum ,         0:   trkAntNum ] = Ucos_Qsin[time_index]* BlAmpMatrix(trkAntNum)
    Pre[   trkBlNum :(2*trkBlNum),         0:   trkAntNum ] = (1.0 - Qcos_Usin[time_index])* DxMatrix(trkAntNum)
    Pre[   trkBlNum :(2*trkBlNum), trkAntNum:(2*trkAntNum)] = (1.0 + Qcos_Usin[time_index])* DyMatrix(trkAntNum)
    Pre[(2*trkBlNum):(3*trkBlNum),         0:   trkAntNum ] = (1.0 - Qcos_Usin[time_index])* DyMatrix(trkAntNum)
    Pre[(2*trkBlNum):(3*trkBlNum), trkAntNum:(2*trkAntNum)] = (1.0 + Qcos_Usin[time_index])* DxMatrix(trkAntNum)
    Pre[(3*trkBlNum):(4*trkBlNum), trkAntNum:(2*trkAntNum)] = Ucos_Qsin[time_index]* BlAmpMatrix(trkAntNum)
    #---- Imag Part, Dim of trkant is definately 0
    Pim[           0:   trkBlNum ,      0:   phsNum ] = -Ucos_Qsin[time_index]* BlPhaseMatrix(trkAntNum)
    Pim[   trkBlNum :(2*trkBlNum),      0:   phsNum ] = (1.0 - Qcos_Usin[time_index])* DxImMatrix(trkAntNum)
    Pim[   trkBlNum :(2*trkBlNum), phsNum:(2*phsNum)] = (1.0 - Qcos_Usin[time_index])* DyImMatrix(trkAntNum)
    Pim[(2*trkBlNum):(3*trkBlNum),      0:   phsNum ] = (1.0 + Qcos_Usin[time_index])* DyImMatrix(trkAntNum)
    Pim[(2*trkBlNum):(3*trkBlNum), phsNum:(2*phsNum)] = (1.0 + Qcos_Usin[time_index])* DxImMatrix(trkAntNum)
    Pim[(3*trkBlNum):(4*trkBlNum), phsNum:(2*phsNum)] = Ucos_Qsin[time_index]* BlPhaseMatrix(trkAntNum)
    #
    PreWPre = np.dot( Pre.T, np.dot(W, Pre))
    PimWPim = np.dot( Pim.T, np.dot(W, Pim))
    DetPre = np.linalg.det( PreWPre )
    DetPim = np.linalg.det( PimWPim )
    ReXXYY = np.dot(W, np.r_[ DxDx[:,time_index].real, DxDy[:,time_index].real, DyDx[:,time_index].real, DyDy[:,time_index].real ])
    ImXXYY = np.dot(W, np.r_[ DxDx[:,time_index].imag, DxDy[:,time_index].imag, DyDx[:,time_index].imag, DyDy[:,time_index].imag ])
    # print( `time_index` + '  det(PTP_real) = ' + `DetPre` + '  det(PTP_imag) = ' + `DetPim` )
    ReD = np.dot( np.linalg.inv( PreWPre ), np.dot( Pre.T, ReXXYY) )
    ImD = np.dot( np.linalg.inv( PimWPim ), np.dot( Pim.T, ImXXYY) )
    Dx[trkAnt, time_index] = ReD[0:trkAntNum]             + 1.0j* np.r_[0, ImD[0:phsNum]]
    Dy[trkAnt, time_index] = ReD[trkAntNum:(2*trkAntNum)] + 1.0j* np.r_[0, ImD[phsNum:(2*phsNum)]]
#
TrkDx = np.mean(Dx[trkAnt], axis=1)
TrkDy = np.mean(Dy[trkAnt], axis=1)
#-------- Stokes Parameters measured by tracking antennas
StokesVis = np.zeros([trkBlNum, timeNum, 4], dtype=complex) 
for time_index in range(timeNum):
    Pinv = InvPAMatrix( PA[time_index] )
    for bl_index in range(trkBlNum):
        ants = Bl2Ant(bl_index)
        Minv = InvMullerMatrix( TrkDx[ants[1]], TrkDy[ants[1]], TrkDx[ants[0]], TrkDy[ants[0]])
        StokesVis[bl_index, time_index] = np.dot(Minv, np.array( [VisXX[trkBlIndex[bl_index], time_index], VisXY[trkBlIndex[bl_index], time_index], VisYX[trkBlIndex[bl_index], time_index], VisYY[trkBlIndex[bl_index], time_index]]))
        StokesVis[bl_index, time_index] = np.dot(Pinv, StokesVis[bl_index, time_index])
    #
#
TrkI = np.mean( StokesVis[:,:,0], axis=(0,1) ).real
TrkQ = np.mean( StokesVis[:,:,1], axis=(0,1) ).real
TrkU = np.mean( StokesVis[:,:,2], axis=(0,1) ).real
TrkV = np.mean( StokesVis[:,:,3], axis=(0,1) ).real
text_sd = 'TrkMeas / Model: I=%6.4f / %6.4f  Q=%6.4f / %6.4f U=%6.4f / %6.4f V=%6.4f / %6.4f' % (TrkI, 1.0, TrkQ, CalQ, TrkU, CalU, TrkV, 0.0)
print text_sd
#-------- Determination of D-terms in scanning antennas
print('-------- Determining Antenna-based D-terms (scan ants) ----')
for ant_index in range(scnAntNum):
    antID = trkAntNum + ant_index
    blWithScnAnt = np.array(ScTrBlIndex)[np.where( ANT0[ScTrBlIndex] == antID )[0].tolist() + np.where( ANT1[ScTrBlIndex] == antID )[0].tolist()].tolist()
    trkAnt_index = range(trkAntNum)
    Dx[scnAnt[ant_index]] = (np.mean( VisYX[blWithScnAnt] - Ucos_Qsin - np.outer(TrkDy, (1.0 + Qcos_Usin)), axis=0 ) / (1.0 - Qcos_Usin)).conjugate()
    Dy[scnAnt[ant_index]] = (np.mean( VisXY[blWithScnAnt] - Ucos_Qsin + np.outer(TrkDx, (1.0 - Qcos_Usin)), axis=0 ) / (1.0 + Qcos_Usin)).conjugate()
#
#-------- Determination of D-terms in scanning antennas
print('-------- Plot D-term Maps for scan ants ----')
for ant_index in range(scnAntNum):
    antID = scnAnt[ant_index]
    #-------- Plot
    fig = plt.figure( figsize = (10,10))
    fig.suptitle(prefix + ' ' + antList[antID] + ' SPW=' + `spw[0]` + ' Scan=' + `scan[0]`)
    fig.text(0.45, 0.05, 'Az Offset [arcsec]')
    fig.text(0.05, 0.45, 'El Offset [arcsec]', rotation=90)
    #
    xi, yi = np.mgrid[ -floor(2.0*FWHM[antID]):floor(2.0*FWHM[antID]):128j, -floor(2.0*FWHM[antID]):floor(2.0*FWHM[antID]):128j]
    IndexCenter = np.where( ScanAz**2 + ScanEl**2 < 0.005* FWHM[antID]**2 )[0]
    Index3dB = np.where( ScanAz**2 + ScanEl**2 < 0.25* FWHM[antID]**2 )[0]
    Index6dB = np.where( ScanAz**2 + ScanEl**2 < 0.5* FWHM[antID]**2 )[0]
    ReDxmap = GridData( Dx[antID].real, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
    ImDxmap = GridData( Dx[antID].imag, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
    ReDymap = GridData( Dy[antID].real, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
    ImDymap = GridData( Dy[antID].imag, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
    #---- plot Re(Dx)
    plt.subplot( 2, 2, 1, aspect=1); plt.contourf(xi, yi, ReDxmap, np.linspace(-0.1, 0.1, 21)); plt.colorbar(); plt.title('Re(Dx)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Re(Dx) at Center = %5.3f' % ( np.mean(Dx[antID, IndexCenter].real) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(Dx[antID, Index3dB].real), np.min(Dx[antID, Index3dB].real) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(Dx[antID, Index6dB].real), np.min(Dx[antID, Index6dB].real) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    #---- plot Im(Dx)
    plt.subplot( 2, 2, 2, aspect=1); plt.contourf(xi, yi, ImDxmap, np.linspace(-0.1, 0.1, 21)); plt.colorbar(); plt.title('Im(Dx)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Im(Dx) at Center = %5.3f' % ( np.mean(Dx[antID, IndexCenter].imag) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(Dx[antID, Index3dB].imag), np.min(Dx[antID, Index3dB].imag) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(Dx[antID, Index6dB].imag), np.min(Dx[antID, Index6dB].imag) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    #---- plot Re(Dy)
    plt.subplot( 2, 2, 3, aspect=1); plt.contourf(xi, yi, ReDymap, np.linspace(-0.1, 0.1, 21)); plt.colorbar(); plt.title('Re(Dy)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Re(Dy) at Center = %5.3f' % ( np.mean(Dy[antID, IndexCenter].real) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(Dy[antID, Index3dB].real), np.min(Dy[antID, Index3dB].real) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(Dy[antID, Index6dB].real), np.min(Dy[antID, Index6dB].real) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    #---- plot Im(Dy)
    plt.subplot( 2, 2, 4, aspect=1); plt.contourf(xi, yi, ImDymap, np.linspace(-0.1, 0.1, 21)); plt.colorbar(); plt.title('Im(Dy)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Im(Dy) at Center = %5.3f' % ( np.mean(Dy[antID, IndexCenter].imag) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(Dy[antID, Index3dB].imag), np.min(Dy[antID, Index3dB].imag) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(Dy[antID, Index6dB].imag), np.min(Dy[ant_index, Index6dB].imag) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    plt.plot( ScanAz, ScanEl, '.', color='k', alpha=0.1)
    plt.axis([-2.0*FWHM[antID], 2.0*FWHM[antID], -2.0*FWHM[antID], 2.0*FWHM[antID]])
    plt.savefig( prefix + '-' + antList[antID] + '-SPW' + `spw[0]` + '-DtermMap.pdf', form='pdf'); plt.close()
    plt.close()
#
#-------- D-term-corrected Stokes parameters --------
print('-------- D-term-corrected Stokes parameters ----')
StokesVis = np.zeros([ScScBlNum, timeNum, 4], dtype=complex) 
for time_index in range(timeNum):
    Pinv = InvPAMatrix( PA[time_index] )
    for bl_index in range(ScScBlNum):
        BlID = ScScBlIndex[bl_index]
        ants = Bl2Ant(BlID)
        Ant0 = scnAnt[ants[0] - trkAntNum]
        Ant1 = scnAnt[ants[1] - trkAntNum]
        Minv = InvMullerMatrix( Dx[Ant1, time_index], Dy[Ant1, time_index], Dx[Ant0, time_index], Dy[Ant0, time_index])
        StokesVis[bl_index, time_index] = np.dot(Minv, np.array( [VisXX[BlID, time_index], VisXY[BlID, time_index], VisYX[BlID, time_index], VisYY[BlID, time_index]]))
        StokesVis[bl_index, time_index] = np.dot(Pinv, StokesVis[bl_index, time_index])
    #
#
ScnI = np.mean( StokesVis[:,:,0], axis=0 ).real
ScnQ = np.mean( StokesVis[:,:,1], axis=0 ).real
ScnU = np.mean( StokesVis[:,:,2], axis=0 ).real
ScnV = np.mean( StokesVis[:,:,3], axis=0 ).real
Qerr = ScnQ - TrkQ
Uerr = ScnU - TrkU
Perr = sqrt(ScnQ**2 + ScnU**2) - sqrt(TrkQ**2 + TrkU**2)
Aerr = np.arccos( (ScnQ* TrkQ + ScnU* TrkU) / sqrt( (ScnQ**2 + ScnU**2)* (TrkQ**2 + TrkU**2)))*90.0/pi
#-------- Plot Stokes Beam Map
fig = plt.figure( figsize = (10,10))
FWHM = np.max(FWHM)
xi, yi = np.mgrid[ -floor(FWHM):floor(FWHM):128j, -floor(FWHM):floor(FWHM):128j]
Imap = GridData( ScnI, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Qmap = GridData( ScnQ, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Umap = GridData( ScnU, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Vmap = GridData( ScnV, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
QEmap = GridData(Qerr, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
UEmap = GridData(Uerr, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
PEmap = GridData(Perr, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
AEmap = GridData(Aerr, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
IndexCenter = np.where( ScanAz**2 + ScanEl**2 < 0.005* FWHM**2 )[0]
Index3dB = np.where( ScanAz**2 + ScanEl**2 < 0.25* FWHM**2 )[0]
Index6dB = np.where( ScanAz**2 + ScanEl**2 < 0.5* FWHM**2 )[0]
#---- Plot I map
plt.subplot(2, 2, 1, aspect=1); plt.contourf(xi, yi, Imap, np.linspace(0.9, 1.1, 21)); plt.colorbar(); plt.title(prefix + '-StokesI')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'I at Center   = %5.3f' % ( np.mean(ScnI[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnI[Index3dB]), np.min(ScnI[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnI[Index6dB]), np.min(ScnI[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot Q map
plt.subplot(2, 2, 2, aspect=1); plt.contourf(xi, yi, Qmap, np.linspace(0.01*(floor(TrkQ* 100)-5), 0.01*(floor(TrkQ* 100)+5), 21)); plt.colorbar(); plt.title(prefix + '-StokesQ')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Q at Center   = %5.3f' % ( np.mean(ScnQ[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnQ[Index3dB]), np.min(ScnQ[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnQ[Index6dB]), np.min(ScnQ[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot U map
plt.subplot(2, 2, 3, aspect=1); plt.contourf(xi, yi, Umap, np.linspace(0.01*(floor(TrkU* 100)-5), 0.01*(floor(TrkU* 100)+5), 21)); plt.colorbar(); plt.title(prefix + '-StokesU')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'U at Center   = %5.3f' % ( np.mean(ScnU[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnU[Index3dB]), np.min(ScnU[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnU[Index6dB]), np.min(ScnU[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot V map
plt.subplot(2, 2, 4, aspect=1); plt.contourf(xi, yi, Vmap, np.linspace(-0.05, 0.05, 21)); plt.colorbar(); plt.title(prefix + '-StokesV')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.plot( ScanAz, ScanEl, '.', color='k', alpha=0.1)
text_sd = 'V at Center   = %5.3f' % ( np.mean(ScnV[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnV[Index3dB]), np.min(ScnV[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnV[Index6dB]), np.min(ScnV[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.savefig( prefix + '-StokesMap.pdf', form='pdf'); plt.close()
plt.close()
#
#---- Plot Q err map
plt.subplot(2, 2, 1, aspect=1); plt.contourf(xi, yi, 100.0*QEmap, np.linspace(-2.0, 2.0, 41)); plt.colorbar(); plt.title(prefix + ' Q_err %')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Qerr at Center       = %4.1f %%' % ( 100.0* np.mean(Qerr[IndexCenter])); plt.text(-0.8*FWHM, -0.7*FWHM, text_sd, size='x-small')
text_sd = 'Max Qerr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(Qerr[Index3dB])) ); plt.text(-0.8*FWHM, -0.8*FWHM, text_sd, size='x-small')
text_sd = 'Max Qerr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(Qerr[Index6dB])) ); plt.text(-0.8*FWHM, -0.9*FWHM, text_sd, size='x-small')
#----Plot U err map
plt.subplot(2, 2, 2, aspect=1); plt.contourf(xi, yi, 100.0*UEmap, np.linspace(-2.0, 2.0, 41)); plt.colorbar(); plt.title(prefix + ' U_err %')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Uerr at Center       = %4.1f %%' % ( 100.0* np.mean(Uerr[IndexCenter])); plt.text(-0.8*FWHM, -0.7*FWHM, text_sd, size='x-small')
text_sd = 'Max Uerr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(Uerr[Index3dB])) ); plt.text(-0.8*FWHM, -0.8*FWHM, text_sd, size='x-small')
text_sd = 'Max Uerr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(Uerr[Index6dB])) ); plt.text(-0.8*FWHM, -0.9*FWHM, text_sd, size='x-small')
#----Plot P error
plt.subplot( 2, 2, 3, aspect=1); plt.contourf(xi, yi, 100.0*PEmap, np.linspace(-2.0, 2.0, 41)); plt.colorbar(); plt.title(prefix + ' P_err %')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Perr at Center       = %4.1f %%' % ( 100.0* np.mean(Perr[IndexCenter])); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = 'Max Perr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(Perr[Index3dB])) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = 'Max Perr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(Perr[Index6dB])) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#----Plot EVPA error
plt.subplot( 2, 2, 4, aspect=1); plt.contourf(xi, yi, AEmap, np.linspace(-5, 5, 21)); plt.colorbar(); plt.title(prefix + ' EVPA_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'EVPAerr at Center       = %4.1f deg' % ( np.mean(Aerr[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = 'Max EVPAerr (-3dB beam) = %4.1f deg' % ( max(abs(Aerr[Index3dB])) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = 'Max EVPAerr (-6dB beam) = %4.1f deg' % ( max(abs(Aerr[Index6dB])) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.savefig( prefix + '-StokesErr.pdf', form='pdf'); plt.close()
plt.close()
#--------- Gridding for 11x11 sampling points



"""
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
AZ = np.zeros([timeNum])
EL = np.zeros([timeNum])
PA = np.zeros([timeNum])
scanTime, AntID, az, el = GetAzEl(msfile)
index = np.where( AntID == refAnt[0])[0].tolist()
azel = np.r_[az[index], el[index]].reshape(2, len(index))
for time_index in range(timeNum):
    AZ[time_index], EL[time_index] = AzElMatch( timeStamp[time_index], scanTime[index], np.median(interval), azel)
    PA[time_index] = AzEl2PA(AZ[time_index], EL[time_index], ALMA_lat)    #
#
#
#-------- Load Visibilities
BlMap  = range(blNum)
BlInv = [False]* blNum      # True -> inverted baseline
blWeight = np.ones([blNum])
for bl_index in range(blNum):
    ants = Bl2Ant(bl_index)
    BlMap[bl_index], BlInv[bl_index] = Ant2BlD(ants[0], ants[1])
    blWeight[bl_index] = antWeight[ants[0]]* antWeight[ants[1]]
#
refBlIndex  = np.where(blWeight == 1.0)[0].tolist(); refBlNum  = len(refBlIndex)        # Ref-Ref baselines
ScRfBlIndex = np.where(blWeight == 0.5)[0].tolist(); ScRfBlNum = len(ScRfBlIndex)       # Ref-Scan baselines
ScScBlIndex = np.where(blWeight == 0.25)[0].tolist(); ScScBlNum = len(ScScBlIndex)      # Scan-Scan baselines
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], scan[0])
#-------- Channel-Averaged Visibilities
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
        XY[bl_index] = YX[bl_index].conjugate()     # Baseline inversion -> swap(X, Y)
        YX[bl_index] = XY[bl_index].conjugate()     # Baseline inversion -> swap(X, Y)
        YY[bl_index] = YY[bl_index].conjugate()
    #
#
GainX, GainY = polariGain(XX, YY, PA, solution[0], solution[1])
Ucos_Qsin = solution[1]* np.cos(2.0*PA) - solution[0]* np.sin(2.0*PA)
Qcos_Usin = solution[0]* np.cos(2.0*PA) + solution[1]* np.sin(2.0*PA)
#
#-------- RefAnt D-term
RefGainX = GainX[refAnt]
RefGainY = GainY[refAnt]
blGain = np.ones([refBlNum])
VisXX = gainCalVis( XX[refBlIndex], RefGainX, RefGainX )
VisYY = gainCalVis( YY[refBlIndex], RefGainY, RefGainY )
VisXY = gainCalVis( XY[refBlIndex], RefGainX, RefGainY )* np.exp(-solution[2]*1.0j)
VisYX = gainCalVis( YX[refBlIndex], RefGainY, RefGainX )* np.exp( solution[2]*1.0j)
DxDx = VisXX - 1.0 - Qcos_Usin
DxDy = VisXY - Ucos_Qsin
DyDx = VisYX - Ucos_Qsin
DyDy = VisYY - 1.0 + Qcos_Usin
print('-------- Determining Antenna-based D-terms (refants) ----')
phsNum = refAntNum - 1     # Number of phase solutions (excluding refant)
Pre = np.zeros([4*refBlNum, 2*refAntNum])  # [ReXX, ReXY, ReYX, ReYY] x [ReDx, ReDy]
Pim = np.zeros([4*refBlNum, 2*(phsNum)])   # [ReXX, ReXY, ReYX, ReYY] x [ReDx, ReDy], no refant component
refDx = np.zeros([refAntNum, timeNum], dtype=complex)
refDy = np.zeros([refAntNum, timeNum], dtype=complex)
for time_index in range(timeNum):
    for bl_index in range(refBlNum):
        ants = Bl2Ant(bl_index)
        blGain[bl_index] = abs(RefGainX[ants[1], time_index])* abs(RefGainX[ants[0], time_index]) + abs(RefGainY[ants[1], time_index])* abs(RefGainY[ants[0], time_index])
    #
    W  = np.diag(np.r_[blGain, blGain, blGain, blGain])
    #---- Real Part
    Pre[           0:   refBlNum ,         0:   refAntNum ] = Ucos_Qsin[time_index]* BlAmpMatrix(refAntNum)
    Pre[   refBlNum :(2*refBlNum),         0:   refAntNum ] = (1.0 - Qcos_Usin[time_index])* DxMatrix(refAntNum)
    Pre[   refBlNum :(2*refBlNum), refAntNum:(2*refAntNum)] = (1.0 + Qcos_Usin[time_index])* DyMatrix(refAntNum)
    Pre[(2*refBlNum):(3*refBlNum),         0:   refAntNum ] = (1.0 - Qcos_Usin[time_index])* DyMatrix(refAntNum)
    Pre[(2*refBlNum):(3*refBlNum), refAntNum:(2*refAntNum)] = (1.0 + Qcos_Usin[time_index])* DxMatrix(refAntNum)
    Pre[(3*refBlNum):(4*refBlNum), refAntNum:(2*refAntNum)] = Ucos_Qsin[time_index]* BlAmpMatrix(refAntNum)
    #---- Imag Part, Dim of refant is definately 0
    Pim[           0:   refBlNum ,      0:   phsNum ] = -Ucos_Qsin[time_index]* BlPhaseMatrix(refAntNum)
    Pim[   refBlNum :(2*refBlNum),      0:   phsNum ] = (1.0 - Qcos_Usin[time_index])* DxImMatrix(refAntNum)
    Pim[   refBlNum :(2*refBlNum), phsNum:(2*phsNum)] = (1.0 - Qcos_Usin[time_index])* DyImMatrix(refAntNum)
    Pim[(2*refBlNum):(3*refBlNum),      0:   phsNum ] = (1.0 + Qcos_Usin[time_index])* DyImMatrix(refAntNum)
    Pim[(2*refBlNum):(3*refBlNum), phsNum:(2*phsNum)] = (1.0 + Qcos_Usin[time_index])* DxImMatrix(refAntNum)
    Pim[(3*refBlNum):(4*refBlNum), phsNum:(2*phsNum)] = Ucos_Qsin[time_index]* BlPhaseMatrix(refAntNum)
    #
    PreWPre = np.dot( Pre.T, np.dot(W, Pre))
    PimWPim = np.dot( Pim.T, np.dot(W, Pim))
    DetPre = np.linalg.det( PreWPre )
    DetPim = np.linalg.det( PimWPim )
    ReXXYY = np.dot(W, np.r_[ DxDx[:,time_index].real, DxDy[:,time_index].real, DyDx[:,time_index].real, DyDy[:,time_index].real ])
    ImXXYY = np.dot(W, np.r_[ DxDx[:,time_index].imag, DxDy[:,time_index].imag, DyDx[:,time_index].imag, DyDy[:,time_index].imag ])
    # print( `time_index` + '  det(PTP_real) = ' + `DetPre` + '  det(PTP_imag) = ' + `DetPim` )
    ReD = np.dot( np.linalg.inv( PreWPre ), np.dot( Pre.T, ReXXYY) )
    ImD = np.dot( np.linalg.inv( PimWPim ), np.dot( Pim.T, ImXXYY) )
    refDx[:, time_index] = ReD[0:refAntNum]             + 1.0j* np.r_[0, ImD[0:phsNum]]
    refDy[:, time_index] = ReD[refAntNum:(2*refAntNum)] + 1.0j* np.r_[0, ImD[phsNum:(2*phsNum)]]
#
RefDx = np.mean(refDx, axis=1)
RefDy = np.mean(refDy, axis=1)
StokesVis = np.zeros([refBlNum, timeNum, 4], dtype=complex) 
for time_index in range(timeNum):
    Pinv = InvPAMatrix( PA[time_index] )
    for bl_index in range(refBlNum):
        ants = Bl2Ant(bl_index)
        Minv = InvMullerMatrix( RefDx[ants[1]], RefDy[ants[1]], RefDx[ants[0]], RefDy[ants[0]])
        StokesVis[bl_index, time_index] = np.dot(Minv, np.array( [VisXX[bl_index, time_index], VisXY[bl_index, time_index], VisYX[bl_index, time_index], VisYY[bl_index, time_index]]))
        StokesVis[bl_index, time_index] = np.dot(Pinv, StokesVis[bl_index, time_index])
    #
#
RefI = np.mean( StokesVis[:,:,0], axis=(0,1) ).real
RefQ = np.mean( StokesVis[:,:,1], axis=(0,1) ).real
RefU = np.mean( StokesVis[:,:,2], axis=(0,1) ).real
RefV = np.mean( StokesVis[:,:,3], axis=(0,1) ).real
text_sd = 'RefMeas / Model: I=%6.4f / %6.4f  Q=%6.4f / %6.4f U=%6.4f / %6.4f V=%6.4f / %6.4f' % (RefI, 1.0, RefQ, CalQ, RefU, CalU, RefV, 0.0)
print text_sd
#-------- Determination of D-terms in scanning antennas
#
VisXX = np.zeros([refAntNum, timeNum], dtype=complex)
VisXY = np.zeros([refAntNum, timeNum], dtype=complex)
VisYX = np.zeros([refAntNum, timeNum], dtype=complex)
VisYY = np.zeros([refAntNum, timeNum], dtype=complex)
scnDx = np.zeros([scanAntNum, timeNum], dtype=complex)
scnDy = np.zeros([scanAntNum, timeNum], dtype=complex)
print('-------- Determining Antenna-based D-terms (scan ants) ----')
for ant_index in range(scanAntNum):
    antID = scanAnt[ant_index]
    blWithScanAnt = np.array(ScRfBlIndex)[np.where( ANT0[ScRfBlIndex] == antID )[0].tolist() + np.where( ANT1[ScRfBlIndex] == antID )[0].tolist()].tolist()
    for bl_index in range(refAntNum):
        blID = blWithScanAnt[bl_index]
        ants = Bl2Ant(blID)
        if ants[0] == antID:        # Baseline Inversed
            VisXX[bl_index] = XX[blID].conjugate() / (GainX[ants[1]]* GainX[antID].conjugate())
            VisXY[bl_index] = YX[blID].conjugate() / (GainY[ants[1]]* GainX[antID].conjugate()) # Baseline inversion -> swap(X,Y)
            VisYX[bl_index] = XY[blID].conjugate() / (GainX[ants[1]]* GainY[antID].conjugate()) # Baseline inversion -> swap(X,Y)
            VisYY[bl_index] = YY[blID].conjugate() / (GainY[ants[1]]* GainY[antID].conjugate())
        else:
            VisXX[bl_index] = XX[blID] / (GainX[ants[0]]* GainX[antID].conjugate())
            VisXY[bl_index] = XY[blID] / (GainY[ants[0]]* GainX[antID].conjugate())
            VisYX[bl_index] = YX[blID] / (GainX[ants[0]]* GainY[antID].conjugate())
            VisYY[bl_index] = YY[blID] / (GainY[ants[0]]* GainY[antID].conjugate())
        #
        VisXX[bl_index] = 2.0* VisXX[bl_index] - 1.0 - Qcos_Usin + RefDx[bl_index]* Ucos_Qsin
        VisXY[bl_index] = 2.0* VisXY[bl_index] - Ucos_Qsin + RefDx[bl_index]* (1.0 - Qcos_Usin)
        VisYX[bl_index] = 2.0* VisYX[bl_index] - Ucos_Qsin - RefDy[bl_index]* (1.0 + Qcos_Usin)
        VisYY[bl_index] = 2.0* VisYY[bl_index] - 1.0 + Qcos_Usin - RefDy[bl_index]* Ucos_Qsin
        #
    #
    scnDx[ant_index] = (np.mean( VisYX, axis=0 ) / (1.0 - Qcos_Usin)).conjugate()
    scnDy[ant_index] = (np.mean( VisXY, axis=0 ) / (1.0 + Qcos_Usin)).conjugate()
#
print('-------- Plot D-term Maps for scan ants ----')
for ant_index in range(scanAntNum):
    antID = scanAnt[ant_index]
    #-------- Plot
    fig = plt.figure( figsize = (10,10))
    fig.suptitle(prefix + ' ' + antList[antID] + ' SPW=' + `spw[0]` + ' Scan=' + `scan[0]`)
    fig.text(0.45, 0.05, 'Az Offset [arcsec]')
    fig.text(0.05, 0.45, 'El Offset [arcsec]', rotation=90)
    #
    xi, yi = np.mgrid[ -floor(2.0*FWHM[antID]):floor(2.0*FWHM[antID]):128j, -floor(2.0*FWHM[antID]):floor(2.0*FWHM[antID]):128j]
    IndexCenter = np.where( ScanAz**2 + ScanEl**2 < 0.005* FWHM[antID]**2 )[0]
    Index3dB = np.where( ScanAz**2 + ScanEl**2 < 0.25* FWHM[antID]**2 )[0]
    Index6dB = np.where( ScanAz**2 + ScanEl**2 < 0.5* FWHM[antID]**2 )[0]
    ReDxmap = GridData( scnDx[ant_index].real, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
    ImDxmap = GridData( scnDx[ant_index].imag, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
    ReDymap = GridData( scnDy[ant_index].real, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
    ImDymap = GridData( scnDy[ant_index].imag, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
    #---- plot Re(Dx)
    plt.subplot( 2, 2, 1, aspect=1); plt.contourf(xi, yi, ReDxmap, np.linspace(-0.1, 0.1, 21)); plt.colorbar(); plt.title('Re(Dx)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Re(Dx) at Center = %5.3f' % ( np.mean(scnDx[ant_index, IndexCenter].real) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(scnDx[ant_index, Index3dB].real), np.min(scnDx[ant_index, Index3dB].real) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(scnDx[ant_index, Index6dB].real), np.min(scnDx[ant_index, Index6dB].real) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    #---- plot Im(Dx)
    plt.subplot( 2, 2, 2, aspect=1); plt.contourf(xi, yi, ImDxmap, np.linspace(-0.1, 0.1, 21)); plt.colorbar(); plt.title('Im(Dx)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Im(Dx) at Center = %5.3f' % ( np.mean(scnDx[ant_index, IndexCenter].imag) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(scnDx[ant_index, Index3dB].imag), np.min(scnDx[ant_index, Index3dB].imag) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(scnDx[ant_index, Index6dB].imag), np.min(scnDx[ant_index, Index6dB].imag) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    #---- plot Re(Dy)
    plt.subplot( 2, 2, 3, aspect=1); plt.contourf(xi, yi, ReDymap, np.linspace(-0.1, 0.1, 21)); plt.colorbar(); plt.title('Re(Dy)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Re(Dy) at Center = %5.3f' % ( np.mean(scnDy[ant_index, IndexCenter].real) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(scnDy[ant_index, Index3dB].real), np.min(scnDy[ant_index, Index3dB].real) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(scnDy[ant_index, Index6dB].real), np.min(scnDy[ant_index, Index6dB].real) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    #---- plot Im(Dy)
    plt.subplot( 2, 2, 4, aspect=1); plt.contourf(xi, yi, ImDymap, np.linspace(-0.1, 0.1, 21)); plt.colorbar(); plt.title('Im(Dy)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Im(Dy) at Center = %5.3f' % ( np.mean(scnDy[ant_index, IndexCenter].imag) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(scnDy[ant_index, Index3dB].imag), np.min(scnDy[ant_index, Index3dB].imag) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(scnDy[ant_index, Index6dB].imag), np.min(scnDy[ant_index, Index6dB].imag) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    plt.plot( ScanAz, ScanEl, '.', color='k', alpha=0.1)
    plt.axis([-2.0*FWHM[antID], 2.0*FWHM[antID], -2.0*FWHM[antID], 2.0*FWHM[antID]])
    plt.savefig( prefix + '-' + antList[antID] + '-SPW' + `spw[0]` + '-DtermMap.pdf', form='pdf'); plt.close()
    plt.close()
#
#
print('-------- D-term-corrected Stokes parameters ----')
ScGainX = GainX[scanAnt]
ScGainY = GainY[scanAnt]
ScScXX = gainCalVis( XX[ScScBlIndex], ScGainX, ScGainX )
ScScYY = gainCalVis( YY[ScScBlIndex], ScGainY, ScGainY )
ScScXY = gainCalVis( XY[ScScBlIndex], ScGainX, ScGainY ) * np.exp(-solution[2]*1.0j)
ScScYX = gainCalVis( YX[ScScBlIndex], ScGainY, ScGainX ) * np.exp( solution[2]*1.0j)
StokesVis = np.zeros([ScScBlNum, timeNum, 4], dtype=complex) 
for time_index in range(timeNum):
    Pinv = InvPAMatrix( PA[time_index] )
    for bl_index in range(ScScBlNum):
        ants = Bl2Ant(bl_index)
        Minv = InvMullerMatrix( scnDx[ants[1], time_index], scnDy[ants[1], time_index], scnDx[ants[0], time_index], scnDy[ants[0], time_index])
        StokesVis[bl_index, time_index] = np.dot(Minv, np.array( [ScScXX[bl_index, time_index], ScScXY[bl_index, time_index], ScScYX[bl_index, time_index], ScScYY[bl_index, time_index]]))
        StokesVis[bl_index, time_index] = np.dot(Pinv, StokesVis[bl_index, time_index])
    #
#
ScnI = np.mean( StokesVis[:,:,0], axis=0 ).real
ScnQ = np.mean( StokesVis[:,:,1], axis=0 ).real
ScnU = np.mean( StokesVis[:,:,2], axis=0 ).real
ScnV = np.mean( StokesVis[:,:,3], axis=0 ).real
#-------- Plot
fig = plt.figure( figsize = (10,8))
FWHM = np.max(FWHM)
xi, yi = np.mgrid[ -floor(FWHM):floor(FWHM):128j, -floor(FWHM):floor(FWHM):128j]
Imap = GridData( ScnI, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Qmap = GridData( ScnQ, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Umap = GridData( ScnU, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Vmap = GridData( ScnV, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
IndexCenter = np.where( ScanAz**2 + ScanEl**2 < 0.005* FWHM**2 )[0]
Index3dB = np.where( ScanAz**2 + ScanEl**2 < 0.25* FWHM**2 )[0]
Index6dB = np.where( ScanAz**2 + ScanEl**2 < 0.5* FWHM**2 )[0]
#---- Plot I map
plt.subplot( 2, 2, 1); plt.contourf(xi, yi, Imap, np.linspace(0.9, 1.1, 21), col='b'); plt.colorbar(); plt.title(prefix + '-StokesI')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'I at Center   = %5.3f' % ( np.mean(ScnI[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnI[Index3dB]), np.min(ScnI[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnI[Index6dB]), np.min(ScnI[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot Q map
plt.subplot( 2, 2, 2); plt.contourf(xi, yi, Qmap, np.linspace(0.01*(floor(RefQ* 100)+5), 0.01*(floor(RefQ* 100)-5), 21), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesQ')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Q at Center   = %5.3f' % ( np.mean(ScnQ[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnQ[Index3dB]), np.min(ScnQ[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnQ[Index6dB]), np.min(ScnQ[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot U map
plt.subplot( 2, 2, 3); plt.contourf(xi, yi, Umap, np.linspace(0.01*(floor(RefU* 100)+5), 0.01*(floor(RefU* 100)-5), 21), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesU')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'U at Center   = %5.3f' % ( np.mean(ScnU[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnU[Index3dB]), np.min(ScnU[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnU[Index6dB]), np.min(ScnU[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot V map
plt.subplot( 2, 2, 4); plt.contourf(xi, yi, Vmap, np.linspace(-0.05, 0.05, 21), cmap=plt.cm.jet); plt.colorbar(); plt.title(prefix + '-StokesV')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.plot( ScanAz, ScanEl, '.', color='k', alpha=0.1)
text_sd = 'V at Center   = %5.3f' % ( np.mean(ScnV[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnV[Index3dB]), np.min(ScnV[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnV[Index6dB]), np.min(ScnV[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.savefig( prefix + '-StokesMap.pdf', form='pdf'); plt.close()
plt.close()
#
#--------- Gridding for 11x11 sampling points
az, el = readGrid(gridFile)
logfile = open(prefix + '_StokesGrid.log', 'w')
text_sd = 'No. dAz      dEl        I        Q        U         V'
logfile.write(text_sd + '\n')
xi, yi = np.mgrid[ min(az):max(az):11j, max(el):min(el):11j]
Imap = GridData( ScnI, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Qmap = GridData( ScnQ, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Umap = GridData( ScnU, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Vmap = GridData( ScnV, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
for x_index in range(11):
    for y_index in range(11):
        text_sd = '%d %f %f %f %f %f %f' % (x_index*11 + y_index + 1, xi[x_index, y_index], yi[x_index, y_index], Imap[x_index, y_index], Qmap[x_index, y_index], Umap[x_index, y_index], Vmap[x_index, y_index])
        logfile.write(text_sd + '\n')
    #
#
logfile.close()
#-------- Plot within 2xFWHM
timeIndex = np.where( ScanAz**2 + ScanEl**2 < FWHM**2 )[0]
IndexCenter = np.where( ScanAz**2 + ScanEl**2 < 0.005* FWHM**2 )[0]
Index3dB = np.where( ScanAz**2 + ScanEl**2 < 0.25* FWHM**2 )[0]
Index6dB = np.where( ScanAz**2 + ScanEl**2 < 0.5* FWHM**2 )[0]
xi, yi = np.mgrid[ -floor(FWHM):floor(FWHM):128j, -floor(FWHM):floor(FWHM):128j]
fig = plt.figure( figsize = (10,8))
QerrScn = ScnQ - CalQ
UerrScn = ScnU - CalU
PerrScn = sqrt(ScnQ**2 + ScnU**2) - sqrt(CalQ**2 + CalU**2)
AerrScn = np.arccos( (ScnQ* CalQ + ScnU* CalU) / sqrt( (ScnQ**2 + ScnU**2)* (CalQ**2 + CalU**2)))*90.0/pi
#
Qerr = GridData( QerrScn[timeIndex], ScanAz[timeIndex], ScanEl[timeIndex], xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Uerr = GridData( UerrScn[timeIndex], ScanAz[timeIndex], ScanEl[timeIndex], xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
Perr = GridData( PerrScn[timeIndex], ScanAz[timeIndex], ScanEl[timeIndex], xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
EVPAerr = GridData( AerrScn[timeIndex], ScanAz[timeIndex], ScanEl[timeIndex], xi.reshape(xi.size), yi.reshape(xi.size), 5 ).reshape(len(xi), len(xi))
#
#----Plot Q error
plt.subplot( 2, 2, 1); plt.contourf(xi, yi, 100.0*Qerr, np.linspace(-2.0, 2.0, 41)); plt.colorbar(); plt.title(prefix + ' Q_err %')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Qerr at Center       = %4.1f %%' % ( 100.0* np.mean(QerrScn[IndexCenter]) ); plt.text(-0.8*FWHM, -0.7*FWHM, text_sd, size='x-small')
text_sd = 'Max Qerr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(QerrScn[Index3dB])) ); plt.text(-0.8*FWHM, -0.8*FWHM, text_sd, size='x-small')
text_sd = 'Max Qerr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(QerrScn[Index6dB])) ); plt.text(-0.8*FWHM, -0.9*FWHM, text_sd, size='x-small')
#----Plot U error
plt.subplot( 2, 2, 2); plt.contourf(xi, yi, 100.0*Uerr, np.linspace(-2.0, 2.0, 41)); plt.colorbar(); plt.title(prefix + ' U_err %')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Uerr at Center       = %4.1f %%' % ( 100.0* np.mean(UerrScn[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = 'Max Uerr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(UerrScn[Index3dB])) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = 'Max Uerr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(UerrScn[Index6dB])) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#----Plot P error
plt.subplot( 2, 2, 3); plt.contourf(xi, yi, 100.0*Perr, np.linspace(-2.0, 2.0, 41)); plt.colorbar(); plt.title(prefix + ' P_err %')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Perr at Center       = %4.1f %%' % ( 100.0* np.mean(PerrScn[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = 'Max Perr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(PerrScn[Index3dB])) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = 'Max Perr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(PerrScn[Index6dB])) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#----Plot EVPA error
plt.subplot( 2, 2, 4); plt.contourf(xi, yi, EVPAerr, np.linspace(-5, 5, 21)); plt.colorbar(); plt.title(prefix + ' EVPA_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'EVPAerr at Center       = %4.1f deg' % ( np.mean(AerrScn[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = 'Max EVPAerr (-3dB beam) = %4.1f deg' % ( max(abs(AerrScn[Index3dB])) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = 'Max EVPAerr (-6dB beam) = %4.1f deg' % ( max(abs(AerrScn[Index6dB])) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
plt.savefig( prefix + '-StokesErr.pdf', form='pdf'); plt.close()
plt.close()
"""
