execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
from scipy.interpolate import UnivariateSpline
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
BP_ant = np.load(wd + prefix + '.BPant.npy')
XYdelay = np.load(wd + prefix + '.XYdelay.npy')
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
#-- Load all-baseline visibilities in CalScan
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], Dscan)
chNum = Xspec.shape[1]
#-------- Bunch/Sel
chNum = Xspec.shape[1]
timeNum = len(timeStamp)
if chNum == 1:
    temp = Xspec[:,0, BlMap]
else:
    temp = np.mean(Xspec[:,chRange], axis=1)
#
#-- Complex Conjugate and swap XY-YX in inverted baselines
Ximag = temp.transpose(0,2,1).imag * (-2.0* np.array(BlInv) + 1.0)      # For inverted baselines
temp.imag = Ximag.transpose(0,2,1)                                      # visibilities become complex conjugate
XX = temp[0].copy(); XY = temp[1].copy(); YX = temp[2].copy(); YY = temp[3].copy()
invIndex = np.where(BlInv)[0].tolist()          # For inverted baselines
XY[invIndex] = temp[2,invIndex].copy()          # XY and YX are swapped
YX[invIndex] = temp[1,invIndex].copy()          #
#-- Gain solutions
AZ = np.zeros([timeNum]); EL = np.zeros([timeNum]); PA = np.zeros([timeNum])
scanTime, AntID, az, el = GetAzEl(msfile)
index = np.where( AntID == refAntIndex[0])[0].tolist()
azel = np.r_[az[index], el[index]].reshape(2, len(index))
for time_index in range(timeNum):
    AZ[time_index], EL[time_index] = AzElMatch( timeStamp[time_index], scanTime[index], AntID, refAntIndex, np.median(diff(timeStamp)), azel)
    PA[time_index] = AzEl2PA(AZ[time_index], EL[time_index], ALMA_lat) - BANDPA   #
#
UCmQS = solution[1]* np.cos(2.0*PA) - solution[0]* np.sin(2.0*PA)   # U cos - Q sin
QCpUS = solution[0]* np.cos(2.0*PA) + solution[1]* np.sin(2.0*PA)   # Q cos + U sin
#-- Gain solutions using parallel-hand correlations
RGainX, RGainY = polariGain(XX, YY, PA, solution[0], solution[1])
RGainY = RGainY* exp(solution[2]* 1.0j)                                   # XY-phase correction
#-- Gain Correction
VisXX = gainCalVis( XX, RGainX, RGainX )
VisXY = gainCalVis( XY, RGainX, RGainY )
VisYX = gainCalVis( YX, RGainY, RGainX )
VisYY = gainCalVis( YY, RGainY, RGainY )
#
#-------- RefAnt D-term
print('-------- Determining Antenna-based D-terms (refants) ----')
Dx = np.zeros([antNum], dtype=complex)
Dy = np.zeros([antNum], dtype=complex)
PS = np.dot(PAMatrix(np.mean(PA)), np.array([1.0, solution[0], solution[1], 0.0])).real
VisTime = np.r_[np.mean(VisXX, axis=1), np.mean(VisXY, axis=1), np.mean(VisYX, axis=1), np.mean(VisYY, axis=1)]
Dx, Dy = Vis2solveDD( VisTime, PS )
#
#-------- Stokes Parameters measured by tracking antennas
StokesVis = np.zeros([blNum, timeNum, 4], dtype=complex)
for time_index in range(timeNum):
    Pinv = InvPAMatrix( PA[time_index] )
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        Minv = InvMullerMatrix( Dx[ants[1]], Dy[ants[1]], Dx[ants[0]], Dy[ants[0]])
        StokesVis[bl_index, time_index] = np.dot(Minv, np.array( [VisXX[bl_index, time_index], VisXY[bl_index, time_index], VisYX[bl_index, time_index], VisYY[bl_index, time_index]]))
        StokesVis[bl_index, time_index] = np.dot(Pinv, StokesVis[bl_index, time_index])
    #
#
TrkI = np.mean( StokesVis[:,:,0], axis=(0,1) ).real
TrkQ = np.mean( StokesVis[:,:,1], axis=(0,1) ).real
TrkU = np.mean( StokesVis[:,:,2], axis=(0,1) ).real
TrkV = np.mean( StokesVis[:,:,3], axis=(0,1) ).real
text_sd = 'TrkMeas / Model: I=%6.4f / %6.4f  Q=%6.4f / %6.4f U=%6.4f / %6.4f V=%6.4f / %6.4f' % (TrkI, 1.0, TrkQ, CalQ, TrkU, CalU, TrkV, 0.0); print text_sd
#-------- Stokes Parameters for scanning antennas
#-- Load all-baseline visibilities
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], scan)
#-------- Bunch/Sel
chNum = Xspec.shape[1]
if chNum == 1:
    temp = Xspec[:,0, BlMap]
else:
    # BPcal, Bunch/Sel process comes here!
    temp = np.mean(Xspec[:,chRange], axis=1)
#
Ximag = temp.transpose(0,2,1).imag * (-2.0* np.array(BlInv) + 1.0)      # For inverted baselines
temp.imag = Ximag.transpose(0,2,1)                                      # visibilities become complex conjugate
XX = temp[0].copy(); XY = temp[1].copy(); YX = temp[2].copy(); YY = temp[3].copy()
invIndex = np.where(BlInv)[0].tolist()          # For inverted baselines
XY[invIndex] = temp[2,invIndex].copy()          # XY and YX are swapped
YX[invIndex] = temp[1,invIndex].copy()          #
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
    PA[time_index] = AzEl2PA(AZ[time_index], EL[time_index], ALMA_lat) - BANDPA   #
#
UCmQS = solution[1]* np.cos(2.0*PA) - solution[0]* np.sin(2.0*PA)   # U cos - Q sin
QCpUS = solution[0]* np.cos(2.0*PA) + solution[1]* np.sin(2.0*PA)   # Q cos + U sin
#-- Gain solutions using parallel-hand correlations
SGainX, SGainY = polariGain(XX, YY, PA, solution[0], solution[1])
SGainY = SGainY* exp(solution[2]* 1.0j)                                   # XY-phase correction
#
GainX = ((SGainX / abs(SGainX)).T * np.max(abs(SGainX), axis=1)).T
GainY = ((SGainY / abs(SGainY)).T * np.max(abs(SGainY), axis=1)).T
#-- Gain Correction
VisXX = gainCalVis( XX, GainX, GainX )
VisXY = gainCalVis( XY, GainX, GainY )
VisYX = gainCalVis( YX, GainY, GainX )
VisYY = gainCalVis( YY, GainY, GainY )
#-------- Stokes Parameters after D-term correction
StokesVis = np.zeros([blNum, timeNum, 4], dtype=complex) 
for time_index in range(timeNum):
    Pinv = InvPAMatrix( PA[time_index] )
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        Minv = InvMullerMatrix( Dx[ants[1]], Dy[ants[1]], Dx[ants[0]], Dy[ants[0]])
        StokesVis[bl_index, time_index] = np.dot(Minv, np.array( [VisXX[bl_index, time_index], VisXY[bl_index, time_index], VisYX[bl_index, time_index], VisYY[bl_index, time_index]]))
        StokesVis[bl_index, time_index] = np.dot(Pinv, StokesVis[bl_index, time_index])
    #
#
#-------- D-term-corrected Stokes parameters --------
ScnI = np.mean( StokesVis[ScScBlIndex,:,0], axis=0 ).real; ScnI = ScnI / max(ScnI)
ScnQ = np.mean( StokesVis[ScScBlIndex,:,1], axis=0 ).real; ScnQ = ScnQ / max(ScnQ)
ScnU = np.mean( StokesVis[ScScBlIndex,:,2], axis=0 ).real; ScnU = ScnU / max(ScnU)
ScnV = np.mean( StokesVis[ScScBlIndex,:,3], axis=0 ).real; ScnV = ScnV / max(ScnI)
Qerr = ScnQ - TrkQ
Uerr = ScnU - TrkU
Perr = sqrt(ScnQ**2 + ScnU**2) - sqrt(TrkQ**2 + TrkU**2)
Aerr = np.arctan( (ScnU* TrkQ - ScnQ* TrkU) / (ScnQ* TrkQ + ScnU* TrkU) )* 90.0/math.pi
#-------- Plot Stokes Beam Map
fig = plt.figure( figsize = (10,10))
fig.text(0.45, 0.05, 'Az Offset [arcsec]')
fig.text(0.05, 0.45, 'El Offset [arcsec]', rotation=90)
fig.text(0.45, 0.95, prefix)
FWHM = np.max(FWHM)
xi, yi = np.mgrid[ -floor(FWHM):floor(FWHM):128j, -floor(FWHM):floor(FWHM):128j]
Imap = GridData( ScnI, -ScanAz, -ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
Qmap = GridData( ScnQ, -ScanAz, -ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
Umap = GridData( ScnU, -ScanAz, -ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
Vmap = GridData( ScnV, -ScanAz, -ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
QEmap = GridData(Qerr, -ScanAz, -ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
UEmap = GridData(Uerr, -ScanAz, -ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
PEmap = GridData(Perr, -ScanAz, -ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
AEmap = GridData(Aerr, -ScanAz, -ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
IndexCenter = np.where( ScanAz**2 + ScanEl**2 < 0.005* FWHM**2 )[0]
Index3dB = np.where( ScanAz**2 + ScanEl**2 < 0.25* FWHM**2 )[0]
Index6dB = np.where( ScanAz**2 + ScanEl**2 < 0.5* FWHM**2 )[0]
#---- Plot I map
plt.subplot(2, 2, 1, aspect=1); plt.contourf(xi, yi, Imap, np.linspace(0.0, 1.0, 11)); plt.colorbar(); plt.title('Stokes I')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'I at Center   = %5.3f' % ( np.mean(ScnI[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnI[Index3dB]), np.min(ScnI[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnI[Index6dB]), np.min(ScnI[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot Q map
plt.subplot(2, 2, 2, aspect=1); plt.contourf(xi, yi, Qmap, np.linspace(0.0, 1.0, 11)); plt.colorbar(); plt.title('Stokes Q')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Q at Center   = %5.3f' % ( np.mean(ScnQ[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnQ[Index3dB]), np.min(ScnQ[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnQ[Index6dB]), np.min(ScnQ[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot U map
plt.subplot(2, 2, 3, aspect=1); plt.contourf(xi, yi, Umap, np.linspace(0.0, 1.0, 11)); plt.colorbar(); plt.title('Stokes U')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'U at Center   = %5.3f' % ( np.mean(ScnU[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnU[Index3dB]), np.min(ScnU[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnU[Index6dB]), np.min(ScnU[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot V map
plt.subplot(2, 2, 4, aspect=1); plt.contourf(xi, yi, Vmap, np.linspace(-0.02, 0.02, 41)); plt.colorbar(); plt.title('Stokes V')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
plt.plot( ScanAz, ScanEl, '.', color='k', alpha=0.1)
text_sd = 'V at Center   = %5.3f' % ( np.mean(ScnV[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnV[Index3dB]), np.min(ScnV[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnV[Index6dB]), np.min(ScnV[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.savefig( prefix + '-SPW' + `spw[0]` + '-StokesMap.pdf', form='pdf'); plt.close()
plt.close()
#
"""
logfile = open(prefix + '-SPW' + `spw[0]` + '-QUerr.log', 'w')
text_sd = 'Qerr Qerr3max Qerr6max Uerr Uerr3max Uerr6max Perr Perr3max Perr6max EVPAerr EVPAerr3max EVPAerr6max' 
logfile.write(text_sd + '\n')
#---- Plot Q err map
fig = plt.figure( figsize = (10,10))
fig.text(0.45, 0.05, 'Az Offset [arcsec]')
fig.text(0.05, 0.45, 'El Offset [arcsec]', rotation=90)
fig.text(0.45, 0.95, prefix)
plt.subplot(2, 2, 1, aspect=1); plt.contourf(xi, yi, 100.0*QEmap, np.linspace(-5.0, 5.0, 11)); plt.colorbar(); plt.title('Q_err')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Qerr at Center       = %4.1f %%' % ( 100.0* np.mean(Qerr[IndexCenter])); plt.text(-0.8*FWHM, -0.7*FWHM, text_sd, size='x-small')
text_sd = 'Max Qerr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(Qerr[Index3dB])) ); plt.text(-0.8*FWHM, -0.8*FWHM, text_sd, size='x-small')
text_sd = 'Max Qerr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(Qerr[Index6dB])) ); plt.text(-0.8*FWHM, -0.9*FWHM, text_sd, size='x-small')
#----Plot U err map
plt.subplot(2, 2, 2, aspect=1); plt.contourf(xi, yi, 100.0*UEmap, np.linspace(-5.0, 5.0, 11)); plt.colorbar(); plt.title('U_err')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Uerr at Center       = %4.1f %%' % ( 100.0* np.mean(Uerr[IndexCenter])); plt.text(-0.8*FWHM, -0.7*FWHM, text_sd, size='x-small')
text_sd = 'Max Uerr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(Uerr[Index3dB])) ); plt.text(-0.8*FWHM, -0.8*FWHM, text_sd, size='x-small')
text_sd = 'Max Uerr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(Uerr[Index6dB])) ); plt.text(-0.8*FWHM, -0.9*FWHM, text_sd, size='x-small')
#----Plot P error
plt.subplot( 2, 2, 3, aspect=1); plt.contourf(xi, yi, 100.0*PEmap, np.linspace(-5.0, 5.0, 11)); plt.colorbar(); plt.title('P_err')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Perr at Center       = %4.1f %%' % ( 100.0* np.mean(Perr[IndexCenter])); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = 'Max Perr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(Perr[Index3dB])) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = 'Max Perr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(Perr[Index6dB])) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#----Plot EVPA error
plt.subplot( 2, 2, 4, aspect=1); plt.contourf(xi, yi, AEmap, np.linspace(-5.0, 5.0, 11)); plt.colorbar(); plt.title('EVPA_error')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'EVPAerr at Center       = %4.1f deg' % ( np.mean(Aerr[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = 'Max EVPAerr (-3dB beam) = %4.1f deg' % ( max(abs(Aerr[Index3dB])) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = 'Max EVPAerr (-6dB beam) = %4.1f deg' % ( max(abs(Aerr[Index6dB])) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#
plt.axis([-FWHM, FWHM, -FWHM, FWHM])
plt.savefig( prefix + '-SPW' + `spw[0]` + '-StokesErr.pdf', form='pdf'); plt.close()
plt.close()
text_sd = '%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f' % (np.mean(Qerr[IndexCenter]), max(abs(Qerr[Index3dB])), max(abs(Qerr[Index6dB])), np.mean(Uerr[IndexCenter]), max(abs(Uerr[Index3dB])), max(abs(Uerr[Index6dB])), np.mean(Perr[IndexCenter]), max(abs(Perr[Index3dB])), max(abs(Perr[Index6dB])), np.mean(Aerr[IndexCenter]), max(abs(Aerr[Index3dB])), max(abs(Aerr[Index6dB])))
logfile.write(text_sd + '\n')
logfile.close()
#--------- Gridding for 11x11 sampling points
az, el = readGrid(gridFile)
logfile = open(prefix + '-SPW' + `spw[0]` + '-StokesGrid.log', 'w')
text_sd = 'No. dAz      dEl        I        Q        U         V'
logfile.write(text_sd + '\n')
xi, yi = np.mgrid[ min(az):max(az):11j, max(el):min(el):11j]
Imap = GridData( ScnI, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
Qmap = GridData( ScnQ, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
Umap = GridData( ScnU, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
Vmap = GridData( ScnV, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
for x_index in range(11):
    for y_index in range(11):
        text_sd = '%d %f %f %f %f %f %f' % (x_index*11 + y_index + 1, xi[x_index, y_index], yi[x_index, y_index], Imap[x_index, y_index], Qmap[x_index, y_index], Umap[x_index, y_index], Vmap[x_index, y_index])
        logfile.write(text_sd + '\n')
    #
#
logfile.close()
"""
