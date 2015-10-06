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
#BP_ant = np.load(wd + prefix + '-SPW' + `spw[0]` + '-BPant.npy')
#XYdelay = np.load(wd + prefix + '-SPW' + `spw[0]` + '-XYdelay.npy')
#solution = np.load(wd + QUXY + '.QUXY.npy')
BP_ant = np.load(wd + BPfile)
XYdelay = np.load(wd + XYdelayfile)
solution = np.load(wd + QUXYfile)
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
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
for bl_index in range(blNum):
    BlMap[bl_index], BlInv[bl_index] = Ant2BlD(AntIndex[ant0[bl_index]], AntIndex[ant1[bl_index]])
    blWeight[bl_index] = antWeight[AntIndex[ant0[bl_index]]]* antWeight[AntIndex[ant1[bl_index]]]
#
trkBlIndex  = np.where(blWeight == 1.0)[0].tolist(); trkBlNum  = len(trkBlIndex)        # Ref-Ref baselines
ScTrBlIndex = np.where(blWeight == 0.5)[0].tolist(); ScTrBlNum = len(ScTrBlIndex)       # Ref-Scan baselines
ScScBlIndex = np.where(blWeight == 0.25)[0].tolist(); ScScBlNum = len(ScScBlIndex)      # Scan-Scan baselines
#-- Load all-baseline visibilities
print('Loading Visibilities ....')
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[0], scan[0])
#-------- Bunch/Sel
print('Applying Bandpass Calibration....')
chNum = Xspec.shape[1]
if chNum > 1:
    chRange = range( int(0.1*chNum), int(0.9*chNum))
#
if chNum == 1:
    temp = Xspec[:,0, BlMap]
else:
    Xspec[0] = (Xspec[0].transpose(2,1,0) / (BP_ant[ant1,0].conjugate()* BP_ant[ant0,0])).transpose(2, 1, 0)
    Xspec[1] = (Xspec[1].transpose(2,1,0) / (BP_ant[ant1,0].conjugate()* BP_ant[ant0,1])).transpose(2, 1, 0)
    Xspec[2] = (Xspec[2].transpose(2,1,0) / (BP_ant[ant1,1].conjugate()* BP_ant[ant0,0])).transpose(2, 1, 0)
    Xspec[3] = (Xspec[3].transpose(2,1,0) / (BP_ant[ant1,1].conjugate()* BP_ant[ant0,1])).transpose(2, 1, 0)
    #
    XYdlSpec = delay_cal( np.ones([chNum], dtype=complex), XYdelay )
    Xspec[1] = (Xspec[1].transpose(1,2,0) * XYdlSpec).transpose(2,0,1)
    Xspec[2] = (Xspec[2].transpose(1,2,0) / XYdlSpec).transpose(2,0,1)
    Xspec = Xspec[:,:,BlMap]
#
Ximag = Xspec.transpose(0, 1, 3, 2).imag * (-2.0* np.array(BlInv) + 1.0)      # For inverted baselines
Xspec.imag = Ximag.transpose(0, 1, 3, 2)                                      # visibilities become complex conjugate
XX = Xspec[0].copy(); XY = Xspec[1].copy(); YX = Xspec[2].copy(); YY = Xspec[3].copy()
invIndex = np.where(BlInv)[0].tolist()          # For inverted baselines
XY[:,invIndex] = Xspec[2][:,invIndex].copy()     # XY and YX are swapped
YX[:,invIndex] = Xspec[1][:,invIndex].copy()     #
#-- BL within antenna groups
TrkXX  = XX[:,trkBlIndex];  TrkXY  = XY[:,trkBlIndex];  TrkYX  = YX[:,trkBlIndex];  TrkYY  = YY[:,trkBlIndex]  # Ref-Ref 
ScTrXX = XX[:,ScTrBlIndex]; ScTrXY = XY[:,ScTrBlIndex]; ScTrYX = YX[:,ScTrBlIndex]; ScTrYY = YY[:,ScTrBlIndex] # Ref-Sc 
ScScXX = XX[:,ScScBlIndex]; ScScXY = XY[:,ScScBlIndex]; ScScYX = YX[:,ScScBlIndex]; ScScYY = YY[:,ScScBlIndex] # Ref-Sc 
#-------- Frequency and Wavelength
chNum, chWid, Freq = GetChNum(msfile, spw[0])
FWHM = GetFWHM(msfile, spw[0], AntD)
Freq *= 1.0e-9  # [GHz]
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
#-------- Loop for each channel
print('Applying Gain Correction....')
VisXX = np.zeros([blNum, timeNum, chNum], dtype=complex)
VisXY = np.zeros([blNum, timeNum, chNum], dtype=complex)
VisYX = np.zeros([blNum, timeNum, chNum], dtype=complex)
VisYY = np.zeros([blNum, timeNum, chNum], dtype=complex)
#-- Gain solutions using parallel-hand correlations
GainX, GainY = polariGain( np.mean(XX[chRange], axis=0), np.mean(YY[chRange], axis=0), PA, solution[0], solution[1])
GainY = GainY* exp(solution[2]* 1.0j)                                   # XY-phase correction
#-- Gain Correction
for ch_index in range(chNum):
    VisXX[:,:,ch_index] = gainCalVis( XX[ch_index], GainX, GainX )
    VisXY[:,:,ch_index] = gainCalVis( XY[ch_index], GainX, GainY )
    VisYX[:,:,ch_index] = gainCalVis( YX[ch_index], GainY, GainX )
    VisYY[:,:,ch_index] = gainCalVis( YY[ch_index], GainY, GainY )
#
print('-------- Determining Antenna-based D-terms (refants) ----')
Dx = np.zeros([antNum, timeNum, chNum], dtype=complex)
Dy = np.zeros([antNum, timeNum, chNum], dtype=complex)
for ch_index in range(chNum):
    print `ch_index` + ' / ' + `chNum` + '\r'
    #-------- RefAnt D-term
    # segmentation of PA
    PAwidth = 0.005
    PAsegNum = int((max(PA) - min(PA))/PAwidth)
    for seg_index in range(PAsegNum):
        timeIndexRange = range( (seg_index* timeNum/PAsegNum), ((seg_index + 1)* timeNum/PAsegNum) )
        PS = np.dot(PAMatrix(np.mean(PA[timeIndexRange])), np.array([1.0, solution[0], solution[1], 0.0])).real
        VisTime = np.r_[np.mean(VisXX[trkBlIndex][:,timeIndexRange, ch_index], axis=1), np.mean(VisXY[trkBlIndex][:,timeIndexRange, ch_index], axis=1), np.mean(VisYX[trkBlIndex][:,timeIndexRange, ch_index], axis=1), np.mean(VisYY[trkBlIndex][:,timeIndexRange, ch_index], axis=1)]
        TrkDx, TrkDy = Vis2solveDD( VisTime, PS )
        for time_index in timeIndexRange:
            Dx[range(trkAntNum), time_index, ch_index] = TrkDx
            Dy[range(trkAntNum), time_index, ch_index] = TrkDy
        #
    #
#
#-------- Record on-axis D-term spectrum
logfile = open(prefix + '-SPW' + `spw[0]` + '-TrkDtermSpec.log', 'w')
text_sd = 'ant ch ReDx ImDx ReDy ImDy'
logfile.write(text_sd + '\n')
DxMean = np.mean(Dx, axis=1)
DyMean = np.mean(Dy, axis=1)
for time_index in range(timeNum):
    Dx[:,time_index] = DxMean
    Dy[:,time_index] = DyMean
#
for ant_index in range(trkAntNum):
    for ch_index in range(chNum):
        text_sd = '%s %d %8.6f %8.6f %8.6f %8.6f' % (antList[trkAnt[ant_index]], ch_index, Dx[ant_index, 0, ch_index].real, Dx[ant_index, 0, ch_index].imag, Dy[ant_index, 0, ch_index].real, Dy[ant_index, 0, ch_index].imag)
        logfile.write(text_sd + '\n')
    #
#
logfile.close()
#-------- Stokes Parameters measured by tracking antennas
StokesVis = np.zeros([trkBlNum, PAsegNum, chNum, 4], dtype=complex) 
for seg_index in range(PAsegNum):
    timeIndexRange = range( (seg_index* timeNum/PAsegNum), ((seg_index + 1)* timeNum/PAsegNum) )
    Pinv = InvPAMatrix( np.mean(PA[timeIndexRange] ))
    for bl_index in range(trkBlNum):
        for ch_index in range(chNum):
            Minv = InvMullerMatrix(
                np.mean(Dx[ant1[bl_index], timeIndexRange, ch_index]),
                np.mean(Dy[ant1[bl_index], timeIndexRange, ch_index]),
                np.mean(Dx[ant0[bl_index], timeIndexRange, ch_index]),
                np.mean(Dy[ant0[bl_index], timeIndexRange, ch_index]))
            StokesVis[bl_index, seg_index, ch_index] = np.dot(Minv, np.array([
                np.mean(VisXX[trkBlIndex[bl_index], timeIndexRange, ch_index]),
                np.mean(VisXY[trkBlIndex[bl_index], timeIndexRange, ch_index]),
                np.mean(VisYX[trkBlIndex[bl_index], timeIndexRange, ch_index]),
                np.mean(VisYY[trkBlIndex[bl_index], timeIndexRange, ch_index])]))
            StokesVis[bl_index, seg_index, ch_index] = np.dot(Pinv, StokesVis[bl_index, seg_index, ch_index])
        #
    #
#
TrkI = np.mean( StokesVis[:,:,chRange,0]).real
TrkQ = np.mean( StokesVis[:,:,chRange,1]).real
TrkU = np.mean( StokesVis[:,:,chRange,2]).real
TrkV = np.mean( StokesVis[:,:,chRange,3]).real
text_sd = 'TrkMeas / Model: I=%6.4f / %6.4f  Q=%6.4f / %6.4f U=%6.4f / %6.4f V=%6.4f / %6.4f' % (TrkI, 1.0, TrkQ, CalQ, TrkU, CalU, TrkV, 0.0); print text_sd
UCmQS = TrkU* np.cos(2.0*PA) - TrkQ* np.sin(2.0*PA)   # U cos - Q sin
QCpUS = TrkQ* np.cos(2.0*PA) + TrkU* np.sin(2.0*PA)   # Q cos + U sin
logfile = open(prefix + '-SPW' + `spw[0]` + '-trkStokes.log', 'w')
text_sd = 'CH I Q U V'; logfile.write(text_sd + '\n')
for ch_index in range(chNum):
    text_sd = '%d %8.6f %8.6f %8.6f %8.6f' % (ch_index, np.mean(StokesVis[:,:,ch_index,0]), np.mean(StokesVis[:,:,ch_index,1]), np.mean(StokesVis[:,:,ch_index,2]), np.mean(StokesVis[:,:,ch_index,3]))
    logfile.write(text_sd + '\n')
#
logfile.close()
#-------- Determination of D-terms in scanning antennas
print('-------- Determining Antenna-based D-terms (scan ants) ----')
for ant_index in range(scnAntNum):
    antID = trkAntNum + ant_index
    print 'Determining D-term of ' + antList[antID]
    TrkScnBL = range(antID* (antID - 1) / 2, antID* (antID - 1) / 2 + trkAntNum)
    for time_index in range(timeNum):
        PS = np.dot(PAMatrix(PA[time_index]), np.array([1.0, TrkQ, TrkU, 0.0])).real
        for ch_index in range(chNum):
            VisTime = np.r_[VisXX[TrkScnBL, time_index, ch_index], VisXY[TrkScnBL, time_index, ch_index], VisYX[TrkScnBL, time_index, ch_index], VisYY[TrkScnBL, time_index, ch_index]]
            Dx[antID, time_index, ch_index], Dy[antID, time_index, ch_index] = Vis2solveD(VisTime, Dx[0:trkAntNum, time_index, ch_index], Dy[0:trkAntNum, time_index, ch_index], PS )
        #
    #
#
#-------- Plot D-term spectrum for beam position
logfile = open(prefix + '-SPW' + `spw[0]` + '-DtermSpec.log', 'w')
text_sd = 'ant beamoff branch ch ReDx ImDx ReDy ImDy'
logfile.write(text_sd + '\n')
thresh = np.r_[0.0, np.linspace( min(np.sqrt(ScanAz**2 + ScanEl**2)), max(np.sqrt(ScanAz**2 + ScanEl**2)), num=16) + min(np.sqrt(ScanAz**2 + ScanEl**2))]
Dist2 = ScanAz**2 + ScanEl**2
xrange, yrange = [min(Freq[chRange]), max(Freq[chRange])], [-0.1, 0.1]
for ant_index in range(scnAntNum):
    antID = scnAnt[ant_index]
    DantID = trkAntNum + ant_index
    for thresh_index in range(6):
        time_index = list(set(np.where(Dist2 > thresh[thresh_index]**2 )[0]) & set(np.where(Dist2 < thresh[thresh_index + 1]**2 )[0]))
        fig = plt.figure(thresh_index, figsize = (8,11))
        fig.text(0.45, 0.05, 'Frequency [GHz]')
        fig.text(0.05, 0.45, 'D-term', rotation=90)
        plt.suptitle(prefix + ' ' + antList[antID] + ' D-term@' + `round(np.median(np.sqrt(Dist2[time_index])),1)` + ' arcsec')
        #-------- Plot beam map
        plt.subplot2grid( (9,6), (0,5), aspect=1)
        plt.plot(ScanAz, ScanEl, ',', color='k', alpha=0.1)
        plt.plot( ScanAz[time_index], ScanEl[time_index], '.', color='r')
        circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y, color='green' )
        circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y, color='blue' )
        plt.axis([-FWHM[0], FWHM[0], -FWHM[0], FWHM[0]], fontsize=3)
        plt.tick_params(labelsize = 6)
        #
        #-------- Plot Mean D-term
        plt.subplot2grid( (9,6), (0,0), colspan=4)
        fill_color='g'
        plt.fill([xrange[0], xrange[0], xrange[1], xrange[1]], [yrange[1], yrange[0], yrange[0], yrange[1]], fill_color, alpha=0.1)
        plt.plot( Freq, np.mean(Dx[DantID, time_index], axis=0).real, color='cyan',     label = 'ReDx', ls='steps-mid')
        plt.plot( Freq, np.mean(Dx[DantID, time_index], axis=0).imag, color='darkblue', label = 'ImDx', ls='steps-mid')
        plt.plot( Freq, np.mean(Dy[DantID, time_index], axis=0).real, color='magenta',  label = 'ReDy', ls='steps-mid')
        plt.plot( Freq, np.mean(Dy[DantID, time_index], axis=0).imag, color='darkred',  label = 'ImDy', ls='steps-mid')
        plt.axis([min(Freq), max(Freq), -0.2,0.2], fontsize=3)
        plt.legend(loc = 'best', prop={'size' :6}, numpoints = 1)
        plt.tick_params(labelsize = 6)
        for index in range(48):
            fill_color='g'
            if max(abs(Dx[DantID, time_index[index], chRange].real)) > yrange[1]: fill_color = 'r'
            if max(abs(Dx[DantID, time_index[index], chRange].imag)) > yrange[1]: fill_color = 'r'
            if max(abs(Dy[DantID, time_index[index], chRange].real)) > yrange[1]: fill_color = 'r'
            if max(abs(Dy[DantID, time_index[index], chRange].imag)) > yrange[1]: fill_color = 'r'
            plt.subplot2grid( (9,6), (int(index/6)+1, index%6))
            plt.fill([xrange[0], xrange[0], xrange[1], xrange[1]], [yrange[1], yrange[0], yrange[0], yrange[1]], fill_color, alpha=0.1)
            plt.plot( Freq, Dx[DantID, time_index[index]].real, ls='steps-mid', color='cyan',     label = 'ReDx')
            plt.plot( Freq, Dx[DantID, time_index[index]].imag, ls='steps-mid', color='darkblue', label = 'ImDx')
            plt.plot( Freq, Dy[DantID, time_index[index]].real, ls='steps-mid', color='magenta',  label = 'ReDy')
            plt.plot( Freq, Dy[DantID, time_index[index]].imag, ls='steps-mid', color='darkred',  label = 'ImDy')
            plt.axis([min(Freq), max(Freq), -0.2,0.2], fontsize=3)
            plt.tick_params(labelsize = 6)
            for ch_index in range(chNum):
                text_sd = '%s %4.1f %d %d %8.6f %8.6f %8.6f %8.6f' % (antList[antID], np.median(np.sqrt(Dist2[time_index])), index, ch_index, Dx[DantID, time_index[index], ch_index].real, Dx[DantID, time_index[index], ch_index].imag, Dy[DantID, time_index[index], ch_index].real, Dy[DantID, time_index[index], ch_index].imag)
                logfile.write(text_sd + '\n')
            #
        #
        plt.savefig( prefix + '-' + antList[antID] + '-SPW' + `spw[0]` + '-OFF' + `round(np.median(np.sqrt(Dist2[time_index])),1)` + '-DtermSpec.pdf', form='pdf'); plt.close()
    #
#
logfile.close()
#-------- Plot channel-averaged D-terms of scanning antennas
print('-------- Plot D-term Maps for scan ants ----')
for ant_index in range(scnAntNum):
    antID = scnAnt[ant_index]
    DantID = trkAntNum + ant_index
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
    ReDxmap = GridData( np.mean(Dx[DantID][:, chRange], axis=1).real, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM[antID]/16).reshape(len(xi), len(xi))
    ImDxmap = GridData( np.mean(Dx[DantID][:, chRange], axis=1).imag, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM[antID]/16).reshape(len(xi), len(xi))
    ReDymap = GridData( np.mean(Dy[DantID][:, chRange], axis=1).real, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM[antID]/16).reshape(len(xi), len(xi))
    ImDymap = GridData( np.mean(Dy[DantID][:, chRange], axis=1).imag, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM[antID]/16).reshape(len(xi), len(xi))
    #---- plot Re(Dx)
    plt.subplot( 2, 2, 1, aspect=1); plt.contourf(xi, yi, ReDxmap, np.linspace(-0.10, 0.10, 11)); plt.colorbar(); plt.title('Re(Dx)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Re(Dx) at Center = %5.3f' % ( np.mean(Dx[DantID, IndexCenter].real) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(Dx[DantID, Index3dB].real), np.min(Dx[DantID, Index3dB].real) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(Dx[DantID, Index6dB].real), np.min(Dx[DantID, Index6dB].real) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    #---- plot Im(Dx)
    plt.subplot( 2, 2, 2, aspect=1); plt.contourf(xi, yi, ImDxmap, np.linspace(-0.10, 0.10, 11)); plt.colorbar(); plt.title('Im(Dx)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Im(Dx) at Center = %5.3f' % ( np.mean(Dx[DantID, IndexCenter].imag) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(Dx[DantID, Index3dB].imag), np.min(Dx[DantID, Index3dB].imag) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(Dx[DantID, Index6dB].imag), np.min(Dx[DantID, Index6dB].imag) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    #---- plot Re(Dy)
    plt.subplot( 2, 2, 3, aspect=1); plt.contourf(xi, yi, ReDymap, np.linspace(-0.10, 0.10, 11)); plt.colorbar(); plt.title('Re(Dy)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Re(Dy) at Center = %5.3f' % ( np.mean(Dy[DantID, IndexCenter].real) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(Dy[DantID, Index3dB].real), np.min(Dy[DantID, Index3dB].real) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(Dy[DantID, Index6dB].real), np.min(Dy[DantID, Index6dB].real) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    #---- plot Im(Dy)
    plt.subplot( 2, 2, 4, aspect=1); plt.contourf(xi, yi, ImDymap, np.linspace(-0.10, 0.10, 11)); plt.colorbar(); plt.title('Im(Dy)')
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, FWHM[antID]/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Im(Dy) at Center = %5.3f' % ( np.mean(Dy[DantID, IndexCenter].imag) ); plt.text(-1.6*FWHM[antID], -1.5*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(Dy[DantID, Index3dB].imag), np.min(Dy[DantID, Index3dB].imag) ); plt.text(-1.6*FWHM[antID], -1.7*FWHM[antID], text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(Dy[DantID, Index6dB].imag), np.min(Dy[DantID, Index6dB].imag) ); plt.text(-1.6*FWHM[antID], -1.9*FWHM[antID], text_sd, size='x-small')
    plt.plot( ScanAz, ScanEl, '.', color='k', alpha=0.1)
    plt.axis([-2.0*FWHM[antID], 2.0*FWHM[antID], -2.0*FWHM[antID], 2.0*FWHM[antID]])
    plt.savefig( prefix + '-' + antList[antID] + '-SPW' + `spw[0]` + '-DtermMap.pdf', form='pdf'); plt.close()
    plt.close()
#
#-------- D-term-corrected Stokes parameters --------
print('-------- D-term-corrected Stokes parameters ----')
StokesVis = np.zeros([ScScBlNum, timeNum, chNum, 4], dtype=complex) 
for time_index in range(timeNum):
    print `time_index` + ' / ' + `timeNum` + '\r'
    Pinv = InvPAMatrix( PA[time_index] )
    for bl_index in range(ScScBlNum):
        BlID = ScScBlIndex[bl_index]
        ants = Bl2Ant(BlID)
        for ch_index in range(chNum):
            Minv = InvMullerMatrix( Dx[ants[1], time_index, ch_index], Dy[ants[1], time_index, ch_index], Dx[ants[0], time_index, ch_index], Dy[ants[0], time_index, ch_index])
            StokesVis[bl_index, time_index, ch_index] = np.dot(Pinv, np.dot(Minv, np.array( [
                VisXX[BlID, time_index, ch_index],
                VisXY[BlID, time_index, ch_index],
                VisYX[BlID, time_index, ch_index],
                VisYY[BlID, time_index, ch_index]])))
        #
    #
#
ScnIspec = np.mean( StokesVis[:,:,:,0], axis=0 ).real
ScnQspec = np.mean( StokesVis[:,:,:,1], axis=0 ).real
ScnUspec = np.mean( StokesVis[:,:,:,2], axis=0 ).real
ScnVspec = np.mean( StokesVis[:,:,:,3], axis=0 ).real
ScnI = np.mean( ScnIspec[:,chRange], axis=1 )
ScnQ = np.mean( ScnQspec[:,chRange], axis=1 )
ScnU = np.mean( ScnUspec[:,chRange], axis=1 )
ScnV = np.mean( ScnVspec[:,chRange], axis=1 )
Qerr = ScnQ - TrkQ
Uerr = ScnU - TrkU
Perr = sqrt(ScnQ**2 + ScnU**2) - sqrt(TrkQ**2 + TrkU**2)
Aerr = np.arctan( (ScnU* TrkQ - ScnQ* TrkU) / (ScnQ* TrkQ + ScnU* TrkU) )* 90.0/math.pi
#-------- Plot Stokes Beam Map
logfile = open(prefix + '-SPW' + `spw[0]` + '-Stokes.log', 'w')
text_sd = 'I I3max I3min I6max I6min Q Q3max Q3min Q6max Q6min U U3max U3min U6max U6min V V3max V3min V6max V6min'
logfile.write(text_sd + '\n')
fig = plt.figure( figsize = (10,10))
fig.text(0.45, 0.05, 'Az Offset [arcsec]')
fig.text(0.05, 0.45, 'El Offset [arcsec]', rotation=90)
fig.text(0.45, 0.95, prefix)
FWHM = np.max(FWHM)
xi, yi = np.mgrid[ -floor(FWHM):floor(FWHM):128j, -floor(FWHM):floor(FWHM):128j]
Imap = GridData( ScnI, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
Qmap = GridData( ScnQ, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
Umap = GridData( ScnU, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
Vmap = GridData( ScnV, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
QEmap = GridData(Qerr, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
UEmap = GridData(Uerr, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
PEmap = GridData(Perr, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
AEmap = GridData(Aerr, ScanAz, ScanEl, xi.reshape(xi.size), yi.reshape(xi.size), FWHM/16).reshape(len(xi), len(xi))
IndexCenter = np.where( ScanAz**2 + ScanEl**2 < 0.005* FWHM**2 )[0]
Index3dB = np.where( ScanAz**2 + ScanEl**2 < 0.25* FWHM**2 )[0]
Index6dB = np.where( ScanAz**2 + ScanEl**2 < 0.5* FWHM**2 )[0]
#---- Plot I map
plt.subplot(2, 2, 1, aspect=1); plt.contourf(xi, yi, Imap, np.linspace(0.9, 1.1, 21)); plt.colorbar(); plt.title('Stokes I')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'I at Center   = %5.3f' % ( np.mean(ScnI[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnI[Index3dB]), np.min(ScnI[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnI[Index6dB]), np.min(ScnI[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot Q map
plt.subplot(2, 2, 2, aspect=1); plt.contourf(xi, yi, Qmap, np.linspace(0.01*(floor(TrkQ* 100)-5), 0.01*(floor(TrkQ* 100)+5), 21)); plt.colorbar(); plt.title('Stokes Q')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Q at Center   = %5.3f' % ( np.mean(ScnQ[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnQ[Index3dB]), np.min(ScnQ[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnQ[Index6dB]), np.min(ScnQ[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot U map
plt.subplot(2, 2, 3, aspect=1); plt.contourf(xi, yi, Umap, np.linspace(0.01*(floor(TrkU* 100)-5), 0.01*(floor(TrkU* 100)+5), 21)); plt.colorbar(); plt.title('Stokes U')
circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'U at Center   = %5.3f' % ( np.mean(ScnU[IndexCenter]) ); plt.text(-0.8*FWHM, -0.75*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnU[Index3dB]), np.min(ScnU[Index3dB]) ); plt.text(-0.8*FWHM, -0.85*FWHM, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnU[Index6dB]), np.min(ScnU[Index6dB]) ); plt.text(-0.8*FWHM, -0.95*FWHM, text_sd, size='x-small')
#---- Plot V map
plt.subplot(2, 2, 4, aspect=1); plt.contourf(xi, yi, Vmap, np.linspace(-0.05, 0.05, 21)); plt.colorbar(); plt.title('Stokes V')
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
text_sd = '%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f ' % (np.mean(ScnI[IndexCenter]), np.max(ScnI[Index3dB]), np.min(ScnI[Index3dB]), np.max(ScnI[Index6dB]), np.min(ScnI[Index6dB]), np.mean(ScnQ[IndexCenter]), np.max(ScnQ[Index3dB]), np.min(ScnQ[Index3dB]), np.max(ScnQ[Index6dB]), np.min(ScnQ[Index6dB]), np.mean(ScnU[IndexCenter]), np.max(ScnU[Index3dB]), np.min(ScnU[Index3dB]), np.max(ScnU[Index6dB]), np.min(ScnU[Index6dB]), np.mean(ScnV[IndexCenter]), np.max(ScnV[Index3dB]), np.min(ScnV[Index3dB]), np.max(ScnV[Index6dB]), np.min(ScnV[Index6dB]))
logfile.write(text_sd + '\n')
logfile.close()
#
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
#-------- Plot Stokes Spectra at Beam Position
logfile = open(prefix + '-SPW' + `spw[0]` + '-StokesSpec.log', 'w')
text_sd = 'beamoff branch CH I Q U V ';  logfile.write(text_sd + '\n')
text_sd = '%4.1f %d %8.6f %8.6f %8.6f %8.6f' % (0.0,  0, TrkI, TrkQ, TrkU, TrkV); logfile.write(text_sd + '\n')
xrange, yrange = [min(Freq[chRange]), max(Freq[chRange])], [-0.01, 0.01]
for thresh_index in range(6):
    time_index = list(set(np.where(Dist2 > thresh[thresh_index]**2 )[0]) & set(np.where(Dist2 < thresh[thresh_index + 1]**2 )[0]))
    fig = plt.figure(thresh_index, figsize = (8,11))
    fig.text(0.45, 0.05, 'Frequency [GHz]')
    fig.text(0.05, 0.45, 'Stokes Residual [scaled by Stokes I]', rotation=90)
    plt.suptitle(prefix + ' ' + ' Stokes Residuals@' + `round(np.median(np.sqrt(Dist2[time_index])),1)` + ' arcsec')
    #-------- Plot beam map
    plt.subplot2grid( (9,6), (0,5), aspect=1)
    plt.plot(ScanAz, ScanEl, ',', color='k', alpha=0.1)
    plt.plot( ScanAz[time_index], ScanEl[time_index], '.', color='r')
    circle_x, circle_y = circlePoints(0, 0, FWHM/2); plt.plot( circle_x, circle_y, color='green' )
    circle_x, circle_y = circlePoints(0, 0, FWHM/sqrt(2)); plt.plot( circle_x, circle_y, color='blue' )
    plt.axis([-FWHM, FWHM, -FWHM, FWHM])
    plt.tick_params(labelsize = 6)
    #
    #-------- Plot Mean IQUV
    plt.subplot2grid( (9,6), (0,0), colspan=4)
    fill_color='g'
    plt.fill([xrange[0], xrange[0], xrange[1], xrange[1]], [yrange[1], yrange[0], yrange[0], yrange[1]], fill_color, alpha=0.1)
    plt.plot( Freq, np.mean(ScnIspec[time_index] - TrkI, axis=0), color='k',     label = 'I', ls='steps-mid')
    plt.plot( Freq, np.mean(ScnQspec[time_index] - TrkQ, axis=0), color='r', label = 'Q', ls='steps-mid')
    plt.plot( Freq, np.mean(ScnUspec[time_index] - TrkU, axis=0), color='g',  label = 'U', ls='steps-mid')
    plt.plot( Freq, np.mean(ScnVspec[time_index] - TrkV, axis=0), color='b',  label = 'V', ls='steps-mid')
    plt.axis([min(Freq), max(Freq), -0.02,0.02], fontsize=3)
    plt.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    plt.tick_params(labelsize = 6)
    #
    for index in range(48):
        fill_color='g'
        if max(abs(ScnQspec[time_index[index], chRange] - TrkQ)) > yrange[1]: fill_color = 'r'
        if max(abs(ScnUspec[time_index[index], chRange] - TrkU)) > yrange[1]: fill_color = 'r'
        if max(abs(ScnVspec[time_index[index], chRange] - TrkV)) > yrange[1]: fill_color = 'r'
        plt.subplot2grid( (9,6), (int(index/6)+1, index%6))
        plt.fill([xrange[0], xrange[0], xrange[1], xrange[1]], [yrange[1], yrange[0], yrange[0], yrange[1]], fill_color, alpha=0.1)
        plt.plot( Freq, ScnIspec[time_index[index]] - TrkI, ls='steps-mid', color='k',     label = 'I')
        plt.plot( Freq, ScnQspec[time_index[index]] - TrkQ, ls='steps-mid', color='r', label = 'Q')
        plt.plot( Freq, ScnUspec[time_index[index]] - TrkU, ls='steps-mid', color='g',  label = 'U')
        plt.plot( Freq, ScnVspec[time_index[index]] - TrkV, ls='steps-mid', color='b',  label = 'V')
        plt.axis([min(Freq), max(Freq), -0.05,0.05], fontsize=3)
        plt.tick_params(labelsize = 6)
        for ch_index in range(chNum):
            text_sd = '%4.1f %d %d %8.6f %8.6f %8.6f %8.6f' % (np.median(np.sqrt(Dist2[time_index])), index, ch_index, ScnIspec[time_index[index], ch_index], ScnQspec[time_index[index], ch_index], ScnUspec[time_index[index], ch_index], ScnVspec[time_index[index], ch_index]); logfile.write(text_sd + '\n')
        #
    #
    plt.savefig( prefix + '-' + '-SPW' + `spw[0]` + '-OFF' + `round(np.median(np.sqrt(Dist2[time_index])),1)` + '-StokesSpec.pdf', form='pdf'); plt.close()
#
logfile.close()