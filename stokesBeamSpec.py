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
#-------- Scanning Offset < threshold
def scanThresh( msfile, ant_index, thresh ):
    Time, AntID, Offset = GetAzOffset(msfile)
    time_index = np.where( AntID == ant_index )[0]
    onAxisIndex = np.where( Offset[0, time_index]**2 + Offset[1, time_index]**2 < thresh**2 )[0]
    return time_index[onAxisIndex].tolist()
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
    for index in range(gridNum): results[index] = GridPoint( value, samp_x, samp_y, grid_x[index], grid_y[index], kernel)
    return results
#
#----------------------------------------- Antenna Mapping
# antList : (eg.) array(['DA41', 'DA42', 'DV01', ... ])     : antenna name list ordered in MS
# antID   : [0,1,2,3,4,5,6,...]                             : one-by-one number on antList 
# fragAnt : (eg.) [2,15]                                    : subset of antID to flag out
# antMap  : (eg.) [32,0,1,4, ... ]                          : subset of antID by canonical order, fragged antennas are not included
# trkAntMap : (eg.) [32,0,1,4, ...]                         : subset of antID for tracking antennas by canonical order
# refantID :  (eg.) 32                                      : antID for the reference antenna
# scnAnt  : (eg.) [33,3,38,...]                             : subset of antID for scanning antennas
#----------------------------------------- Procedures
msfile = wd + prefix + '.ms'
BPantList, BP_ant, XYdelay, solution = np.load(BPprefix + '-REF' + refantName + '.Ant.npy'), np.load(wd + BPfile), np.load(wd + XYdelayfile), np.load(wd + QUXYfile)
GYtwiddle = np.exp( (0.0 + 1.0j)* solution[2])
if QUmodel: CalQ, CalU = solution[0], solution[1]
mjdSec, Az, El, dAz, dEl = np.ones([0]), np.ones([0]), np.ones([0]), np.ones([0]), np.ones([0])
chNum, chWid, Freq = GetChNum(msfile, spw); Freq = Freq* 1.0e-9
##-------- Antenna List
antList = GetAntName(msfile)
flagAnt = indexList(antFlag, antList)
interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spw, scan); timeNum = len(timeStamp)
trkAnt, scnAnt, scanTime, AzElOffset = antRefScan(msfile, [min(timeStamp), max(timeStamp)])
trkAnt, scnAnt  = list(set(trkAnt) - set(flagAnt)), list(set(scnAnt) - set(flagAnt))
if refantName not in antList[trkAnt]:
    print refantName + ' does not exist in this MS.'
    sys.exit()
#
refantID = np.where( antList == refantName)[0][0]
trkAntMap = [refantID] + np.array(trkAnt)[range(trkAnt.index(refantID))].tolist() + np.array(trkAnt)[range(trkAnt.index(refantID)+1, len(trkAnt))].tolist()
antMap = trkAntMap + scnAnt
antNum, trkAntNum, scnAntNum = len(antMap), len(trkAntMap), len(scnAnt)
blNum, trkBlNum = antNum* (antNum - 1) / 2, trkAntNum* (trkAntNum - 1) / 2
blMap, blInv = range(blNum), [False]* blNum
#-------- Baseline Indexing
ant0, ant1, antWeight = ANT0[0:blNum], ANT1[0:blNum], np.ones([antNum])
antWeight[range(trkAntNum, antNum)] = 0.5
for bl_index in range(blNum):    blMap[bl_index], blInv[bl_index] = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
blWeight = antWeight[ant0]* antWeight[ant1]
#-------- BP table
BP_ant = BP_ant[indexList(antList[antMap], BPantList)]
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
#-------- scan pattern
AntD = np.zeros([antNum])
for ant_index in range(antNum): AntD[ant_index] = GetAntD(antList[antMap[ant_index]])
scanTime, AntID, AZ, EL = GetAzEl(msfile)
azelTime_index = np.where( AntID == refantID )[0].tolist()
timeThresh = np.median( np.diff( scanTime[azelTime_index]))
FWHM = GetFWHM(msfile, spw, AntD)      # FWHM in arcsec
centerIndex = scanThresh(msfile, scnAnt[0], FWHM[scnAnt[0]]/10.0); centerTime = scanTime[centerIndex]
matchNum  = np.zeros([timeNum])
for time_index in range(timeNum):
    matchNum[time_index] = timeMatch(timeStamp[time_index], centerTime, np.median(interval))
    scanAz, scanEl = AzElMatch(timeStamp[time_index], scanTime, AntID, refantID, timeThresh, AZ, EL)
    diffAz, diffEl = AzElMatch(timeStamp[time_index], scanTime, AntID, scnAnt[0], timeThresh, AzElOffset[0], AzElOffset[1])
    Az, El = np.append(Az, scanAz), np.append(El, scanEl)
    dAz, dEl = np.append(dAz, diffAz), np.append(dEl, diffEl) 
#
onAxisIndex = np.where( matchNum > 0 )[0].tolist()
#-------- Visibility self-cal
print '-- Loading visibility data %s SPW=%d ...' % (prefix, spw)
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
print '---- Bandpass cal'
tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3, 2, 0, 1)
Xspec = (tempSpec / (BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate())).transpose(2,3,1,0)
print '---- XY delay cal'
XYdlSpec = delay_cal( np.ones([chNum], dtype=complex), XYdelay )
Xspec[1] = (Xspec[1].transpose(1,2,0)* XYdlSpec).transpose(2,0,1)
Xspec[2] = (Xspec[2].transpose(1,2,0)* XYdlSpec.conjugate()).transpose(2,0,1)
print '---- Antenna-based gain correction'
chAvgVis = np.mean(Xspec[:,chRange], axis=1)
PA = AzEl2PA(Az, El, ALMA_lat) - BANDPA
PA = np.arctan2( np.sin(PA), np.cos(PA))
GainX, GainY = polariGain(chAvgVis[0], chAvgVis[3], PA, CalQ, CalU)
Gain = np.array([GainX, GYtwiddle* GainY])
CaledXspec = (Xspec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())).transpose(1,0,2,3)
#-------- D-term of tracking antennas
Dx, Dy = np.zeros([antNum, timeNum, chNum], dtype=complex), np.zeros([antNum, timeNum, chNum], dtype=complex)
print('-------- Determining Antenna-based D-terms (refants) ----')
PAwidth = 0.01; PAsegNum = int((max(PA) - min(PA))/PAwidth)
if max(np.diff(PA)) > np.pi: PAsegNum = int((max(np.sin(PA)) - min(np.sin(PA)))/PAwidth) # for +/-pi gap
if PAsegNum > timeNum/2: PAsegNum = timeNum/2
#
trkDx, trkDy = np.zeros([trkAntNum, PAsegNum, chNum], dtype=complex), np.zeros([trkAntNum, PAsegNum, chNum], dtype=complex)
for seg_index in range(PAsegNum):
    progress = (seg_index + 1.0) / PAsegNum
    sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
    timeIndexRange = range( (seg_index* timeNum/PAsegNum), ((seg_index + 1)* timeNum/PAsegNum) )
    PS = np.dot(PAMatrix(np.mean(PA[timeIndexRange])), np.array([1.0, CalQ, CalU, 0.0])).real
    for ch_index in range(chNum):
        VisTime = np.mean(CaledXspec[:, ch_index][:,range(trkBlNum)][:,:,timeIndexRange], axis=2)
        trkDx[:,seg_index,ch_index], trkDy[:,seg_index,ch_index]  = Vis2solveDD(VisTime, PS)
    #
#
sys.stderr.write('\n'); sys.stderr.flush()
DxMean, DyMean = np.mean(trkDx, axis=1), np.mean(trkDy, axis=1)
Dx[range(trkAntNum)] = (Dx[range(trkAntNum)].transpose(1,0,2) + DxMean).transpose(1,0,2)
Dy[range(trkAntNum)] = (Dy[range(trkAntNum)].transpose(1,0,2) + DyMean).transpose(1,0,2)
#-------- Record on-axis D-term spectrum
logfile = open(prefix + '-SPW' + `spw` + '-TrkDtermSpec.log', 'w')
text_sd = 'ant ch ReDx ImDx ReDy ImDy'
logfile.write(text_sd + '\n')
for ant_index in range(trkAntNum):
    for ch_index in range(chNum):
        text_sd = '%s %d %8.6f %8.6f %8.6f %8.6f' % (antList[trkAntMap[ant_index]], ch_index, DxMean[ant_index, ch_index].real, DxMean[ant_index, ch_index].imag, DyMean[ant_index, ch_index].real, DyMean[ant_index, ch_index].imag)
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
            Minv = InvMullerMatrix(DxMean[ant1[bl_index], ch_index], DyMean[ant1[bl_index], ch_index], DxMean[ant0[bl_index], ch_index], DyMean[ant0[bl_index], ch_index])
            StokesVis[bl_index, seg_index, ch_index] = np.dot(Pinv, np.dot(Minv, np.mean(CaledXspec[:,ch_index][:,bl_index][:,timeIndexRange], axis=1)))
        #
    #
#
#-------- Record Stokes Fluxes measured by tracking antennas
StokesFlux = np.mean(StokesVis[:, :, chRange], axis=(0,1,2)).real
text_sd = 'TrkMeas / Model: I=%6.4f / %6.4f  Q=%6.4f / %6.4f U=%6.4f / %6.4f V=%6.4f / %6.4f' % (StokesFlux[0], 1.0, StokesFlux[1], CalQ, StokesFlux[2], CalU, StokesFlux[3], 0.0); print text_sd
UCmQS = StokesFlux[2]* np.cos(2.0*PA) - StokesFlux[1]* np.sin(2.0*PA)   # U cos - Q sin
QCpUS = StokesFlux[1]* np.cos(2.0*PA) + StokesFlux[2]* np.sin(2.0*PA)   # Q cos + U sin
logfile = open(prefix + '-SPW' + `spw` + '-trkStokes.log', 'w')
text_sd = 'CH I Q U V'; logfile.write(text_sd + '\n')
for ch_index in range(chNum):
    text_sd = '%d %8.6f %8.6f %8.6f %8.6f' % (ch_index, np.mean(StokesVis[:,:,ch_index,0]).real, np.mean(StokesVis[:,:,ch_index,1]).real, np.mean(StokesVis[:,:,ch_index,2]).real, np.mean(StokesVis[:,:,ch_index,3]).real)
    logfile.write(text_sd + '\n')
#
logfile.close()
#-------- Determination of D-terms in scanning antennas
print('-------- Determining Antenna-based D-terms (scan ants) ----')
for ant_index in range(scnAntNum):
    antID = trkAntNum + ant_index
    print 'Determining D-term of ' + antList[antMap[antID]]
    TrkScnBL = range(antID* (antID - 1) / 2, antID* (antID - 1) / 2 + trkAntNum)
    for time_index in range(timeNum):
        PS = np.dot(PAMatrix(PA[time_index]), np.array([1.0, StokesFlux[1], StokesFlux[2], 0.0])).real
        for ch_index in range(chNum):
            progress = (time_index* chNum + ch_index + 1.0) / (timeNum* chNum)
            sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
            VisTime = CaledXspec[:, ch_index, TrkScnBL, time_index]
            Dx[antID, time_index, ch_index], Dy[antID, time_index, ch_index] = Vis2solveD(VisTime, Dx[0:trkAntNum, time_index, ch_index], Dy[0:trkAntNum, time_index, ch_index], PS )
        #
    #
    sys.stderr.write('\n'); sys.stderr.flush()
#
trkBlIndex  = np.where(blWeight == 1.0)[0].tolist(); trkBlNum  = len(trkBlIndex)        # Ref-Ref baselines
ScTrBlIndex = np.where(blWeight == 0.5)[0].tolist(); ScTrBlNum = len(ScTrBlIndex)       # Ref-Scan baselines
ScScBlIndex = np.where(blWeight == 0.25)[0].tolist(); ScScBlNum = len(ScScBlIndex)      # Scan-Scan baselines
#
#-------- Plot D-term spectrum for beam position
logfile = open(prefix + '-SPW' + `spw` + '-DtermSpec.log', 'w')
text_sd = 'ant beamoff branch ch ReDx ImDx ReDy ImDy'
logfile.write(text_sd + '\n')
OffBeam = np.sqrt(dAz**2 + dEl**2)
ScanInterval = median( np.diff( dAz[0:50] ))
SortOffBeam = np.sort( OffBeam )
BreakIndex = np.where( np.diff(SortOffBeam) > 0.5* ScanInterval)[0]
thresh = 0.5*( SortOffBeam[BreakIndex] + SortOffBeam[BreakIndex + 1])
#thresh = np.r_[0.0, np.linspace( min(np.sqrt(ScanAz**2 + ScanEl**2)), max(np.sqrt(ScanAz**2 + ScanEl**2)), num=16) + min(np.sqrt(ScanAz**2 + ScanEl**2))]
Dist2 = dAz**2 + dEl**2
xrange, yrange = [min(Freq[chRange]), max(Freq[chRange])], [-0.1, 0.1]
for ant_index in range(scnAntNum):
    antID = scnAnt[ant_index]
    DantID = trkAntNum + ant_index
    fwhm = FWHM[DantID]
    for thresh_index in range(6):
        time_index = list(set(np.where(Dist2 > thresh[thresh_index]**2 )[0]) & set(np.where(Dist2 < thresh[thresh_index + 1]**2 )[0]))
        fig = plt.figure(thresh_index, figsize = (8,11))
        fig.text(0.45, 0.05, 'Frequency [GHz]')
        fig.text(0.05, 0.45, 'D-term', rotation=90)
        plt.suptitle(prefix + ' ' + antList[antID] + ' D-term@' + `round(np.median(np.sqrt(Dist2[time_index])),1)` + ' arcsec')
        #-------- Plot beam map
        plt.subplot2grid( (9,6), (0,5), aspect=1)
        plt.plot(dAz, dEl, ',', color='k', alpha=0.1)
        plt.plot( dAz[time_index], dEl[time_index], '.', color='r')
        circle_x, circle_y = circlePoints(0, 0, fwhm/2); plt.plot( circle_x, circle_y, color='green' )
        circle_x, circle_y = circlePoints(0, 0, fwhm/sqrt(2)); plt.plot( circle_x, circle_y, color='blue' )
        plt.axis([-fwhm, fwhm, -fwhm, fwhm], fontsize=3)
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
        plt.savefig( prefix + '-' + antList[antID] + '-SPW' + `spw` + '-OFF' + `round(np.median(np.sqrt(Dist2[time_index])),1)` + '-DtermSpec.pdf', form='pdf'); plt.close()
    #
#
logfile.close()
#-------- Plot channel-averaged D-terms of scanning antennas
print('-------- Plot D-term Maps for scan ants ----')
chAvgDx, chAvgDy = np.mean(Dx[:,:,chRange], axis=2), np.mean(Dy[:,:,chRange], axis=2)
for ant_index in range(scnAntNum):
    antID = scnAnt[ant_index]
    DantID = trkAntNum + ant_index
    fwhm = FWHM[DantID]
    #-------- Plot
    fig = plt.figure( figsize = (10,10))
    fig.suptitle(prefix + ' ' + antList[antID] + ' SPW=' + `spw` + ' Scan=' + `scan`)
    fig.text(0.45, 0.05, 'Az Offset [arcsec]')
    fig.text(0.05, 0.45, 'El Offset [arcsec]', rotation=90)
    #
    xi, yi = np.mgrid[ -floor(2.0*fwhm):floor(2.0*fwhm):128j, -floor(2.0*fwhm):floor(2.0*fwhm):128j]
    IndexCenter = np.where( dAz**2 + dEl**2 < 0.005* fwhm**2 )[0]
    Index3dB = np.where( dAz**2 + dEl**2 < 0.25* fwhm**2 )[0]
    Index6dB = np.where( dAz**2 + dEl**2 < 0.5* fwhm**2 )[0]
    #-------- Save Gain and D-term to logfile
    logfile = open(prefix + '-' + antList[antID] + '-SPW' + `spw` + '-beamGainD.log', 'w')
    GXA, GXP = smoothGain(timeStamp[IndexCenter], Gain[0, DantID, IndexCenter])
    GYA, GYP = smoothGain(timeStamp[IndexCenter], Gain[1, DantID, IndexCenter])
    NormGX = Gain[0, DantID] / GXA(timeStamp) * np.exp( (0.0 - 1.0j)* GXP(timeStamp))
    NormGY = Gain[1, DantID] / GYA(timeStamp) * np.exp( (0.0 - 1.0j)* GYP(timeStamp))
    text_sd = '#Antenna-based Complex Gain and D-term : %s %s %9.4f GHz' % (prefix, antList[antID], np.median(Freq)) ; logfile.write(text_sd + '\n')
    text_sd = '#dAZ   dEL    ReGX    ImGX     ReGY     ImGY     ReDX     ImDX      RdDY     ImDY   '; logfile.write(text_sd + '\n')
    text_sd = '#-----------------------------------------------------------------------------------'; logfile.write(text_sd + '\n')
    for time_index in range(timeNum):
        text_sd = '%4.1f %4.1f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f' % (dAz[time_index], dEl[time_index], NormGX[time_index].real, NormGX[time_index].imag, NormGY[time_index].real, NormGY[time_index].imag, chAvgDx[DantID,time_index].real, chAvgDx[DantID,time_index].imag, chAvgDy[DantID,time_index].real, chAvgDy[DantID,time_index].imag)
        logfile.write(text_sd + '\n')
    #
    logfile.close()
    Dx3dB, Dy3dB = chAvgDx[DantID][Index3dB], chAvgDy[DantID][Index3dB]
    Dx6dB, Dy6dB = chAvgDx[DantID][Index6dB], chAvgDy[DantID][Index6dB]
    ReDxmap = GridData( chAvgDx[DantID].real, dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), fwhm/16).reshape(len(xi), len(xi))
    ImDxmap = GridData( chAvgDx[DantID].imag, dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), fwhm/16).reshape(len(xi), len(xi))
    ReDymap = GridData( chAvgDy[DantID].real, dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), fwhm/16).reshape(len(xi), len(xi))
    ImDymap = GridData( chAvgDy[DantID].imag, dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), fwhm/16).reshape(len(xi), len(xi))
    #---- plot Re(Dx)
    plt.subplot( 2, 2, 1, aspect=1); plt.contourf(xi, yi, ReDxmap, np.linspace(-0.10, 0.10, 11)); plt.colorbar(); plt.title('Re(Dx)')
    circle_x, circle_y = circlePoints(0, 0, fwhm/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, fwhm/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Re(Dx) at Center = %5.3f' % ( np.mean(Dx[DantID, IndexCenter].real) ); plt.text(-1.6*fwhm, -1.5*fwhm, text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( max(Dx3dB.real), min(Dx3dB.real) ); plt.text(-1.6*fwhm, -1.7*fwhm, text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( max(Dx6dB.real), min(Dx6dB.real) ); plt.text(-1.6*fwhm, -1.9*fwhm, text_sd, size='x-small')
    #---- plot Im(Dx)
    plt.subplot( 2, 2, 2, aspect=1); plt.contourf(xi, yi, ImDxmap, np.linspace(-0.10, 0.10, 11)); plt.colorbar(); plt.title('Im(Dx)')
    circle_x, circle_y = circlePoints(0, 0, fwhm/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, fwhm/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Im(Dx) at Center = %5.3f' % ( np.mean(Dx[DantID, IndexCenter].imag) ); plt.text(-1.6*fwhm, -1.5*fwhm, text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( max(Dx3dB.imag), min(Dx3dB.imag) ); plt.text(-1.6*fwhm, -1.7*fwhm, text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( max(Dx6dB.imag), min(Dx6dB.imag) ); plt.text(-1.6*fwhm, -1.9*fwhm, text_sd, size='x-small')
    #---- plot Re(Dy)
    plt.subplot( 2, 2, 3, aspect=1); plt.contourf(xi, yi, ReDymap, np.linspace(-0.10, 0.10, 11)); plt.colorbar(); plt.title('Re(Dy)')
    circle_x, circle_y = circlePoints(0, 0, fwhm/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, fwhm/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Re(Dy) at Center = %5.3f' % ( np.mean(Dy[DantID, IndexCenter].real) ); plt.text(-1.6*fwhm, -1.5*fwhm, text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( max(Dy3dB.real), min(Dy3dB.real) ); plt.text(-1.6*fwhm, -1.7*fwhm, text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( max(Dy6dB.real), min(Dy6dB.real) ); plt.text(-1.6*fwhm, -1.9*fwhm, text_sd, size='x-small')
    #---- plot Im(Dy)
    plt.subplot( 2, 2, 4, aspect=1); plt.contourf(xi, yi, ImDymap, np.linspace(-0.10, 0.10, 11)); plt.colorbar(); plt.title('Im(Dy)')
    circle_x, circle_y = circlePoints(0, 0, fwhm/2); plt.plot( circle_x, circle_y )
    circle_x, circle_y = circlePoints(0, 0, fwhm/sqrt(2)); plt.plot( circle_x, circle_y )
    text_sd = 'Im(Dy) at Center = %5.3f' % ( np.mean(Dy[DantID, IndexCenter].imag) ); plt.text(-1.6*fwhm, -1.5*fwhm, text_sd, size='x-small')
    text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( max(Dy3dB.imag), min(Dy3dB.imag) ); plt.text(-1.6*fwhm, -1.7*fwhm, text_sd, size='x-small')
    text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( max(Dy6dB.imag), min(Dy6dB.imag) ); plt.text(-1.6*fwhm, -1.9*fwhm, text_sd, size='x-small')
    plt.plot( dAz, dEl, '.', color='k', alpha=0.1)
    plt.axis([-2.0*fwhm, 2.0*fwhm, -2.0*fwhm, 2.0*fwhm])
    plt.savefig( prefix + '-' + antList[antID] + '-SPW' + `spw` + '-DtermMap.pdf', form='pdf'); plt.close()
    plt.close()
#
#-------- D-term-corrected Stokes parameters --------
print('-------- D-term-corrected Stokes parameters ----')
StokesVis = np.zeros([ScScBlNum, timeNum, chNum, 4], dtype=complex) 
for time_index in range(timeNum):
    progress = (time_index + 1.0) / timeNum
    sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
    Pinv = InvPAMatrix( PA[time_index] )
    for bl_index in range(ScScBlNum):
        BlID = ScScBlIndex[bl_index]
        for ch_index in range(chNum):
            Minv = InvMullerMatrix( Dx[ant1[BlID], time_index, ch_index], Dy[ant1[BlID], time_index, ch_index], Dx[ant0[BlID], time_index, ch_index], Dy[ant0[BlID], time_index, ch_index])
            StokesVis[bl_index, time_index, ch_index] = np.dot(Pinv, np.dot(Minv, CaledXspec[:,ch_index, BlID, time_index]))
        #
    #
#
sys.stderr.write('\n'); sys.stderr.flush()
ScnStokesSpec = np.mean(StokesVis, axis=0).real
ScnStokes = np.mean(ScnStokesSpec[:,chRange], axis=1)
Qerr = ScnStokes[:,1] - StokesFlux[1]
Uerr = ScnStokes[:,2] - StokesFlux[2]
Perr = sqrt(ScnStokes[:,1]**2 + ScnStokes[:,2]**2) - sqrt(StokesFlux[1]**2 + StokesFlux[2]**2)
Aerr = np.arctan( (ScnStokes[:,2]* StokesFlux[1] - ScnStokes[:,1]* StokesFlux[2]) / (ScnStokes[:,1]* StokesFlux[1] + ScnStokes[:,2]* StokesFlux[2]) )* 90.0/math.pi

#-------- Plot Stokes Beam Map
logfile = open(prefix + '-SPW' + `spw` + '-Stokes.log', 'w')
text_sd = 'I I3max I3min I6max I6min Q Q3max Q3min Q6max Q6min U U3max U3min U6max U6min V V3max V3min V6max V6min'
logfile.write(text_sd + '\n')
fig = plt.figure( figsize = (10,10))
fig.text(0.45, 0.05, 'Az Offset [arcsec]')
fig.text(0.05, 0.45, 'El Offset [arcsec]', rotation=90)
fig.text(0.45, 0.95, prefix)
mapW = np.max(FWHM)
xi, yi = np.mgrid[ -floor(mapW):floor(mapW):128j, -floor(mapW):floor(mapW):128j]
Imap = GridData( ScnStokes[:,0], dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), mapW/16).reshape(len(xi), len(xi))
Qmap = GridData( ScnStokes[:,1], dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), mapW/16).reshape(len(xi), len(xi))
Umap = GridData( ScnStokes[:,2], dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), mapW/16).reshape(len(xi), len(xi))
Vmap = GridData( ScnStokes[:,3], dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), mapW/16).reshape(len(xi), len(xi))
QEmap = GridData(Qerr, dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), mapW/16).reshape(len(xi), len(xi))
UEmap = GridData(Uerr, dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), mapW/16).reshape(len(xi), len(xi))
PEmap = GridData(Perr, dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), mapW/16).reshape(len(xi), len(xi))
AEmap = GridData(Aerr, dAz, dEl, xi.reshape(xi.size), yi.reshape(xi.size), mapW/16).reshape(len(xi), len(xi))
IndexCenter = np.where( dAz**2 + dEl**2 < 0.005* mapW**2 )[0]
Index3dB = np.where( dAz**2 + dEl**2 < 0.25* mapW**2 )[0]
Index6dB = np.where( dAz**2 + dEl**2 < 0.5* mapW**2 )[0]
#---- Plot I map
plt.subplot(2, 2, 1, aspect=1); plt.contourf(xi, yi, Imap, np.linspace(0.9, 1.1, 21)); plt.colorbar(); plt.title('Stokes I')
circle_x, circle_y = circlePoints(0, 0, mapW/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, mapW/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'I at Center   = %5.3f' % ( np.mean(ScnStokes[IndexCenter, 0]) ); plt.text(-0.8*mapW, -0.75*mapW, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnStokes[Index3dB, 0]), np.min(ScnStokes[Index3dB, 0]) ); plt.text(-0.8*mapW, -0.85*mapW, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnStokes[Index6dB, 0]), np.min(ScnStokes[Index6dB, 0]) ); plt.text(-0.8*mapW, -0.95*mapW, text_sd, size='x-small')
#---- Plot Q map
plt.subplot(2, 2, 2, aspect=1); plt.contourf(xi, yi, Qmap, np.linspace(0.01*(floor(StokesFlux[1]* 100)-5), 0.01*(floor(StokesFlux[1]* 100)+5), 21)); plt.colorbar(); plt.title('Stokes Q')
circle_x, circle_y = circlePoints(0, 0, mapW/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, mapW/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Q at Center   = %5.3f' % ( np.mean(ScnStokes[IndexCenter, 1]) ); plt.text(-0.8*mapW, -0.75*mapW, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnStokes[Index3dB, 1]), np.min(ScnStokes[Index3dB, 1]) ); plt.text(-0.8*mapW, -0.85*mapW, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnStokes[Index6dB, 1]), np.min(ScnStokes[Index6dB, 1]) ); plt.text(-0.8*mapW, -0.95*mapW, text_sd, size='x-small')
#---- Plot U map
plt.subplot(2, 2, 3, aspect=1); plt.contourf(xi, yi, Umap, np.linspace(0.01*(floor(StokesFlux[2]* 100)-5), 0.01*(floor(StokesFlux[2]* 100)+5), 21)); plt.colorbar(); plt.title('Stokes U')
circle_x, circle_y = circlePoints(0, 0, mapW/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, mapW/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'U at Center   = %5.3f' % ( np.mean(ScnStokes[IndexCenter, 2]) ); plt.text(-0.8*mapW, -0.75*mapW, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnStokes[Index3dB, 2]), np.min(ScnStokes[Index3dB, 2]) ); plt.text(-0.8*mapW, -0.85*mapW, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnStokes[Index6dB, 2]), np.min(ScnStokes[Index6dB, 2]) ); plt.text(-0.8*mapW, -0.95*mapW, text_sd, size='x-small')
#---- Plot V map
plt.subplot(2, 2, 4, aspect=1); plt.contourf(xi, yi, Vmap, np.linspace(-0.05, 0.05, 21)); plt.colorbar(); plt.title('Stokes V')
circle_x, circle_y = circlePoints(0, 0, mapW/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, mapW/sqrt(2));   plt.plot( circle_x, circle_y )
plt.plot( dAz, dEl, '.', color='k', alpha=0.1)
text_sd = 'V at Center   = %5.3f' % ( np.mean(ScnStokes[IndexCenter, 3]) ); plt.text(-0.8*mapW, -0.75*mapW, text_sd, size='x-small')
text_sd = '(max,min)_3dB = (%5.3f %5.3f) ' % ( np.max(ScnStokes[Index3dB, 3]), np.min(ScnStokes[Index3dB, 3]) ); plt.text(-0.8*mapW, -0.85*mapW, text_sd, size='x-small')
text_sd = '(max,min)_6dB = (%5.3f %5.3f) ' % ( np.max(ScnStokes[Index6dB, 3]), np.min(ScnStokes[Index6dB, 3]) ); plt.text(-0.8*mapW, -0.95*mapW, text_sd, size='x-small')
#
plt.axis([-mapW, mapW, -mapW, mapW])
plt.savefig( prefix + '-SPW' + `spw` + '-StokesMap.pdf', form='pdf'); plt.close()
plt.close()
text_sd = '%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f ' % (np.mean(ScnStokes[IndexCenter, 0]), np.max(ScnStokes[Index3dB, 0]), np.min(ScnStokes[Index3dB, 0]), np.max(ScnStokes[Index6dB, 0]), np.min(ScnStokes[Index6dB, 0]), np.mean(ScnStokes[IndexCenter, 1]), np.max(ScnStokes[Index3dB, 1]), np.min(ScnStokes[Index3dB, 1]), np.max(ScnStokes[Index6dB, 1]), np.min(ScnStokes[Index6dB, 1]), np.mean(ScnStokes[IndexCenter, 2]), np.max(ScnStokes[Index3dB, 2]), np.min(ScnStokes[Index3dB, 2]), np.max(ScnStokes[Index6dB, 2]), np.min(ScnStokes[Index6dB, 2]), np.mean(ScnStokes[IndexCenter, 3]), np.max(ScnStokes[Index3dB, 3]), np.min(ScnStokes[Index3dB, 3]), np.max(ScnStokes[Index6dB, 3]), np.min(ScnStokes[Index6dB, 3]))
logfile.write(text_sd + '\n')
logfile.close()
#
logfile = open(prefix + '-SPW' + `spw` + '-QUerr.log', 'w')
text_sd = 'Qerr Qerr3max Qerr6max Uerr Uerr3max Uerr6max Perr Perr3max Perr6max EVPAerr EVPAerr3max EVPAerr6max' 
logfile.write(text_sd + '\n')
#---- Plot Q err map
fig = plt.figure( figsize = (10,10))
fig.text(0.45, 0.05, 'Az Offset [arcsec]')
fig.text(0.05, 0.45, 'El Offset [arcsec]', rotation=90)
fig.text(0.45, 0.95, prefix)
plt.subplot(2, 2, 1, aspect=1); plt.contourf(xi, yi, 100.0*QEmap, np.linspace(-5.0, 5.0, 11)); plt.colorbar(); plt.title('Q_err')
circle_x, circle_y = circlePoints(0, 0, mapW/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, mapW/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Qerr at Center       = %4.1f %%' % ( 100.0* np.mean(Qerr[IndexCenter])); plt.text(-0.8*mapW, -0.7*mapW, text_sd, size='x-small')
text_sd = 'Max Qerr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(Qerr[Index3dB])) ); plt.text(-0.8*mapW, -0.8*mapW, text_sd, size='x-small')
text_sd = 'Max Qerr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(Qerr[Index6dB])) ); plt.text(-0.8*mapW, -0.9*mapW, text_sd, size='x-small')
#----Plot U err map
plt.subplot(2, 2, 2, aspect=1); plt.contourf(xi, yi, 100.0*UEmap, np.linspace(-5.0, 5.0, 11)); plt.colorbar(); plt.title('U_err')
circle_x, circle_y = circlePoints(0, 0, mapW/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, mapW/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Uerr at Center       = %4.1f %%' % ( 100.0* np.mean(Uerr[IndexCenter])); plt.text(-0.8*mapW, -0.7*mapW, text_sd, size='x-small')
text_sd = 'Max Uerr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(Uerr[Index3dB])) ); plt.text(-0.8*mapW, -0.8*mapW, text_sd, size='x-small')
text_sd = 'Max Uerr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(Uerr[Index6dB])) ); plt.text(-0.8*mapW, -0.9*mapW, text_sd, size='x-small')
#----Plot P error
plt.subplot( 2, 2, 3, aspect=1); plt.contourf(xi, yi, 100.0*PEmap, np.linspace(-5.0, 5.0, 11)); plt.colorbar(); plt.title('P_err')
circle_x, circle_y = circlePoints(0, 0, mapW/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, mapW/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'Perr at Center       = %4.1f %%' % ( 100.0* np.mean(Perr[IndexCenter])); plt.text(-0.8*mapW, -0.75*mapW, text_sd, size='x-small')
text_sd = 'Max Perr (-3dB beam) = %4.1f %%' % ( 100.0* max(abs(Perr[Index3dB])) ); plt.text(-0.8*mapW, -0.85*mapW, text_sd, size='x-small')
text_sd = 'Max Perr (-6dB beam) = %4.1f %%' % ( 100.0* max(abs(Perr[Index6dB])) ); plt.text(-0.8*mapW, -0.95*mapW, text_sd, size='x-small')
#----Plot EVPA error
plt.subplot( 2, 2, 4, aspect=1); plt.contourf(xi, yi, AEmap, np.linspace(-5.0, 5.0, 11)); plt.colorbar(); plt.title('EVPA_error')
circle_x, circle_y = circlePoints(0, 0, mapW/2); plt.plot( circle_x, circle_y )
circle_x, circle_y = circlePoints(0, 0, mapW/sqrt(2));   plt.plot( circle_x, circle_y )
text_sd = 'EVPAerr at Center       = %4.1f deg' % ( np.mean(Aerr[IndexCenter]) ); plt.text(-0.8*mapW, -0.75*mapW, text_sd, size='x-small')
text_sd = 'Max EVPAerr (-3dB beam) = %4.1f deg' % ( max(abs(Aerr[Index3dB])) ); plt.text(-0.8*mapW, -0.85*mapW, text_sd, size='x-small')
text_sd = 'Max EVPAerr (-6dB beam) = %4.1f deg' % ( max(abs(Aerr[Index6dB])) ); plt.text(-0.8*mapW, -0.95*mapW, text_sd, size='x-small')
#
plt.axis([-mapW, mapW, -mapW, mapW])
plt.savefig( prefix + '-SPW' + `spw` + '-StokesErr.pdf', form='pdf'); plt.close()
plt.close()
text_sd = '%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f' % (np.mean(Qerr[IndexCenter]), max(abs(Qerr[Index3dB])), max(abs(Qerr[Index6dB])), np.mean(Uerr[IndexCenter]), max(abs(Uerr[Index3dB])), max(abs(Uerr[Index6dB])), np.mean(Perr[IndexCenter]), max(abs(Perr[Index3dB])), max(abs(Perr[Index6dB])), np.mean(Aerr[IndexCenter]), max(abs(Aerr[Index3dB])), max(abs(Aerr[Index6dB])))
logfile.write(text_sd + '\n')
logfile.close()
#-------- Plot Stokes Spectra at Beam Position
logfile = open(prefix + '-SPW' + `spw` + '-StokesSpec.log', 'w')
text_sd = 'beamoff branch CH I Q U V ';  logfile.write(text_sd + '\n')
text_sd = '%4.1f %d %8.6f %8.6f %8.6f %8.6f' % (0.0,  0, StokesFlux[0], StokesFlux[1], StokesFlux[2], StokesFlux[3]); logfile.write(text_sd + '\n')
xrange, yrange = [min(Freq[chRange]), max(Freq[chRange])], [-0.01, 0.01]
for thresh_index in range(6):
    time_index = list(set(np.where(Dist2 > thresh[thresh_index]**2 )[0]) & set(np.where(Dist2 < thresh[thresh_index + 1]**2 )[0]))
    fig = plt.figure(thresh_index, figsize = (8,11))
    fig.text(0.45, 0.05, 'Frequency [GHz]')
    fig.text(0.05, 0.45, 'Stokes Residual [scaled by Stokes I]', rotation=90)
    plt.suptitle(prefix + ' ' + ' Stokes Residuals@' + `round(np.median(np.sqrt(Dist2[time_index])),1)` + ' arcsec')
    #-------- Plot beam map
    plt.subplot2grid( (9,6), (0,5), aspect=1)
    plt.plot(dAz, dEl, ',', color='k', alpha=0.1)
    plt.plot( dAz[time_index], dEl[time_index], '.', color='r')
    circle_x, circle_y = circlePoints(0, 0, mapW/2); plt.plot( circle_x, circle_y, color='green' )
    circle_x, circle_y = circlePoints(0, 0, mapW/sqrt(2)); plt.plot( circle_x, circle_y, color='blue' )
    plt.axis([-mapW, mapW, -mapW, mapW])
    plt.tick_params(labelsize = 6)
    #
    #-------- Plot Mean IQUV
    plt.subplot2grid( (9,6), (0,0), colspan=4)
    fill_color='g'
    plt.fill([xrange[0], xrange[0], xrange[1], xrange[1]], [yrange[1], yrange[0], yrange[0], yrange[1]], fill_color, alpha=0.1)
    plt.plot( Freq, np.mean(ScnStokesSpec[time_index, :, 0] - StokesFlux[0], axis=0), color='k', label = 'I', ls='steps-mid')
    plt.plot( Freq, np.mean(ScnStokesSpec[time_index, :, 1] - StokesFlux[1], axis=0), color='r', label = 'Q', ls='steps-mid')
    plt.plot( Freq, np.mean(ScnStokesSpec[time_index, :, 2] - StokesFlux[2], axis=0), color='g', label = 'U', ls='steps-mid')
    plt.plot( Freq, np.mean(ScnStokesSpec[time_index, :, 3] - StokesFlux[3], axis=0), color='b', label = 'V', ls='steps-mid')
    plt.axis([min(Freq), max(Freq), -0.02,0.02], fontsize=3)
    plt.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    plt.tick_params(labelsize = 6)
    #
    for index in range(48):
        fill_color='g'
        if max(abs(ScnStokesSpec[time_index[index], chRange, 1] - StokesFlux[1])) > yrange[1]: fill_color = 'r'
        if max(abs(ScnStokesSpec[time_index[index], chRange, 2] - StokesFlux[2])) > yrange[1]: fill_color = 'r'
        if max(abs(ScnStokesSpec[time_index[index], chRange, 3] - StokesFlux[3])) > yrange[1]: fill_color = 'r'
        plt.subplot2grid( (9,6), (int(index/6)+1, index%6))
        plt.fill([xrange[0], xrange[0], xrange[1], xrange[1]], [yrange[1], yrange[0], yrange[0], yrange[1]], fill_color, alpha=0.1)
        plt.plot( Freq, ScnStokesSpec[time_index[index], :, 0] - StokesFlux[0], ls='steps-mid', color='k',     label = 'I')
        plt.plot( Freq, ScnStokesSpec[time_index[index], :, 1] - StokesFlux[1], ls='steps-mid', color='r', label = 'Q')
        plt.plot( Freq, ScnStokesSpec[time_index[index], :, 2] - StokesFlux[2], ls='steps-mid', color='g',  label = 'U')
        plt.plot( Freq, ScnStokesSpec[time_index[index], :, 3] - StokesFlux[3], ls='steps-mid', color='b',  label = 'V')
        plt.axis([min(Freq), max(Freq), -0.05,0.05], fontsize=3)
        plt.tick_params(labelsize = 6)
        for ch_index in range(chNum):
            text_sd = '%4.1f %d %d %8.6f %8.6f %8.6f %8.6f' % (np.median(np.sqrt(Dist2[time_index])), index, ch_index, ScnStokesSpec[time_index[index], ch_index, 0], ScnStokesSpec[time_index[index], ch_index, 1], ScnStokesSpec[time_index[index], ch_index, 2], ScnStokesSpec[time_index[index], ch_index, 3]); logfile.write(text_sd + '\n')
        #
    #
    plt.savefig( prefix + '-' + '-SPW' + `spw` + '-OFF' + `round(np.median(np.sqrt(Dist2[time_index])),1)` + '-StokesSpec.pdf', form='pdf'); plt.close()
#
logfile.close()
