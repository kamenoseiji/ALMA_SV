execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#
#-------- Scanning Offset < threshold
def scanThresh( msfile, ant_index, thresh ):
    Time, AntID, Offset = GetAzOffset(msfile)
    time_index = np.where( AntID == ant_index )[0]
    onAxisIndex = np.where( Offset[0, time_index]**2 + Offset[1, time_index]**2 < thresh**2 )[0]
    return time_index[onAxisIndex].tolist()
    #return onAxisIndex.tolist()
#
#-------- Time-based matching between time tags in visibilities and in scan pattern 
#def AzElMatch( refTime, scanTime, thresh, Az, El ):
#    index = np.where( abs(scanTime - refTime) < thresh)[0]
#    return np.median(Az[index]), np.median(El[index])
#
def timeMatch( refTime, scanTime, thresh):
    match = np.where( abs(scanTime - refTime) < thresh)[0].tolist()
    return len(match)
#
#----------------------------------------- Procedures
#mjdSec, Az, El, PA = [], [], []. []
msfile = wd + prefix + '.ms'
msmd.open(msfile)
#-------- Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
antDia = np.zeros(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
#-------- Check Scans and SPWs
BPScan, ASHScan = msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE")[0], msmd.scansforintent("MAP_ANTENNA_SURFACE#ON_SOURCE")[0]
spw = list(set(msmd.tdmspws()) & set(msmd.spwsforintent("MAP_ANTENNA_SURFACE#ON_SOURCE"))); spw.sort(); spwNum = len(spw)
#-------- Scanning and Tracking antennas
interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spw[0], ASHScan)
timeNum = len(timeStamp)
trkAnt, scnAnt, scanTime, AzElOffset = antRefScan(msfile, [min(timeStamp), max(timeStamp)])
trkAntNum, scnAntNum = len(trkAnt), len(scnAnt)
trkBlNum, blNum = trkAntNum* (trkAntNum - 1) / 2, antNum* (antNum - 1)/2
trkBlMap, blMap, trkBlInv, blInv = range(trkBlNum), range(blNum), [False]* trkBlNum, [False]* blNum
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
for bl_index in range(trkBlNum): trkBlMap[bl_index] = Ant2Bl(trkAnt[ant0[bl_index]], trkAnt[ant1[bl_index]])
#-------- Choose Refant from Tracking antennas
timeBPScan, UVW = GetUVW(msfile, spw[0], BPScan)
uvw = np.mean(UVW, axis=2)[:,trkBlMap]; uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
refantID = trkAnt[bestRefant(uvDist)]
trkAnt.remove(refantID); trkAnt.insert(0, refantID)
print 'Tracking antennas ',; print antList[trkAnt]
print '-- Select %s as Refant' % (antList[refantID])
print 'Scanning antennas ',; print antList[scnAnt]
antMap = trkAnt + scnAnt
#-------- Indexing Time Stamps
scanTime, AntID, az, el = GetAzEl(msfile)
FWHM, Freq = [], []
for spw_index in range(spwNum):
    chNum, chWid, freq = GetChNum(msfile, spw[spw_index])
    Freq = Freq + [np.median(freq)]
    FWHM = FWHM + [GetFWHM(msfile, spw[spw_index], antDia[antMap])]      # FWHM in arcsec
#
FWHM = np.array(FWHM).transpose(1,0)    # [ant, spw]
centerIndex = scanThresh(msfile, scnAnt[0], FWHM[scnAnt[0], 0]/10); centerTime = scanTime[centerIndex]
matchNum  = np.zeros([timeNum])
for time_index in range(timeNum):
    matchNum[time_index] = timeMatch(timeStamp[time_index], centerTime, np.median(interval))
#
onAxisIndex = np.where( matchNum > 0 )[0].tolist()
#-------- Baseline Indexing
for bl_index in range(trkBlNum): trkBlMap[bl_index], trkBlInv[bl_index] = Ant2BlD(trkAnt[ant0[bl_index]], trkAnt[ant1[bl_index]])
for bl_index in range(blNum): blMap[bl_index], blInv[bl_index] = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#-------- Bandpass Table
print '---Generating antenna-based bandpass table'
BPList, XYdelayList = [], []
for spw_index in spw:
    BP_ant, XYdelay = BPtable(msfile, spw_index, BPScan, blMap, blInv)
    BPList = BPList + [BP_ant]
    XYdelayList = XYdelayList + [XYdelay]
#
#-------- Scanning Visibilities
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
GainList = []
for spw_index in range(spwNum):
#for spw_index in range(1):
    print '-- Loading visibility data of SPW=%d ...' % (spw[spw_index])
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], ASHScan)  # Xspec[POL, CH, BL, TIME]
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    print '---- Bandpass cal'
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3, 2, 0, 1)
    Xspec = (tempSpec / (BPList[spw_index][ant0][:,polXindex]* BPList[spw_index][ant1][:,polYindex].conjugate())).transpose(2,3,1,0) # [:,:,SAblMap]
    #-------- XY delay cal
    print '---- XY delay cal'
    XYdlSpec = delay_cal( np.ones([chNum], dtype=complex), XYdelayList[spw_index] )
    Xspec[1] = (Xspec[1].transpose(1,2,0)* XYdlSpec).transpose(2,0,1)
    Xspec[2] = (Xspec[2].transpose(1,2,0)* XYdlSpec.conjugate()).transpose(2,0,1)
    #-------- Antenna-based Gain
    print '---- Gain cal using tracking antennasl'
    chAvgVis = np.mean(Xspec[:,chRange], axis=1)
    timeNum = chAvgVis.shape[2]
    trkVis = chAvgVis[:,0:trkBlNum]
    Gain  = np.ones([2, antNum, timeNum], dtype=complex)
    Gain[0, 0:trkAntNum] = np.apply_along_axis( gainComplex, 0, trkVis[0])
    Gain[1, 0:trkAntNum] = np.apply_along_axis( gainComplex, 0, trkVis[3])
    #-------- Antenna-based Gain Calibration
    caledVis   = chAvgVis / (Gain[polXindex][:,ant0]* Gain[polYindex][:,ant1].conjugate())
    #-------- Antenna-based Gain Cal at the beam center
    print '---- Gain cal with smoothed beam-center solutions'
    Gain[0][trkAntNum:antNum][:,onAxisIndex] = np.apply_along_axis(gainComplex, 0, caledVis[0][:,onAxisIndex])[trkAntNum:antNum]
    Gain[1][trkAntNum:antNum][:,onAxisIndex] = np.apply_along_axis(gainComplex, 0, caledVis[3][:,onAxisIndex])[trkAntNum:antNum]
    for ant_index in range(trkAntNum,antNum):
        for pol_index in range(2):
            GR, GI = smoothGain(timeStamp[onAxisIndex], Gain[pol_index, ant_index, onAxisIndex])
            Gain[pol_index, ant_index] = np.median(abs(Gain[pol_index, ant_index, onAxisIndex]))* np.exp((0.0 + 1j)* np.arctan2(GI(timeStamp), GR(timeStamp)))
        #
    #
    caledVis = chAvgVis / (Gain[polXindex][:,ant0]* Gain[polYindex][:,ant1].conjugate())
    #-------- Antenna-based Gain Cal at every position
    print '---- Gain pattern of scanning antennas'
    Gain[0][trkAntNum:antNum] = np.apply_along_axis(gainComplex, 0, caledVis[0])[trkAntNum:antNum]
    Gain[1][trkAntNum:antNum] = np.apply_along_axis(gainComplex, 0, caledVis[3])[trkAntNum:antNum]
    GainList = GainList + [Gain]
#
scanGain = np.array(GainList).transpose(2,0,1,3)[trkAntNum:antNum]
#-------- Beam pattern for each antenna
AZEL = []
for ant_index in range(trkAntNum,antNum):
    index = np.where( AntID == antMap[ant_index])[0].tolist()
    scanAZEL = np.zeros([2, timeNum])
    for time_index in range(timeNum):
        scanAZEL[0,time_index], scanAZEL[1,time_index] = AzElMatch(timeStamp[time_index], scanTime[index], np.median(interval), AzElOffset[0][index], AzElOffset[1][index])
    #
    AZEL = AZEL + [scanAZEL]
#
AZEL = np.array(AZEL)

np.save(prefix +  '-antList.npy', antList[scnAnt])
np.save(prefix +  '-scnAZEL.npy', AZEL)
np.save(prefix +  '-scnGain.npy', scanGain)
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
#-------- For Plot
def circlePoints( x, y, radius ):
    angle = np.arange(-pi, (130/128)*pi, pi/128)
    return x + radius* np.cos(angle), y + radius* np.sin(angle)
#
#-------- Plot D-terms of scanning antennas
if BEAMPLOT:
    print('-------- Plot Beam Maps for scan ants ----')
    for spw_index in range(spwNum):
        for ant_index in range(scnAntNum):
            antID = scnAnt[ant_index]
            #-------- Plot
            fig = plt.figure( figsize = (10,10))
            fig.suptitle('%s %s %.2f GHz' % (prefix, antList[antID], Freq[spw_index]*1.0e-9))
            fig.text(0.45, 0.05, 'Az Offset [arcsec]')
            fig.text(0.05, 0.45, 'El Offset [arcsec]', rotation=90)
            #
            xi, yi = np.mgrid[ min(AZEL[ant_index,0]):max(AZEL[ant_index,0]):128j, min(AZEL[ant_index,1]):max(AZEL[ant_index,1]):128j]
            ampBeamX = GridData( 10.0*np.log10(abs(scanGain[ant_index, spw_index, 0])), AZEL[ant_index, 0], AZEL[ant_index,1], xi.reshape(xi.size), yi.reshape(xi.size), FWHM[antID, spw_index]/16).reshape(len(xi), len(xi))
            phsBeamX = GridData( np.angle(scanGain[ant_index, spw_index, 0]), AZEL[ant_index, 0], AZEL[ant_index,1], xi.reshape(xi.size), yi.reshape(xi.size), FWHM[antID, spw_index]/16).reshape(len(xi), len(xi))
            ampBeamY = GridData( 10.0*np.log10(abs(scanGain[ant_index, spw_index, 1])), AZEL[ant_index, 0], AZEL[ant_index,1], xi.reshape(xi.size), yi.reshape(xi.size), FWHM[antID, spw_index]/16).reshape(len(xi), len(xi))
            phsBeamY = GridData( np.angle(scanGain[ant_index, spw_index, 1]), AZEL[ant_index, 0], AZEL[ant_index,1], xi.reshape(xi.size), yi.reshape(xi.size), FWHM[antID, spw_index]/16).reshape(len(xi), len(xi))
            #---- plot BeamX
            plt.subplot( 2, 2, 1, aspect=1); plt.contourf(xi, yi, ampBeamX, np.linspace(-25.0, 1.0, 27)); plt.colorbar(); plt.title('Amp(BeamX)')
            circle_x, circle_y = circlePoints(0, 0, FWHM[antID, spw_index]/2); plt.plot( circle_x, circle_y )
            circle_x, circle_y = circlePoints(0, 0, FWHM[antID, spw_index]); plt.plot( circle_x, circle_y )
            plt.subplot( 2, 2, 2, aspect=1); plt.contourf(xi, yi, 180.0*phsBeamX/np.pi, np.linspace(-180, 180, 37)); plt.colorbar(); plt.title('Phs(BeamX)')
            circle_x, circle_y = circlePoints(0, 0, FWHM[antID, spw_index]/2); plt.plot( circle_x, circle_y )
            circle_x, circle_y = circlePoints(0, 0, FWHM[antID, spw_index]); plt.plot( circle_x, circle_y )
            plt.subplot( 2, 2, 3, aspect=1); plt.contourf(xi, yi, ampBeamY, np.linspace(-25.0, 1.0, 27)); plt.colorbar(); plt.title('Amp(BeamX)')
            circle_x, circle_y = circlePoints(0, 0, FWHM[antID, spw_index]/2); plt.plot( circle_x, circle_y )
            circle_x, circle_y = circlePoints(0, 0, FWHM[antID, spw_index]); plt.plot( circle_x, circle_y )
            plt.subplot( 2, 2, 4, aspect=1); plt.contourf(xi, yi, 180.0*phsBeamY/np.pi, np.linspace(-180, 180, 37)); plt.colorbar(); plt.title('Phs(BeamX)')
            circle_x, circle_y = circlePoints(0, 0, FWHM[antID, spw_index]/2); plt.plot( circle_x, circle_y )
            circle_x, circle_y = circlePoints(0, 0, FWHM[antID, spw_index]); plt.plot( circle_x, circle_y )
            plt.savefig( prefix + '-' + antList[antID] + '-SPW' + `spw[spw_index]` + '-Beam.pdf', form='pdf'); plt.close()
        #
    #
#
