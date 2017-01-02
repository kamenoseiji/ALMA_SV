execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import griddata
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#----------------------------------------- Procedures
fileNum, spwNum = len(prefixList), len(spwList)
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
#
trkAntSet = set(range(64))
scansFile = []
for file_index in range(fileNum):
    prefix = prefixList[file_index]
    msfile = wd + prefix + '.ms'
    print '-- Checking %s ' % (msfile)
    antList = GetAntName(msfile)
    refAntID  = indexList([refantName], antList)
    flagAntID = indexList(antFlag, antList)
    if len(refAntID) < 1:
        print 'Antenna %s didn not participate in this file.' % refantName
        continue
    else:
        refAntID = refAntID[0]
    #
    msmd.open(msfile)
    scanList = msmd.scannumbers().tolist()
    for scan in scanList:
        timeStamp = msmd.timesforscans(scan).tolist()
        trkAnt, scanAnt, Time, Offset = antRefScan( msfile, [timeStamp[0], timeStamp[-1]] )
        trkAnt = list(set(trkAnt) - set(flagAntID))
        if refAntID in trkAnt:
            trkAntSet = set(trkAnt) & trkAntSet
        else:
            scanList = list( set(scanList) - set([scan]) )
        #
        print '---- Scan %d : %d tracking antennas' % (scan, len(trkAnt))
    #
    scansFile.append(scanList)
#
antMap = [refAntID] + list(trkAntSet - set([refAntID]))
antNum = len(antMap); blNum = antNum * (antNum - 1)/2
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
blMap, blInv= range(blNum), [False]* blNum
for bl_index in range(blNum):
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#
#-------- Loop for SPW
for spw_index in range(spwNum):
    spw = spwList[spw_index]
    caledVis = np.ones([4,blNum, 0], dtype=complex)
    mjdSec, Az, El = np.ones([0]), np.ones([0]), np.ones([0])
    if BPprefix != '':  # Bandpass file
        BPantList, BP_ant, XYspec, XYdelay = np.load(BPprefix + '-REF' + refantName + '.Ant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-BPant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYspec.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYdelay.npy')
        BP_ant = BP_ant[indexList(antList[antMap], BPantList)]
        BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()
    #
#
    chAvgVis = np.zeros([4, blNum, 0], dtype=complex)
    for file_index in range(fileNum):
        prefix = prefixList[file_index]
        msfile = wd + prefix + '.ms'
        #-------- AZ, EL, PA
        azelTime, AntID, AZ, EL = GetAzEl(msfile)
        azelTime_index = np.where( AntID == refAntID )[0].tolist()
        timeThresh = np.median( np.diff( azelTime[azelTime_index]))
        for scan in scansFile[file_index]:
            #-------- Load Visibilities
            print '-- Loading visibility data %s SPW=%d SCAN=%d...' % (prefix, spw, scan)
            timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)  # Xspec[POL, CH, BL, TIME]
            chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
            if chNum == 1:
                print '  -- Channel-averaged data: no BP and delay cal'
                chAvgVis= np.c_[chAvgVis, CrossPolBL(Xspec[:,:,blMap], blInv)[:,0]]
            else:
                print '  -- Apply bandpass cal'
                tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3, 2, 0, 1)
                BPCaledXspec = (tempSpec / BP_bl).transpose(2,3,1,0) 
                #-------- XY delay cal
                print '  -- XY delay cal'
                BPCaledXspec[1] = (BPCaledXspec[1].transpose(1,2,0)* XYspec.conjugate()).transpose(2,0,1)
                BPCaledXspec[2] = (BPCaledXspec[2].transpose(1,2,0)* XYspec).transpose(2,0,1)
                #-------- Antenna-based Gain
                print '  -- Channel-averaging'
                chAvgVis = np.c_[chAvgVis, np.mean(BPCaledXspec[:,chRange], axis=1)]
            #
            #-------- Time index at on axis
            timeNum = len(timeStamp)
            mjdSec = np.append(mjdSec, timeStamp)
            for time_index in range(timeNum):
                scanAz, scanEl = AzElMatch(timeStamp[time_index], azelTime, AntID, refAntID, timeThresh, AZ, EL)
                Az, El = np.append(Az, scanAz), np.append(El, scanEl)
            #
        #
    #
    print '---- Antenna-based gain solution using tracking antennas'
    Gain  = np.ones([2, antNum, len(mjdSec)], dtype=complex)
    Gain[0, 0:antNum] = gainComplexVec(chAvgVis[0])
    Gain[1, 0:antNum] = gainComplexVec(chAvgVis[3])
    Gamp = np.sqrt(np.mean(abs(Gain)**2, axis=0))
    Gain[0, 0:antNum] = Gain[0, 0:antNum] * Gamp / abs(Gain[0, 0:antNum])
    Gain[1, 0:antNum] = Gain[1, 0:antNum] * Gamp / abs(Gain[1, 0:antNum])
    #-------- Gain-calibrated visibilities
    print '  -- Apply gain calibration'
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    PA = AzEl2PA(Az, El, ALMA_lat) + BANDPA     # Apply feed orientation
    PA = np.arctan2( np.sin(PA), np.cos(PA))    # to set in [-pi, pi]
    Vis    = np.mean(caledVis, axis=1)
    #-------- Solve for Stokes Parameters and XY phase
    print '  -- Solution for Q and U'
    solution = np.r_[XXYY2QU(PA, Vis[[0,3]]), np.zeros(5)]              # XX*, YY* to estimate Q, U
    solution[2] = XY2Phase(PA, solution[0], solution[1], Vis[[1,2]])    # XY*, YX* to estimate X-Y phase
    solution, solerr = XY2Stokes(PA, Vis[[1,2]], solution)
    GainX, GainY = polariGain(caledVis[0], caledVis[3], PA, solution[0], solution[1]); Gain = np.array([GainX, GainY])
    Vis    = np.mean(caledVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate()), axis=1)
    solution, solerr = XY2Stokes(PA, Vis[[1,2]], solution)
    text_sd = '  Q/I= %6.3f+-%6.4f  U/I= %6.3f+-%6.4f  X-Y phase= %6.3f+-%6.4f rad EVPA = %6.2f deg' % (solution[0], solerr[0], solution[1], solerr[1], solution[2], solerr[2], np.arctan2(solution[1],solution[0])*90.0/pi); print text_sd
    #
    #-------- Plot
    if np.mean(np.cos(PA)) < 0.0: PA = np.arctan2(-np.sin(PA), -np.cos(PA)) +  np.pi
    PArange = np.arange(min(PA), max(PA), 0.01)
    plt.plot(PArange,  np.cos(2.0*PArange)* solution[0] + np.sin(2.0* PArange)* solution[1], '-', color='green')
    plt.plot(PArange,  np.cos(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[3], '-', color='cyan')
    plt.plot(PArange,  np.sin(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[4], '-', color='darkblue')
    plt.plot(PArange,  np.cos(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[5], '-', color='magenta')
    plt.plot(PArange, -np.sin(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[6], '-', color='darkred')
    plt.plot(PArange, -np.cos(2.0*PArange)* solution[0] - np.sin(2.0* PArange)* solution[1], '-', color='orange')
    plt.plot(PA, Vis[0].real - 1.0, '.', label = 'XX* - 1.0',   color='green')
    plt.plot(PA, Vis[1].real, '.', label = 'ReXY*', color='cyan')
    plt.plot(PA, Vis[1].imag, '.', label = 'ImXY*', color='darkblue')
    plt.plot(PA, Vis[2].real, '.', label = 'ReYX*', color='magenta')
    plt.plot(PA, Vis[2].imag, '.', label = 'ImYX*', color='darkred')
    plt.plot(PA, Vis[3].real - 1.0, '.', label = 'YY* - 1.0',   color='orange')
    plt.xlabel('X-Feed Position Angle [rad]'); plt.ylabel('Normalized cross correlations')
    plt.ylim([-0.15,0.15])
    plt.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    text_sd = '(Q, U)/I = (%6.3f+-%6.3f, %6.3f+-%6.3f) X-Y phase=%6.3f+-%6.3f rad (Ref:%s)' % (solution[0], solerr[0], solution[1], solerr[1], solution[2], solerr[2], antList[refAntID]); plt.text(min(PA), 0.16, text_sd, size='x-small')
    plt.savefig(prefixList[0] + '-SPW' + `spw` + '-' + refantName + 'QUXY.pdf', form='pdf')
    #-------- Save Results
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Ant.npy', antList[antMap])
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Azel.npy', np.array([mjdSec, Az, El, PA]))
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.QUXY.npy', solution )
    plt.close('all')
#
