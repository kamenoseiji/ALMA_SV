execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import griddata
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#----------------------------------------- Procedures
fileNum, spwNum = len(prefixList), len(spwList)
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
antNum = len(trkAnt); blNum = antNum * (antNum - 1)/2
#
for spw_index in range(spwNum):
    spw = spwList[spw_index]
    caledVis = np.ones([4,blNum, 0], dtype=complex)
    mjdSec, Az, El = np.ones([0]), np.ones([0]), np.ones([0])
    for file_index in range(fileNum):
        prefix, scan = prefixList[file_index], scanList[file_index]
        msfile = wd + prefix + '.ms'
        #BPantList = np.load(wd + prefix + '.Ant.npy') 
        #BP_ant = np.load(wd + BPfile[file_index])       # 
        #XYdelay = np.load(wd + XYdelayfile[file_index])
        interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spw, scan); timeNum = len(timeStamp)
        #-------- Subarray with Tracking antennas
        antList = GetAntName(msfile)
        interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spw, scan)
        if refantName not in antList[trkAnt]:
            print refantName + ' does not exist in this MS.'
            sys.exit()
        #
        refantID = np.where( antList == refantName)[0][0]
        antMap = [refantID] + np.array(trkAnt)[range(trkAnt.index(refantID))].tolist() + np.array(trkAnt)[range(trkAnt.index(refantID)+1, len(trkAnt))].tolist()
        #-------- AZ, EL, PA
        azelTime, AntID, AZ, EL = GetAzEl(msfile)
        azelTime_index = np.where( AntID == refantID )[0].tolist()
        timeThresh = np.median( np.diff( azelTime[azelTime_index]))
        #-------- Time index at on axis
        for time_index in range(timeNum):
            scanAz, scanEl = AzElMatch(timeStamp[time_index], azelTime, timeThresh, AZ, EL)
            Az, El = np.append(Az, scanAz), np.append(El, scanEl)
        #
        mjdSec = np.append(mjdSec, timeStamp)
        #-------- Baseline mapping
        ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
        blMap, blInv= range(blNum), [False]* blNum
        for bl_index in range(blNum):
            blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
        #
        #-------- Bandpass table
        print '-- Bandpass table %s ...' % (prefix)
        BP_ant, XYdelay = BPtable(msfile, spw, scan, blMap, blInv)
        #-------- Load Visibilities
        print '-- Loading visibility data of SPW=%d ...' % (spw)
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)  # Xspec[POL, CH, BL, TIME]
        chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
        #-------- Bandpass Calibration
        print '---- Bandpass cal'
        tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3, 2, 0, 1)
        #Xspec = (tempSpec / (BP_ant[ant0][:,polXindex]* BP_ant[ant1][:,polYindex].conjugate())).transpose(2,3,1,0) # 
        Xspec = (tempSpec / (BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # 
        #-------- XY delay cal
        print '---- XY delay cal'
        XYdlSpec = delay_cal( np.ones([chNum], dtype=complex), XYdelay )
        Xspec[1] = (Xspec[1].transpose(1,2,0)* XYdlSpec).transpose(2,0,1)
        Xspec[2] = (Xspec[2].transpose(1,2,0)* XYdlSpec.conjugate()).transpose(2,0,1)
        #-------- Antenna-based Gain
        print '---- Gain cal using tracking antennasl'
        chAvgVis = np.mean(Xspec[:,chRange], axis=1)
        timeNum = chAvgVis.shape[2]
        Gain  = np.ones([2, antNum, timeNum], dtype=complex)
        Gain[0, 0:antNum] = np.apply_along_axis( gainComplex, 0, chAvgVis[0])
        Gain[1, 0:antNum] = np.apply_along_axis( gainComplex, 0, chAvgVis[3])
        #-------- Gain-calibrated visibilities
        print '---- Calibrated visibilities'
        caledVis = concatenate([caledVis, chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())], axis=2)
    #
    PA = AzEl2PA(Az, El, ALMA_lat) - BANDPA
    solution = np.zeros([7])
    for iter_index in range(3):
        print '---- Iteration ' + `iter_index` + ' for Stokes (Q, U) and Gain ----'
        GainX, GainY = polariGain(caledVis[0], caledVis[3], PA, solution[0], solution[1]); Gain = np.array([GainX, GainY])
        Vis    = np.mean(caledVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate()), axis=1)
        solution, solerr = XY2Stokes(PA, Vis[1], Vis[2])
        text_sd = 'Q/I= %6.3f+-%6.4f  U/I= %6.3f+-%6.4f  XYphase= %6.3f+-%6.4f rad EVPA = %6.2f deg' % (solution[0], solerr[0], solution[1], solerr[1], solution[2], solerr[2], np.arctan(solution[1]/solution[0])*90.0/pi)
        print text_sd
    #
    plt.plot(PA, Vis[1].real, '.', label = 'ReXY', color='cyan')
    plt.plot(PA, Vis[1].imag, '.', label = 'ImXY', color='darkblue')
    plt.plot(PA, Vis[2].real, '.', label = 'ReYX', color='magenta')
    plt.plot(PA, Vis[2].imag, '.', label = 'ImYX', color='darkred')
    PArange = np.arange(min(PA), max(PA), 0.01)
    plt.plot(PArange,  np.cos(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[3], '-', color='cyan')
    plt.plot(PArange,  np.sin(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[4], '-', color='darkblue')
    plt.plot(PArange,  np.cos(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[5], '-', color='magenta')
    plt.plot(PArange, -np.sin(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[6], '-', color='darkred')
    plt.xlabel('PA [rad]'); plt.ylabel('XY, YX (real and imaginary)')
    plt.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    text_sd = 'Q/I=%6.3f+-%6.3f U/I=%6.3f+-%6.3f XYphase=%6.3f+-%6.3f rad (RefAnt:%s)' % (solution[0], solerr[0], solution[1], solerr[1], solution[2], solerr[2], antList[refantID]); plt.text(min(PA), min(Vis[1].real), text_sd, size='x-small')
    plt.savefig(prefixList[0] + '-SPW' + `spw` + '-' + refantName + 'QUXY.pdf', form='pdf')
    #-------- Save Results
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Ant.npy', antList[antMap])
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Azel.npy', np.array([mjdSec, Az, El, PA]))
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.QUXY.npy', solution )
    plt.close('all')
#
