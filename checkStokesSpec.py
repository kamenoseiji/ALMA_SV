execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Grid.py')
execfile(SCR_DIR + 'Plotters.py')
import pickle
from matplotlib.backends.backend_pdf import PdfPages
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#----------------------------------------- Procedures
spwNum = len(spwList)
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
#
trkAntSet = set(range(64))
scansFile = []
pattern = r'RB_..'
timeNum = 0
sourceList = []
msfile = wd + prefix + '.ms'
sources, posList = GetSourceList(msfile); sourceList = sourceList + sourceRename(sources)
sourceList = unique(sourceList).tolist()
sourceScan = []
scanDic   = dict(zip(sourceList, [[]]*len(sourceList))) # Scan list index for each source
timeDic   = dict(zip(sourceList, [[]]*len(sourceList))) # Time index list for each source
#StokesDic = dict(zip(sourceList, [[]]*len(sourceList))) # Stokes parameters for each source
fileDic   = open(SourceDicFile, mode='rb')
StokesDic = pickle.load(fileDic)
fileDic.close()
scanIndex = 0
msfile = wd + prefix + '.ms'
msmd.open(msfile)
print '-- Checking %s ' % (msfile)
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList)
antList = GetAntName(msfile)
refAntID  = indexList([refantName], antList)
flagAntID = indexList(antFlag, antList)
if len(refAntID) < 1: print 'Antenna %s didn not participate in this file.' % refantName
else: refAntID = refAntID[0]
#
#
if 'scanList' in locals(): scanLS = scanList
else: scanLS = msmd.scannumbers().tolist()
#
spwName = msmd.namesforspws(spwList[0])[0]; BandName = re.findall(pattern, spwName)[0]; bandID = int(BandName[3:5])
BandPA = (BANDPA[bandID] + 90.0)*pi/180.0
for scan_index in range(len(scanLS)):
    scan = scanLS[scan_index]
    interval, timeStamp = GetTimerecord(msfile, 0, 1, spwList[0], scan)
    timeNum += len(timeStamp)
    trkAnt, scanAnt, Time, Offset = antRefScan( msfile, [timeStamp[0], timeStamp[-1]], antFlag )
    trkAnt = list(set(trkAnt) - set(flagAntID))
    if refAntID in trkAnt: trkAntSet = set(trkAnt) & trkAntSet
    #
'''
    sourceName = sourceList[msmd.sourceidforfield(msmd.fieldsforscan(scan)[0])]
    sourceScan = sourceScan + [sourceName]
    #scanDic[sourceName] = scanDic[sourceName] + [scanIndex]
    if AprioriDic[sourceName] :
        StokesDic[sourceName] = AprioriDic[sourceName]
    else: 
        IQU = GetPolQuery(sourceName, timeStamp[0], BANDFQ[bandID], SCR_DIR, R_DIR)
        if len(IQU[0]) > 0:
            StokesDic[sourceName] = [IQU[0][sourceName], IQU[1][sourceName], IQU[2][sourceName], 0.0]
        else:
            StokesDic[sourceName] = [0.01, 0.0, 0.0, 0.0]
        #
    print '---- Scan%3d : %d tracking antennas : %s, %d records, expected I=%.1f p=%.1f%%' % (scan, len(trkAnt), sourceName, len(timeStamp), StokesDic[sourceName][0], 100.0*sqrt(StokesDic[sourceName][1]**2 + StokesDic[sourceName][2]**2)/StokesDic[sourceName][0])
    scanIndex += 1
#
'''
scansFile.append(scanLS)
BPsrc = sourceList[msmd.sourceidforfield(msmd.fieldsforscan(BPscan)[0])]
#
if not 'scanList' in locals(): scanList = scansFile
#-------- Check source list and Stokes Parameters
antMap = [refAntID] + list(trkAntSet - set([refAntID]))
antNum = len(antMap); blNum = antNum * (antNum - 1)/2
#-------- Load Aeff file
antDia = np.ones(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[antMap[ant_index]])['value']
etaA = GetAeff(TBL_DIR, antList[antMap], bandID, np.mean(timeStamp)).T
Ae = 0.0025* np.pi* etaA* antDia[antMap]**2
if 'BLCORR' in locals():
    if BLCORR: Ae = 1.15* Ae         # BL Correlator correction factor
#-------- Antenna-Baseline mapping
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
blMap, blInv= range(blNum), [False]* blNum
for bl_index in range(blNum):
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#
if not 'bunchNum' in locals(): bunchNum = 1
def bunchVecCH(spec): return bunchVec(spec, bunchNum)
print 'Total %d integration periods.' % (timeNum)
msmd.done()
#-------- Loop for SPW
DxList, DyList, FreqList = [], [], []
for spw_index in range(spwNum):
    spw = spwList[spw_index]
    chNum, chWid, Freq = GetChNum(msfile, spw); Freq *= 1.0e-9; chRange = range(int(0.05*chNum/bunchNum), int(0.95*chNum/bunchNum)); FreqList = FreqList + [bunchVecCH(Freq) ]
    #-------- Load D-term files
    DxSpec, DySpec = np.zeros([antNum, chNum/bunchNum], dtype=complex), np.zeros([antNum, chNum/bunchNum], dtype=complex)
    for antName in antList[antMap]:
        Dfile = DPATH + '-SPW' + `spwList[0]` + '-' + antName + '.DSpec.npy'
        if os.path.exists(Dfile):
            Dterm = np.load(Dfile)
            Dfreq = Dterm[0]
            if len(Dfreq) != chNum/bunchNum:
                SB = int(np.sign( np.median(np.diff(Dfreq))))   # LSB -> -1, UWB -> +1
                DX_real, DX_imag = scipy.interpolate.splrep(Dfreq[::SB], Dterm[1][::SB], s=0.5), scipy.interpolate.splrep(Dfreq[::SB], Dterm[2][::SB], s=0.5)
                DY_real, DY_imag = scipy.interpolate.splrep(Dfreq[::SB], Dterm[3][::SB], s=0.5), scipy.interpolate.splrep(Dfreq[::SB], Dterm[4][::SB], s=0.5)
                DxList = DxList + [scipy.interpolate.splev(FreqList[0], DX_real)+ (0.0 + 1.0j)* scipy.interpolate.splev(FreqList[0], DX_imag)]
                DyList = DyList + [scipy.interpolate.splev(FreqList[0], DY_real)+ (0.0 + 1.0j)* scipy.interpolate.splev(FreqList[0], DY_imag)]
            else:
                DxList = DxList + [Dterm[1] + (0.0 + 1.0j)* Dterm[2]]
                DyList = DyList + [Dterm[3] + (0.0 + 1.0j)* Dterm[4]]
            #
        else:
            print 'D-term file [%s] is not found.' % (Dfile)
            continue
        #
    #
    DxSpec, DySpec = np.array(DxList).reshape([antNum, chNum/bunchNum]), np.array(DyList).reshape([antNum, chNum/bunchNum])
    M = InvMullerVector(DxSpec[ant1], DySpec[ant1], DxSpec[ant0], DySpec[ant0], np.ones([blNum,chNum/bunchNum])).transpose(2,3,0,1)
    D = MullerVector(DxSpec[ant1], DySpec[ant1], DxSpec[ant0], DySpec[ant0],np.ones([blNum, chNum/bunchNum])).transpose(2,3,0,1)
    #-------- Load Aeff file
    '''
    msmd.open(msfile)
    antDia = np.ones(antNum)
    AeX, AeY, etaX, etaY = [], [], [], []
    for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[antMap[ant_index]])['value']
    msmd.done()
    Afile = open(SCR_DIR + 'AeB' + `bandID` + '.data')
    Alines = Afile.readlines()
    Afile.close()
    for ant_index in range(antNum):
        for Aline in Alines:
            if antList[antMap[ant_index]] in Aline:
                etaX = etaX + [float(Aline.split()[1])]
                etaY = etaY + [float(Aline.split()[2])]
                AeX = AeX + [(0.0025* np.pi* float(Aline.split()[1]))* antDia[ant_index]**2]
                AeY = AeY + [(0.0025* np.pi* float(Aline.split()[2]))* antDia[ant_index]**2]
                # print '%s  : etaX = %.2f  etaY = %.2f' % (antList[antMap[ant_index]], float(Aline.split()[1]), float(Aline.split()[2]))
            #
        #
    #
    Ae = np.array([AeX, AeY])
    '''
    #-------- Load Flag Table
    caledVis = np.ones([4,blNum, 0], dtype=complex)
    if 'FGprefix' in locals():  # Flag table
        FG = np.load(FGprefix + '-SPW' + `spw` + '.FG.npy'); FG = np.min(FG, axis=0)
        TS = np.load(FGprefix + '-SPW' + `spw` + '.TS.npy')
    #
    #-------- Load Bandpass Table
    if 'BPprefix' in locals():  # Bandpass file
        BPantList, BP_ant = np.load(BPprefix + '-REF' + refantName + '.Ant.npy'), np.load(BPprefix + '-REF' + refantName + '-SC' + `BPscan` + '-SPW' + `spw` + '-BPant.npy')
        BP_ant = BP_ant[indexList(antList[antMap], BPantList)]      # BP antenna mapping
    #-------- Load XY phase spectral table
    if 'XYprefix' in locals():
    #-------- Load XY phase file
        XYspec = np.load(XYprefix + '-REF' + refantName + '-SC' + `BPscan` + '-SPW' + `spw` + '-XYspec.npy')
        print 'Apply XY phase into Y-pol Bandpass.'; BP_ant[:,1] *= XYspec  # XY phase cal
        XYPH = np.load(XYprefix + '-SPW' + `spwList[0]` + '-' + refantName + '.XYPH.npy')
        XYTS = np.load(XYprefix + '-SPW' + `spwList[0]` + '-' + refantName + '.TS.npy')
        XYPHSP = scipy.interpolate.splrep(XYTS, XYPH, k=3, task=0, s=0.1)
    #
    BP_ant = np.apply_along_axis(bunchVecCH, 2, BP_ant)
    #
    #-------- Load Tsys table
    TrxSpec  = np.ones([antNum, 2, chNum/bunchNum])
    TrxFreq  = np.load(prefix +  '-' + BandName + '-SPW' + `TsysSPW[spw_index]` + '.TrxFreq.npy') # TrxFreq[ch]
    TrxAnts  = np.load(prefix +  '-' + BandName + '.TrxAnt.npy') # TrxAnts[ant]
    Tau0spec = np.load(prefix +  '-' + BandName + '-SPW' + `TsysSPW[spw_index]` + '.Tau0.npy') # Tau0spec[ch]
    Trxspec  = np.load(prefix +  '-' + BandName + '-SPW' + `TsysSPW[spw_index]` + '.Trx.npy')  # Trxspec[spw, pol, ch, ant, scan]
    Tau0E    = np.load(prefix +  '-' + BandName + '.TauE.npy')[spw_index] # Tau0E[atmScan]
    atmTimeRef = np.load(prefix +  '-' + BandName + '.atmTime.npy') # atmTimeRef[atmScan]
    TrxMap = indexList(TrxAnts, antList[antMap]); TrxFlag = np.zeros([antNum]); TrxFlag[TrxMap] = 1.0
    TrxMed = np.median(Trxspec, axis=(1,3))
    #-------- Tsys channel interpolation
    if len(TrxFreq) != chNum/bunchNum:
        SB = int(np.sign( np.median(np.diff(TrxFreq))))   # LSB -> -1, UWB -> +1
        TAU0SP = scipy.interpolate.splrep(TrxFreq[::SB], Tau0spec[::SB], k=3)
        Tau0Spec = scipy.interpolate.splev(FreqList[0], TAU0SP)
        for ant_index in range(len(TrxAnts)):
            for pol_index in range(2):
                TRXSP = scipy.interpolate.splrep(TrxFreq[::SB], np.median(Trxspec, axis=3)[pol_index, :, ant_index][::SB])
                TrxSpec[ant_index, pol_index] = scipy.interpolate.splev(FreqList[0], TRXSP)
            #
        #
    else:
        TrxSpec, Tau0Spec = np.median(Trxspec, axis=3), Tau0spec
    #
    tempAtm = GetTemp(msfile)       # atmosphere temperature
    #-------- Tau0 time interoplation
    TAU0ESP = scipy.interpolate.splrep(atmTimeRef, Tau0E, k=3, task=0, s=0.1*np.std(Tau0E))
    #-------- Flagging Trx outliers
    for pol_index in range(2):
        TrxFlag[np.where(TrxMed[pol_index] - np.median(TrxMed[pol_index]) > 1.5* np.median(TrxMed[pol_index]))[0].tolist()] *= 0.0
        TrxFlag[np.where(TrxMed[pol_index] < 0.25* np.median(TrxMed[pol_index]))[0].tolist()] *= 0.0
    #
    #-------- AZ, EL, PA 
    azelTime, AntID, AZ, EL = GetAzEl(msfile)
    azelTime_index = np.where( AntID == refAntID )[0].tolist()
    if len(azelTime_index) == 0: azelTime_index = np.where(AntID == 0)[0].tolist()
    timeThresh = np.median( np.diff( azelTime[azelTime_index]))
    #-------- Bandpass Scan
    IQUV = np.array(StokesDic[BPsrc])
    BPinterval, BPtimeStamp = GetTimerecord(msfile, 0, 1, spwList[spw_index], BPscan)
    BPAz, BPEl = AzElMatch(BPtimeStamp, azelTime, AntID, refAntID, AZ, EL)
    exp_Tau = np.exp(-(Tau0Spec + np.mean(scipy.interpolate.splev(BPtimeStamp, TAU0ESP))) / np.median(np.sin(BPEl)))
    atmCorrect = 1.0 / exp_Tau
    TsysBPScan = atmCorrect* (TrxSpec.transpose(2,0,1)[TrxMap] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))    # [ant, pol, ch]
    TsysBPShape = (TsysBPScan.transpose(2,0,1) / np.median(TsysBPScan, axis=2)).transpose(1,2,0)    # [ant, pol, ch]
    SEFD_BP = 2.0* kb* (TsysBPScan.transpose(2,1,0) / Ae).transpose(2,1,0)  # [ant, pol, ch]
    #-------- Equalization using BP scan
    PA = AzEl2PA(BPAz, BPEl) + BandPA; PAnum = len(PA)
    PS = PAVector(PA, np.ones(PAnum)).real.transpose(2,0,1).dot(IQUV)    # PS[time, pol] only real part, assuming Stokes V = 0
    XYsign = np.sign(np.mean(PS[:,1]))
    BP_ant[:,1] *= XYsign 
    BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()    # Baseline-based bandpass table
    print '-- Loading visibility data %s SPW=%d SCAN=%d...' % (prefix, spw, BPscan)
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan, -1, True)  # Xspec[POL, CH, BL, TIME]
    del Pspec
    if bunchNum > 1: Xspec = np.apply_along_axis(bunchVecCH, 1, Xspec)
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3, 2, 0, 1)    # tempSpec[time, BL, pol, ch]
    print '  -- Normalizing by source model'
    for time_index in range(PAnum): tempSpec[time_index] /= D.dot(PS[time_index]).transpose(0,2,1)
    print '  -- Apply bandpass cal'
    VisSpec = (tempSpec[:,:,[0,3]] * (np.sqrt(SEFD_BP[ant0]* SEFD_BP[ant1])/BP_bl[:,[0,3]])).transpose(2,3,1,0) # VisSpec[pol, ch, BL, time]
    print '  -- Dtermination of antenna gain'
    chAvgVis = np.mean(VisSpec[:,chRange], axis=1) # chAvgVis[pol, bl, time]
    Gain = np.array([ gainComplexVec(chAvgVis[0]), gainComplexVec(chAvgVis[1]) ])
    gainCaledVis = chAvgVis / (Gain[:,ant0]* Gain[:,ant1].conjugate())
    Gamp = abs(np.mean(Gain, axis=2))
    Ae *= (Gamp**2) # Equalized aperture area
    #-------- For each scan
    StokesVis = np.zeros([4, timeNum, chNum/bunchNum], dtype=complex)
    mjdSec = np.zeros(timeNum)
    timeIndex = 0
    for scan in scanList:
        print '-- Loading visibility data %s SPW=%d SCAN=%d...' % (prefix, spw, scan)
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan, -1, True)  # Xspec[POL, CH, BL, TIME]
        del Pspec
        if bunchNum > 1: Xspec = np.apply_along_axis(bunchVecCH, 1, Xspec)
        timeIndexRange = range(timeIndex, timeIndex + len(timeStamp))
        mjdSec[timeIndexRange] = timeStamp
        print '  -- Apply bandpass cal'
        tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3, 2, 0, 1)    # tempSpec[time, BL, pol, ch]
        VisSpec = (tempSpec / BP_bl).transpose(2,3,1,0)     # VisSpec[pol, ch, BL, time]
        print '  -- Apply parallel-pol phase cal'
        chAvgVis = np.mean(VisSpec[[0,3]][:,chRange], axis=1) # chAvgVis[pol, bl, time]
        GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
        twiddle = np.exp((1.0j)* scipy.interpolate.splev(timeStamp, XYPHSP))
        GainP[1] *= twiddle     # XY phase correction
        pCalVis = (VisSpec.transpose(1,0,2,3) / (GainP[polYindex][:,ant0]* GainP[polXindex][:,ant1].conjugate()))
        print '  -- Apply amplitude cal'
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refAntID, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
        exp_Tau = np.exp(-(Tau0Spec + np.mean(scipy.interpolate.splev(timeStamp, TAU0ESP))) / np.median(np.sin(ElScan)))
        atmCorrect = 1.0 / exp_Tau
        TsysScan = atmCorrect* (TrxSpec.transpose(2,0,1)[TrxMap] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))    # [ant, pol, ch]
        SEFD = 2.0* kb* ((TsysBPScan / TsysBPShape).transpose(2,1,0) / Ae).transpose(2,1,0)  # [ant, pol, ch]
        GainCaledVis = pCalVis.transpose(3,2,1,0)* np.sqrt(SEFD[ant0][:,polYindex]*SEFD[ant1][:,polXindex])
        print '  -- Equalizing antenna-based gain variation'
        chAvgVis = np.mean(GainCaledVis[:,:,[0,3]][:,:,:,chRange], axis=(0,3))  # [pol, bl]
        indivRelGain = abs(gainComplexVec(chAvgVis))
        indivRelGain /= np.median(indivRelGain, axis=0)
        SEFD = (SEFD.transpose(2,0,1) / indivRelGain**2).transpose(1,2,0)
        GainCaledVis = pCalVis.transpose(3,2,1,0)* np.sqrt(SEFD[ant0][:,polYindex]*SEFD[ant1][:,polXindex])
        print '  -- Stokes spectra'
        for time_index in range(PAnum):
            for bl_index in range(blNum):
                for ch_index in range(chNum/bunchNum):
                    StokesVis[:,(time_index + timeIndex), ch_index] += PS[:, :, time_index].dot(M[bl_index, ch_index].dot(GainCaledVis[time_index, bl_index, :, ch_index]))
                #
            #
        #
        timeIndex += len(timeStamp)
    #
    StokesVis /= blNum
    text_sd = '(I,Q,U,V) = %.4f %.4f %.4f %.4f' % (np.mean(StokesVis[0]).real, np.mean(StokesVis[1]).real, np.mean(StokesVis[2]).real, np.mean(StokesVis[3]).real)
    print text_sd
    np.save(prefix + '.Stokes.npy', StokesVis)
    np.save(prefix + '.TS.npy', mjdSec)

    """
    VisSpec  = np.zeros([4, chNum/bunchNum, blNum, timeNum], dtype=complex)
    chAvgVis = np.zeros([4, blNum, timeNum], dtype=complex)
    mjdSec, scanST, scanET, Az, El, PA = [], [], [], [], [], []
    timeIndex, scanIndex  = 0, 0
    for scan in scanList:
        #-------- Load Visibilities
        print '-- Loading visibility data %s SPW=%d SCAN=%d...' % (prefix, spw, scan)
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan, -1, True)  # Xspec[POL, CH, BL, TIME]
        del Pspec
        if bunchNum > 1: Xspec = np.apply_along_axis(bunchVecCH, 1, Xspec)
        if 'FG' in locals():
            flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
        else:
            flagIndex = range(len(timeStamp))
        if chNum == 1:
            print '  -- Channel-averaged data: no BP and delay cal'
            chAvgVis[:, :, timeIndex:timeIndex + len(timeStamp)] =  CrossPolBL(Xspec[:,:,blMap], blInv)[:,0]
        else:
            print '  -- Apply bandpass cal'
            tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3, 2, 0, 1)[flagIndex]
            del Xspec
            VisSpec[:,:,:,timeIndex:timeIndex + len(timeStamp)] = (tempSpec / BP_bl).transpose(2,3,1,0) 
            del tempSpec
        #
        scanST = scanST + [timeIndex]
        timeIndex += len(timeStamp)
        scanET = scanET + [timeIndex]
        #-------- Expected polarization responses
        scanAz, scanEl = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refAntID, AZ, EL)
        mjdSec, scanPA = mjdSec + timeStamp[flagIndex].tolist(), AzEl2PA(scanAz, scanEl, ALMA_lat) + BandPA
        Az, El, PA = Az + scanAz.tolist(), El + scanEl.tolist(), PA + scanPA.tolist()
        scanIndex += 1
    #
    scanNum = scanIndex
    if chNum > 1: chAvgVis = np.mean(VisSpec[:,chRange], axis=1)
    mjdSec, Az, El = np.array(mjdSec), np.array(Az), np.array(El)
    print '---- Antenna-based gain solution using tracking antennas'
    Gain = np.array([ gainComplexVec(chAvgVis[0]), gainComplexVec(chAvgVis[3]) ])   # Parallel-pol gain
    Gamp = np.sqrt(np.mean(abs(Gain)**2, axis=0))
    Gain = Gamp* Gain/abs(Gain)
    #-------- Gain-calibrated visibilities
    print '  -- Apply parallel-hand gain calibration'
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    PAnum = len(PA); PA = np.array(PA)
    QCpUS, UCmQS =  np.zeros(PAnum), np.zeros(PAnum)
    print '  -- Solution for Q and U'
    """
    # break point
    '''
    #-------- Coarse estimation of Q and U using XX and YY
    if 'QUmodel' not in locals(): QUmodel = False
    CS, SN = np.cos(2.0* PA), np.sin(2.0* PA)
    for sourceName in sourceList:
        scanLS = scanDic[sourceName]
        if len(scanLS) < 1 : continue
        timeIndex = []
        for scanIndex in scanLS: timeIndex = timeIndex + range(scanST[scanIndex], scanET[scanIndex])
        timeDic[sourceName] = timeIndex
        if QUmodel:
            QUsol = np.array(StokesDic[sourceName])[[1,2]]/StokesDic[sourceName][0]
        else:
            QUsol   = XXYY2QU(PA[timeIndex], Vis[[0,3]][:,timeIndex])             # XX*, YY* to estimate Q, U
            text_sd = '[XX,YY] %s: Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, QUsol[0], QUsol[1], 100.0* np.sqrt(QUsol[0]**2 + QUsol[1]**2), np.arctan2(QUsol[1],QUsol[0])*90.0/pi); print text_sd
        #
        QCpUS[timeIndex] = QUsol[0]* CS[timeIndex] + QUsol[1]* SN[timeIndex]
        UCmQS[timeIndex] = QUsol[1]* CS[timeIndex] - QUsol[0]* SN[timeIndex]
    #
    ##-------- XY phase determination
    print '  -- Degenerating pi-ambiguity in XY phase'
    XYphase = XY2Phase(UCmQS, Vis[[1,2]])    # XY*, YX* to estimate X-Y phase
    XYsign = np.sign(np.cos(XYphase))
    text_sd = '  XY Phase = %6.2f [deg]  sign = %3.0f' % (XYphase* 180.0 / pi, XYsign); print text_sd
    #-------- Gain adjustment
    print '  -- Polarized gain calibration'
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY* XYsign])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    #-------- Fine estimation of Q and U using XY and YX
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        Qsol, Usol = XY2Stokes(PA[timeIndex], Vis[[1,2]][:,timeIndex])
        text_sd = '[XY,YX] %s:  Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, Qsol, Usol, 100.0* np.sqrt(Qsol**2 + Usol**2), np.arctan2(Usol,Qsol)*90.0/pi); print text_sd
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
    #
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis    = np.mean(caledVis, axis=1)
    #-------- XY phase correction
    XYphase, XYproduct = XY2PhaseVec(mjdSec - np.median(mjdSec), UCmQS, Vis[[1,2]])
    twiddle = np.exp((1.0j)* XYphase)
    caledVis[1] /= twiddle
    caledVis[2] *= twiddle
    #-------- Fine estimation of Q and U using XY and YX
    print '  -- XY phase correction'
    Vis    = np.mean(caledVis, axis=1)
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        Qsol, Usol = XY2Stokes(PA[timeIndex], Vis[[1,2]][:,timeIndex])
        text_sd = '[XY,YX] %s:  Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, Qsol, Usol, 100.0* np.sqrt(Qsol**2 + Usol**2), np.arctan2(Usol,Qsol)*90.0/pi); print text_sd
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
    #
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    # Vis    = np.mean(caledVis, axis=1)
    GainCaledVisSpec = VisSpec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    del VisSpec
    #-------- Antenna-based on-axis D-term (chAvg)
    StokesI = np.ones(PAnum)
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        caledVis[:,:,timeIndex] *= StokesDic[sourceName][0]
        QCpUS[timeIndex] *= StokesDic[sourceName][0]
        UCmQS[timeIndex] *= StokesDic[sourceName][0]
        StokesI[timeIndex] *= StokesDic[sourceName][0]
    #
    Dx, Dy = VisMuiti_solveD(caledVis, QCpUS, UCmQS, [], [], StokesI)
    #-------- D-term-corrected Stokes parameters
    Minv = InvMullerVector(Dx[ant1], Dy[ant1], Dx[ant0], Dy[ant0], np.ones(blNum, dtype=complex))
    print '  -- D-term-corrected visibilities'
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        srcTimeNum = len(timeIndex)
        if srcTimeNum < 1 : continue
        PS = InvPAVector(PA[timeIndex], np.ones(srcTimeNum))
        StokesVis = PS.reshape(4, 4*srcTimeNum).dot(Minv.reshape(4, 4*blNum).dot(caledVis[:,:,timeIndex].reshape(4*blNum, srcTimeNum)).reshape(4*srcTimeNum)) / (srcTimeNum* blNum)
        Isol, Qsol, Usol = StokesVis[0].real, StokesVis[1].real, StokesVis[2].real
        Ierr, Qerr, Uerr = abs(StokesVis[0].imag), abs(StokesVis[1].imag), abs(StokesVis[2].imag)
        text_sd = '%s: I= %6.3f  Q= %6.3f+-%6.4f  U= %6.3f+-%6.4f EVPA = %6.2f deg' % (sourceName, Isol, Qsol, Qerr, Usol, Uerr, np.arctan2(Usol,Qsol)*90.0/pi); print text_sd
        StokesDic[sourceName] = (np.array([Isol, Qsol, Usol, 0.0])).tolist()
        StokesI[timeIndex] = Isol
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
        for index in timeIndex: GainCaledVisSpec[:,:,:,index] *= Isol
    #
    #-------- get D-term spectra
    print '  -- Determining D-term spectra'
    for ch_index in range(int(chNum/bunchNum)):
        DxSpec[:,ch_index], DySpec[:,ch_index] = VisMuiti_solveD(GainCaledVisSpec[ch_index], QCpUS, UCmQS, Dx, Dy, StokesI)
        progress = (ch_index + 1.0) / (chNum / bunchNum)
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
    #
    sys.stderr.write('\n'); sys.stderr.flush()
    if 'Dsmooth' in locals():
        node_index = range(3, chNum/bunchNum, Dsmooth)
        for ant_index in range(antNum):
            DX_real, DX_imag = scipy.interpolate.splrep(Freq, DxSpec[ant_index].real, k=3, t=Freq[node_index]), scipy.interpolate.splrep(Freq, DxSpec[ant_index].imag, k=3, t=Freq[node_index])
            DY_real, DY_imag = scipy.interpolate.splrep(Freq, DySpec[ant_index].real, k=3, t=Freq[node_index]), scipy.interpolate.splrep(Freq, DySpec[ant_index].imag, k=3, t=Freq[node_index])
            DxSpec[ant_index] =scipy.interpolate.splev(Freq, DX_real) + (0.0 + 1.0j)* scipy.interpolate.splev(Freq, DX_imag)
            DySpec[ant_index] =scipy.interpolate.splev(Freq, DY_real) + (0.0 + 1.0j)* scipy.interpolate.splev(Freq, DY_imag)
        #
    #
    '''
    #-------- D-term-corrected visibilities (invD dot Vis = PS)
    #del chAvgVis, StokesVis
    #print '  -- Applying D-term spectral correction'
    #M  = InvMullerVector(DxSpec[ant0], DySpec[ant0], DxSpec[ant1], DySpec[ant1], np.ones([blNum,chNum/bunchNum])).transpose(0,3,1,2)
    '''
    StokesVis = np.zeros([4, chNum/bunchNum, PAnum], dtype=complex )
    for time_index in range(PAnum): StokesVis[:, :, time_index] = 4.0* np.mean(M* GainCaledVisSpec[:,:,:,time_index], axis=(2,3))
    chAvgVis = np.mean(StokesVis[:,chRange], axis=1)
    PS = InvPAVector(PA, np.ones(PAnum))
    for ch_index in range(chNum/bunchNum): StokesVis[:,ch_index] = np.sum(PS* StokesVis[:,ch_index], axis=1)
    maxP = 0.0
    for sourceName in sourceList:
        colorIndex = lineCmap(sourceList.index(sourceName) / 8.0)
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        QUsol = np.array([StokesDic[sourceName][1], StokesDic[sourceName][2]])
        maxP = max(maxP, sqrt(QUsol.dot(QUsol)))
        EVPA = 0.5* np.arctan2(QUsol[1], QUsol[0])
        ThetaPlot = PA[timeIndex] - EVPA; ThetaPlot = np.arctan(np.tan(ThetaPlot))
        ThetaMin, ThetaMax = min(ThetaPlot), max(ThetaPlot)
        PArange = np.arange(ThetaMin + EVPA, ThetaMax + EVPA, 0.01)
        ThetaRange = np.arange(ThetaMin, ThetaMax, 0.01)
        CSrange, SNrange = np.cos(2.0*PArange), np.sin(2.0*PArange)
        UCmQS, QCpUS = QUsol[1]*CSrange - QUsol[0]* SNrange, QUsol[0]*CSrange + QUsol[1]* SNrange
        ThetaRange[ThetaRange >  1.56] = np.inf
        ThetaRange[ThetaRange < -1.56] = -np.inf
        plt.plot(RADDEG* ThetaRange,  QCpUS, '-', color=colorIndex, linestyle='dashed', label=sourceName + ' XX* - I')     # XX* - 1.0
        plt.plot(RADDEG* ThetaRange, -QCpUS, '-', color=colorIndex, linestyle='dashdot', label=sourceName + ' YY* - I')     # YY* - 1.0
        plt.plot(RADDEG* ThetaRange,  UCmQS, '-', color=colorIndex, linestyle='solid', label=sourceName + ' ReXY*')
        plt.plot(RADDEG* ThetaRange,  np.zeros(len(ThetaRange)), '-', color=colorIndex, linestyle='dotted', label=sourceName + ' ImXY*')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[0][timeIndex].real - StokesDic[sourceName][0], ',', color=colorIndex)
        plt.plot(RADDEG* ThetaPlot, chAvgVis[1][timeIndex].real, ',', color=colorIndex)
        plt.plot(RADDEG* ThetaPlot, chAvgVis[1][timeIndex].imag, ',', color=colorIndex)
        plt.plot(RADDEG* ThetaPlot, chAvgVis[2][timeIndex].real, ',', color=colorIndex)
        plt.plot(RADDEG* ThetaPlot, chAvgVis[2][timeIndex].imag, ',', color=colorIndex)
        plt.plot(RADDEG* ThetaPlot, chAvgVis[3][timeIndex].real - StokesDic[sourceName][0], ',', color=colorIndex)
    #
    plt.xlabel('Linear polarization angle w.r.t. X-Feed [deg]'); plt.ylabel('Cross correlations [Jy]')
    plt.xlim([-90.0, 90.0])
    plt.ylim([-1.5* maxP, 1.5*maxP])
    plt.legend(loc = 'best', prop={'size' :6}, numpoints = 1)
    plt.savefig(prefixList[0] + '-SPW' + `spw` + '-' + refantName + 'QUXY.pdf', form='pdf')
    #-------- Save Results
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Ant.npy', antList[antMap])
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.Azel.npy', np.array([mjdSec, Az, El, PA]))
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.TS.npy', mjdSec )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.GA.npy', Gain )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.XYPH.npy', XYphase )
    np.save(prefixList[0] + '-SPW' + `spw` + '-' + refantName + '.XYC.npy', chAvgVis[[1,2]])
    for ant_index in range(antNum):
        DtermFile = np.array([FreqList[spw_index], DxSpec[ant_index].real, DxSpec[ant_index].imag, DySpec[ant_index].real, DySpec[ant_index].imag])
        np.save(prefixList[0] + '-SPW' + `spw` + '-' + antList[antMap[ant_index]] + '.DSpec.npy', DtermFile)
    #
    plt.close('all')
    #-------- Plot Stokes spectra
    polLabel, Pcolor = ['I', 'Q', 'U', 'V'], ['black', 'blue', 'red', 'green']
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        pp = PdfPages('SP_' + prefixList[0] + '-REF' + refantName + '-' + sourceName + '-SPW' + `spw` + '.pdf')
        figSP = plt.figure(figsize = (11, 8))
        figSP.suptitle(prefixList[0] + ' ' + sourceName)
        figSP.text(0.45, 0.05, 'Frequency [GHz]')
        figSP.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90)
        #
        StokesI_SP = figSP.add_subplot( 2, 1, 1 )
        StokesP_SP = figSP.add_subplot( 2, 1, 2 )
        StokesSpec, StokesErr =  np.mean(StokesVis[:,:,timeIndex], axis=2).real, abs(np.mean(StokesVis[:,:,timeIndex], axis=2).imag)
        np.save(prefixList[0] + '-REF' + refantName + '-' + sourceName + '-SPW' + `spw` + '.StokesSpec.npy', StokesSpec)
        np.save(prefixList[0] + '-REF' + refantName + '-' + sourceName + '-SPW' + `spw` + '.StokesErr.npy', StokesErr)
        np.save(prefixList[0] + '-REF' + refantName + '-' + sourceName + '-SPW' + `spw` + '.Freq.npy', Freq)
        #
        IMax = np.max(StokesSpec[0])
        StokesI_SP.plot(Freq[chRange], StokesSpec[0][chRange], ls='steps-mid', label=polLabel[0], color=Pcolor[0])
        StokesP_SP.plot(Freq[chRange], StokesSpec[1][chRange], ls='steps-mid', label=polLabel[1], color=Pcolor[1])
        StokesP_SP.plot(Freq[chRange], StokesSpec[2][chRange], ls='steps-mid', label=polLabel[2], color=Pcolor[2])
        StokesP_SP.plot(Freq[chRange], StokesSpec[3][chRange], ls='steps-mid', label=polLabel[3], color=Pcolor[3])
        StokesI_SP.tick_params(axis='both', labelsize=6)
        StokesP_SP.tick_params(axis='both', labelsize=6)
        StokesI_SP.axis([np.min(Freq[chRange]), max(Freq[chRange]), 0.0, 1.25*IMax])
        StokesP_SP.axis([np.min(Freq[chRange]), max(Freq[chRange]), -0.1*IMax, 0.1*IMax])
        StokesI_SP.text(min(Freq[chRange]), IMax*1.35, sourceName)
        StokesI_SP.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        StokesP_SP.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        plt.show()
        figSP.savefig(pp, format='pdf')
        pp.close()
    #
    '''
#
