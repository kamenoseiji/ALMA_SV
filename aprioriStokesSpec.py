import sys
from scipy import stats
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
from matplotlib.backends.backend_pdf import PdfPages
RADperHzMeterArcsec = 2.0* pi / 299792458 / (180*3600/pi)
#-------- Configure Array
msmd.open(msfile)
print '---Checking array configulation'
antDia = np.ones(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
print '  -- usable antenna checking for BP scan : '
spwList = scnspw
gainFlag = np.ones([antNum])
for spw_index in range(spwNum):
    #-------- Checking usable baselines and antennas
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPScan)
    timeNum, chNum, blNum = Xspec.shape[3], Xspec.shape[1], Xspec.shape[2]; chRange, timeRange = range(int(0.05*chNum), int(0.95*chNum)), range(timeNum-4, timeNum-1)
    for polID in pPol:
        blD, blA = np.apply_along_axis(delay_search, 0, np.mean(Xspec[polID][chRange][:,:,timeRange], axis=2))
        blA = blA / np.sqrt(antDia[ANT0[0:blNum]]* antDia[ANT1[0:blNum]])
        errD, errA = np.where(abs(blD - np.median(blD)) > 4.0)[0].tolist(), np.where(abs(blA - np.median(blA)) > 0.5* np.median(blA))[0].tolist()
        errCount = np.zeros(antNum)
        for bl in set(errD) or set(errA): errCount[ list(Bl2Ant(bl)) ] += 1
        gainFlag[np.where(errCount > 2.5 )[0].tolist()] *= 0.0
    #
#-------- Check D-term files
Dloaded = False
if not 'DPATH' in locals(): DPATH = SCR_DIR
print '---Checking D-term files in ' + DPATH
DantList, noDlist, Dflag = [], [], np.ones([antNum])
Dpath = DPATH + 'DtermB' + `int(UniqBands[band_index][3:5])` + '/'
for ant_index in range(antNum):
    if flagAnt[ant_index] == 1:
        Dfile = Dpath + 'B' + `int(UniqBands[band_index][3:5])` + '-SPW' + `spwList[0]` + '-' + antList[ant_index] + '.DSpec.npy'
        if os.path.exists(Dfile): DantList += [ant_index]
        else: noDlist += [ant_index]; Dflag[ant_index] *= 0.0
    #
#
DantNum, noDantNum = len(DantList), len(noDlist)
print 'Antennas with D-term file (%d):' % (DantNum),
for ant_index in DantList: print '%s ' % antList[ant_index],
print ''
if noDantNum > 0:
    print 'Antennas without D-term file (%d) : ' % (noDantNum),
    for ant_index in noDlist: print '%s ' % antList[ant_index],
    sys.exit(' Run DtermTransfer first!!')
#
#-------- Load Tsys table
TrxList = []
for spw in spwList: TrxList = TrxList + [np.median(np.load(prefix +  '-' + UniqBands[band_index] + '-SPW' + `spw` + '.Trx.npy'), axis=3)]  # TrxList[spw][pol, ch, ant]
TrxFreq  = np.load(Tsysprefix +  '-' + UniqBands[band_index] + '.TrxFreq.npy') # TrxFreq[spw][ch]
TrxAnts  = np.load(Tsysprefix +  '-' + UniqBands[band_index] + '.TrxAnt.npy') # TrxAnts[ant]
Tau0spec = np.load(Tsysprefix +  '-' + UniqBands[band_index] + '.Tau0.npy') # Tau0spec[spw][ch]
Tau0E    = np.load(Tsysprefix +  '-' + UniqBands[band_index] + '.TauE.npy') # Tau0E[spw, atmScan]
atmTimeRef = np.load(Tsysprefix +  '-' + UniqBands[band_index] + '.atmTime.npy') # atmTimeRef[atmScan]
TrxMap = indexList(TrxAnts, antList); TrxFlag = np.zeros([antNum]); TrxFlag[TrxMap] = 1.0
Tau0E = np.nanmedian(Tau0E, axis=0); Tau0E[np.isnan(Tau0E)] = np.nanmedian(Tau0E); Tau0E[np.isnan(Tau0E)] = 0.0
for spw_index in range(spwNum):
    TrxMed = np.median(TrxList[spw_index], axis=1)  # TrxMed[pol, ant]
    for pol_index in range(2): TrxFlag[np.where(abs(TrxMed[pol_index] - np.median(TrxMed[pol_index])) > 0.8* np.median(TrxMed[pol_index]))[0].tolist()] *= 0.0
if np.min(np.median(Tau0spec[:,chRange], axis=1)) < 0.0: TrxFlag *= 0.0    # Negative Tau(zenith)
tempAtm = GetTemp(msfile)
'''
#-------- Tsys channel interpolation
chNum, chWid, Freq = GetChNum(msfile, spwList[0])
if TrxFreq.shape[1] != chNum:
    tmpTAU0, tmpTRX = np.zeros([spwNum, chNum]), np.zeros([spwNum, antNum, 2, chNum])
    for spw_index in range(spwNum):
        chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq *= 1.0e-9
        TAU0 = interpolate.interp1d(TrxFreq[spw_index], Tau0spec[spw_index])
        tmpTAU0[spw_index] = TAU0(Freq)
        for ant_index in range(len(TrxAnts)):
            for pol_index in range(2):
                TRX = interpolate.interp1d(TrxFreq[spw_index],TrxList[spw_index][ant_index, pol_index])
                tmpTRX[spw_index, ant_index, pol_index] = TRX(Freq)
        #
    #
    Tau0spec = tmpTAU0
    Trxspec  = tmpTRX
#
for spw_index in range(spwNum):
    for pol_index in range(2): TrxFlag[np.where(TrxMed[spw_index][:,pol_index] - np.median(TrxMed[spw_index][:,pol_index]) > 1.5* np.median(TrxMed[spw_index][:,pol_index]))[0].tolist()] *= 0.0
    for pol_index in range(2): TrxFlag[np.where(TrxMed[spw_index][:,pol_index] < 0.3* np.median(TrxMed[spw_index][:,pol_index]))[0].tolist()] *= 0.0
if np.min(np.median(Tau0spec[:,chRange], axis=1)) < 0.0: TrxFlag *= 0.0    # Negative Tau(zenith) 
'''
#
#-------- Time-variable excess zenith opacity
if len(atmTimeRef) > 5:
    exTauSP  = UnivariateSpline(atmTimeRef, Tau0E, np.ones(len(atmTimeRef)), s=0.1*np.std(Tau0E))
else:
    tempTime = np.arange(np.min(atmTimeRef) - 3600.0,  np.max(atmTimeRef) + 3600.0, 300.0)
    tempTauE = np.repeat(np.median(Tau0E), len(tempTime))
    exTauSP = UnivariateSpline(tempTime, tempTauE, np.ones(len(tempTime)), s=0.1)
#
print 'Ant:',
for ant_index in range(antNum): print antList[ant_index],
print; print 'givn',
for ant_index in range(antNum): print '   %.0f' % (flagAnt[ant_index]),
print; print 'Trx ',
for ant_index in range(antNum): print '   %.0f' % (TrxFlag[ant_index]),
print; print 'gain',
for ant_index in range(antNum): print '   %.0f' % (gainFlag[ant_index]),
print; print 'Dtrm',
for ant_index in range(antNum): print '   %.0f' % (Dflag[ant_index]),
print
flagAnt = flagAnt* TrxFlag* gainFlag* Dflag
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
if len(UseAnt) < 5: sys.exit('Too few usable antennas. Reduction failed.')
#-------- Check Scans for atmCal
ingestFile = open(prefix + '-' + UniqBands[band_index] + '-Ingest.log', 'w')
text_sd = '#source,   RA,eRA,dec,edec,frequency,  flux,eflux,degree,edeg,EVPA,eEVPA,uvmin,uvmax,date\n'; ingestFile.write(text_sd)
logfile = open(prefix + '-' + UniqBands[band_index] + '-Flux.log', 'w')
logfile.write(BPcalText + '\n'); logfile.write(EQcalText + '\n')
scanList = onsourceScans
#-------- Array Configuration
if 'refant' in locals():
    refantID = np.where(antList == refant)[0][0]
else:
    print '---Determining refant'
    timeStamp, UVW = GetUVW(msfile, spwList[0], msmd.scansforspw(spwList[0])[0])
    uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = bestRefant(uvDist, UseAnt)
print '  Use ' + antList[refantID] + ' as the refant.'
#
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
#-------- Load D-term file
DxList, DyList = [], []
for ant_index in range(UseAntNum):
    for spw_index in range(spwNum):
        Dfile = Dpath + 'B' + `int(UniqBands[band_index][3:5])` + '-SPW' + `spwList[spw_index]` + '-' + antList[UseAnt[ant_index]] + '.DSpec.npy'
        # print 'Loading %s' % (Dfile)
        Dterm = np.load(Dfile)
        DxList = DxList + [Dterm[1] + (0.0 + 1.0j)* Dterm[2]]
        DyList = DyList + [Dterm[3] + (0.0 + 1.0j)* Dterm[4]]
    #
#
chNum = np.array(DxList).shape[1]
DxSpec, DySpec = np.array(DxList).reshape([UseAntNum, spwNum, chNum]), np.array(DyList).reshape([UseAntNum, spwNum, chNum])
Dloaded = True
if 'gainRef' in locals(): flagRef = np.zeros([UseAntNum]); refIndex = indexList(gainRef, antList[antMap]); flagRef[refIndex] = 1.0; del(gainRef)
if 'refIndex' not in locals(): refIndex = range(UseAntNum)
if len(refIndex) == 0: refIndex = range(UseAntNum)
#-------- Load Aeff file
Afile = open(SCR_DIR + 'AeB' + `int(UniqBands[band_index][3:5])` + '.data')
Alines = Afile.readlines()
Afile.close()
AeX, AeY = 0.25* np.pi* antDia**2, 0.25* np.pi* antDia**2       # antenna collecting area (100% efficiency)
AeX, AeY, etaX, etaY = [], [], [], []
for ant_index in antMap:
    for Aline in Alines:
        if antList[ant_index] in Aline:
            etaX = etaX + [float(Aline.split()[1])]
            etaY = etaY + [float(Aline.split()[2])]
            AeX = AeX + [(0.0025* np.pi* float(Aline.split()[1]))* antDia[ant_index]**2]
            AeY = AeY + [(0.0025* np.pi* float(Aline.split()[2]))* antDia[ant_index]**2]
            #print '%s  : etaX = %.2f  etaY = %.2f' % (antList[ant_index], float(Aline.split()[1]), float(Aline.split()[2]))
        #
    #
#
Ae = np.array([AeX, AeY])
Trx2antMap = indexList( antList[antMap], antList[TrxMap] )
#AeX, AeY = np.array(AeX), np.array(AeY) # in antMap order 
#-------- Flag table
if 'FGprefix' in locals():
    print '---Checking Flag File'
    FGList = []
    for spw_index in range(spwNum): FG = np.load(FGprefix + '-SPW' + `spwList[spw_index]` + '.FG.npy'); FGList = FGList + [np.min(FG, axis=0)]
    FG = np.min( np.array(FGList), axis=0)
    TS = np.load(FGprefix + '-SPW' + `spwList[spw_index]` + '.TS.npy')
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
    flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
else :
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
    flagIndex = range(timeNum)
#
BPList, XYList = [], []
if 'BPprefix' in locals():
    for spw_index in range(spwNum):
        BP_ant = np.load(BPprefix + '-REF' + refant + '-SPW' + `spwList[spw_index]` + '-BPant.npy')
        exp_Tau = np.exp(-Tau0spec[spw_index] / np.sin(BPEL))
        atmCorrect = 1.0 / exp_Tau
        TsysBPScan = atmCorrect* (TrxList[spw_index].transpose(2,0,1)[Trx2antMap] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau)) # [antMap, pol, ch]
        TsysBPShape = (TsysBPScan.transpose(2,0,1) / np.median(TsysBPScan, axis=2)).transpose(1,2,0)
        BPList = BPList + [BP_ant* np.sqrt(TsysBPShape)]
    #
else:
    print '---Generating antenna-based bandpass table'
    for spw_index in range(spwNum):
        BP_ant, XY_BP, XYdelay, Gain = BPtable(msfile, spwList[spw_index], BPScan, blMap, blInv)
        BP_ant[:,1] *= XY_BP
        exp_Tau = np.exp(-Tau0spec[spw_index] / np.sin(BPEL))
        atmCorrect = 1.0 / exp_Tau
        TsysBPScan = atmCorrect* (TrxList[spw_index].transpose(2,0,1)[Trx2antMap] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau)) # [antMap, pol, ch]
        TsysBPShape = (TsysBPScan.transpose(2,0,1) / np.median(TsysBPScan, axis=2)).transpose(1,2,0)
        BPList = BPList + [BP_ant* np.sqrt(TsysBPShape)]
    #
#
if PLOTBP:
    pp = PdfPages('BP_' + prefix + '_REF' + antList[UseAnt[refantID]] + '_Scan' + `BPScan` + '.pdf')
    plotBP(pp, prefix, antList[antMap], spwList, BPScan, BPList)
#
BPDone = True
#
##-------- Gain solutions using polcal scans
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
if 'PolCalScans' in locals():
    print 'Gain Correction using Polarization Calibrator '
    #-------- Check scan time
    mjdSec, Az, El = [], [], []
    for scan in PolCalScans:
        interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], scan)
        scanAz, scanEl = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
        mjdSec = mjdSec + timeStamp.tolist()
        Az = Az + scanAz.tolist()
        El = El + scanEl.tolist()
    #
    mjdSec, Az, El = np.array(mjdSec), np.array(Az), np.array(El)
    PA = AzEl2PA(Az, El) + BandPA[band_index]
    timeNum = len(mjdSec)
    SP_GR, SP_XYPH = [], []
    #-------- Load QU solutions for each SPW
    for spw_index in range(spwNum):
        StokesModel = np.array([1.0] + np.load(QUprefix + '-SPW' + `spwList[spw_index]` + '-' + refant + '.QUXY.npy').tolist() + [0.0])
        BP_bl = BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate()  # BP_bl[BL, corr, ch]
        #-------- Load visibilities of polCal scans
        VisSpec  = np.zeros([4, chNum, UseBlNum, timeNum], dtype=complex)
        for scan in PolCalScans:
            timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], scan, -1, True)  # Xspec[POL, CH, BL, TIME]
            tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3, 2, 0, 1)
            timeIndex = indexList(timeStamp, mjdSec)
            exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(El[timeIndex])))
            TsysScan = (TrxList[spw_index].transpose(2,0,1) + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))[Trx2antMap] # [antMap, pol, ch]
            aprioriSEFD = 2.0* kb* ((TsysScan / exp_Tau).transpose(2,1,0) / Ae).transpose(2,1,0)
            VisSpec[:,:,:,timeIndex] = (tempSpec / BP_bl * np.sqrt(aprioriSEFD[ant0][:,polYindex]* aprioriSEFD[ant1][:,polXindex])).transpose(2, 3, 1, 0)
        #
        StokesModel = StokesModel* np.median(np.mean(abs(np.mean(VisSpec[[0,3]][:,chRange], axis=1)), axis=0))
        print 'SPW ' + `spwList[spw_index]` + ' : Stokes Model [I, Q, U, V] = ' + `StokesModel.tolist()`
        #-------- D-term-corrected visibilities
        PS = PAVector(PA, np.ones(timeNum)).transpose(2,0,1).dot(StokesModel)
        M  = MullerVector(DxSpec[ant0][:,spw_index], DySpec[ant0][:,spw_index], DxSpec[ant1][:,spw_index], DySpec[ant1][:,spw_index], np.ones([UseBlNum,chNum])).transpose(0,3,2,1)
        DcaledChAvg = np.mean((VisSpec / M.dot(PS.T))[:,chRange], axis=1)
        Gain = np.array([gainComplexVec(DcaledChAvg[0]), gainComplexVec(DcaledChAvg[3])])
        GainCaledVis = DcaledChAvg / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
        XYvis = np.mean(GainCaledVis[[1,2]], axis=1)
        XYphase = np.angle(XYvis[0] + XYvis[1].conjugate())
        #SP_XYPH = SP_XYPH + [UnivariateSpline(mjdSec, XYphase, s=0.1)]
        SP_XYPH = SP_XYPH + [UnivariateSpline(mjdSec, XYphase)]
        for ant_index in range(UseAntNum):
            #SP_GR   = SP_GR   + [UnivariateSpline(mjdSec, abs(Gain[1,ant_index]/Gain[0,ant_index]), s=0.1)]
            SP_GR   = SP_GR   + [UnivariateSpline(mjdSec, abs(Gain[1,ant_index]/Gain[0,ant_index]))]
        #
        timePlot = arange(mjdSec[0], mjdSec[-1], 1)
        plt.plot(timePlot, SP_XYPH[spw_index](timePlot), ls='steps-mid')
        plt.plot(mjdSec, XYphase, '.')
    #
#
FreqList = []
for spw_index in range(spwNum):
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index])
    FreqList = FreqList + [Freq* 1.0e-9]
#
print '---Flux densities of sources ---'
pp, polLabel, Pcolor = PdfPages('SP_' + prefix + '_' + UniqBands[band_index] + '.pdf'), ['I', 'Q', 'U', 'V'], ['black', 'blue', 'red', 'green']
figSP = plt.figure(figsize = (11, 8))
figSP.suptitle(prefix + ' ' + UniqBands[band_index])
figSP.text(0.45, 0.05, 'Frequency [GHz]')
figSP.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90)
#-------- Stokes Spectrum for each scan
for scan_index in range(scanNum):
    scan = scanList[scan_index]
    if scan_index > 0:
        for SP in IList: figSP.delaxes(SP)
        for SP in PList: figSP.delaxes(SP)
    #
    timeStamp, UVW = GetUVW(msfile, spwList[0], scan);  timeNum = len(timeStamp)
    uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    #-------- Flagging
    if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
    else: flagIndex = range(timeNum)
    #-------- Parallactic Angle
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    PA = (AzEl2PA(AzScan, ElScan) + BandPA[band_index])[flagIndex]; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
    BPCaledXspec, IList, PList = [], [], []
    text_time = qa.time('%fs' % np.median(timeStamp), form='ymd')[0]
    text_src  = ' %02d %010s EL=%4.1f deg PA=%6.1f deg' % (scanList[scan_index], sourceList[sourceIDscan[scan_index]], 180.0* np.median(ElScan)/pi, 180.0* np.median(AzEl2PA(AzScan, ElScan))/pi); logfile.write(text_src + ' ' + text_time + '\n'); print text_src + ' ' + text_time
    #-------- Subarray formation
    SAantMap, SAblMap, SAblInv, SAant0, SAant1 = antMap, blMap, blInv, ant0, ant1
    bpAntMap = indexList(antList[SAantMap],antList[antMap])
    Trx2antMap = indexList( antList[SAantMap], antList[TrxMap] )
    #-------- Baseline-based cross power spectra
    IMax = 0.0
    for spw_index in range(spwNum):
        if 'Field' in locals():
            timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], scan, Field)
        else:
            timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], scan)
        timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); UseChNum = len(chRange)
        if np.max(abs(Xspec)) < 1.0e-9: continue
        XYtwiddle = np.exp((1.0j)* SP_XYPH[spw_index](timeStamp))
        #-------- Position offset phase correction
        if 'offAxis' in locals():
            lm = np.array(offAxis[scan])
            Twiddle =  np.exp((0.0 + 1.0j)* np.outer(FreqList[spw_index]*1.0e9, uvw[0:2].transpose(1,2,0).dot(lm)).reshape([chNum, UseBlNum, timeNum])* RADperHzMeterArcsec)
            tempSpec = CrossPolBL(Xspec[:,:,SAblMap]*Twiddle, SAblInv)
            tempSpec[1] /= XYtwiddle
            tempSpec[2] *= XYtwiddle
            tempSpec = tempSpec.transpose(3,2,0,1)[flagIndex]      # Cross Polarization Baseline Mapping
        else:
            tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv)
            tempSpec[1] /= XYtwiddle 
            tempSpec[2] *= XYtwiddle 
            tempSpec = tempSpec.transpose(3,2,0,1)[flagIndex]
        #-------- Bandpass Calibration
        BPCaledXspec = (tempSpec / (BPList[spw_index][SAant0][:,polYindex]* BPList[spw_index][SAant1][:,polXindex].conjugate())).transpose(2,3,1,0)
        chAvgVis = np.mean(BPCaledXspec[[0,3]][:,chRange], axis=1) # chAvgVis[pol, bl, time]
        GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
        for ant_index in range(UseAntNum):
            GainRatio = np.sqrt(SP_GR[spw_index* UseAntNum + ant_index](timeStamp))     # sqrt(GainY / GainX)
            GainP[0,ant_index] /= np.sqrt(GainRatio)
            GainP[1,ant_index] *= np.sqrt(GainRatio)
        #
        pCalVis = (BPCaledXspec.transpose(1,0,2,3) / (GainP[polYindex][:,SAant0]* GainP[polXindex][:,SAant1].conjugate()))[chRange]
        #-------- SEFD amplitude calibration
        StokesI_SP = figSP.add_subplot( 2, spwNum, spw_index + 1 )
        StokesP_SP = figSP.add_subplot( 2, spwNum, spwNum + spw_index + 1 )
        IList = IList + [StokesI_SP]
        PList = PList + [StokesP_SP]
        exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
        atmCorrect = 1.0 / exp_Tau
        TsysSPW = (TrxList[spw_index].transpose(2,0,1) + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))[Trx2antMap]    # [ant, pol, ch]
        SEFD = 2.0* kb* (TsysSPW * atmCorrect).transpose(2,1,0) / Ae[:,bpAntMap]                # SEFD[ch,pol,antMap]
        SAantNum = len(SAantMap); SAblNum = len(SAblMap)
        #AmpCalVis = np.mean(np.mean(pCalVis, axis=3)[:,[0,3]] * np.sqrt(SEFD[chRange][:,:,ant0[0:SAblNum]]* SEFD[chRange][:,:,ant1[0:SAblNum]]), axis=0)
        #indivRelGain = abs(gainComplexVec(AmpCalVis.T)); indivRelGain /= np.percentile(indivRelGain, 75, axis=0)
        #indivRelGain = abs(gainComplexVec(AmpCalVis.T)); indivRelGain /= np.median(indivRelGain, axis=0)
        #SEFD /= (indivRelGain**2).T
        AmpCalVis = (pCalVis.transpose(3,0,1,2)* np.sqrt(SEFD[chRange][:,polYindex][:,:,SAant0]* SEFD[chRange][:,polXindex][:,:,SAant1])).transpose(3,2,1,0)
        #-------- Stokes Spectra
        Stokes = np.zeros([4,blNum, UseChNum], dtype=complex)  # Stokes[stokes, bl, ch]
        for bl_index in range(SAblNum):
            Minv = InvMullerVector(DxSpec[SAant1[bl_index], spw_index][chRange], DySpec[SAant1[bl_index], spw_index][chRange], DxSpec[SAant0[bl_index], spw_index][chRange], DySpec[SAant0[bl_index], spw_index][chRange], np.ones(UseChNum))
            for time_index in range(timeNum):
                for ch_index in range(UseChNum):
                    Stokes[:,bl_index,ch_index] += PS[:,:,time_index].dot(Minv[:,:,ch_index].dot(AmpCalVis[bl_index, :, ch_index, time_index]))
                #
        #
        Stokes = np.mean(Stokes, axis=1) / timeNum
        StokesSpec, StokesErr = Stokes.real, abs(Stokes.imag)
        IMax = max(IMax, np.max(StokesSpec[0]))
        np.save(prefix + '-' + UniqBands[band_index] + '-SPW' + `spwList[spw_index]` + '-Scan' + `scan` + '.StokesSpec.npy', StokesSpec)
        np.save(prefix + '-' + UniqBands[band_index] + '-SPW' + `spwList[spw_index]` + '-Scan' + `scan` + '.StokesErr.npy', StokesErr)
        np.save(prefix + '-' + UniqBands[band_index] + '-SPW' + `spwList[spw_index]` + '-Scan' + `scan` + '.Freq.npy', np.array(FreqList)[:,chRange])
        StokesI_SP.plot(FreqList[spw_index][chRange], StokesSpec[0], ls='steps-mid', label=polLabel[0], color=Pcolor[0])
        # StokesI_SP.errorbar(FreqList[spw_index][chRange], StokesSpec[0], yerr=StokesErr[0], fmt='', ecolor=Pcolor[0], capsize=0)
        StokesP_SP.plot(FreqList[spw_index][chRange], StokesSpec[1], ls='steps-mid', label=polLabel[1], color=Pcolor[1])
        StokesP_SP.plot(FreqList[spw_index][chRange], StokesSpec[2], ls='steps-mid', label=polLabel[2], color=Pcolor[2])
        StokesP_SP.plot(FreqList[spw_index][chRange], StokesSpec[3], ls='steps-mid', label=polLabel[3], color=Pcolor[3])
        # StokesP_SP.errorbar(FreqList[spw_index][chRange], StokesSpec[1], yerr=StokesErr[1], fmt='', ecolor=Pcolor[1], capsize=0)
        # StokesP_SP.errorbar(FreqList[spw_index][chRange], StokesSpec[2], yerr=StokesErr[2], fmt='', ecolor=Pcolor[2], capsize=0)
        # StokesP_SP.errorbar(FreqList[spw_index][chRange], StokesSpec[3], yerr=StokesErr[3], fmt='', ecolor=Pcolor[3], capsize=0)
    #
    for spw_index in range(spwNum):
        StokesI_SP, StokesP_SP = IList[spw_index], PList[spw_index]
        StokesI_SP.tick_params(axis='both', labelsize=6)
        StokesP_SP.tick_params(axis='both', labelsize=6)
        StokesI_SP.axis([min(FreqList[spw_index]), max(FreqList[spw_index]), 0.0, 1.25*IMax])
        StokesP_SP.axis([min(FreqList[spw_index]), max(FreqList[spw_index]), -0.25*IMax, 0.25*IMax])
    #
    IList[0].text(min(FreqList[0]), IMax*1.35, text_src)
    IList[-1].text(max(FreqList[-1]), IMax*1.35, text_time, ha='right')
    StokesI_SP.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    StokesP_SP.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    plt.show()
    figSP.savefig(pp, format='pdf')
#
logfile.close()
plt.close('all')
pp.close()
