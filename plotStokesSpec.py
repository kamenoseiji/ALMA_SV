import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
from matplotlib.backends.backend_pdf import PdfPages
#RADperHzMeterArcsec = 2.0* pi / 299792458 / (180*3600/pi)
#-------- Configure Array
msfile = prefix + '.ms'; msmd.open(msfile)
"""
azelTime, AntID, AZ, EL = GetAzEl(msfile)
antList = GetAntName(msfile)
antNum = len(antList)
blNum = antNum* (antNum - 1) / 2
spwNum = len(spwList)
spwName = msmd.namesforspws(spwList[0])[0]
BandName = re.findall(r'RB_..', spwName)[0]; BandID = int(BandName[3:5])
BandPA = (BANDPA[int(BandName[3:5])] + 90.0)*pi/180.0
#-------- Configure Array
print '---Checking array configulation'
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
if 'antFlag' in locals(): flagAnt[indexList(antFlag, antList)] = 0.0; del(antFlag)
Tatm_OFS  = 5.0     # Ambient-load temperature - Atmosphere temperature
kb        = 1.38064852e3
msmd.close()
msmd.done()
#-------- Load Tsys table
logfile = open(prefix + '-' + BandName + '-Flux.log', 'w')
if 'QUfile' in locals(): QUmodel = np.load(QUfile + '.QUXY.npy')
Tau0spec = np.load(prefix +  '-' + BandName + '.Tau0.npy') # Tau0spec[spw][ch]
Trxspec  = np.load(prefix +  '-' + BandName + '.Trx.npy')  # Trxspec[ant*spw][pol, ch]
OnEL = np.load(prefix +  '-' + BandName + '.OnEL.npy')
AtmEL = np.load(prefix +  '-' + BandName+ '.AtmEL.npy')
OnEL = np.median(OnEL, axis=0)
UseAntNum = AtmEL.shape[0]; TrxSPWnum = Trxspec.shape[0]/UseAntNum
chAvgTrx = np.median(Trxspec, axis=2).reshape([UseAntNum, TrxSPWnum, 2])
msmd.open(msfile)
#-------- Array Configuration
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt)
antDia = np.ones(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
flagList = np.where(np.median(chAvgTrx.reshape(UseAntNum, 2* spwNum), axis=1) > 2.0* np.percentile(chAvgTrx, 75))[0].tolist()  # Avoid too-high Trx
flagList = unique(flagList + np.where(np.min(chAvgTrx.reshape(UseAntNum, 2* spwNum), axis=1) < 1.0 )[0].tolist()).tolist()     # Avoid too-low Trx
if len(flagList) >0 :
    for index in flagList: del UseAnt[index]
#
UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
text_sd = '  Usable antennas: '
for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
logfile.write(text_sd + '\n'); print text_sd
text_sd = '  Flagged by Trx:  '
for ants in antList[flagList].tolist(): text_sd = text_sd + ants + ' '
logfile.write(text_sd + '\n'); print text_sd
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
timeStamp, UVW = GetUVW(msfile, spwList[0], msmd.scansforspw(spwList[0])[0])
uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
try:
    refantID = np.where(antList[UseAnt] == refant )[0][0]
except:
    refantID = bestRefant(uvDist)
#
refantName = antList[UseAnt[refantID]]
print '  Use ' + refantName + ' as the refant.'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
if 'gainRef' in locals():
    flagRef = np.zeros([UseAntNum])
    refIndex = indexList(gainRef, antList[antMap])
    flagRef[refIndex] = 1.0; del(gainRef)
else:
    refIndex = range(UseAntNum)
#
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
#-------- Load Aeff file
Afile = open(SCR_DIR + 'AeB' + `int(BandName[3:5])` + '.data')
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
AeX, AeY = np.array(AeX), np.array(AeY) # in antMap order 
#-------- Check D-term files
if not 'DPATH' in locals(): DPATH = SCR_DIR
print '---Checking D-term files in ' + DPATH
DantList, noDlist = [], []
Dpath = DPATH + 'DtermB' + `int(BandName[3:5])` + '/'
for ant_index in UseAnt:
    Dfile = Dpath + 'B' + `int(BandName[3:5])` + '-SPW' + `spwList[0]` +'-' + antList[ant_index] + '.DSpec.npy'
    if os.path.exists(Dfile): DantList += [ant_index]
    else: noDlist += [ant_index]
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
#-------- Load D-term file
DxList, DyList = [], []
print '---Loading D-term table'
for ant_index in range(UseAntNum):
    # Dpath = SCR_DIR + 'DtermB' + `int(BandName[3:5])` + '/'
    for spw in spwList:
        Dfile = Dpath + 'B' + `int(BandName[3:5])` + '-SPW' + `spw` + '-' + antList[antMap[ant_index]] + '.DSpec.npy'
        Dterm = np.load(Dfile)
        DxList = DxList + [Dterm[1] + (0.0 + 1.0j)* Dterm[2]]
        DyList = DyList + [Dterm[3] + (0.0 + 1.0j)* Dterm[4]]
    #
#
chNum = np.array(DxList).shape[1]
DxSpec, DySpec = np.array(DxList).reshape([UseAntNum, spwNum, chNum]), np.array(DyList).reshape([UseAntNum, spwNum, chNum])
#
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
#-------- Bandpass Table
FreqList = []
for spw_index in range(spwNum):
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* Freq  # GHz
    FreqList = FreqList + [Freq]
#
BPList = []
if 'BPprefix' in locals():      # Load Bandpass table
    print '---Loding bandpass table...'
    for spw_index in spwList:
        BPantList, BP_ant, XYspec = np.load(BPprefix + '-REF' + refantName + '.Ant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw_index` + '-BPant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw_index` + '-XYspec.npy')
        BP_ant[:,1] *= XYspec
        BPList = BPList + [BP_ant]
    #
else:
    print '---Generating antenna-based bandpass table'
    for spw_index in spwList:
        BP_ant, XY_BP, XYdelay, Gain = BPtable(msfile, spw_index, BPScan, blMap, blInv)
        BP_ant[:,1] *= XY_BP
        BPList = BPList + [BP_ant]
    #
#
BPDone = True
if 'PLOTBP' not in locals(): PLOTBP = False
if PLOTBP:
    pp = PdfPages('BP_' + prefix + '_REF' + antList[UseAnt[refantID]] + '_Scan' + `BPScan` + '.pdf')
    plotMax = 1.5* np.median(abs(BP_ant))
    #-------- Prepare Plots
    for ant_index in range(UseAntNum):
        figAnt = plt.figure(ant_index, figsize = (11, 8))
        figAnt.suptitle(prefix + ' ' + antList[antMap[ant_index]] + ' Scan = ' + `BPScan`)
        figAnt.text(0.45, 0.05, 'Frequency [GHz]')
        figAnt.text(0.03, 0.45, 'Bandpass Amplitude and Phase', rotation=90)
    #
    #-------- Plot BP
    for ant_index in range(UseAntNum):
        figAnt = plt.figure(ant_index)
        for spw_index in range(spwNum):
            # chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* Freq  # GHz
            Freq = FreqList[spw_index]
            BPampPL = figAnt.add_subplot( 2, spwNum, spw_index + 1 )
            BPphsPL = figAnt.add_subplot( 2, spwNum, spw_index + spwNum + 1 )
            for pol_index in range(ppolNum):
                plotBP = BPList[spw_index][ant_index, pol_index]
                BPampPL.plot( Freq, abs(plotBP), ls='steps-mid', label = 'Pol=' + PolList[pol_index])
                BPampPL.axis([np.min(Freq), np.max(Freq), 0.0, 1.25* plotMax])
                BPampPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
                BPampPL.yaxis.offsetText.set_fontsize(10)
                BPampPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
                BPphsPL.plot( Freq, np.angle(plotBP), '.', label = 'Pol=' + PolList[pol_index])
                BPphsPL.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
            #
            BPampPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            BPphsPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
            BPampPL.text( np.min(Freq), 1.1* plotMax, 'SPW=' + `spwList[spw_index]` + ' Amp')
            BPphsPL.text( np.min(Freq), 2.5, 'SPW=' + `spwList[spw_index]` + ' Phase')
        #
        figAnt.savefig(pp, format='pdf')
    #
    plt.close('all')
    pp.close()
#
##-------- Equalization using EQ scan
relGain = np.ones([spwNum, 2, UseAntNum])
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA; PA = np.arctan2( np.sin(PA), np.cos(PA))
QCpUS = QUmodel[0]* np.cos(2.0* PA) + QUmodel[1]* np.sin(2.0* PA)
useAntMapRev = indexList(np.array(antMap), np.array(UseAnt))
exp_Tau = np.exp(-Tau0spec / np.mean(np.sin(ElScan)))
tempAtm = GetTemp(msfile)
TsysEQScan = (np.mean(Trxspec[:,:,chRange],axis=2).reshape([UseAntNum, spwNum, 2]).transpose(2,0,1) + Tcmb* np.mean(exp_Tau, axis=1) + tempAtm* (1.0 - np.mean(exp_Tau[:,chRange], axis=1))).transpose(1,2,0)[useAntMapRev]
for spw_index in range(spwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)[flagIndex]       # Cross Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)[pPol]
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
    pCaledVis = np.array([chAvgVis[0] / (GainP[0,ant0]* GainP[0,ant1].conjugate()), chAvgVis[1]/(GainP[1,ant0]* GainP[1,ant1].conjugate())])
    aprioriSEFD = 2.0* kb* TsysEQScan[:,spw_index].T / np.array([AeX, AeY])
    #
    aprioriVisX = np.mean(pCaledVis[0] / (1.0 + QCpUS), axis=1) * np.sqrt(aprioriSEFD[0, ant0]* aprioriSEFD[0, ant1])
    aprioriVisY = np.mean(pCaledVis[1] / (1.0 - QCpUS), axis=1) * np.sqrt(aprioriSEFD[1, ant0]* aprioriSEFD[1, ant1])
    #-------- Determine Antenna-based Gain
    relGain[spw_index, 0] = abs(gainComplex(aprioriVisX)); relGain[spw_index, 0] /= np.median( abs(relGain[spw_index, 0, refIndex]) ) # X-pol delta gain
    relGain[spw_index, 1] = abs(gainComplex(aprioriVisY)); relGain[spw_index, 1] /= np.median( abs(relGain[spw_index, 1, refIndex]) ) # Y-pol delta gain
#
print '---Equalized aperture efficiencies (Pol-X, Pol-Y) in %'
antRelGain = np.median(relGain, axis=0)
for ant_index in range(UseAntNum):
    print '%s  %.2f  %.2f' % (antList[antMap[ant_index]], etaX[ant_index]* antRelGain[0,ant_index]**2, etaY[ant_index]* antRelGain[1,ant_index]**2)
#
##-------- Iteration for Equalization using EQ scan
print '---XY phase determination in bandpass scan'
#-------- XY phase using BP scan
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], BPScan); timeNum = len(timeStamp)
if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
else: flagIndex = range(timeNum)
#
AzScan, ElScan = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA; PA = np.arctan2( np.sin(PA), np.cos(PA))
GainP, XYphase, caledVis = [], [], []
TsysShape = []
for spw_index in range(spwNum):
    atmCorrect = np.exp(Tau0spec[spw_index] / np.mean(np.sin(ElScan)))
    exp_Tau = 1.0 / atmCorrect
    TsysSPW =  (Trxspec[spw_index::spwNum].transpose(1,0,2) + Tcmb* exp_Tau + tempAtm * (1.0 - exp_Tau))[:,useAntMapRev]
    TsysShape = TsysShape + [(TsysSPW.transpose(2,0,1) / np.mean(TsysSPW, axis=2)).transpose(1,2,0)]    # Normalized Tsys spectrum
    TsysBL  = np.sqrt( TsysSPW[polYindex][:,ant0]* TsysSPW[polXindex][:,ant1] ).transpose(1,0,2)
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPScan); timeNum = len(timeStamp)
    if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
    else : flagIndex = range(timeNum)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)[flagIndex]      # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec * TsysBL/ (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainP = GainP + [np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])]
    SEFD = 2.0* kb / (np.array([AeX, AeY]) * (relGain[spw_index]**2))
    caledVis.append(np.mean((chAvgVis / (GainP[spw_index][polYindex][:,ant0]* GainP[spw_index][polXindex][:,ant1].conjugate())).transpose(2, 0, 1)* np.sqrt(SEFD[polYindex][:,ant0]* SEFD[polXindex][:,ant1]), axis=2).T)
#
caledVis = np.array(caledVis)
#-------- SPW-specific phase using BP scan
GainP = np.array(GainP) # GainP[spw, pol, ant, time]
spwPhase = [0.0]* 2* spwNum
for ant_index in range(1,UseAntNum):
    for pol_index in range(2):
        spwPhase = spwPhase + [0.0]
        for spw_index in range(1,spwNum): spwPhase = spwPhase + [np.angle(GainP[spw_index, pol_index, ant_index].dot(GainP[0, pol_index, ant_index].conjugate()))]
    #
#
spwPhase = np.array(spwPhase).reshape([UseAntNum, 2, spwNum]); spwTwiddle = exp(1.0j *spwPhase)
#-------- XY phase cal in Bandpass table
XYsign = np.ones(spwNum)
for spw_index in range(spwNum):
    XYphase = XY2Phase(PA, QUmodel[0], QUmodel[1], caledVis[spw_index][[1,2]])
    XYsign[spw_index] = np.sign(np.cos(XYphase))
    BPList[spw_index] = (BPList[spw_index].transpose(2,0,1)* spwTwiddle[:,:,spw_index]).transpose(1,2,0)
    BPList[spw_index][:,1] *= XYsign[spw_index]
    print 'SPW[%d] : XYphase = %6.1f [deg] sign = %3.0f' % (spwList[spw_index], 180.0*XYphase/pi, XYsign[spw_index])
#
XYD, XYC = [], []      # XY delay and correlation
for scan in scanList:
    timeStamp, UVW = GetUVW(msfile, spwList[0], scan);  timeNum = len(timeStamp)
    SAantennas, SAbl, SAblFlag, SAant0, SAant1 =  range(UseAntNum), range(UseBlNum), np.ones([blNum]), ant0, ant1
    SAantMap, SAblMap, SAblInv = antMap, blMap, blInv
    SAantNum = len(SAantennas); SAblNum = len(SAblMap)
    #---- AZ, EL, PA
    BPCaledXspec = []
    if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
    else: flagIndex = range(timeNum)
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    PA = (AzEl2PA(AzScan, ElScan) + BandPA)[flagIndex]; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
    #-------- Flagging
    for spw_index in range(spwNum):
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], scan)
        timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); UseChNum = len(chRange)
        if np.max(abs(Xspec)) < 1.0e-9: continue
        #-------- Position offset phase correction
        if 'offAxis' in locals():
            lm = np.array(offAxis[scan])
            Twiddle =  np.exp((0.0 + 1.0j)* np.outer(FreqList[spw_index], uvw[0:2].transpose(1,2,0).dot(lm)).reshape([chNum, UseBlNum, timeNum])* RADperHzMeterArcsec)
            tempSpec = CrossPolBL(Xspec[:,:,SAblMap]*Twiddle, SAblInv).transpose(3,2,0,1)[flagIndex]      # Cross Polarization Baseline Mapping
        else:
            tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)[flagIndex]      # Cross Polarization Baseline Mapping
        #-------- Bandpass Calibration
        BPCaledXspec = BPCaledXspec + [(tempSpec / (BPList[spw_index][SAant0][:,polYindex]* BPList[spw_index][SAant1][:,polXindex].conjugate())).transpose(2,3,1,0)]
    #
    #-------- Antenna-based Phase Solution
    BPCaledXspec = np.array(BPCaledXspec)   # BPCaledXspec[spw, pol, ch, bl, time]
    chAvgVis = np.mean(np.array(BPCaledXspec)[:,:,chRange], axis=(0,2)) # chAvgVis[pol, bl, time]
    if 'timeBunch' in locals():
        useTimeNum = timeNum / timeBunch * timeBunch
        leapNum = timeNum - useTimeNum
        timeAvgVis = np.mean(chAvgVis[:,:,range(useTimeNum)].reshape(polNum, UseBlNum, timeNum / timeBunch, timeBunch), axis=3)
        GainP = np.array([np.apply_along_axis(clphase_solve, 0, timeAvgVis[0]), np.apply_along_axis(clphase_solve, 0, timeAvgVis[3])]).repeat(timeBunch, axis=2)
        if leapNum > 0: GainP = np.append(GainP,  GainP[:,:,(useTimeNum-1):(useTimeNum)].repeat(leapNum, axis=2), axis=2)
    else:
        GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    #
    pCalVis = (BPCaledXspec.transpose(0,2,1,3,4) / (GainP[polYindex][:,SAant0]* GainP[polXindex][:,SAant1].conjugate()))[:,chRange]
    #-------- XY phase spectra
    for spw_index in range(spwNum):
        delayFact = (chNum + 0.0)/len(chRange)
        XYspec = np.mean(pCalVis[spw_index, :, 1:3, :], axis=(2,3)).T
        XYdelay, XYamp = delay_search(XYspec[:,0]); YXdelay, YXamp = delay_search(XYspec[:,1])
        XYD = XYD + [XYdelay* delayFact, YXdelay* delayFact]
        XYC = XYC + [np.mean(delay_cal(XYspec[:,0], XYdelay)), np.mean(delay_cal(XYspec[:,1], YXdelay))]
    #-------- Full-Stokes parameters
    for spw_index in range(spwNum):
        atmCorrect = np.exp(Tau0spec[spw_index] / np.mean(np.sin(ElScan)))
        exp_Tau = 1.0 / atmCorrect
        TsysSPW =  (Trxspec[spw_index::spwNum].transpose(1,0,2) + Tcmb* exp_Tau + tempAtm * (1.0 - exp_Tau))[:,useAntMapRev]
        TsysSPW = TsysSPW / TsysShape[spw_index]
        #---- Flagged by Tsys
        tsysFlagAntIndex = unique(np.where(TsysSPW <0.0)[1]).tolist()
        if len(tsysFlagAntIndex) > 0:
            for ant_index in tsysFlagAntIndex: TsysSPW[:,ant_index] = Trxspec[spw_index::spwNum][ant_index] + tempAtm* (1.0 - exp_Tau)
        #
        SEFD = 2.0* kb* (TsysSPW * atmCorrect).transpose(2,0,1) /  (np.array([AeX[SAantennas], AeY[SAantennas]]) / relGain[spw_index][:,SAantennas]**2)
        SEFD = np.mean(SEFD[chRange], axis=0)
        #AmpCalVis = np.mean(np.mean(pCalVis[spw_index], axis=3)[:,[0,3]] * np.sqrt(SEFD[chRange][:,:,ant0[0:SAblNum]]* SEFD[chRange][:,:,ant1[0:SAblNum]]), axis=0)
        AmpCalVis = np.mean(np.mean(pCalVis[spw_index], axis=3)[:,[0,3]] * np.sqrt(SEFD[:,ant0[0:SAblNum]]* SEFD[:,ant1[0:SAblNum]]), axis=0)
        indivRelGain = abs(gainComplexVec(AmpCalVis.T)); indivRelGain /= np.percentile(indivRelGain, 75, axis=0)
        SEFD /= (indivRelGain**2).T
        #AmpCalVis = (pCalVis[spw_index].transpose(3,0,1,2)* np.sqrt(SEFD[chRange][:,polYindex][:,:,ant0[0:SAblNum]]* SEFD[chRange][:,polXindex][:,:,ant1[0:SAblNum]])).transpose(3,2,1,0)
        AmpCalVis = (pCalVis[spw_index].transpose(3,0,1,2)* np.sqrt(SEFD[polYindex][:,ant0[0:SAblNum]]* SEFD[polXindex][:,ant1[0:SAblNum]])).transpose(3,2,1,0)
        Stokes = np.zeros([4,blNum, UseChNum], dtype=complex)  # Stokes[stokes, bl, ch]
        for bl_index in range(SAblNum):
            Minv = InvMullerVector(DxSpec[SAant1[bl_index], spw_index][chRange], DySpec[SAant1[bl_index], spw_index][chRange], DxSpec[SAant0[bl_index], spw_index][chRange], DySpec[SAant0[bl_index], spw_index][chRange], np.ones(UseChNum))
            for time_index in range(timeNum):
                for ch_index in range(UseChNum):
                    Stokes[:,bl_index,ch_index] += PS[:,:,time_index].dot(Minv[:,:,ch_index].dot(AmpCalVis[bl_index, :, ch_index, time_index]))
                #
        #
        Stokes = np.mean(Stokes, axis=1) / timeNum
        StokesSpec, StokesErr = Stokes.real, Stokes.imag
        np.save(prefix + '-' + BandName + '-SPW' + `spwList[spw_index]` + '-Scan' + `scan` + '.StokesSpec.npy', StokesSpec)
        np.save(prefix + '-' + BandName + '-SPW' + `spwList[spw_index]` + '-Scan' + `scan` + '.StokesErr.npy', StokesErr)
    #
#
logfile.close()
msmd.done()
msmd.close()
"""
chRange = range(int(0.05*chNum), int(0.95*chNum))
#-------- Plot
print '---Plot Stokes Spectrum densities of sources ---'
pp, polLabel, Pcolor = PdfPages('SP_' + prefix + '_' + BandName + '.pdf'), ['I', 'LP', 'CP', 'EVPA'], ['black', 'red', 'green', 'black']
page_index = 0
figScan = plt.figure(page_index, figsize = (11,8))
#figScan.suptitle(prefix + ' ' + BandName + ' Scan' + `scan`)
figScan.suptitle(prefix + ' ' + BandName )
figScan.text(0.45, 0.05, 'Frequency [GHz]') 
# figScan.text(0.03, 0.45, 'Stokes Parameters [Jy]', rotation=90) 
figScan.text(0.03, 0.45, 'EVPA [deg] and (polarized) Flux Desity [Jy]', rotation=90) 
#-------- Ploting I and (Q, U, V) spectra
for spw_index in range(spwNum):
    Freq = FreqList[spw_index]
    StokesI_PL = figScan.add_subplot( 2, spwNum, spw_index+1)
    StokesP_PL = figScan.add_subplot( 2, spwNum, spw_index+spwNum+1)
    StokesList, StokesErrList = [], []
    for scan in scanList:
        StokesSpec = np.load(prefix + '-' + BandName + '-SPW' + `spwList[spw_index]` + '-Scan' + `scan` + '.StokesSpec.npy')
        StokesErr  = np.load(prefix + '-' + BandName + '-SPW' + `spwList[spw_index]` + '-Scan' + `scan` + '.StokesErr.npy')
        StokesList = StokesList + [StokesSpec]
        StokesErrList = StokesErrList + [StokesErr]
    #
    StokesSpec = np.array(StokesList)
    weight = 1.0 / np.mean(np.array(StokesErrList), axis=2)[:,3]
    StokesSpec = StokesSpec.transpose(1,2,0).dot(weight) / np.sum(weight)
    StokesI_PL.plot( Freq[chRange], StokesSpec[0], ls='steps-mid', label=polLabel[0], color=Pcolor[0])
    """
    StokesP_PL.plot( np.array([min(Freq[chRange]), max(Freq[chRange])]), np.array([0.0,0.0]), '-', color='gray')
    StokesP_PL.plot( Freq[chRange], StokesSpec[1], ls='steps-mid', label=polLabel[1], color=Pcolor[1])
    StokesP_PL.plot( Freq[chRange], StokesSpec[2], ls='steps-mid', label=polLabel[2], color=Pcolor[2])
    StokesP_PL.plot( Freq[chRange], StokesSpec[3], ls='steps-mid', label=polLabel[3], color=Pcolor[3])
    """
    #plotMax = max(StokesSpec[0])
    plotMax = 6.0
    #StokesP_PL.axis([min(Freq[chRange]), max(Freq[chRange]), -0.15*plotMax, 0.15*plotMax ])
    LPflux = np.sqrt(StokesSpec[1]**2 + StokesSpec[2]**2)
    StokesI_PL.plot( Freq[chRange], LPflux, ls='steps-mid', label=polLabel[1], color=Pcolor[1])
    StokesI_PL.plot( Freq[chRange], StokesSpec[3], ls='steps-mid', label=polLabel[2], color=Pcolor[2])
    StokesI_PL.axis([min(Freq[chRange]), max(Freq[chRange]), -0.1*plotMax, 1.2* plotMax])
    EVPA = 90.0* np.arctan2(StokesSpec[2], StokesSpec[1]) / pi
    StokesP_PL.plot( Freq[chRange], EVPA, '.', label=polLabel[3], color=Pcolor[3])
    #StokesP_PL.axis([min(Freq[chRange]), max(Freq[chRange]), min(EVPA), max(EVPA) ])
    StokesP_PL.axis([min(Freq[chRange]), max(Freq[chRange]), -70.0, -55.0])
    #
    StokesI_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    StokesP_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    #StokesP_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    #StokesP_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
#
figScan.savefig(pp, format='pdf')
page_index += 1
#
plt.close('all')
pp.close()
