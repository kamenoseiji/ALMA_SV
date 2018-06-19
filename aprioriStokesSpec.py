import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import analysisUtils as au
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
from matplotlib.backends.backend_pdf import PdfPages
#-------- Initial Settings
msfile = prefix + '.ms'; msmd.open(msfile)
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
#-------- Load Tsys File
if 'QUfile' in locals(): QUmodel = np.load(QUfile + '.QUXY.npy')
#Tau0spec = np.load(prefix +  '-' + UniqBands[band_index] + '.Tau0.npy') # Tau0spec[spw][ch]
#Trxspec  = np.load(prefix +  '-' + UniqBands[band_index] + '.Trx.npy')  # Trxspec[ant*spw][pol, ch]
#OnEL = np.load(prefix +  '-' + UniqBands[band_index] + '.OnEL.npy')
#AtmEL = np.load(prefix +  '-' + UniqBands[band_index] + '.AtmEL.npy')
Tau0spec = np.load(TSprefix + '.Tau0.npy') # Tau0spec[spw][ch]
Trxspec  = np.load(TSprefix + '.Trx.npy')  # Trxspec[ant*spw][pol, ch]
OnEL = np.load(TSprefix + '.OnEL.npy')
AtmEL = np.load(TSprefix + '.AtmEL.npy')
OnEL = np.median(OnEL, axis=0)
UseAntNum = AtmEL.shape[0]; TrxSPWnum = Trxspec.shape[0]/UseAntNum
chAvgTrx = np.median(Trxspec, axis=2).reshape([UseAntNum, TrxSPWnum, 2])
#TRchNum = Trec.shape[2]; TRchRange = range(int(0.05*TRchNum), int(0.95*TRchNum))
#chAvgTrx = np.median(np.mean(Trec[:,:,TRchRange], axis=2), axis=2)
#Tau0 = np.median(Tau0)
#-------- Array Configuration
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt)
print '---Checking array configuration'
flagList = np.where(np.max(chAvgTrx, axis=(1,2)) > 2.0* np.percentile(chAvgTrx, 75))[0].tolist()  # Avoid too-high Trx
flagList = unique(flagList + np.where(np.min(chAvgTrx, axis=(1,2)) < Tatm_OFS )[0].tolist()).tolist()     # Avoid too-low Trx
if len(flagList) >0 :
    for index in flagList: del UseAnt[index]
#
UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
text_sd = '  Usable antennas: '
for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
print text_sd
text_sd = '  Flagged by Trx:  '
for ants in antList[flagList].tolist(): text_sd = text_sd + ants + ' '
print text_sd
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
timeStamp, UVW = GetUVW(msfile, spwList[0], msmd.scansforspw(spwList[0])[0])
uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
if 'refant' in locals():    refantID = indexList(np.array([refant]), antList)[0]
else: refantID = bestRefant(uvDist)
print '  Use ' + antList[UseAnt[refantID]] + ' as the refant.'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
antDia = np.ones(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
#-------- Load Aeff file
Afile = open(SCR_DIR + 'AeB' + `BandID` + '.data')
Alines = Afile.readlines()
Afile.close()
AeX, AeY = 0.25* np.pi* antDia**2, 0.25* np.pi* antDia**2       # antenna collecting area (100% efficiency)
AeX, AeY = [], []
for ant_index in antMap:
    for Aline in Alines:
        if antList[ant_index] in Aline:
            AeX = AeX + [(0.0025* np.pi* float(Aline.split()[1]))* antDia[ant_index]**2]
            AeY = AeY + [(0.0025* np.pi* float(Aline.split()[2]))* antDia[ant_index]**2]
        #
    #
#
AeX, AeY = np.array(AeX), np.array(AeY) # in antMap order 
#-------- Check D-term files
if not 'DPATH' in locals(): DPATH = SCR_DIR
print '---Checking D-term files in ' + DPATH
DantList, noDlist = [], []
for ant_index in UseAnt:
    Dfile = DPATH + 'B' + `BandID` + '-SPW' + `spwList[0]` + '-' + antList[ant_index] + '.DSpec.npy'
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
for ant_index in antMap:
    for spw in spwList:
        Dfile = DPATH + 'B' + `BandID` + '-SPW' + `spw` + '-' + antList[ant_index] + '.DSpec.npy'
        Dterm = np.load(Dfile)
        DxList = DxList + [splineComplex(Dterm[0], Dterm[1] + (0.0 + 1.0j)* Dterm[2])]
        DyList = DyList + [splineComplex(Dterm[0], Dterm[3] + (0.0 + 1.0j)* Dterm[4])]
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
BPList = []
print '---Loading bandpass table'
for spw in spwList:
    BP_ant = np.load(BPprefix + '-SPW' + `spw` + '-BPant.npy')
    XY_BP  = splineComplex(Dterm[0], np.load(BPprefix + '-SPW' + `spw` + '-XYspec.npy'))
    BP_ant[:,1] *= XY_BP
    if 'XYsign' in locals(): BP_ant[:,1] *= XYsign
    BPList = BPList + [BP_ant]
#
##-------- Equalization using EQ scan
relGain = np.ones([2, UseAntNum])
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA; PA = np.arctan2( np.sin(PA), np.cos(PA))
useAntMapRev = indexList(np.array(antMap), np.array(UseAnt))
exp_Tau = np.exp(-Tau0spec / np.sin(np.mean(ElScan)))
TsysEQScan = (np.mean(Trxspec[:,:,chRange],axis=2).reshape([UseAntNum, spwNum, 2]).transpose(2,0,1) + Tcmb* np.mean(exp_Tau, axis=1) + tempAtm* (1.0 - np.mean(exp_Tau[:,chRange], axis=1))).transpose(1,2,0)[useAntMapRev]
#
#for spw_index in range(spwNum):
spw_index = 0
#-------- Baseline-based cross power spectra
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], EQScan)
AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + (BANDPA[BandID] + 90.0)*pi/180.0
QCpUS = QUmodel[0]* np.cos(2.0* PA) + QUmodel[1]* np.sin(2.0* PA)         # Q cos + U sin
    # TauObs = Tau0 / np.sin(np.mean(ElScan)); atmCorrect = np.exp(-TauObs)
    # chAvgTsys = chAvgTrx + 285.0* (1.0 - np.exp(-TauObs))
chNum = Xspec.shape[1]
tempSpec = CrossPolBL(Xspec[:,chRange][:,:,blMap], blInv).transpose(3,2,0,1)       # Cross Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex][:,:,chRange]* BPList[spw_index][ant1][:,polXindex][:,:,chRange].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
#-------- Antenna-based Gain correction
chAvgVis = np.mean(BPCaledXspec, axis=1)[pPol]
GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
#aprioriSEFD = 2.0* kb* chAvgTsys[antMap].T / np.array([AeX, AeY])
aprioriSEFD = 2.0* kb* TsysEQScan[:,spw_index].T / np.array([AeX, AeY])
aprioriVisX = np.mean(chAvgVis[0] / ((1.0 + QCpUS)* GainP[0,ant0]* GainP[0,ant1].conjugate()), axis=1) * np.sqrt(aprioriSEFD[0, ant0]* aprioriSEFD[0, ant1])
aprioriVisY = np.mean(chAvgVis[1] / ((1.0 - QCpUS)* GainP[1,ant0]* GainP[1,ant1].conjugate()), axis=1) * np.sqrt(aprioriSEFD[1, ant0]* aprioriSEFD[1, ant1])
#-------- Determine Antenna-based Gain
relGain[0] = abs(gainComplex(aprioriVisX))
relGain[1] = abs(gainComplex(aprioriVisY))
medGain = np.median( abs(relGain) )
relGain[0] /= medGain # X-pol delta gain
relGain[1] /= medGain # Y-pol delta gain
#
#-------- Flux Density
spw_index = 0
chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq *= 1.0e-9
print '---Stokes Spectrum densities of sources ---'
pp, polLabel, Pcolor = PdfPages('SP_' + prefix + '_' + BandName + '.pdf'), ['I', 'Q', 'U', 'V'], ['black', 'blue', 'red', 'green']
page_index = 0
for scan in scanList:
    figScan = plt.figure(page_index, figsize = (8,11))
    figScan.suptitle(prefix + ' ' + BandName + ' Scan' + `scan`)
    figScan.text(0.45, 0.05, 'Frequency [GHz]') 
    figScan.text(0.03, 0.45, 'Stokes Parameters [Jy]', rotation=90) 
    BPCaledXspec = []
    #-------- UV distance
    timeStamp, UVW = GetUVW(msfile, spwList[0], scan);  uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    figScan.text(0.75, 0.95, qa.time('%fs' % timeStamp[0], form='ymd')[0]) 
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    PA = AzEl2PA(AzScan, ElScan) + (BANDPA[BandID] + 90.0)*pi/180.0; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
    text_sd = ' AZ=%4.1f EL=%4.1f X-feed-PA=%4.1f' % (np.median(AzScan)*180.0/pi, np.median(ElScan)*180.0/pi, np.median(PA)*180.0/pi)
    print text_sd
    #-------- Tsys calibration
    atmCorrect = np.exp(Tau0spec[spw_index] / np.sin(np.mean(ElScan)))
    exp_Tau = 1.0 / atmCorrect
    TsysSPW =  (Trxspec[spw_index::spwNum].transpose(1,0,2) + Tcmb* exp_Tau + tempAtm * (1.0 - exp_Tau))[:,useAntMapRev]
    #---- Flagged by Tsys
    tsysFlagAntIndex = unique(np.where(TsysSPW <0.0)[1]).tolist()
    if len(tsysFlagAntIndex) > 0:
        for ant_index in tsysFlagAntIndex: TsysSPW[:,ant_index] = Trxspec[spw_index::spwNum][ant_index] + tempAtm* (1.0 - np.exp(-Tau0spec[spw_index] / np.sin(OnEL[scan_index])))
    #
    SEFD = 2.0* kb* (TsysSPW * atmCorrect).transpose(2,0,1) /  (np.array([AeX[0:UseAntNum], AeY[0:UseAntNum]]) / relGain[:,0:UseAntNum]**2)
    #TauObs = Tau0 / np.sin(np.mean(ElScan)); atmCorrect = np.exp(-TauObs)
    #chAvgTsys = chAvgTrx + 285.0* (1.0 - np.exp(-TauObs))
    figScan.text(0.05, 0.95, text_sd) 
    #-------- Baseline-based cross power spectra
    if 'FIELDID' in locals():
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], scan, FIELDID)
    else :
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], scan)
    #
    timeNum, chNum = Xspec.shape[3], Xspec.shape[1]
    UseChNum = len(chRange)
    if np.max(abs(Xspec)) < 1.0e-9: continue
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
    #-------- Bandpass Calibration
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    #-------- Antenna-based Phase Solution
    chAvgVis = np.mean(BPCaledXspec[:,PCchRange], axis=1) # chAvgVis[pol, bl, time]
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    pCalVis = (BPCaledXspec.transpose(1,0,2,3) / (GainP[polYindex][:,ant0]* GainP[polXindex][:,ant1].conjugate()))
    # SEFD = 2.0* kb* chAvgTsys[antMap].T / (np.array([AeX, AeY]) * atmCorrect)/ (relGain[spw_index]**2)
    # AmpCalVis = (pCalVis.transpose(0,3,1,2)* np.sqrt(SEFD[polYindex][:,ant0]* SEFD[polXindex][:,ant1])).transpose(3,2,0,1) # AmpCalVis[bl,pol,ch,time]
    AmpCalVis = (pCalVis[chRange].transpose(3,0,1,2)* np.sqrt(SEFD[chRange][:,polYindex][:,:,ant0[0:UseBlNum]]* SEFD[chRange][:,polXindex][:,:,ant1[0:UseBlNum]])).transpose(3,2,1,0)
    StokesI_PL = figScan.add_subplot( 2, 1, 1)
    StokesP_PL = figScan.add_subplot( 2, 1, 2)
    Stokes = np.zeros([4,blNum, UseChNum], dtype=complex)  # Stokes[stokes, bl, ch]
    for bl_index in range(UseBlNum):
        Minv = InvMullerVector(DxSpec[ant1[bl_index],spw_index], DySpec[ant1[bl_index],spw_index], DxSpec[ant0[bl_index],spw_index], DySpec[ant0[bl_index],spw_index], np.ones(chNum))
        for time_index in range(timeNum):
            for ch_index in range(UseChNum):
                Stokes[:,bl_index,ch_index] += PS[:,:,time_index].dot(Minv[:,:,ch_index].dot(AmpCalVis[bl_index, :, ch_index, time_index]))
            #
        #
    #
    Stokes = np.mean(Stokes, axis=1) / timeNum
    StokesSpec, StokesErr = Stokes.real, Stokes.imag
    np.save(prefix + '-' + BandName + '-Scan' + `scan` + '.StokesSpec.npy', StokesSpec)
    np.save(prefix + '-' + BandName + '-Scan' + `scan` + '.StokesErr.npy', StokesErr)
    StokesI_PL.plot( Freq[chRange], StokesSpec[0], ls='steps-mid', label=polLabel[0], color=Pcolor[0])
    StokesP_PL.plot( np.array([min(Freq[chRange]), max(Freq[chRange])]), np.array([0.0,0.0]), '-', color='gray')
    StokesP_PL.plot( Freq[chRange], StokesSpec[1], ls='steps-mid', label=polLabel[1], color=Pcolor[1])
    StokesP_PL.plot( Freq[chRange], StokesSpec[2], ls='steps-mid', label=polLabel[2], color=Pcolor[2])
    StokesP_PL.plot( Freq[chRange], StokesSpec[3], ls='steps-mid', label=polLabel[3], color=Pcolor[3])
    #StokesI_PL.plot( np.array(chRange), StokesSpec[0, chRange], ls='steps-mid', label=polLabel[0], color=Pcolor[0])
    #StokesP_PL.plot( np.array([min(chRange), max(chRange)]), np.array([0.0,0.0]), '-', color='gray')
    #StokesP_PL.plot( np.array(chRange), StokesSpec[1, chRange], ls='steps-mid', label=polLabel[1], color=Pcolor[1])
    #StokesP_PL.plot( np.array(chRange), StokesSpec[2, chRange], ls='steps-mid', label=polLabel[2], color=Pcolor[2])
    #StokesP_PL.plot( np.array(chRange), StokesSpec[3, chRange], ls='steps-mid', label=polLabel[3], color=Pcolor[3])
    plotMax = max(StokesSpec[0])
    StokesI_PL.axis([min(Freq[chRange]), max(Freq[chRange]), 0.0, 1.2* plotMax])
    StokesP_PL.axis([min(Freq[chRange]), max(Freq[chRange]), -0.15*plotMax, 0.15*plotMax ])
    #StokesI_PL.axis([min(chRange), max(chRange), 0.0, 1.2* plotMax])
    #StokesP_PL.axis([min(chRange), max(chRange), -0.10*plotMax, 0.10*plotMax ])
    #
    StokesI_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    StokesP_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    StokesP_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    StokesP_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    figScan.savefig(pp, format='pdf')
    page_index += 1
#
plt.close('all')
pp.close()
msmd.close()
