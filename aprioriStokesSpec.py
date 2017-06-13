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
spwNum = len(spw)
spwName = msmd.namesforspws(spw[0])[0]
BandName = re.findall(r'RB_..', spwName)[0]; BandID = int(BandName[3:5])
#-------- Configure Array
print '---Checking array configulation'
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
Tatm_OFS  = 15.0     # Ambient-load temperature - Atmosphere temperature
kb        = 1.38064852e3
#-------- Load Tsys File
Trec = np.load(TSprefix + '.Trx.npy') 
Tau0 = np.load(TSprefix + '.Tau0.npy') 
TRchNum = Trec.shape[2]; TRchRange = range(int(0.05*TRchNum), int(0.95*TRchNum))
chAvgTrx = np.median(np.mean(Trec[:,:,TRchRange], axis=2), axis=2)
Tau0 = np.median(Tau0)
#-------- Array Configuration
print '---Checking array configuration'
flagList = np.where(np.median(chAvgTrx.reshape(antNum, 2* spwNum), axis=1) > 2.0* np.median(chAvgTrx))[0].tolist()
flagList = unique(flagList + np.where(np.min(chAvgTrx.reshape(antNum, 2* spwNum), axis=1) < 1.0 )[0].tolist()).tolist()
flagAnt[flagList] = 0.0 # Flagging by abnormal Trx
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
text_sd = '  Usable antennas: '
for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
print text_sd
text_sd = '  Flagged by Trx:  '
for ants in antList[flagList].tolist(): text_sd = text_sd + ants + ' '
print text_sd
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
timeStamp, UVW = GetUVW(msfile, spw[0], msmd.scansforspw(spw[0])[0])
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
    Dfile = DPATH + 'B' + `BandID` + '-SPW0-' + antList[ant_index] + '.DSpec.npy'
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
    for spw_index in spw:
        Dfile = DPATH + 'B' + `BandID` + '-SPW' + `spw_index` + '-' + antList[ant_index] + '.DSpec.npy'
        Dterm = np.load(Dfile)
        DxList = DxList + [splineComplex(Dterm[0], Dterm[1] + (0.0 + 1.0j)* Dterm[2])]
        DyList = DyList + [splineComplex(Dterm[0], Dterm[3] + (0.0 + 1.0j)* Dterm[4])]
    #
#
chNum = np.array(DxList).shape[1]
DxSpec, DySpec = np.array(DxList).reshape([UseAntNum, spwNum, chNum]), np.array(DyList).reshape([UseAntNum, spwNum, chNum])
#-------- Bandpass Table
BPList = []
print '---Loading bandpass table'
for spw_index in spw:
    BP_ant = np.load(BPprefix + '-SPW' + `spw_index` + '-BPant.npy')
    XY_BP  = splineComplex(Dterm[0], np.load(BPprefix + '-SPW' + `spw_index` + '-XYspec.npy'))
    BP_ant[:,1] *= XY_BP
    if 'XYsign' in locals(): BP_ant[:,1] *= XYsign
    BPList = BPList + [BP_ant]
#
##-------- Equalization using EQ scan
relGain = np.ones([spwNum, 2, UseAntNum])
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
for spw_index in range(spwNum):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], EQScan)
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    TauObs = Tau0 / np.sin(np.mean(ElScan)); atmCorrect = np.exp(-TauObs)
    chAvgTsys = chAvgTrx + 285.0* (1.0 - np.exp(-TauObs))
    chNum = Xspec.shape[1]
    tempSpec = CrossPolBL(Xspec[:,BPchRange][:,:,blMap], blInv).transpose(3,2,0,1)       # Cross Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex][:,:,BPchRange]* BPList[spw_index][ant1][:,polXindex][:,:,BPchRange].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec, axis=1)[pPol]
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])
    aprioriSEFD = 2.0* kb* chAvgTsys[antMap].T / np.array([AeX, AeY])
    aprioriVisX = np.mean(chAvgVis[0] / (GainP[0,ant0]* GainP[0,ant1].conjugate()), axis=1) * np.sqrt(aprioriSEFD[0, ant0]* aprioriSEFD[0, ant1])
    aprioriVisY = np.mean(chAvgVis[1] / (GainP[1,ant0]* GainP[1,ant1].conjugate()), axis=1) * np.sqrt(aprioriSEFD[1, ant0]* aprioriSEFD[1, ant1])
    #-------- Determine Antenna-based Gain
    relGain[spw_index, 0] = abs(gainComplex(aprioriVisX))
    relGain[spw_index, 1] = abs(gainComplex(aprioriVisY))
    medGain = np.median( abs(relGain) )
    relGain[spw_index, 0] /= medGain # X-pol delta gain
    relGain[spw_index, 1] /= medGain # Y-pol delta gain
#
#-------- Flux Density
spw_index = 0
chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq *= 1.0e-9
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
    timeStamp, UVW = GetUVW(msfile, spw[0], scan);  uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    figScan.text(0.75, 0.95, qa.time('%fs' % timeStamp[0], form='ymd')[0]) 
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    PA = AzEl2PA(AzScan, ElScan) + (BANDPA[BandID] + 90.0)*pi/180.0; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
    text_sd = ' AZ=%4.1f EL=%4.1f X-feed-PA=%4.1f' % (np.median(AzScan)*180.0/pi, np.median(ElScan)*180.0/pi, np.median(PA)*180.0/pi)
    print text_sd
    TauObs = Tau0 / np.sin(np.mean(ElScan)); atmCorrect = np.exp(-TauObs)
    chAvgTsys = chAvgTrx + 285.0* (1.0 - np.exp(-TauObs))
    figScan.text(0.05, 0.95, text_sd) 
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], scan)
    timeNum, chNum = Xspec.shape[3], Xspec.shape[1]
    UseChNum = len(chRange)
    if np.max(abs(Xspec)) < 1.0e-9: continue
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
    #-------- Bandpass Calibration
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    #-------- Antenna-based Phase Solution
    chAvgVis = np.mean(BPCaledXspec[:,chRange], axis=1) # chAvgVis[pol, bl, time]
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    pCalVis = (BPCaledXspec.transpose(1,0,2,3) / (GainP[polYindex][:,ant0]* GainP[polXindex][:,ant1].conjugate()))
    SEFD = 2.0* kb* chAvgTsys[antMap].T / (np.array([AeX, AeY]) * atmCorrect)/ (relGain[spw_index]**2)
    AmpCalVis = (pCalVis.transpose(0,3,1,2)* np.sqrt(SEFD[polYindex][:,ant0]* SEFD[polXindex][:,ant1])).transpose(3,2,0,1) # AmpCalVis[bl,pol,ch,time]
    StokesI_PL = figScan.add_subplot( 2, spwNum, 0* spwNum + spw_index + 1)
    StokesP_PL = figScan.add_subplot( 2, spwNum, 1* spwNum + spw_index + 1)
    Stokes = np.zeros([4,blNum, chNum], dtype=complex)  # Stokes[stokes, bl, ch]
    for bl_index in range(blNum):
        Minv = InvMullerVector(DxSpec[ant1[bl_index],spw_index], DySpec[ant1[bl_index],spw_index], DxSpec[ant0[bl_index],spw_index], DySpec[ant0[bl_index],spw_index], np.ones(chNum))
        for time_index in range(timeNum):
            for ch_index in range(chNum):
                Stokes[:,bl_index,ch_index] += PS[:,:,time_index].dot(Minv[:,:,ch_index].dot(AmpCalVis[bl_index, :, ch_index, time_index]))
            #
        #
    #
    Stokes = np.mean(Stokes, axis=1) / timeNum
    StokesSpec, StokesErr = Stokes.real, Stokes.imag
    StokesI_PL.plot( Freq[chRange], StokesSpec[0, chRange], ls='steps-mid', label=polLabel[0], color=Pcolor[0])
    StokesP_PL.plot( np.array([min(Freq[chRange]), max(Freq[chRange])]), np.array([0.0,0.0]), '-', color='gray')
    StokesP_PL.plot( Freq[chRange], StokesSpec[1, chRange], ls='steps-mid', label=polLabel[1], color=Pcolor[1])
    StokesP_PL.plot( Freq[chRange], StokesSpec[2, chRange], ls='steps-mid', label=polLabel[2], color=Pcolor[2])
    StokesP_PL.plot( Freq[chRange], StokesSpec[3, chRange], ls='steps-mid', label=polLabel[3], color=Pcolor[3])
    plotMax = max(StokesSpec[0, chRange])
    StokesI_PL.axis([min(Freq[chRange]), max(Freq[chRange]), 0.0, 1.2* plotMax])
    StokesP_PL.axis([min(Freq[chRange]), max(Freq[chRange]), -0.12*plotMax, 0.12*plotMax ])
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
"""
np.save(prefix + '-' + UniqBands[band_index] + '.Flux.npy', ScanFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Ferr.npy', ErrFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Source.npy', np.array(sourceList)[sourceIDscan])
np.save(prefix + '-' + UniqBands[band_index] + '.EL.npy', ScanEL)
"""
msmd.close()
