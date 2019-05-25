EQflux = np.ones([2*spwNum])
#-------- Flux density of the equalizer
spwStokesDic = dict(zip(sourceList, [[]]*len(sourceList))) # spwStokesDic[sourceName][pol* spw]
for spw_index in range(spwNum):
    FLX, FLY = [], []
    for sso_index in SSOUseList:
        index = np.where(AeX[:, spw_index, sso_index] > 1.0)[0].tolist()
        FLX = FLX + (AeSeqX[spw_index, index] / AeX[index, spw_index, sso_index]).tolist()
        FLY = FLY + (AeSeqY[spw_index, index] / AeY[index, spw_index, sso_index]).tolist()
    #
    EQflux[spw_index], EQflux[spw_index + spwNum] = np.median(np.array(FLX)), np.median(np.array(FLY))
#
#-------- Power-law fit
P = np.c_[ np.r_[np.log(centerFreqList),np.log(centerFreqList)], np.r_[np.ones(spwNum), np.zeros(spwNum)], np.r_[np.zeros(spwNum), np.ones(spwNum)]]
EQmodel = scipy.linalg.solve(np.dot(P.T, P), np.dot(P.T, np.log(EQflux)))   # alpha, logSx, logSy
EQflux = np.c_[np.exp(EQmodel[0]* np.log(centerFreqList) + EQmodel[1]), np.exp(EQmodel[0]* np.log(centerFreqList) + EQmodel[2])]
##-------- Aperture efficiencies
Ae     = np.c_[AeSeqX.T / EQflux[:,0], AeSeqY.T / EQflux[:,1]].reshape(UseAntNum, ppolNum, spwNum)
AEFF   = (Ae.transpose(1,2,0) / (0.25* pi*antDia[antMap]**2)).transpose(2,0,1)
np.save(prefix + '-' + UniqBands[band_index] + '.AntList.npy', antList[antMap]) 
np.save(prefix + '-' + UniqBands[band_index] + '.Aeff.npy', AEFF) 
text_sd =  ' Aeff:'; logfile.write(text_sd); print text_sd,
for spw_index in range(spwNum):
    for pol_index in range(2):
        text_sd = 'SPW%02d-%s ' % (spwList[spw_index], PolList[pol_index])
        logfile.write(text_sd); print text_sd,
    #
#
logfile.write('\n'); print ''
logjy = open(prefix + '-' + UniqBands[band_index] + '-JyK.log', 'w')
for ant_index in range(UseAntNum):
    text_sd = '%s :' % (antList[antMap[ant_index]]); logfile.write(text_sd); print text_sd,
    for spw_index in range(spwNum):
        for pol_index in range(2):
            #logjy.write( '%s %s %d %' % (prefix, antList[antMap[ant_index]], spwList[spw_index])
            text_sd = '  %4.1f%% ' % (100.0* AEFF[ant_index, pol_index, spw_index]); logfile.write(text_sd); print text_sd,
            text_jy = '%s %s %d %s %f %s %s %fe+9 %e %5.1f' % (prefix, antList[antMap[ant_index]], spwList[spw_index], PolList[pol_index], 2.0* kb / (0.25* AEFF[ant_index, pol_index, spw_index]* pi*antDia[antMap[ant_index]]**2), timeLabel, UniqBands[band_index], centerFreqList[spw_index], len(chRange)* chWid[0], tempAtm); logjy.write(text_jy + '\n')
        #
    #
    logfile.write('\n'); print ''; logjy.write('\n')
##
logfile.write('\n'); print ''
logjy.write('\n'); logjy.close()
#-------- XY phase using BP scan
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], BPScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PA = np.arctan2( np.sin(PA), np.cos(PA))
XYphase, caledVis = [], []
scan_index = onsourceScans.index(BPScan)
Trx2antMap = indexList( antList[antMap], antList[TrxMap] )
for spw_index in range(spwNum):
    exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
    atmCorrect = 1.0 / exp_Tau
    TsysSPW = (TrxList[spw_index].transpose(2,0,1) + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))[Trx2antMap] # [antMap, pol, ch]
    TsysBL  = np.sqrt( TsysSPW[ant0][:,polYindex]* TsysSPW[ant1][:,polXindex])
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPScan); timeNum = len(timeStamp)
    if 'FG' in locals(): flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
    else : flagIndex = range(timeNum)
    chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
    BPCaledXspec = (tempSpec * TsysBL/ (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    SEFD = 2.0* kb / Ae[:,:,spw_index].T
    caledVis.append(np.mean((chAvgVis / (GainP[polYindex][:,ant0]* GainP[polXindex][:,ant1].conjugate())).transpose(2, 0, 1)* np.sqrt(SEFD[polYindex][:,ant0]* SEFD[polXindex][:,ant1]), axis=2).T)
#
caledVis = np.array(caledVis)   # [spw, pol, time]
if 'QUMODEL' in locals():
    StokesBP = np.array(StokesDic[BPcal])
    if QUMODEL: QUsolution = np.array(StokesBP[[1,2]])/StokesBP[0]
else: QUsolution = XXYY2QU(PA, np.mean(caledVis[:,[0,3]], axis=0))
#-------- XY phase cal in Bandpass table
for spw_index in range(spwNum):
    XYphase = XY2Phase(QUsolution[1]* np.cos(2.0* PA) - QUsolution[0]* np.sin(2.0* PA), caledVis[spw_index][[1,2]])
    XYsign = np.sign(np.cos(XYphase))
    BPList[spw_index][:,1] *= XYsign
    print 'SPW[%d] : XY phase = %6.1f [deg] sign = %3.0f' % (spwList[spw_index], 180.0*XYphase/pi, XYsign)
#
#-------- Flux Density
ScanFlux = np.zeros([scanNum, spwNum, 4])
ScanSlope= np.zeros([scanNum, spwNum, 4])
ErrFlux  = np.zeros([scanNum, spwNum, 4])
print '---Flux densities of sources ---'
pp = PdfPages('FL_' + prefix + '_' + UniqBands[band_index] + '.pdf')
polLabel = ['I', 'Q', 'U', 'V']
Pcolor   = ['black', 'blue', 'red', 'green']
XYD, XYC, scanTime = [], [],[]      # XY delay and correlation
scanDic   = dict(zip(sourceList, [[]]*len(sourceList))) # Scan list index for each source
timeDic   = dict(zip(sourceList, [[]]*len(sourceList))) # Time index list for each source
PADic     = dict(zip(sourceList, [[]]*len(sourceList))) # PA list for each source
figFL = plt.figure(figsize = (11, 8))
figFL.suptitle(prefix + ' ' + UniqBands[band_index])
figFL.text(0.45, 0.05, 'Projected baseline [m]')
figFL.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90)
VisSpec = np.zeros([spwNum, 4, chNum, UseBlNum, timeSum], dtype=complex)
timePointer = 0
for scan_index in range(scanNum):
    sourceName = sourceList[sourceIDscan[scan_index]]
    scanDic[sourceName] = scanDic[sourceName] + [scan_index]
    if scan_index > 0:
        for PL in IList: figFL.delaxes(PL)
        for PL in PList: figFL.delaxes(PL)
    #
    #-------- UV distance
    timeStamp, UVW = GetUVW(msfile, spwList[0], onsourceScans[scan_index])
    scanTime = scanTime + [np.median(timeStamp)]
    timeNum = len(timeStamp)
    timeDic[sourceName] = range(timePointer, timePointer + timeNum)
    uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    if min(ElScan) < 20.0 / 180.0* pi: continue
    PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PAnum = len(PA); PS = InvPAVector(PA, np.ones(PAnum))
    PADic[sourceName] = PA.tolist()
    #-------- Prepare plots
    IList, PList = [], []      # XY delay and correlation
    text_time = qa.time('%fs' % np.median(timeStamp), form='ymd')[0]
    text_src  = ' %02d %010s EL=%4.1f deg' % (scanList[scan_index], sourceName, 180.0* OnEL[scan_index]/pi); logfile.write(text_src + ' ' + text_time + '\n'); print text_src + ' ' + text_time
    if(onsourceScans[scan_index] in SSOscanID):
        SSO_flag = True
        SSO_ID = SSOscanID.index(onsourceScans[scan_index])
        text_sd = ' SPW  Frequency    I               Q               U               V             | Model I'; logfile.write(text_sd + '\n'); print text_sd
    else:
        SSO_flag = False
        text_sd = ' SPW  Frequency    I                 Q                 U                 V                 %Pol     EVPA '; logfile.write(text_sd + '\n'); print text_sd
    text_sd = ' --------------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
    BPCaledXspec = []
    #-------- Sub-array formation
    SAantMap, SAblMap, SAblInv, SAant0, SAant1 = antMap, blMap, blInv, ant0, ant1
    if SSO_flag:
        SAantMap, SAblMap, SAblInv = subArrayIndex(uvFlag[SSO_ID], refantID) # antList[SAantMap] lists usable antennas
        SAblIndex = indexList(np.array(SAblMap), np.array(blMap))
        SAant0, SAant1 = np.array(ant0)[SAblIndex].tolist(), np.array(ant1)[SAblIndex].tolist()
    #
    bpAntMap = indexList(antList[SAantMap],antList[antMap])
    Trx2antMap = indexList( antList[SAantMap], antList[TrxMap] )
    SAantNum = len(SAantMap); SAblNum = len(SAblMap)
    if SAblNum < 6:
        text_sd = ' Only %d baselines for short enough sub-array. Skip!' % (SAblNum) ; logfile.write(text_sd + '\n'); print text_sd
        continue
    #
    #-------- Baseline-based cross power spectra
    for spw_index in range(spwNum):
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], onsourceScans[scan_index])
        chNum = Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum)); UseChNum = len(chRange)
        tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        #-------- Bandpass Calibration
        BPCaledXspec = BPCaledXspec + [(tempSpec / (BPList[spw_index][SAant0][:,polYindex]* BPList[spw_index][SAant1][:,polXindex].conjugate())).transpose(2,3,1,0)] # Bandpass Cal
    #
    #-------- Antenna-based Gain
    BPCaledXspec = np.array(BPCaledXspec)   # BPCaledXspec[spw, pol, ch, bl, time]
    chAvgVis = np.mean(BPCaledXspec[:, :, chRange], axis=(0,2))
    if(SSO_flag): chAvgVis =(np.mean(BPCaledXspec[:,:, chRange], axis=(0,2)).transpose(0,2,1) / SSOmodelVis[SSO_ID, spw_index, SAblMap]).transpose(0,2,1)
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    pCalVis = (BPCaledXspec.transpose(0,2,1,3,4) / (GainP[polYindex][:,ant0[0:SAblNum]]* GainP[polXindex][:,ant1[0:SAblNum]].conjugate()))[:,chRange]
    #-------- XY phase spectra
    for spw_index in range(spwNum):
        delayFact = (chNum + 0.0)/len(chRange)
        XYspec = np.mean(pCalVis[spw_index, :, 1:3, :], axis=(2,3)).T
        XYdelay, XYamp = delay_search(XYspec[:,0]); YXdelay, YXamp = delay_search(XYspec[:,1])
        XYD = XYD + [XYdelay* delayFact, YXdelay* delayFact]
        xyc, yxc = np.mean(delay_cal(XYspec[:,0], XYdelay)), np.mean(delay_cal(XYspec[:,1], YXdelay))
        XYC = XYC + [xyc, yxc]
        #text_sd = ' spw%d : Delay (samples) XY= %.3f YX=%.3f  phase (deg) XY= %.2f YX=%.2f' % (spwList[spw_index], XYdelay* delayFact, YXdelay* delayFact, 90.0* np.angle(np.exp((0.0 + 2.0j)* np.angle(xyc))) / pi, 90.0* np.angle(np.exp((0.0 + 2.0j)* np.angle(yxc))) / pi)
        #print text_sd
    #
    for spw_index in range(spwNum):
        StokesI_PL = figFL.add_subplot( 2, spwNum, spw_index + 1 )
        StokesP_PL = figFL.add_subplot( 2, spwNum, spwNum + spw_index + 1 )
        IList = IList + [StokesI_PL]
        PList = PList + [StokesP_PL]
        exp_Tau = np.exp(-(Tau0spec[spw_index] + exTauSP(np.median(timeStamp))) / np.mean(np.sin(ElScan)))
        atmCorrect = 1.0 / exp_Tau
        TsysSPW = (TrxList[spw_index].transpose(2,0,1) + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau))[Trx2antMap]    # [ant, pol, ch]
        if SSO_flag:
            Ta = SSOflux[sso_index, spw_index]* Ae[bpAntMap, :, spw_index]* np.mean(atmCorrect) / (2.0* kb)
            TsysSPW = (TsysSPW.transpose(2,0,1) + Ta).transpose(1,2,0)
        #
        SEFD = 2.0* kb* (TsysSPW * atmCorrect).transpose(2,1,0) / Ae[bpAntMap][:,:,spw_index].T   # SEFD[ch,pol,antMap]
        SAantNum = len(SAantMap); SAblNum = len(SAblMap)
        #-------- Additional equalizaiton
        if not SSO_flag:        # Additional equalization for point sources
            AmpCalVis = np.mean(np.mean(pCalVis[spw_index], axis=3)[:,[0,3]] * np.sqrt(SEFD[chRange][:,:,ant0[0:SAblNum]]* SEFD[chRange][:,:,ant1[0:SAblNum]]), axis=0)
            indivRelGain = abs(gainComplexVec(AmpCalVis.T)); indivRelGain /= np.percentile(indivRelGain, 75, axis=0)
            SEFD /= (indivRelGain**2).T
        #
        AmpCalVis = (pCalVis[spw_index].transpose(3,0,1,2)* np.sqrt(SEFD[chRange][:,polYindex][:,:,ant0[0:SAblNum]]* SEFD[chRange][:,polXindex][:,:,ant1[0:SAblNum]])).transpose(3,2,1,0)
        VisSpec[spw_index][:,chRange,:,timePointer:timePointer+timeNum] = AmpCalVis.transpose(1,2,0,3)
        text_sd = ' SPW%02d %5.1f GHz ' % (spwList[spw_index], centerFreqList[spw_index]); logfile.write(text_sd); print text_sd,
        Stokes = np.zeros([4,SAblNum], dtype=complex)
        for bl_index in range(SAblNum):
            Minv = InvMullerVector(DxSpec[SAant1[bl_index], spw_index][chRange], DySpec[SAant1[bl_index], spw_index][chRange], DxSpec[SAant0[bl_index], spw_index][chRange], DySpec[SAant0[bl_index], spw_index][chRange], np.ones(UseChNum))
            Stokes[:,bl_index] = PS.reshape(4, 4*PAnum).dot(Minv.reshape(4, 4*UseChNum).dot(AmpCalVis[bl_index].reshape(4*UseChNum, PAnum)).reshape(4*PAnum)) / (PAnum* UseChNum)
        #
        StokesVis, StokesErr = Stokes.real, Stokes.imag
        if SSO_flag: StokesVis /= SSOmodelVis[SSO_ID, spw_index][SAblMap]
        percent75 = np.percentile(StokesVis[0], 75); sdvis = np.std(StokesVis[0])
        visFlag = np.where(abs(StokesVis[0] - percent75) < 2.0* sdvis )[0]
        if len(visFlag) < 2 : continue
        weight = np.zeros(SAblNum); weight[visFlag] = 1.0/np.var(StokesVis[0][visFlag])
        if SSO_flag: weight *= (SSOmodelVis[SSO_ID, spw_index][SAblMap]**2)
        P, W = np.c_[np.ones(len(weight)), uvDist[SAblMap]], np.diag(weight)
        PtWP_inv = scipy.linalg.inv(P.T.dot(W.dot(P)))
        solution, solerr = PtWP_inv.dot(P.T.dot(weight* StokesVis[0])),  np.sqrt(np.diag(PtWP_inv)) # solution[0]:intercept, solution[1]:slope
        if abs(solution[1]) < 5.0* solerr[1]: solution[0], solution[1] = np.median(StokesVis[0][visFlag]), 0.0
        ScanFlux[scan_index, spw_index, 0], ScanSlope[scan_index, spw_index, 0], ErrFlux[scan_index, spw_index, 0] = solution[0], solution[1], solerr[0]
        for pol_index in range(1,4):
            ScanSlope[scan_index, spw_index, pol_index] = ScanSlope[scan_index, spw_index, 0] * np.median(StokesVis[pol_index])/ScanFlux[scan_index, spw_index, 0]
            solution[0] = (weight.dot(StokesVis[pol_index]) - ScanSlope[scan_index, spw_index, pol_index]* weight.dot(uvDist[SAblMap]))/(np.sum(weight))
            ScanFlux[scan_index, spw_index, pol_index] = solution[0]
            resid = StokesVis[pol_index] - ScanSlope[scan_index, spw_index, pol_index]* uvDist[SAblMap] - solution[0]; ErrFlux[scan_index, spw_index, pol_index] = np.sqrt(weight.dot(resid**2)/np.sum(weight))
        #
        for pol_index in range(4):
            if len(visFlag) < 6:
                text_sd = 'Only %d vis.    ' % (len(visFlag)) ; logfile.write(text_sd + '\n'); print text_sd,
                continue
            #
            text_sd = '%7.4f (%.4f) ' % (ScanFlux[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index]); logfile.write(text_sd); print text_sd,
        #
        StokesI_PL.plot( uvDist[SAblMap], StokesVis[0], '.', label=polLabel[0], color=Pcolor[0])
        StokesP_PL.plot( uvDist[SAblMap], StokesVis[1], '.', label=polLabel[1], color=Pcolor[1])
        StokesP_PL.plot( uvDist[SAblMap], StokesVis[2], '.', label=polLabel[2], color=Pcolor[2])
        StokesP_PL.plot( uvDist[SAblMap], StokesVis[3], '.', label=polLabel[3], color=Pcolor[3])
        if(SSO_flag):
            text_sd = '| %6.3f ' % (SSOflux0[SSO_ID, spw_index]); logfile.write(text_sd); print text_sd,
            logfile.write('\n'); print ''
        else: 
            text_sd = '%6.3f   %6.1f ' % (100.0* np.sqrt(ScanFlux[scan_index, spw_index, 1]**2 + ScanFlux[scan_index, spw_index, 2]**2)/ScanFlux[scan_index, spw_index, 0], np.arctan2(ScanFlux[scan_index, spw_index, 2],ScanFlux[scan_index, spw_index, 1])*90.0/pi); logfile.write(text_sd); print text_sd,
            logfile.write('\n'); print ''
        #
        #
    #
    uvMin, uvMax, IMax = min(uvDist[blMap]), max(uvDist[blMap]), max(ScanFlux[scan_index,:,0])
    for spw_index in range(spwNum):
        StokesI_PL, StokesP_PL = IList[spw_index], PList[spw_index]
        if spw_index == 0: StokesI_PL.text(0.0, IMax*1.35, text_src)
        if spw_index == spwNum - 1: StokesI_PL.text(0.0, IMax*1.35, text_time)
        StokesI_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 0], ScanFlux[scan_index, spw_index, 0]+ uvMax* ScanSlope[scan_index, spw_index, 0]]), '-', color=Pcolor[0])
        StokesP_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 1], ScanFlux[scan_index, spw_index, 1]+ uvMax* ScanSlope[scan_index, spw_index, 1]]), '-', color=Pcolor[1])
        StokesP_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 2], ScanFlux[scan_index, spw_index, 2]+ uvMax* ScanSlope[scan_index, spw_index, 2]]), '-', color=Pcolor[2])
        StokesP_PL.plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 3], ScanFlux[scan_index, spw_index, 3]+ uvMax* ScanSlope[scan_index, spw_index, 3]]), '-', color=Pcolor[3])
        StokesI_PL.axis([0.0, uvMax, 0.0, 1.25*IMax]); StokesP_PL.axis([0.0, uvMax, -0.25*IMax, 0.25*IMax])
        StokesI_PL.text(0.0, 1.26*IMax, 'SPW%2d %5.1f GHz' % (spwList[spw_index], centerFreqList[spw_index]))
    #
    StokesI_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    StokesP_PL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    plt.show()
    figFL.savefig(pp, format='pdf')
    #
    freqArray = np.array(centerFreqList)[range(spwNum)]; meanFreq = np.mean(freqArray); relFreq = freqArray - meanFreq
    text_sd = ' --------------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print text_sd
    pflux, pfluxerr = np.zeros(4), np.zeros(4)
    text_sd = ' mean  %5.1f GHz' % (meanFreq); logfile.write(text_sd); print text_sd,
    for pol_index in range(4):
        sol, solerr = linearRegression(relFreq, ScanFlux[scan_index, :, pol_index], ErrFlux[scan_index, :, pol_index] ); pflux[pol_index], pfluxerr[pol_index] = sol[0], solerr[0]
        text_sd = '%7.4f (%.4f) ' % (pflux[pol_index], pfluxerr[pol_index]) ; logfile.write(text_sd); print text_sd,
        if SSO_flag and (pol_index > 0):  sol = np.zeros(2)
        spwStokesDic[sourceName] = spwStokesDic[sourceName] + (sol[0] + sol[1]* relFreq).tolist()
    #
    text_sd = '%6.3f   %6.1f ' % (100.0* np.sqrt(pflux[1]**2 + pflux[2]**2)/pflux[0], np.arctan2(pflux[2],pflux[1])*90.0/pi); logfile.write(text_sd); print text_sd,
    logfile.write('\n')
    text_sd = 'UV_min_max  %6.1f  %6.1f ' % (uvMin, uvMax); logfile.write(text_sd); print text_sd,
    logfile.write('\n \n'); print '\n'
    if SSO_flag:
        StokesDic[sourceName] = [pflux[0], 0.0, 0.0, 0.0]
    else:
        StokesDic[sourceName] = pflux.tolist()
        waveLength = 299.792458/meanFreq    # wavelength in mm
        text_sd = '%s, NE, NE, NE, NE, %.2fE+09, %.3f, %.3f, %.3f, %.3f, %.2f, %.2f, %.2f, %.2f, %s\n' % (sourceName, meanFreq, pflux[0], pfluxerr[0], np.sqrt(pflux[1]**2 + pflux[2]**2)/pflux[0], np.sqrt(pfluxerr[1]**2 + pfluxerr[2]**2)/pflux[0], np.arctan2(pflux[2],pflux[1])*90.0/pi, np.sqrt(pfluxerr[1]**2 + pfluxerr[2]**2)/np.sqrt(pflux[1]**2 + pflux[2]**2)*90.0/pi, uvMin/waveLength, uvMax/waveLength, timeLabel[0:10].replace('/','-'))
        ingestFile.write(text_sd)
    #
    if COMPDB & (not SSO_flag):
        print ' -------- Comparison with ALMA Calibrator Catalog --------'
        au.searchFlux(sourcename='%s' % (sourceList[sourceIDscan[scan_index]]), band=int(UniqBands[band_index][3:5]), date=timeLabel[0:10], maxrows=3)
        print '\n'
    #
    timePointer += timeNum
#
for PL in IList: figFL.delaxes(PL)
for PL in PList: figFL.delaxes(PL)
#
ingestFile.close()
logfile.close()
plt.close('all')
pp.close()
np.save(prefix + '-' + UniqBands[band_index] + '.Flux.npy', ScanFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Ferr.npy', ErrFlux)
np.save(prefix + '-' + UniqBands[band_index] + '.Source.npy', np.array(sourceList)[sourceIDscan])
np.save(prefix + '-' + UniqBands[band_index] + '.AZEL.npy', np.array([scanTime, OnAZ, OnEL, OnPA]))
np.save(prefix + '-' + UniqBands[band_index] + '.XYC.npy', np.array(XYC).reshape([len(XYC)/spwNum/2, spwNum, 2]))
np.save(prefix + '-' + UniqBands[band_index] + '.XYD.npy', np.array(XYD).reshape([len(XYC)/spwNum/2, spwNum, 2]))
#----
StokesI, QCpUS, UCmQS = np.zeros(timeSum), np.zeros(timeSum), np.zeros(timeSum)
DxNew, DyNew = np.zeros([UseAntNum, spwNum, UseChNum], dtype=complex), np.zeros([UseAntNum, spwNum, UseChNum], dtype=complex)
for spw_index in range(spwNum):
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = Freq * 1.0e-9
    for sourceName in sourceList:
        scanList = scanDic[sourceName]
        timeIndex = timeDic[sourceName]
        if len(scanList) < 1 : continue
        PA = np.array(PADic[sourceName]); timeNum = len(PA)
        CS, SN = np.cos(2.0* PA), np.sin(2.0* PA)
        Isol, Qsol, Usol = spwStokesDic[sourceName][spw_index], spwStokesDic[sourceName][spwNum + spw_index], spwStokesDic[sourceName][2*spwNum + spw_index]
        QCpUS[timeIndex] = Qsol* CS + Usol* SN
        UCmQS[timeIndex] = Usol* CS - Qsol* SN
        StokesI[timeIndex] = Isol
    #
    #print '  -- Determining D-term spectra for spw ' + `spwList[spw_index]`
    for ch_index in range(UseChNum):
        #DxNew[:,spw_index, ch_index], DyNew[:,spw_index, ch_index] = VisMuiti_solveD(VisSpec[spw_index][:,chRange[ch_index]], QCpUS, UCmQS, DxSpec[:,spw_index, chRange[ch_index]], DySpec[:,spw_index, chRange[ch_index]], StokesI)
        DxNew[:,spw_index, ch_index], DyNew[:,spw_index, ch_index] = VisMuiti_solveD(VisSpec[spw_index][:,chRange[ch_index]], QCpUS, UCmQS, [], [], StokesI)
    #
    DxSpec[:,:,chRange], DySpec[:,:,chRange] = DxNew, DyNew
    #
    for ant_index in range(UseAntNum):
        Dfile = prefix + '-B' + `int(UniqBands[band_index][3:5])` + '-SPW' + `spw_index` + '-' + antList[antMap[ant_index]] + '.DSpec.npy'
        Dterm = np.array([Freq, DxSpec[ant_index, spw_index].real, DxSpec[ant_index, spw_index].imag, DySpec[ant_index, spw_index].real, DySpec[ant_index, spw_index].imag])
        np.save(Dfile, Dterm)
    #
#
#----
msmd.close()
msmd.done()
del flagAnt, TrxFlag, gainFlag, Dflag, AntID, BPCaledXspec, BP_ant, Gain, GainP, Minv, SEFD, TrxList, TsysSPW, TsysBL, azelTime, azelTime_index, chAvgVis, W
