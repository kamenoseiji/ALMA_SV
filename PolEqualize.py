print '---- Gain Equalization between X and Y polarizations'
interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spw[0], EQScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PA = np.arctan2( np.sin(PA), np.cos(PA))
scan_index = onsourceScans.index(EQScan)
caledVis = []
for spw_index in range(spwNum):
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], EQScan)
    timeNum, chNum = Xspec.shape[3], Xspec.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1) 
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    GainP = np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[3])])
    SEFD = 2.0* kb* chAvgTsys[:,spw_index, :,scan_index] / np.array([AeX, AeY]).T
    SEFD /= (relGain[spw_index]**2)
    caledVis.append(np.mean((chAvgVis / (GainP[polYindex][:,ant0]* GainP[polXindex][:,ant1].conjugate())).transpose(2, 0, 1)* np.sqrt(SEFD[ant0][:,polYindex].T* SEFD[ant1][:,polXindex].T), axis=2).T)
#
caledVisAmp = abs(np.mean(np.array(caledVis), axis=(0,2))[[0,3]])
catalogQU = np.array([catalogStokesQ.get(EQcal, 0.0), catalogStokesU.get(EQcal, 0.0)])/catalogStokesI.get(EQcal, 1.0)
CS_SN = np.array([np.mean(np.cos(PA)), np.mean(np.sin(PA))])
GainCorrR = 0.5* np.diff(caledVisAmp) / np.mean(caledVisAmp) - catalogQU.dot(CS_SN)
relGain[:,:,0] /= (1.0 + 0.5*GainCorrR)
relGain[:,:,1] /= (1.0 - 0.5*GainCorrR)
