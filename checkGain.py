#---- Script for Band-3 Astroholograpy Data
from scipy import stats
import matplotlib.pyplot as plt
execfile(SCR_DIR + 'interferometry.py')
#-------- Initial Settings
if 'SNR_THRESH' not in locals(): SNR_THRESH = 3.0
msfile = prefix + '.ms'; msmd.open(msfile)
antList = GetAntName(msfile)
antNum = len(antList)
blNum = antNum* (antNum - 1) / 2
spwName = msmd.namesforspws(spw)[0]
BandName = re.findall(r'RB_..', spwName)[0]; BandID = int(BandName[3:5])
#-------- Array Configuration
print '---Checking array configuration'
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
text_sd = '  Usable antennas: '
for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
print text_sd
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
timeStamp, UVW = GetUVW(msfile, spw, msmd.scansforspw(spw)[0])
uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
if 'refant' in locals():    refantID = indexList(np.array([refant]), antList)[0]
else: refantID = bestRefant(uvDist)
print '  Use ' + antList[UseAnt[refantID]] + ' as the refant.'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
msmd.done(); msmd.close()
#-------- Bandpass Table
if 'BPprefix' in locals():
    BPfileName = BPprefix + '-REF' + antList[UseAnt[refantID]] +'-SPW' + `spw` + '-BPant.npy'
    print '---Loading bandpass table : ' + BPfileName
    BP_ant = np.load(BPfileName)
#
#-------- Loop for Scan
GainAP0, GainAP1, timeList, uvwList, flagList = [], [], [], [], []
for scan in scanList:
    print 'Processing Scan ' + `scan`
    #-------- Baseline-based cross power spectra
    timeUVW, UVW = GetUVW(msfile, spw, scan)
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
    timeNum, polNum, chNum = Xspec.shape[3], Xspec.shape[0], Xspec.shape[1]
    if polNum == 4: polIndex = [0, 3]
    if polNum == 2: polIndex = [0, 1]
    if polNum == 1: polIndex = [0]
    polNum = len(polIndex)
    tempSpec = ParaPolBL(Xspec[polIndex][:,:,blMap], blInv).transpose(3,2,0,1)  # Parallel Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    if 'BP_ant' in locals():
        BPCaledXspec = (tempSpec / (BP_ant[ant0]* BP_ant[ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    else:
        BPCaledXspec = tempSpec.transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    for time_index in range(timeNum):
        gainFlag = np.ones(antNum)
        tempGain, tempErr = gainComplexErr(chAvgVis[0, :, time_index]); GainAP0 = GainAP0 + [tempGain]
        gainFlag[np.where(abs(tempGain) / tempErr  < SNR_THRESH)[0]] = 0.0
        tempGain, tempErr = gainComplexErr(chAvgVis[1, :, time_index]); GainAP1 = GainAP1 + [tempGain]
        gainFlag[np.where(abs(tempGain) / tempErr  < SNR_THRESH)[0]] = 0.0
        flagList = flagList + [gainFlag]
    #
    timeList.extend(timeStamp.tolist())
    uvwList = uvwList + [UVW[:,blMap]]
#
antFlag = np.array(flagList).T                          # [ant, time]
Gain = np.array([GainAP0, GainAP1]).transpose(2,0,1)    # [ant, pol, time]
np.save(prefix + '.Ant.npy', antList[antMap]) 
np.save(prefix + '-SPW' + `spw` + '.TS.npy', np.array(timeList)) 
np.save(prefix + '-SPW' + `spw` + '.GA.npy', Gain) 
np.save(prefix + '-SPW' + `spw` + '.FG.npy', antFlag) 
"""
Gain, uvw = np.array(GainAP) , uvwList[0]
for scan_index in range(1, len(scanList)):
    Gain = np.append(Gain, GainAP[scan_index], axis=2) 
    uvw  = np.append(uvw, uvwList[scan_index], axis=2) 
#
np.save(prefix + '.UVW.npy', uvw[:,KERNEL_BL[0:(antNum-1)]]) 
np.save(prefix + '-SPW' + `spw` + '.TS.npy', np.array(timeList)) 
np.save(prefix + '-SPW' + `spw` + '.GA.npy', Gain) 
"""
