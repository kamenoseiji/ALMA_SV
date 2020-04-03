#---- Script for Band-3 Astroholograpy Data
import sys
from scipy import stats
import matplotlib.pyplot as plt
execfile(SCR_DIR + 'interferometry.py')
#-------- Initial Settings
if 'SNR_THRESH' not in locals(): SNR_THRESH = 3.0
msfile = wd + prefix + '.ms'; msmd.open(msfile)
antList = GetAntName(msfile)
antNum = len(antList)
blNum = antNum* (antNum - 1) / 2
spwName = msmd.namesforspws(spw)[0]
BandName = re.findall(r'RB_..', spwName)[0]; BandID = int(BandName[3:5])
chNum, chWid, Freq = GetChNum(msfile, spw); bandWidth = np.sum(chWid)
chRange = range(chNum)
if 'bandEdge' in locals():
    if bandEdge:
        chRange = range(int(0.05*chNum), int(0.95*chNum))
        bandWidth *= (len(chRange) + 0.0)/(chNum + 0.0)
    #
#
#-------- Array Configuration
print '---Checking array configuration'
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
text_sd = '  Usable antennas (%d) : ' % (len(UseAnt))
for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
print text_sd
blMap, blInv= range(UseBlNum), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
timeStamp, UVW = GetUVW(msfile, spw, scanList[0])
uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
if 'refant' in locals():    refantID = indexList(np.array([refant]), antList[UseAnt])[0]
else: refantID = bestRefant(uvDist)
print '  Use ' + antList[UseAnt[refantID]] + ' as the refant.'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
#-------- For channel bunching
if 'bunchCH' in locals(): bunchCH = 1
def bunchVecCH(Spec):
    return bunchVec(Spec, bunchCH)
#-------- For time bunching
if 'bunchTM' in locals(): bunchTM = 1
def bunchVecTM(Spec):
    return bunchVec(Spec, bunchTM)
#-------- Polarization List
polList = [[],[0],[0,1],[],[0,3]]
#-------- Loop for Scan
if 'scanTime' in locals(): del scanTime
if 'antDelay' in locals(): del antDelay
for scan in scanList:
    field_names = msmd.fieldsforscan(scan, True)
    print 'Loading Visibilities: Scan ' + `scan` + ' : ' + field_names[0]
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
    polNum, chNum, timeNum = Xspec.shape[0], Xspec.shape[1], Xspec.shape[3]
    tempSpec = np.apply_along_axis(bunchVecTM, 3,  np.apply_along_axis(bunchVecCH, 1, Xspec[polList[polNum]][:,chRange][:,:,blMap]))  # tempSpec[pol, ch, bl, time]
    tempTime = bunchVecTM(timeStamp)
    antGainSpec = np.apply_along_axis(gainComplex, 2, tempSpec)   # [pol, ch, ant, time]
    antDelayAmp = np.apply_along_axis(delay_search, 1, antGainSpec)    # [pol, delay-amp, ant, time]
    if 'antDelay' not in locals():
        antDelay = antDelayAmp[:,0].transpose(1,0,2) / (2.0* bandWidth)
    else: 
        antDelay = np.append(antDelay, antDelayAmp[:,0].transpose(1,0,2) / (2.0* bandWidth), axis=2)
    #
    if 'scanTime' not in locals():
        scanTime = tempTime
    else:
        scanTime = np.append(scanTime, tempTime)
#
msmd.done(); msmd.close()
np.save(prefix + '.Ant.npy', antList[antMap]) 
np.save(prefix + '-SPW' + `spw` + '.TS.npy', scanTime) 
np.save(prefix + '-SPW' + `spw` + '.DL.npy',  antDelay)
