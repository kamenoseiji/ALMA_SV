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
if 'bunchNum' in locals():
    def bunchVecN(Spec):
        return bunchVec(Spec, bunchNum)
    #
#-------- Loop for Scan
antDelay, scanTime = [], []
for scan in scanList:
    field_names = msmd.fieldsforscan(scan, True)
    print 'Loading Visibilities: Scan ' + `scan` + ' : ' + field_names[0]
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
    Xspec = np.mean(Xspec, axis=3)                          # Time average
    Xspec = np.apply_along_axis( bunchVecN, 1, Xspec )      # channel bunching     
    polNum, chNum = Xspec.shape[0], Xspec.shape[1]
    if polNum == 4: polIndex = [0, 3]
    if polNum == 2: polIndex = [0, 1]
    if polNum == 1: polIndex = [0]
    polNum = len(polIndex)
    Xspec = Xspec[polIndex]
    if polNum == 2:
        antGainSpec = np.array([np.apply_along_axis(clphase_solve, 1, Xspec[0]).T, np.apply_along_axis(clphase_solve, 1, Xspec[1]).T])
        antDelay = antDelay + [np.mean(np.apply_along_axis(delay_search, 2, antGainSpec), axis=0)[:,0]]
    else:
        antGainSpec = np.apply_along_axis(clphase_solve, 1, Xspec[0]).T
        antDelay = antDelay + [np.apply_along_axis(delay_search, 2, antGainSpec)[0][:,0]]
    #
    scanTime = scanTime + [np.median(timeStamp)]
#
msmd.done(); msmd.close()
np.save(prefix + '.Ant.npy', antList) 
np.save(prefix + '-SPW' + `spw` + '.TS.npy', np.array(scanTime)) 
np.save(prefix + '-SPW' + `spw` + '.DL.npy', np.array(antDelay).T) 
"""
#-------- Plot
fig = plt.figure(figsize = (8,11))
#-- pol loop
for pol_index in range(2):
    plt.subplot(2, 1, pol_index + 1 )
    plt.pcolormesh(antSNR[:,pol_index], cmap='gray', vmin=0.0, vmax=SNR_THRESH)
    #plt.pcolormesh(antFlag)
    plt.colorbar()
    plt.xlim(1, timeNum); plt.ylim(0, UseAntNum)
    #plt.xticks(arange(timeNum), timeStamp)
    plt.yticks(arange(UseAntNum)+0.5, antList[antMap])
#
mjdTick = range(int(min(timeStamp+59.9))/60*60, int(max(timeStamp))/60*60+1, 60)
tickPos = (max(timeStamp) - min(timeStamp)) / timeNum
timeLabel = []
for tickTime in mjdTick: timeLabel = timeLabel + [qa.time('%fs' % (tickTime), form='hm')[0]]


fig = plt.figure(figsize = (11,8))
ax = fig.add_subplot(2, 1, 1 )
ax.imshow(antSNR[:,0], cmap='gray', extent=[min(timeStamp), max(timeStamp), -0.5, UseAntNum-0.5], interpolation='none', clim=(0.0, 10.0), aspect=(max(timeStamp)-min(timeStamp))/(2.0* UseAntNum))
Xticks = ax.set_xticks(labelTime); labels = ax.set_xticklabels(timeLabel, rotation=90, fontsize='small')
Yticks = ax.set_yticks(range(UseAntNum)); labels = ax.set_yticklabels(antList[antMap], fontsize='small')
#
ax = fig.add_subplot(2, 1, 2)
ax.imshow(antSNR[:,1], cmap='gray', extent=[min(timeStamp), max(timeStamp), -0.5, UseAntNum-0.5], interpolation='none', clim=(0.0, 10.0), aspect=(max(timeStamp)-min(timeStamp))/(2.0* UseAntNum))
#
Xticks = ax.set_xticks(labelTime); labels = ax.set_xticklabels(timeLabel, rotation=90, fontsize='small')
Yticks = ax.set_yticks(range(UseAntNum)); labels = ax.set_yticklabels(antList[antMap], fontsize='small')

#plt.subplot(2, 1, 2 ); plt.imshow(antSNR[:,1], cmap='gray', extent=[min(timeStamp), max(timeStamp), 0, antNum], interpolation='none', clim=(0.0, 10.0), aspect=(max(timeStamp)-min(timeStamp))/(2.0* antNum))
#xi, yi = np.mgrid[0:antNum:1, min(timeStamp):max(timeStamp):(1.0j* timeNum)]
#plt.contourf(xi, yi, antSNR[:,0])
#plt.contourf(antSNR[:,0], xi, yi, np.linspace(0.0, 10.0))
#SNRmap = GridData(antSNR[:,0].real, ,timeStamp, xi.reshape(xi.size), yi.reshape(xi.size), 1.0).reshape(len(xi), len(xi))
"""
