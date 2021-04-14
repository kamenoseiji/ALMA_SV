#---- Script for Band-3 Astroholograpy Data
import sys
import subprocess
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
#
#-------- Procedures
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile); antNum = len(antList); blNum = antNum* (antNum - 1)/2
#-------- Configure Array
print '---Checking array configulation'
flagAnt = indexList(antFlag, antList)
UseAnt = list(set(range(antNum)) - set(flagAnt)); UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
blMap, blInv= range(UseBlNum), [False]* UseBlNum
try:
    refantID = np.where(antList[UseAnt] == refant )[0][0]
except:
    ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
    for bl_index in range(UseBlNum): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
    timeStamp, UVW = GetUVW(msfile, spwList[0], BPscan)
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = bestRefant(uvDist)
#
print '  Use ' + antList[UseAnt[refantID]] + ' as the refant.'
#-------- Baseline Mapping
print '---Baseline Mapping'
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
for bl_index in range(UseBlNum):
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
#-------- Bandpass Table
print '---Generating antenna-based bandpass table'
SideBand = ['LSB', 'USB']
BPList, XYList, XYdelayList = [], [], []
spwNum = len(spwList)
if 'bunchNum' not in locals(): bunchNum = 1
for spw_index in range(spwNum):
    if 'FGprefix' in locals():  # Flag table
        try:
            FG = np.load(FGprefix + '-SPW' + `spwList[spw_index]` + '.FG.npy'); FG = np.min(FG, axis=0)
            TS = np.load(FGprefix + '-SPW' + `spwList[spw_index]` + '.TS.npy')
            BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spwList[spw_index], BPscan, blMap, blInv, FG, TS)
        except:
            BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spwList[spw_index], BPscan, blMap, blInv)
        #
    else:
        BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spwList[spw_index], BPscan, blMap, blInv, bunchNum)
    #
    BPList = BPList + [BP_ant]
    XYList = XYList + [XY_BP]
    XYdelayList = XYdelayList + [XYD]
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index])
    BW = chNum* np.median(chWid)    # Bandwidth
    print 'SPW%2d: [%s] XY delay = %+f [ns] : SNR = %f' % (spwList[spw_index], SideBand[int(np.sign(np.median(chWid))+1)/2], 0.5* XYD / (BW * 1.0e-9), XYsnr)
#
BP_ant = np.array(BPList).transpose(1,0,2,3)
XYspec = np.array(XYList)
XYdelay = np.array(XYdelayList)
ppolNum = BP_ant.shape[2]
PolList = ['X', 'Y']
#
#-------- Save CalTables
np.save(prefix + '-REF' + antList[UseAnt[refantID]] + '.Ant.npy', antList[antMap]) 
for spw_index in range(spwNum):
    np.save(prefix + '-REF' + antList[UseAnt[refantID]] + '-SC' + `BPscan` + '-SPW' + `spwList[spw_index]` + '-BPant.npy', BP_ant[:,spw_index]) 
    np.save(prefix + '-REF' + antList[UseAnt[refantID]] + '-SC' + `BPscan` + '-SPW' + `spwList[spw_index]` + '-XYspec.npy', XYspec[spw_index]) 
    np.save(prefix + '-REF' + antList[UseAnt[refantID]] + '-SC' + `BPscan` + '-SPW' + `spwList[spw_index]` + '-XYdelay.npy', XYdelay[spw_index]) 
#
#-------- Plots
if BPPLOT:
    pp = PdfPages('XYP_' + prefix + '_REF' + antList[UseAnt[refantID]] + '_Scan' + `BPscan` + '.pdf')
    plotXYP(pp, prefix, spwList, XYspec, bunchNum) 
    pp = PdfPages('BP_' + prefix + '_REF' + antList[UseAnt[refantID]] + '_Scan' + `BPscan` + '.pdf')
    if 'spurRFLists' in locals():
        plotBP(pp, prefix, antList[antMap], spwList, BPscan, BPList, bunchNum, 1.2, spurRFLists) 
    else:
        plotBP(pp, prefix, antList[antMap], spwList, BPscan, BPList, bunchNum) 
#
