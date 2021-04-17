import sys
import subprocess
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
wd = './'
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'ASDM_XML.py')
BPPLOT   = True
#prefix = 'uid___A002_Xe20b32_X6931'
#prefix = 'uid___A002_Xe37224_Xb141'
msfile = prefix + '.ms'
polName = ['X', 'Y']
spurLog = open(prefix + '-LO2Spur.log', 'w')
#-------- Get LO1 and LO2 frequencies
LO1, LO2List = BBLOfreq(prefix)
BBNum = len(LO2List)
SpuriousRFLists = []
#-------- Predicted RF frequencies of spurious signals
for BBIndex1 in range(BBNum):
    SpuriousRFList = []
    for BBIndex2 in list(set(range(BBNum)) - set([BBIndex1])):
        SpuriousRFList = SpuriousRFList + [ LO1 + LO2List[BBIndex1] - abs(LO2List[BBIndex2] - LO2List[BBIndex1]) ] # for USB
        SpuriousRFList = SpuriousRFList + [ LO1 - LO2List[BBIndex1] + abs(LO2List[BBIndex2] - LO2List[BBIndex1]) ] # for LSB
    #
    SpuriousRFLists = SpuriousRFLists + [SpuriousRFList]
#
#-------- Check Bandpass SPWs and scans
msmd.open(msfile)
BPScanList = msmd.scansforintent("CALIBRATE_BANDPASS*").tolist()
BPspwList  = list(set(msmd.fdmspws()) & set(msmd.spwsforintent("CALIBRATE_BANDPASS*")))
if 'spwFlag' in locals():
    flagIndex = indexList(np.array(spwFlag), np.array(BPspwList))
    for index in flagIndex: del BPspwList[index]
#
spwNum = len(BPspwList)
#-------- Check BB for SPW
BPspwNameList = msmd.namesforspws(BPspwList)
BBspwList = range(spwNum)
for spw_index in range(spwNum):
    BBspwList[spw_index] = int([BBname for BBname in BPspwNameList[spw_index].split('#') if 'BB_' in BBname][0].split('_')[1]) - 1
#
#-------- Check antenna List
antList = GetAntName(msfile)
antNum = len(antList)
antDia = np.ones(antNum)
for ant_index in range(antNum): antDia[ant_index] = msmd.antennadiameter(antList[ant_index])['value']
pPol = [0,1]
polNum = msmd.ncorrforpol(msmd.polidfordatadesc(BPspwList[0]))
if polNum == 4: pPol = [0, 3]
msmd.close()
#-------- Check frequency range of each SPW
spurSPWList, spurRFLists = [], []
for spwIndex in range(spwNum):
    BB_index = BBspwList[spwIndex]
    chNum, chWid, freq = GetChNum(msfile, BPspwList[spwIndex])
    spurFlag = False
    spurRFList = []
    for SpurRF in SpuriousRFLists[BB_index]:
        if SpurRF > min(freq) and SpurRF < max(freq):
            spurFlag = True
            spurRFList = spurRFList + [SpurRF]
            text_sd = 'SPW %d : BW=%.1f - %.1f : LO2@ %.3f GHz' % (BPspwList[spwIndex], min(freq)*1.0e-9, max(freq)*1.0e-9, SpurRF*1.0e-9)
            spurLog.write(text_sd + '\n'); print text_sd
        #
    #
    if spurFlag:
        spurSPWList = spurSPWList + [BPspwList[spwIndex]]
        spurRFLists = spurRFLists + [list(set(spurRFList))]
    #
#-------- check usable antennas
if 'antFlag' not in locals(): antFlag = []
if len(spurSPWList) == 0:
    print 'No LO2 leakages'
else:
    spwList = spurSPWList
    BPscan  = BPScanList[0]
    flagAnt = np.ones([antNum])
    for spw_index in range(len(spwList)):
        #-------- Checking usable baselines and antennas
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPscan)
        timeNum, chNum, blNum = Xspec.shape[3], Xspec.shape[1], Xspec.shape[2]; chRange, timeRange = range(int(0.05*chNum), int(0.95*chNum)), range(timeNum-4, timeNum-3)
        for polID in pPol:
            blD, blA = np.apply_along_axis(delay_search, 0, np.mean(Xspec[polID][chRange][:,:,timeRange], axis=2))
            blA = blA / np.sqrt(antDia[ANT0[0:blNum]]* antDia[ANT1[0:blNum]])
            errD = np.where(abs(blD - np.median(blD)) > 4.0)[0].tolist()
            errA = np.where(blA / np.median(blA) > 2.5)[0].tolist() + np.where(blA / np.median(blA) < 0.4)[0].tolist()
            errCount = np.zeros(antNum)
            for bl in set(errD) or set(errA): errCount[ list(Bl2Ant(bl)) ] += 1
            flagAnt[np.where(errCount > 2.5 )[0].tolist()] *= 0.0
        #
    #
    UseAnt = np.where(flagAnt > 0.0)[0]; UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
    print `UseAntNum` + ' / ' + `antNum` + ' usable antennas'
    flagAntList = np.where(flagAnt < 0.1)[0].tolist()
    antFlag = list(set(antFlag + antList[flagAntList].tolist()))
    #-------- antenna-based Bandpass
    execfile(SCR_DIR + 'checkBP.py')
    #-------- check Spur SNR
    for spw_index in range(len(spwList)):
        chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index])
        for spurRF in spurRFLists[spw_index]:
            spurCH = np.where( abs(Freq - spurRF) < abs(np.median(chWid)))[0].tolist()
            spurBL = range(min(spurCH) - 16, min(spurCH) - 2) + range(max(spurCH) + 2, min(spurCH) + 16)
            for ant_index in range(UseAntNum):
                for pol_index in range(2):
                    BPmean, BPsigma = np.mean(abs(BPList[spw_index][ant_index, pol_index][spurBL])), np.std(abs(BPList[spw_index][ant_index, pol_index][spurBL]))
                    SNR =  abs(abs(BPList[spw_index][ant_index, pol_index][spurCH]) - BPmean) / (BPsigma + 1.0e-8)
                    if max(SNR) > 3.0:
                        text_sd = 'Spur %s SPW=%d POL-%s %.3f GHz (SNR = %.1f)' % (antList[UseAnt[ant_index]], spwList[spw_index], polName[pol_index], 1.0e-9* spurRF, max(SNR))
                        spurLog.write(text_sd + '\n'); print text_sd
                    #
                #
            #
        #
    #
#
spurLog.close()
