import numpy as np
from numpy import *
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import lines
import re
execfile('/Volumes/SSD/ALMA_SV/Scripts/gaussNbit.py')
execfile('/Volumes/SSD/ALMA_SV/Scripts/interferometry.py')
#---- Functions for data handling in Measurement Set
def readBitDist(fname):
    file = open(fname)
    readLines = file.readlines()
    file.close()
    #
    #return readLines
    bitCounts = np.zeros([len(readLines), 8], dtype=uint64)
    mjdSec = np.zeros(len(readLines))
    quadrant = np.zeros(len(readLines), dtype=uint32)
    CAI = np.zeros(len(readLines), dtype=uint32)
    pol = np.zeros(len(readLines), dtype=uint32)
    for index in range(len(readLines)):
        words = re.split(r' +', readLines[index])
        numWords = len(words)
        bitCounts[index] = String2BitCount( words[(numWords-9):(numWords-1)] )
        mjdSec[index] = qa.convert(words[5], 's')['value']
        quadrant[index] = int(words[2])
        CAI[index] = int(words[3])
        pol[index] = ord(words[4]) - ord('X')
    #        
    return mjdSec, quadrant, CAI, pol, bitCounts
#
def String2BitCount( string ):
    bitCount = np.zeros(8, dtype=uint64)
    for index in range(8):
        bitCount[index] = int(string[index])
    #
    return bitCount
#
def BB_filter(timeBB, dataBB):
    index = np.where( abs(dataBB[0]) > 0.0)[0]
    return timeBB[index], abs(dataBB[0,index])
#
def timeMatch( timeMin, timeMax, timeSeries ):
     return np.where( (timeSeries < timeMax) & (timeSeries > timeMin) )[0].tolist() 
#
#
#-------- Procedures
#
"""
bitFile = 'ACA20140930.bitDist'     # This file is generated via the script cf_extr-alma201402.py
prefix = 'uid___A002_X8df33c_X3c'   # 4-bit correction=on
CAIList = [4, 6, 7, 8, 9, 10, 13, 14, 15, 0, 1, 3, 2]
BBList  = [0, 1, 2, 3]
SPWList = [5, 7, 9, 11]
polList = ['X', 'Y']
#
bitFile = 'ACA20140929.bitDist'     # This file is generated via the script cf_extr-alma201402.py
prefix = 'uid___A002_X8dd3fc_X3e'   # 4-bit correction=on
CAIList = [4, 6, 7, 8, 9, 10, 13, 14, 15, 0, 3, 2]
BBList  = [0, 1, 2, 3]
SPWList = [5, 7, 9, 11]
polList = ['X', 'Y']
#
bitFile = 'ACA20140927.bitDist'     # This file is generated via the script cf_extr-alma201402.py
prefix = 'uid___A002_X8da5e6_X2a'   # 4-bit correction=on
CAIList = [4, 6, 7, 8, 9, 10, 13, 14, 15, 0, 3, 2]
BBList  = [0, 1, 2, 3]
SPWList = [5, 7, 9, 11]
polList = ['X', 'Y']
#
bitFile = 'ACA20140924.bitDist'     # This file is generated via the script cf_extr-alma201402.py
prefix = 'uid___A002_X8d6837_X77'   # 4-bit correction=on
#prefix = 'uid___A002_X8d6837_X12b'    # 4-bit correction=off
CAIList = [4, 6, 7, 8, 9, 10, 13, 14, 15, 1, 3, 2]   # This mapping is described in /alma/logs/COJ-CC-1/ACACORR/OBSERVATION_CONTROL/acsStartContainer_cppContainer_2014-09-23_21.57.51.010
#
"""
BBList  = [0, 1, 2, 3]
SPWList = [5, 7, 9, 11]
polList = ['X', 'Y']
bitFile = 'ACA20140926.bitDist'     # This file is generated via the script cf_extr-alma201402.py
#prefix = 'uid___A002_X8d8bc4_X60'   #
prefix = 'uid___A002_X8d8bc4_X91'   #
msfile = prefix + '.ms'
antList = GetAntName(msfile)
antNum  = len(antList)
BBNum   = len(BBList)
polNum  = len(polList)
ACAtimeThresh = 5.0            # [sec], tolerance for time delay between ACA 
BBtimeThresh  = 5.0            # [sec], tolerance for time delay between ACA 
mjdSec, quadrant, CAI, pol, bitCounts = readBitDist(bitFile)
#bitCounts = bitCounts[:,::2] + bitCounts[:,1::2]    # 8-level -> 4-level
#--------- antenna list -> CAI mapping
CAIList = []
CAI_ant = ['DV01', 'PM02', 'PM04', 'PM03', 'CM01', 'CM02', 'CM03', 'CM04', 'CM05', 'CM06', 'CM07', 'CM08', 'CM09', 'CM10', 'CM11', 'CM12']
for ant_index in range(antNum):
    CAIList = CAIList + [where(array(CAI_ant) == antList[ant_index])[0][0]]
    #
#
plotChAvg    = np.zeros([antNum, BBNum, polNum, 12])
plotBitPower = np.zeros([antNum, BBNum, polNum, 12])
plotBBPower  = np.zeros([antNum, BBNum, polNum, 12])
plotThresh   = np.zeros([antNum, BBNum, polNum, 12, 7])
#for ant_index in range(1):
f = open('CSV2311_hist.table', 'w')
f.write('ant  BB pol level PE\n')
for ant_index in range(antNum):
    fig = plt.figure(figsize = (11, 8))
    plt.suptitle(prefix + ' ' + antList[ant_index])
    #-------- Plot Time-Power diagram
    for BB_index in range(BBNum):
        chNum, chWid, freq = GetChNum(msfile, BBList[BB_index]); freq = freq * 1.0e-9    # Hz -> GHz
        chRange = range( int(0.1*chNum), int(0.9*chNum))
        for pol_index in range(polNum):
            ax = plt.subplot(polNum,  BBNum, BBNum* pol_index + BB_index + 1)
            #-------- Plot Autocorr ChAvg
            timeXY, dataXY = GetPSpec(msfile, ant_index, pol_index, SPWList[BB_index])
            chAvg = np.mean(abs(dataXY), axis=0)
            plt.plot(timeXY, chAvg/chAvg[0], '.', color='green')
            #-------- Plot BB detector power
            timeBB, dataBB = GetPSpec(msfile, ant_index, pol_index, BBList[BB_index])
            timeBB,  dataBB = BB_filter( timeBB, dataBB )
            plt.plot(timeBB, 1e3* dataBB, '.', color='blue')
            #
            #-------- Bit power time matching
            bd_index = np.where( (CAI==CAIList[ant_index]) & (pol==pol_index) & (quadrant == BB_index))[0]
            bitPower_index = bd_index[timeMatch( min(timeBB), max(timeBB)+4, mjdSec[bd_index])]
            bitPower = np.zeros( len(bitPower_index) )
            thresh   = np.zeros([len(bitPower_index), 7])   # 8 -level
            refThresh = threshLevel(bitCounts[bitPower_index[0],:])     # Reference Threshold
            refThresh = (refThresh - refThresh[3])* 6.0 / (refThresh[6] - refThresh[0])     # Normalization of thresholds
            for index in range(len(bitPower_index)):
                bitPower[index] = 1.0 / gaussNbitThresh( bitCounts[bitPower_index[index],:], 8, refThresh)[0][0]**2
                #thresh[index] =  threshLevel(bitCounts[bitPower_index[index],:])
                thresh[index] =  threshLevelRef(bitCounts[bitPower_index[index],:], refThresh)
                plotBitPower[ant_index, BB_index, pol_index, index]    = bitPower[index]
                plotThresh[ant_index, BB_index, pol_index, index, :]   = thresh[index, :]
            #
            #-------- Time-Power plot
            plt.plot(mjdSec[bitPower_index], bitPower, 'o', color='red')
            plt.axis([min(timeBB)-4, max(timeBB)+4, 0, 15], fontsize=3)
            if BB_index == 0:
                plt.ylabel('Power (arbitrary unit)')
            if pol_index == 1:
                plt.xlabel('MJD [sec]')
            text_sd = antList[ant_index] + ' BB=' + `BBList[BB_index]` + ' Pol=' + polList[pol_index]
            plt.text( min(timeBB), 13, text_sd, size='x-small')
            plt.legend(('BB (x 1000)', 'Level Hisogram'), 'upper left', prop={'size' :7})
            #-------- Time average for every subscan
            BB_endIndex = np.where( diff(timeBB) > 1.0 )[0]
            BB_startIndex = np.append([0], BB_endIndex + 1)
            BB_endIndex   = np.append(BB_endIndex, [len(timeBB)-1])
            for scan_index in range( len(BB_endIndex) ):
                plotBBPower[ant_index, BB_index, pol_index, scan_index] = np.mean(dataBB[range(BB_startIndex[scan_index], BB_endIndex[scan_index])])
            #
            for scan_index in range(len(bitPower_index)):
                index = timeMatch(mjdSec[bitPower_index[scan_index]]-0.5, mjdSec[bitPower_index[scan_index]]+0.5, timeXY)
                plotChAvg[ant_index, BB_index, pol_index, scan_index] = np.mean(chAvg[index])
            #
        #
    #
    plt.savefig(prefix + '-' + antList[ant_index] + '-' + '3bitPowerTime.pdf')
    plt.close()
    #-------- Plot BB-bitpower diagram
    fig = plt.figure(figsize = (11, 8))
    plt.suptitle(prefix + ' ' + antList[ant_index])
    xlim=[0.0, 4.5]; ylim=[0., 15.]
    for BB_index in range(BBNum):
        for pol_index in range(polNum):
            ax = plt.subplot(polNum,  BBNum, BBNum* pol_index + BB_index + 1)
            xpower = 1e3*plotBBPower[ant_index, BB_index, pol_index]; xpower = xpower[ np.where(xpower != 0.0)[0]]
            ypower = plotBitPower[ant_index, BB_index, pol_index];    ypower = ypower[ np.where(xpower != 0.0)[0]]
            zpower = plotChAvg[ant_index, BB_index, pol_index];       zpower = zpower[ np.where(xpower != 0.0)[0]]
            plt.plot(xpower, ypower, 'o', color='blue')
            plt.plot(xpower, zpower/zpower[0], 'o', color='green')
            plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]], fontsize=3)
            if BB_index == 0:
                plt.ylabel('Level-Histogram Power')
            if pol_index == 1:
                plt.xlabel('BB Power (x 1000)')
            #
            text_sd = antList[ant_index] + ' BB=' + `BBList[BB_index]` + ' Pol=' + polList[pol_index]; plt.text( 0, 0.95*ylim[1], text_sd, size='x-small')
            slope, intercept, r_value, _, _ = stats.linregress(xpower[0:-1], ypower[0:-1])
            text_sd = 'slope=%5.2e icept=%5.2e' % (slope, intercept);  plt.text( 0.0, 0.82*ylim[1], text_sd, size='x-small')
            line = lines.Line2D([0, xlim[1]], [intercept, slope* xlim[1] + intercept], color='black'); ax.add_line(line)
            slope, intercept, r_value, _, _ = stats.linregress(ypower, zpower)
            text_sd = 'slope=%5.2e icept=%5.2e' % (slope, intercept);  plt.text( 0.0, 0.70*ylim[1], text_sd, size='x-small')
            #plt.plot(xpower, ypower - (slope* xpower + intercept), 'o', color='red')
            #line = lines.Line2D([0, xlim[1]], [0.0, 0.0], color='red'); ax.add_line(line)
            #print '%s %02d %s %6.2e %6.2e' % (antList[ant_index], BB_index, polList[pol_index], slope, intercept)
        #
    #
    plt.savefig(prefix + '-' + antList[ant_index] + '-' + '3bitPowerPower.pdf')
    plt.close()
    #-------- Plot Threshold diagram
    fig = plt.figure(figsize = (11, 8))
    plt.suptitle(prefix + ' ' + antList[ant_index])
    xlim=[0.0, 5.0]; ylim=[-4., 4.]
    for BB_index in range(BBNum):
        for pol_index in range(polNum):
            ax = plt.subplot(polNum,  BBNum, BBNum* pol_index + BB_index + 1)
            for level_index in range(7):
                plt.axhline(y=level_index - 3.0, color='gray')
                xpower = 1e3*plotBBPower[ant_index, BB_index, pol_index]; xpower = xpower[ np.where(xpower != 0.0)[0]]
                ypower = plotThresh[ant_index, BB_index, pol_index, :, level_index];    ypower = ypower[ np.where(xpower != 0.0)[0]]
                #plt.plot(1e3*plotBBPower[ant_index, BB_index, pol_index], plotThresh[ant_index, BB_index, pol_index, :, level_index], 'o')
                plt.plot(xpower, ypower, 'o')
                #PE = plotThresh[ant_index, BB_index, pol_index, :, level_index] - level_index + 3.0
                PE = ypower - level_index + 3.0
                text_sd = '%4.1f%%' %(100.0*max(PE)); plt.text(3.5, level_index - 2.9, text_sd);
                text_sd = '%4.1f%%' %(100.0*min(PE)); plt.text(3.5, level_index - 3.3, text_sd);
                for index in range(len(PE)):
                    f.write('%s %02d %s %02d %6.3f\n' % (antList[ant_index], BB_index, polList[pol_index], level_index-3, 100.0*PE[index]) )
            #
            plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]], fontsize=3)
            text_sd = antList[ant_index] + ' BB=' + `BBList[BB_index]` + ' Pol=' + polList[pol_index]; plt.text( 0, 3.6, text_sd, size='x-small')
            if BB_index == 0:
                plt.ylabel('Estimated Threshold Voltage')
            if pol_index == 1:
                plt.xlabel('BB Power (x 1000)')
            #
        #
    #
    plt.savefig(prefix + '-' + antList[ant_index] + '-' + '3bitThresh.pdf')
    plt.close()
#
f.close()
