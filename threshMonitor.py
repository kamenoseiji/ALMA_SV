import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import scipy
from scipy.fftpack import fft, ifft
from scipy.interpolate import interp1d
import scipy.optimize
import time
import datetime
import matplotlib.dates as mdates
import dateutil.parser
execfile(SCR_DIR + 'gaussNbit.py')

CAI_mjdSec = np.load('CIT.' + DATE + '.npy')
CAIMAP     = np.load('CCM.' + DATE + '.npy')
CAIANTNAME = np.load('CAN.' + DATE + '.npy')
#
XC = np.load('XLC.' + DATE + '.npy')
YC = np.load('YLC.' + DATE + '.npy')
mjdSec = np.load('TIM.' + DATE + '.npy')
CAI    = np.load('CAI.' + DATE + '.npy')
BBID   = np.load('BBI.' + DATE + '.npy')
#
CAI_index = np.where( CAIANTNAME == monitorAnt )[0].tolist()


CAIBBindex = np.where( (BBID == 0) & (CAI == CAI_index[0]) )[0]
flag  = np.where( XC[CAIBBindex, 0] > 0.75* np.median(XC[CAIBBindex,0]) )[0]     #  Filter by enough-power input
CAIBBindex = CAIBBindex[flag]
flag  = np.where( XC[CAIBBindex, 0] < 1.5* np.median(XC[CAIBBindex,0]) )[0]     #  Filter by enough-power input
CAIBBindex = CAIBBindex[flag].tolist()
timeNum = len(CAIBBindex)
thresh = np.zeros([2, 4, timeNum, 7])
for BB_index in range(4):
    for time_index in range(timeNum):
        thresh[0, BB_index, time_index] = threshLevel(XC[CAIBBindex[time_index]])
        thresh[1, BB_index, time_index] = threshLevel(YC[CAIBBindex[time_index]])
    #
#
bitPower = np.zeros([2, 4, timeNum])
for BB_index in range(4):
    for time_index in range(timeNum):
        solution, err = gaussNbitThresh( XC[CAIBBindex[time_index]], 8, thresh[0, BB_index, time_index] ); bitPower[0, BB_index, time_index] = 1/solution[0]**2
        solution, err = gaussNbitThresh( YC[CAIBBindex[time_index]], 8, thresh[1, BB_index, time_index] ); bitPower[1, BB_index, time_index] = 1/solution[0]**2
    #
#
flag = np.where( (bitPower[0,0] > 2.9) & (bitPower[0,0] < 3.1) )[0]


#for time_index in range(timeNum):
DateMeasure = np.array([ dateutil.parser.parse(qa.time('%fs' % (mjdSec[CAIBBindex[time_index]]), form='fits', prec=9)[0]) for time_index in range(timeNum)])
#
for levelIndex in range(7):
    plt.plot(DateMeasure[flag] , thresh[0,0,flag,levelIndex], '.')
    meanLevel = np.mean(thresh[0,0,flag,levelIndex]) 
    sdLevel   = np.std(thresh[0,0,flag,levelIndex])
    plt.plot( [np.min(DateMeasure), np.max(DateMeasure)], [meanLevel, meanLevel], color='gray')
    text_sd = 'mean=%5.3f SD=%5.3f' % (meanLevel, sdLevel)
    plt.text( DateMeasure[2000],  meanLevel + 0.05, text_sd, size='x-small')
#
