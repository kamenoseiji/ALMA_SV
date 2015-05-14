import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import scipy
from scipy.fftpack import fft, ifft
from scipy.interpolate import interp1d
import scipy.optimize
import time
import datetime
execfile('/users/skameno/Scripts/gaussNbit.py')
def readACAbit( textfile ):
	file = open( textfile )
	line = True
	mjdSec = array([])
	while(line):
		try:
			line = file.readline()
			mjdSec = append(mjdSec, qa.convert(line.split()[5], 's')['value'])
			bitDist = map(int, line.split()[8:16])
			param, paramErr = gaussNbit(bitDist, 8)
			print param
		except:
			print bitDist
		#
	#
	file.close()
	return mjdSec
#
