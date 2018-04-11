#---- Script to find delay solutions
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
execfile('/users/skameno/Scripts/interferometry.py')

msfile = prefix + '.ms'
antList = GetAntName(msfile)[refant]
antNum = len(refant)
spwNum = len(spwList)
polNum = len(polList)
blNum = antNum* (antNum - 1) / 2
blMap = range(blNum)
blInv = [False]* blNum      # True -> inverted baseline
#-------- Mapping baselines into canonical order
for bl_index in range(blNum):
	ants = Bl2Ant(bl_index)
	blMap[bl_index], blInv[bl_index] = Ant2BlD(refant[ants[0]], refant[ants[1]])
#
#-------- Loop for scans
scanNum = len(scan)
scanPtr = np.zeros(scanNum)		# Cumulative time num
timeStamp = np.zeros(0)
for scan_index in range(scanNum):
	interval, tempTime = GetTimerecord(msfile, 0, 0, spwList[0], scan[scan_index])
	timeStamp = np.r_[timeStamp, tempTime]
	if(scan_index < (scanNum - 1)):
		scanPtr[scan_index + 1] = len(timeStamp)
	#
#
timeNum = len(timeStamp)
delay_ant = np.zeros([antNum, spwNum, polNum, timeNum])
delay_err = np.zeros([antNum, spwNum, polNum, timeNum])
gain_ant   = np.zeros([antNum, spwNum, polNum, timeNum], dtype=complex)
#-------- Determine Antenna-based Delay and Gain
for scan_index in range(scanNum):
	for spw_index in range(spwNum):
		interval, tempTime = GetTimerecord(msfile, 0, 0, spwList[spw_index], scan[scan_index])
		chNum, chWid, freq = GetChNum(msfile, spwList[spw_index])
		tempTime, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], scan[scan_index])
		chRange = range(int(0.1*chNum), int(0.9*chNum))
		Xspec = Xspec[:,:,blMap,:] / np.median(Pspec)
		#------ Inverted Baselines
		for bl_index in range(blNum):
			if blInv[bl_index]:
				Xspec[:,:,bl_index,:] = Xspec[:,:,bl_index,:].conjugate()
			#
		#
		delay_bl = np.zeros(blNum,); amp_bl = np.zeros(blNum);
		for time_index in range(len(tempTime)):
			vis_sigma = 1.0 / sqrt(np.median(abs(chWid))* interval[time_index]* 0.8* chNum)
			print 'SPW=%d Scan=%d: %d / %d' % (spwList[spw_index], scan[scan_index], time_index, len(tempTime))
			for pol_index in range(polNum):
				delay_ant[:, spw_index, pol_index, (scanPtr[scan_index] + time_index)], delay_err[:, spw_index, pol_index, (scanPtr[scan_index] + time_index)], delayCalSpec = delayCalSpec2(Xspec[pol_index, :, :, time_index].T, chRange, vis_sigma)
				visBl = np.mean(delayCalSpec[:, chRange] + 1.0e-30, axis=1)	# add 1e-30 to avoid division by 
				gain_ant[:, spw_index, pol_index, (scanPtr[scan_index] + time_index)] = clcomplex_solve( visBl, vis_sigma/abs(visBl) )
			#
		#
	#
#
np.save(prefix + '_delay.npy', delay_ant)
np.save(prefix + '_gain.npy',  gain_ant)
#
maxGain = np.max(abs(gain_ant))
for ant_index in range(antNum):
	#--------- Plot Delay
	figDelay = plt.figure(figsize = (11, 8))
	figDelay.suptitle( 'Residual Delay (Refant = ' + antList[0] + '): ' + prefix + ' ' + antList[ant_index])
	for spw_index in range(spwNum):
		for pol_index in range(polNum):
			fig_Delay = plt.subplot(polNum, spwNum, spwNum* pol_index + spw_index + 1)
			fig_Delay.errorbar(timeStamp, delay_ant[ant_index, spw_index, pol_index], yerr=delay_err[ant_index, spw_index, pol_index], fmt='.', capsize=0)
			if(spw_index == spwNum - 1):
				plt.xlabel('MJD [sec]')
			#
			text_sd = 'SPW=%d Pol%s' % (spwList[spw_index], polList[pol_index])
			plt.text( min(timeStamp), 0.9*max(delay_ant[ant_index, spw_index, pol_index]) + 0.1*min(delay_ant[ant_index, spw_index, pol_index]), text_sd, size='x-small')
		#
	plt.savefig( prefix + '-' + antList[ant_index] + '-delay.pdf')
	plt.close()
	#--------- Plot Gain
	figGain = plt.figure(figsize = (11, 8))
	figGain.suptitle( 'Antenna Gain: ' + prefix + ' ' + antList[ant_index])
	for spw_index in range(spwNum):
		for pol_index in range(polNum):
			fig_Gain = plt.subplot(polNum, spwNum, spwNum* pol_index + spw_index + 1)
			fig_Gain.plot(timeStamp, abs(gain_ant[ant_index, spw_index, pol_index]), '.')
			if(spw_index == spwNum - 1):
				plt.xlabel('MJD [sec]')
			#
			plt.ylim(0.0, 1.05* maxGain)
			text_sd = 'SPW=%d Pol%s' % (spwList[spw_index], polList[pol_index])
			plt.text( min(timeStamp), 0.9*maxGain, text_sd, size='x-small')
		#
	plt.savefig( prefix + '-' + antList[ant_index] + '-gain.pdf')
	plt.close()
#
figGain = plt.figure(figsize = (8, 11))
figGain.suptitle( 'Antenna Gain: ' + prefix)
for ant_index in range(antNum):
	fig_Gain = plt.subplot(antNum, 1, ant_index + 1)
	fig_Gain.plot(timeStamp, abs(gain_ant[ant_index, 0, 0]), '.')
	if(ant_index == antNum - 1):
		plt.xlabel('MJD [sec]')
	#
	plt.ylim(0.0, 1.05* maxGain)
	text_sd = '%s' % (antList[ant_index])
	plt.text( min(timeStamp), 0.9*maxGain, text_sd, size='x-small')
#
plt.savefig( prefix + '-gainSummary.pdf')
plt.close('all')
