#---- Script for Band-3 Astroholograpy Data
from scipy import stats
import matplotlib.pyplot as plt
execfile('/users/skameno/Scripts/interferometry.py')

#Real and Imaginary Solution
def antBP(Xspec):
	blNum, chNum = Xspec.shape[0], Xspec.shape[1]
	antNum =  Bl2Ant(blNum)[0]
	#
	BP_ant = np.zeros([antNum, chNum], dtype=complex)
	bl_err = np.std(np.angle(Xspec[:, 1:(chNum-1)]), axis=1) 
	for ch_index in range(chNum):
		BP_ant[:, ch_index] = clcomplex_solve(Xspec[:, ch_index], bl_err)
	#
	return BP_ant
#-------- Definitions
antNum = len(refant)
blNum = antNum* (antNum - 1) / 2 
spwNum  = len(spw)
polNum  = len(pol)
#
#-------- Procedures
#-------- Loop for MS file (scan)
scanNum = len(scan)
timeRange = np.zeros([scanNum, 2])	# Start and End time period 
BP_ant = np.zeros([scanNum, polNum, spwNum, antNum, chNum/chBunch], dtype=complex)	# Antenna-based complex bandpass
Vis    = np.zeros([scanNum, polNum, spwNum, chNum/chBunch], dtype=complex)			# Baseline-averaged visibility
msfile = wd + prefix[0] + '.ms'
antList = GetAntName(msfile)[refant]
chRange = range(int(round(chNum/chBunch * 0.05)), int(round(chNum/chBunch * 0.95)))
blMap = range(blNum)
blInv = [False]* blNum		# True -> inverted baseline
for bl_index in range(blNum):
	ants = Bl2Ant(bl_index)
	blMap[bl_index], blInv[bl_index]  = Ant2BlD(refant[ants[0]], refant[ants[1]])
#

figAmp = plt.figure(figsize = (8, 11))
figPhs = plt.figure(figsize = (8, 11))
text_sd = '%s SPW=%d Gain (Amp)' % (prefix[0], spw[0])
figAmp.text(0.45, 0.95, text_sd)
figAmp.text(0.05, 0.45, 'Amplitude of Antenna-based Gain', rotation=90)
figAmp.text(0.45, 0.05, 'Time')
text_sd = '%s SPW=%d Gain (Phase)' % (prefix[0], spw[0])
figPhs.text(0.45, 0.95, text_sd)
figPhs.text(0.05, 0.45, 'Phase of Antenna-based Gain [rad]', rotation=90)
figPhs.text(0.45, 0.05, 'Time')
for scan_index in range(scanNum):
#
	#-------- Procedures
	interval, timeStamp = GetTimerecord(msfile, 0, 0, pol[0], spw[0], scan[scan_index])
	timeNum = len(timeStamp)
	BunchXspec =  np.zeros([polNum, chNum/chBunch, blNum, timeNum/timeBunch], dtype=complex)
	#
	#---- antenna-based bandpass
	bpCalSpec =  np.zeros([spwNum, polNum, blNum, chNum/chBunch, timeNum], dtype=complex)
	print 'Processing Scan = ' + `scan_index` + ' / ' + `scanNum`
	for spw_index in range(spwNum):
		print ' Loading SPW = ' + `spw[spw_index]`
		chNum, chWid, freq = GetChNum(msfile, spw[spw_index])
		timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], scan[scan_index]); Xspec = Xspec[:,:,blMap,:]
		for bl_index in range(blNum):
			if blInv[bl_index]:
				Xspec[:,:,bl_index,:] = Xspec[:,:,bl_index,:].conjugate()
			#
		#
		XAvg = np.mean(Xspec[:, chRange, :, :], axis=1)
		for pol_index in range(polNum):
			#Gain_ant = np.zeros([antNum, len(timeStamp)], dtype=complex)
			GainAnt  = np.zeros([antNum, timeNum], dtype=complex)
			#-------- Time and Spectral Bunch
			#for bl_index in range(blNum):
			#	BunchXspec[pol_index, :, bl_index, :] = specBunch(Xspec[pol_index, :, bl_index], chBunch, timeBunch )
			#
			#-------- Gain Calibration
			print ' Gain Cal in pol ' + `pol_index` + '...'
			for time_index in range(len(timeStamp)):
				vis_sigma = 1.0 / sqrt(np.median(abs(chWid))* interval[time_index]* 0.8* chNum)
				GainAnt[:,time_index] = clcomplex_solve(XAvg[pol_index, :, time_index], vis_sigma/(abs(XAvg[pol_index, :, time_index]) + 1.0e-9*vis_sigma))
			#
			"""
			for ant_index in range(antNum):
				interpRe = interp1d(timeStamp, Gain_ant[ant_index].real)
				interpIm = interp1d(timeStamp, Gain_ant[ant_index].imag)
				GainAnt[ant_index] = interpRe(timeStamp) + (1j)*interpIm(timeStamp)
				#ampPlot = figAmp.add_subplot(antNum, polNum, ant_index* polNum + pol_index + 1)
				#ampPlot.plot( timeStamp, abs(GainAnt[ant_index]), '.')
				#maxAmp = max(np.max(abs(GainAnt[ant_index])), maxAmp)
				#
				#phsPlot = figPhs.add_subplot(antNum, polNum, ant_index* polNum + pol_index + 1)
				#phsPlot.plot( timeStamp, np.angle(GainAnt[ant_index]), '.')
			#
			"""
			#temp = np.mean(gainCal(np.transpose(BunchXspec[pol_index], (1,0,2)), GainAnt ), axis=2)
			#temp = np.mean(gainCal(np.transpose(Xspec[pol_index], (1,0,2)), GainAnt ), axis=2)
			temp = np.mean(phaseCal(np.transpose(Xspec[pol_index], (1,0,2)), GainAnt ), axis=2)
			#-------- Delay Calibration
			if DELAYCAL:
				delay_ant, delay_err, delayCalXspec = delayCalSpec2( temp, chRange, vis_sigma )
			else:
				delayCalXspec = temp
			#-------- Visibility-averated spectrum
			Vis[scan_index, pol_index, spw_index] = np.mean(delayCalXspec, axis=0)
			#-------- Time Average and BL -> Ant solution
			print ' Baseline to Antenna in pol ' + `pol_index` + '...'
			BP_ant[scan_index, pol_index, spw_index] = antBP(delayCalXspec)
		#
	#
	#plt.plot(timeStamp, np.angle(GainAnt[1]), 'o')
	timeRange[scan_index] = min(timeStamp), max(timeStamp)
#
#figAmp.savefig('Amp.pdf')
#figPhs.savefig('Phs.pdf')
np.save(prefix[0] + '.BPant.npy', BP_ant) 
np.save(prefix[0] + '.Vis.npy', Vis) 
np.save(prefix[0] + '.Time.npy', timeRange) 
np.save(prefix[0] + '.Ant.npy', antList) 
np.save(prefix[0] + '.Scan.npy', scan) 
np.save(prefix[0] + '.Freq.npy', freq) 
plt.close('all')
