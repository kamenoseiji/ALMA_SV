#---- Script for Band-3 Astroholograpy Data
from scipy import stats
import matplotlib.pyplot as plt
execfile('/users/skameno/Scripts/interferometry.py')

def delayCal(Xspec, chBunch, timeBunch, startCH, stopCH):	# Xspec[chNum, blNum, timeNum]
	chNum, blNum, timeNum = Xspec.shape[0], Xspec.shape[1], Xspec.shape[2]
	antNum =  Bl2Ant(blNum)[0];  chNum /= chBunch; timeNum /= timeBunch
	XX = np.zeros([blNum, chNum, timeNum], dtype=complex)
	delayCalXX = np.zeros([blNum, chNum, timeNum], dtype=complex)
	for bl_index in range(blNum):
		ants = Bl2Ant(bl_index)
		#bl_RevLex = Ant2Bl_RevLex(ants[1], ants[0], antNum)
		#XX[bl_index] = specBunch(Xspec[:, bl_RevLex], chBunch, timeBunch)
		XX[bl_index] = specBunch(Xspec[:, bl_index], chBunch, timeBunch)
	#	
	delay_bl, amp_bl, delay_ant = np.zeros(blNum), np.zeros(blNum), np.zeros([antNum, timeNum])
	for time_index in range(timeNum):
	 	for bl_index in range(blNum):
 			delay_bl[bl_index], amp_bl[bl_index] = delay_search(XX[bl_index, startCH:stopCH, time_index])
		delay_ant[:,time_index] = cldelay_solve(delay_bl, 1.0/amp_bl)[0]
 		for bl_index in range(blNum):
			ants = Bl2Ant(bl_index); ant1, ant2 = ants[1], ants[0]
			delayCalXX[bl_index, :, time_index] = delay_cal(XX[bl_index, :, time_index], delay_ant[ant2, time_index] - delay_ant[ant1, time_index])
	#	
	return delay_ant, delayCalXX
#
#-------- Delay Cal for single-time cross-power spectra
def delayCalSpec( Xspec, chRange ):		# chRange = [startCH:stopCH] specifies channels to determine delay 
	chNum, blNum = Xspec.shape[0], Xspec.shape[1]
	delay_bl, amp_bl = np.zeros(blNum), np.zeros(blNum)
	delayCalXspec = np.zeros([blNum, chNum], dtype=complex)
	for bl_index in range(blNum):
		delay_bl[bl_index], amp_bl[bl_index] = delay_search(Xspec[bl_index, chRange])
	#
	delay_ant = cldelay_solve(delay_bl, 1.0/amp_bl)[0]
 	for bl_index in range(blNum):
		ants = Bl2Ant(bl_index); ant1, ant2 = ants[1], ants[0]
		delayCalXspec[bl_index] = delay_cal(Xspec[bl_index], delay_ant[ant2] - delay_ant[ant1])
	#	
	return delay_ant, delayCalXspec
#
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
"""
#Amplitude and Phase
def antBP(Xspec):
	blNum, chNum = Xspec.shape[0], Xspec.shape[1]
	antNum =  Bl2Ant(blNum)[0]
	#
	BPAmp_ant = np.zeros([antNum, chNum]); BPPhs_ant = np.zeros([antNum, chNum]); BP_ant = np.zeros([antNum, chNum], dtype=complex)
	for ch_index in range(chNum):
		BPAmp_ant[:,ch_index] = clamp_solve(abs(Xspec[:,ch_index]), 1.0e-3/abs(np.mean(Xspec, axis=1)) )[0]
		BPPhs_ant[:,ch_index] = clphase_solve(np.angle(Xspec[:,ch_index]), 1.0e-3/abs(np.mean(Xspec, axis=1)))[0]
	#
#-------- Segment Spectra
	BP_ant = BPAmp_ant* np.exp(1j* BPPhs_ant)
	return BP_ant
"""
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
antList = GetAntName(msfile)
startCH = int(round(chNum * 0.10)); stopCH  = int(round(chNum * 0.90))
blMap = range(blNum)
for bl_index in range(blNum):
	ants = Bl2Ant(bl_index)
	blMap[bl_index] = Ant2Bl(refant[ants[0]], refant[ants[1]])
#
for scan_index in range(scanNum):
#
	#-------- Procedures
	interval, timeStamp = GetTimerecord(msfile, 0, 0, pol[0], spw[0], scan[scan_index])
	timeNum = len(timeStamp)
	Gain_ant = np.zeros([antNum, timeNum], dtype=complex)
	gainCalSpec =  np.zeros([timeNum, polNum, spwNum, blNum, chNum/chBunch], dtype=complex)
	#
	#---- antenna-based bandpass
	bpCalSpec =  np.zeros([spwNum, polNum, blNum, chNum/chBunch, timeNum], dtype=complex)
	print 'Processing Scan = ' + `scan_index` + ' / ' + `scanNum`
	for spw_index in range(spwNum):
		print ' Loading SPW = ' + `spw[spw_index]`
		timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], scan[scan_index])
		BPac = np.ones(Xspec.shape)
		for bl_index in range(blNum):
			ants = Bl2Ant(blMap[bl_index])
			BPac[:,:,bl_index,:] = np.sqrt( Pspec[:,:,ants[0],:] * Pspec[:,:,ants[1],:] )
		#
		Xspec = Xspec[:,:,blMap,:] * BPac[:,:,blMap,:]
		for pol_index in range(polNum):
			print ' Delay Cal in pol ' + `pol_index` + '...'
			delay_ant, delayCalSpec = delayCal(Xspec[pol[pol_index]], chBunch, timeBunch, startCH, stopCH)
			#-------- Bandpass Phase Correction
			print ' Bandpass Phase Cal in pol ' + `pol_index` + '...'
			BPant = bpPhsAnt(delayCalSpec)
			bpCalSpec[spw_index, pol_index, :, :, :] = bpCal(delayCalSpec, BPant)
			#-------- Gain Calibration
			print ' Gain Cal in pol ' + `pol_index` + '...'
			for time_index in range(timeNum):
				Gain_ant[:,time_index] = gainAnt( np.mean(bpCalSpec[spw_index, pol_index, :, startCH:stopCH, time_index], axis=1), 1e-6/abs( np.mean(bpCalSpec[spw_index, pol_index, : ,startCH:stopCH, time_index], axis=1)))
			if( (scan_index == 0) & (spw_index == 0) & (pol_index == 0)):
				for ant_index in range(antNum):
					plt.plot( abs(Gain_ant[ant_index]), '.'); plt.plot( abs(Gain_ant[ant_index]), ls='steps-mid')
			temp = gainCal(bpCalSpec[spw_index, pol_index], Gain_ant )
			gainCalSpec[:, pol_index, spw_index, :, :] = np.transpose(temp, (2, 0, 1))
			#-------- Time Average and BL -> Ant solution
			print ' Baseline to Antenna in pol ' + `pol_index` + '...'
			BP_ant[scan_index, pol_index, spw_index] = antBP(np.mean(gainCalSpec, axis=0)[pol_index, spw_index])
		#
	#
	Vis[scan_index] = np.mean(np.mean(gainCalSpec, axis=3), axis=0)
	timeRange[scan_index] = min(timeStamp), max(timeStamp)
#
np.save(prefix[0] + '.BPant.npy', BP_ant) 
np.save(prefix[0] + '.Vis.npy', Vis) 
np.save(prefix[0] + '.Time.npy', timeRange) 
np.save(prefix[0] + '.Ant.npy', antList[refant]) 
#sd_spec, pe_spec, sd_phas, pe_phas = plotBPstat(BP_ant, antList[refant], ['XX','YY'], spw)
#histPlot(pe_spec, antList[refant], 1e-3)
#histPlot(pe_phas, antList[refant], 1e-3)
