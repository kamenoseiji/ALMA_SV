#---- Load modules
import numpy as np
import scipy
from scipy import stats
import matplotlib.pyplot as plt

#-------- Definition of functions
def GetAntName(msfile):
	tb.open(msfile+'/'+'ANTENNA')
	antList = tb.getcol("NAME")
	tb.close()
	return antList

def GetTimerecord(msfile, ant1, ant2, spwID, scanID):
	Out='ANTENNA1 == '+`ant1`+' && ANTENNA2 == '+`ant2` + ' && DATA_DESC_ID == '+`spwID`  + ' && SCAN_NUMBER == ' + `scanID`
	tb.open(msfile)
	antXantYspw = tb.query(Out)
	timeXY = antXantYspw.getcol('TIME')
	tb.close()
	return timeXY

def specBunch( spec, chBunch, timeBunch ):
	chNum, timeNum = spec.shape[0], spec.shape[1]
	chNum_b   = chNum/chBunch*chBunch
	timeNum_b = timeNum/timeBunch*timeBunch
	tmp = np.mean( spec[range(chNum_b),:].reshape(chNum/chBunch, chBunch, timeNum), axis=1 )
	return np.mean(tmp[:,range(timeNum_b)].reshape(chNum/chBunch, timeNum/timeBunch, timeBunch), axis=2)

def Ant2Bl(ant1, ant2):     # Antenna -> baseline index (without autocorr)
	antenna1 = max(ant1, ant2); antenna2 = min(ant1, ant2)
	return antenna1* (antenna1 - 1)/2 + antenna2

def Bl2Ant(baseline):
	ant1 = int(sqrt(2.0* baseline))
	while( (ant1* (ant1 + 1)/2 ) > baseline):
		ant1 -= 1
	return [(ant1+1), int(baseline - ant1*(ant1 + 1)/2)]

def Ant2Bla_RevLex(ant0, ant1, antNum):     # Reverse Lexical, with autcorr
	antenna0 = min(ant0, ant1); antenna1 = max(ant0, ant1)
	kernel = antNum* antenna0 - antenna0* (antenna0 - 1)/2
	return kernel + antenna1 - antenna0

def BlPhaseMatrix(num_ant):
	num_bl = num_ant* (num_ant - 1) / 2
	blphs_matrix = np.zeros((num_bl, (num_ant - 1)))
	for bl_index in range(num_bl):
		ants = Bl2Ant(bl_index)
		blphs_matrix[bl_index, (ants[0] - 1)] = 1
		if(ants[1] > 0):
			blphs_matrix[bl_index, (ants[1] - 1)] = -1
		#
	return blphs_matrix

def BlAmpMatrix(num_ant):
	num_bl = num_ant* (num_ant - 1) / 2
	blamp_matrix = np.zeros((num_bl, num_ant))
	for bl_index in range(num_bl):
		ants = Bl2Ant(bl_index)
		blamp_matrix[bl_index, ants[0]] = 1
		blamp_matrix[bl_index, ants[1]] = 1
	return blamp_matrix

def GetAzEl(msfile):
	Out = msfile + '/' + 'POINTING'
	tb.open(Out)
	Direction=tb.getcol("DIRECTION")
	tb.close()
	return Direction[0,0], Direction[1,0]

def GetAzOffset(msfile):
	Out = msfile + '/' + 'POINTING'
	tb.open(Out)
	Offset = tb.getcol("POINTING_OFFSET")
	Time   = tb.getcol("TIME")
	AntID  = tb.getcol("ANTENNA_ID")
	tb.close()
	return Time, AntID, Offset[:,0]*180*3600/pi      # Output in [arcsec]

def GetVisAllBL(msfile, spwID, scanID):
	antNum = len(GetAntName(msfile))
	corrNum= antNum* (antNum + 1)/2     # Number of correlations (with autocorr)
	blNum  = corrNum - antNum
	#
	Out='DATA_DESC_ID == '+`spwID` + ' && SCAN_NUMBER == ' + `scanID`
	tb.open(msfile)
	antXantYspw = tb.query(Out, sortlist='noduplicates ANTENNA1, ANTENNA2, TIME')
	timeXY = antXantYspw.getcol('TIME')
	timeNum = len(timeXY) / corrNum
	dataXY = antXantYspw.getcol('DATA')     # dataXY in array[pol, ch, baselinextime]
	tb.close()
	polNum, chNum = dataXY.shape[0], dataXY.shape[1]
	acorr_index = range(antNum)
	xcorr_index = range(blNum)
	for ant_index in range(antNum):
		acorr_index[ant_index] = Ant2Bla_RevLex(ant_index, ant_index, antNum)
	for bl_index in range(blNum):
		ant1, ant0 = Bl2Ant(bl_index)
		xcorr_index[bl_index] = Ant2Bla_RevLex(ant0, ant1, antNum)
	Pspec = dataXY.reshape(polNum, chNum, corrNum, timeNum)[:,:,acorr_index,:]
	Xspec = dataXY.reshape(polNum, chNum, corrNum, timeNum)[:,:,xcorr_index,:]
	return Pspec, Xspec

def clphase_solve(bl_phase, bl_error):
	blnum  =  len(bl_phase)
	antnum =  Bl2Ant(blnum)[0]
	weight =  np.divide(1.0, np.multiply(bl_error, bl_error))
	resid  =  np.zeros(blnum)
	niter  = 0
	correction = np.ones(antnum - 1)
	solution   = np.zeros(antnum - 1)
    #---- Initial solution
	for ant_index in range(1,antnum) :
		solution[ant_index - 1] = bl_phase[Ant2Bl(0, ant_index )]
	#
	while (np.dot(correction, correction) > 1e-6) and (niter < 5): # Loop for iterations
		#---- Residual Vector
		for bl_index in range(blnum):
			ants = Bl2Ant(bl_index)
			phs_diff  =  bl_phase[bl_index] - solution[ants[0] - 1]
			if  ants[1] != 0 :
				phs_diff = phs_diff + solution[ants[1]-1]   # Baselines without refant
			resid[bl_index] = atan2( sin(phs_diff), cos(phs_diff) )* weight[bl_index]   # 2pi ambibuity removal
		#
		#---- Partial Matrix
		p_matrix   = BlPhaseMatrix(antnum)
		ptwp       = np.dot(p_matrix.T, np.dot(np.diag(weight), p_matrix))
		ptwp_inv   = scipy.linalg.inv(ptwp)
		correction = np.dot(ptwp_inv,  np.dot(p_matrix.T, resid))
		solution   = np.add(solution, correction)
		niter      =  niter + 1
	#
	return np.append(0, solution), np.append(0, np.sqrt(np.diag(ptwp_inv)))

def clamp_solve(bl_amp, bl_error):
	blnum  =  len(bl_amp)
	antnum =  Bl2Ant(blnum)[0]      # Number of baselines and antennas
	weight =  np.divide(1.0, np.multiply(bl_error, bl_error))
	log_bl =  np.log(bl_amp)    # weights from standard error
	p_matrix = BlAmpMatrix(antnum)                                  # Elements of the matrix
	ptwp     = np.dot(p_matrix.T, np.dot(np.diag(weight), p_matrix))
	ptwp_inv = scipy.linalg.inv(ptwp)
	solution = np.exp(np.dot(ptwp_inv,  np.dot(p_matrix.T,  np.dot(np.diag(weight), log_bl))))
	return solution**2, 2.0* np.sqrt(solution* np.diag(ptwp_inv))

def phaseCal(vis, antPhase):
	blnum = len(vis)
	phaseCalVis = np.zeros(blnum, dtype=complex)
	for bl_index in range(blnum):
		ants = Bl2Ant(bl_index)
		phaseCalVis
		twidle = np.exp((0-1j)* (antPhase[ants[0]] - antPhase[ants[1]]))
		phaseCalVis[bl_index] = vis[bl_index] * twidle
	return phaseCalVis

def AzElExtract(antNum, AntID, timeXY, scanTime, Offset):
	timeIndex = range(len(timeXY))
	for index in range(len(timeXY)):
		timeIndex[index] = np.where( (abs(scanTime - timeXY[index]) < 0.001) & (AntID == 0) )[0] 
	return Offset[0, timeIndex], Offset[1, timeIndex]

def residGauss(param, x, y, yerr):
	return (y - param[0]*np.exp( -0.5* ((x - param[1])/ param[2])**2 )) / yerr

def residPoly(param, x, y):
	return (y - (param[0]* (x - param[1])**2 + param[2]))

def initGauss(x, y):
	Y = np.log(y)
	param = [-1.0, 0.0, 1.0]
	result = scipy.optimize.leastsq(residPoly, param, args=(x, Y))
	return result[0]

def fitGauss( x, y, yerr ):
	param = [np.max(y), 0.0, np.max(x)]
	result = scipy.optimize.leastsq( residGauss, param, args=(x, y, yerr), full_output=True)
	return result[0], np.sqrt([result[1][0,0], result[1][1,1], result[1][2,2]])

def circlePoints( x, y, radius ):
	angle = np.arange(-pi, (130/128)*pi, pi/128)
	return x + radius* np.cos(angle), y + radius* np.sin(angle)

def IFPointing( prefix, polID, spwID, scanID):
#-------- Procedures
	msfile = prefix + '.ms'
	antList = GetAntName(msfile)
	antnum = len(antList)
	antennas= range(antnum)
	blNum = antnum* (antnum - 1) / 2 
	spwNum = len(spwID)
	polNum = len(polID)
	#-------- Time Recores for the scan
	timeXY = GetTimerecord(msfile, 0, 1, spwID[0], scanID)
	timeNum = len(timeXY)
	#-------- Scan pattern
	scanTime, AntID, Offset = GetAzOffset(msfile)
	AzOff, ElOff = AzElExtract(antnum, AntID, timeXY, scanTime, Offset)
	posEnd = 1+ np.append(np.where( abs(AzOff[1:timeNum] - AzOff[0:(timeNum-1)]) + abs(ElOff[1:timeNum] - ElOff[0:(timeNum-1)]) > 0)[0], timeNum-1)
	posStart = np.append(0, np.where( abs(AzOff[1:timeNum] - AzOff[0:(timeNum-1)]) + abs(ElOff[1:timeNum] - ElOff[0:(timeNum-1)]) > 0)[0] + 1)
	posNum = len(posEnd)
	az_off = AzOff[posStart,0]; el_off = ElOff[posStart,0]
	#-------- Visibilities
	Vis    = np.zeros([polNum, spwNum, blNum, posNum], dtype=complex)
	Viserr = np.zeros([polNum, spwNum, blNum, posNum])

	phaseCalVis  = np.zeros([polNum, spwNum, blNum, posNum], dtype=complex)
	phaseAnt = np.zeros([polNum, spwNum, antNum, posNum])

	phaseAnterr = np.zeros([antnum, polNum, spwNum, posNum])
	Gain_ant = np.ones([antnum, posNum])
	Gain_err = np.ones([antnum, posNum])
	#-------- Coherent averaging of visibilities
	for spw_index in range(spwNum):
		Pspec, Xspec = GetVisAllBL(msfile, spwID[spw_index], scanID)
		Vis[:, spw_index, :]    = np.mean(np.mean(Xspec, axis=1).reshape(polNum, blNum, posNum, Xspec.shape[3]/posNum), axis=3)
		Viserr[:, spw_index, :] = np.std( np.angle(np.mean(Xspec, axis=1)).reshape(polNum, blNum, posNum, Xspec.shape[3]/posNum), axis=3)
		for pol_index in range(polNum):
			for pos_index in range(posNum):
				#-------- Antenna-based phase correction
				phaseAnt[pol_index, spw_index, :, pos_index], phaseAnterr[:, pol_index, spw_index, pos_index] = clphase_solve(np.angle(Vis[pol_index, spw_index, :, pos_index]), Viserr[pol_index, spw_index, :, pos_index] )
				phaseCalVis[pol_index, spw_index, :, pos_index] = phaseCal(Vis[pol_index, spw_index, :, pos_index], phaseAnt[pol_index, spw_index, :, pos_index])
	#-------- Coherent averaging acroww SPW and Polarization
	avgVis = np.mean(np.mean(phaseCalVis, axis=0), axis=0)
	avgViserr = np.mean(np.mean(Viserr, axis=0), axis=0)/sqrt(spwNum* polNum)
	#-------- Antenna-based Gain solution
	for pos_index in range(posNum):
		Gain_ant[:,pos_index], Gain_err[:,pos_index] = clamp_solve(abs(avgVis[:,pos_index]), avgViserr[:, pos_index])
	#
	"""
	"""
	offset = np.arange( 2.5* min(az_off), 2.5* max(az_off), 1 )
	fig = plt.figure(figsize = (8, 11))
	fig.text(0.4, 0.92, prefix + ' Scan' + `scanID`, fontsize=18)
	fig.text(0.2, 0.905, 'Az Scan', fontsize=10)
	fig.text(0.8, 0.905, 'El Scan', fontsize=10)
	"""
	"""
	#-------- Pointing determination
	Az_result = np.zeros([3, antnum]); Az_err = np.zeros([3, antnum])
	El_result = np.zeros([3, antnum]); El_err = np.zeros([3, antnum])
	for ant_index in range(antnum):
		#-------- Gaussian fit for Az scan
		azScan = np.where(el_off == 0)[0]
		Az_result[:,ant_index], Az_err[:,ant_index] = fitGauss(az_off[azScan], Gain_ant[ant_index,azScan], Gain_err[ant_index,azScan])
		#-------- Gaussian fit for El scan
		elScan = np.where(az_off == 0)[0]
		El_result[:,ant_index], El_err[:,ant_index] = fitGauss(el_off[elScan], Gain_ant[ant_index,elScan], Gain_err[ant_index,elScan])
		ylim = [0.0, 1.2* max( Gain_ant[ant_index] )]
		#-------- Az Plot
		plt.subplot(antnum, 2, 2* ant_index + 1)
		plt.plot(az_off[azScan], Gain_ant[ant_index, azScan], 'ro')
		plt.plot(offset, Az_result[0,ant_index]* np.exp(-0.5* ((offset - Az_result[1,ant_index])/Az_result[2,ant_index])**2))
		legend = '%s: %5.2f (%5.3f)' % (antList[ant_index], Az_result[1,ant_index], Az_err[1, ant_index]) 
		plt.text(0.95*min(offset), 0.75*ylim[1], legend, fontsize=6)
		plt.xticks(fontsize=6); plt.yticks(fontsize=6)
		plt.axis([min(offset),  max(offset), ylim[0], ylim[1]], fontsize=6)
		#
		#-------- El Plot
		plt.subplot(antnum, 2, 2* ant_index + 2)
		plt.plot(el_off[elScan], Gain_ant[ant_index, elScan], 'ro')
		plt.plot(offset, El_result[0,ant_index]* np.exp(-0.5* ((offset - El_result[1,ant_index])/El_result[2,ant_index])**2))
		legend = '%s: %5.2f (%5.3f)' % (antList[ant_index], El_result[1,ant_index], El_err[1, ant_index]) 
		plt.text(0.95*min(offset), 0.75*ylim[1], legend, fontsize=6)
		plt.xticks(fontsize=6); plt.yticks(fontsize=6)
		plt.axis([min(offset),  max(offset), ylim[0], ylim[1]], fontsize=6)
	#
	plt.savefig(prefix + '.' + `scanID` + '.GaussFit.pdf', form='pdf')
	plt.close()
	#return antList, Gain_ant, Gain_err
	return antList, Az_result, Az_err, El_result, El_err
