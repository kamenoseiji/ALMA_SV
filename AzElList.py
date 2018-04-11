#---- Script for Band-3 Astroholograpy Data
execfile(SCR_DIR + 'interferometry.py')

scanNum = len(scan)
msfile = wd + prefix + '.ms'
timeAzEl, AntID, AZ, EL = GetAzEl(msfile)

AZEL = np.zeros([scanNum, 2])	# (AZ, EL) x scan num
#-------- Azimuth Plot
for scan_index in range(scanNum):
	interval, keyTime = GetTimerecord( msfile, 0, 0, spw[0], scan[scan_index]) 
	rec_index = TimeExtract(AntID, timeAzEl, keyTime)
	AZEL[scan_index, 0], AZEL[scan_index,1] = np.median(AZ[rec_index]), np.median(EL[rec_index])
#
np.save(prefix[0] + '.AzEl.npy', AZEL)
