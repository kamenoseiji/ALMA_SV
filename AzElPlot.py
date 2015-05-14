#---- Script for Band-3 Astroholograpy Data
import matplotlib.pyplot as plt
execfile('/Users/kameno/ALMA_SV/Scripts/interferometry.py')
plotColor = ['black', 'red', 'blue', 'red']

msfile = prefix + '.ms'
timeAzEl, AntID, AZ, EL = GetAzEl(msfile)
antList = GetAntName(msfile)
antNum  = len(antList)
timeOffset = int(timeAzEl[0]) / 86400 * 86400

fig = plt.figure(figsize = (8,11))
fig.text(0.4, 0.92, prefix, fontsize=20)
#-------- Azimuth Plot
#plt.subplot( 211 )
xlim = [1e20, -1e20]
for sourceIndex in range(len(scan)):
	az_src = np.zeros(len(scan[sourceIndex])); el_src = np.zeros(len(scan[sourceIndex])); secz_src = np.zeros(len(scan[sourceIndex]))
	for index in range(len(scan[sourceIndex])):
		scanID = scan[sourceIndex][index]
		interval, timeXY = GetTimerecord( msfile, refant, refant, polID, spwID, scanID) 
		timeIndex = TimeExtract(antNum, refant, timeXY, timeAzEl)
		xlim = [min(np.append(timeAzEl[timeIndex]-timeOffset, xlim[0])), max(np.append(timeAzEl[timeIndex] - timeOffset, xlim[1]))]

		#---- Azimuth
		plt.subplot( 211 )
		ylim = [-180, 180]
		plt.plot(timeAzEl[timeIndex] - timeOffset, AZ[timeIndex]* 180 / pi, 'o', color=plotColor[sourceIndex])
		plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]])
		plt.xlabel('Time (second of day)'); plt.ylabel('Azimuth [deg]')
		az_src[index] = np.mean(AZ[timeIndex])*180/pi
		el_src[index] = np.mean(EL[timeIndex])*180/pi
		secz_src[index] = np.mean(1.0 / np.sin(EL[timeIndex]))
		#print "%d  %f  %f %f" % (scanID, az_src[index], el_src[index], secz_src[index])

		#---- Elevation
		plt.subplot( 212 )
		ylim = [20, 85]
		plt.plot(timeAzEl[timeIndex] - timeOffset, EL[timeIndex]* 180 / pi, 'o', color= plotColor[sourceIndex])
		plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]])
		plt.xlabel('Time (second of day)'); plt.ylabel('Elevation [deg]')
	#
	if sourceIndex == 0:
		refAz = np.median(az_src)
		refEl = np.median(el_src)
		refSecZ = np.median(secz_src)
	#
	print "Source#%d  AZ= %6.2f  EL= %6.2f SECZ= %8.6f dAZ=%6.2f dEL=%6.2f dSECZ=%8.6f" % (sourceIndex, np.median(az_src), np.median(el_src), np.median(secz_src), abs(np.median(az_src) - refAz), abs(np.median(el_src) - refEl), abs(np.median(secz_src) - refSecZ))
plt.savefig(prefix + 'AzEl.pdf')
plt.close()

