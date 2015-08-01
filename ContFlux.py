#-------- Script to compare ACA power and BB power
execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
from scipy.interpolate import griddata
ALMA_lat = -23.029/180.0*pi
ALMA_long= -67.755/180.0*pi
#-------- Script to compare ACA power and BB power
def BB_filter(timeBB, dataBB):
    index = np.where( abs(dataBB[0]) > 0.0)[0]
    return timeBB[index], abs(dataBB[0,index])
#
#-------- Scan Time in ACA data
def timeMatch( timeBB, timeACA ):
    return np.where( (timeACA < np.max(timeBB)) & (timeACA > np.min(timeBB)) )[0] 
#

def timeRange(refTime, timeBB, timeWidth):
    return np.where( abs(timeBB - refTime) < 0.5*timeWidth)[0]
#
def antIndex(prefix, antList):
    msfile = prefix + '.ms'
    antListInMS = GetAntName(msfile)
    antIndexInMS = []
    for ant_index in range(len(antList)):
        antIndexInMS.append( antListInMS.tolist().index(antList[ant_index]))
    #
    return antIndexInMS
#
def scanPattern(scanTime, scanGap):
    gap = np.where( diff(scanTime) > scanGap )[0]
    scanNum = len(gap) + 1
    ST_index = append(0, gap+1)
    ED_index = append(gap, len(scanTime)-1)
    return ST_index, ED_index
#
def tsysSpec(prefix, TsysScan, TsysSPW):
    msfile = prefix + '.ms'
    antList = GetAntName(msfile)
    antNum  = len(antList)
    polNum  = len(pol)
    spwNum  = len(TsysSPW)
    TrxList = []
    TsysList = []
    for spw_index in range(spwNum):
        chNum, chWid, freq = GetChNum(msfile, TsysSPW[spw_index]); Tsysfreq = freq* 1.0e-9 # GHz
        TrxSpec = np.zeros([antNum, polNum, chNum]); TsysSpec = np.zeros([antNum, polNum, chNum])
        #
        #-------- Get Physical Temperature of loads
        for ant_index in range(antNum):
            tempAmb, tempHot = GetLoadTemp(msfile, ant_index, TsysSPW[spw_index])
            #
            for pol_index in range(polNum):
                timeXY, dataXY = GetVisibility(msfile, ant_index, ant_index, pol_index, TsysSPW[spw_index], TsysScan)
                #
                #-------- Time range of Sky/Amb/Hot
                edge = np.where( diff(timeXY) > 1.0 )[0]
                skyRange = range(0, edge[0])
                ambRange = range(edge[0]+1, edge[1])
                hotRange = range(edge[1]+1, len(timeXY))
                #
                #-------- Calc. Tsys Spectrum
                Psky, Pamb, Phot = np.mean(dataXY[:,skyRange].real, 1), np.mean(dataXY[:,ambRange].real, 1), np.mean(dataXY[:,hotRange].real, 1)
                TrxSpec[ant_index, pol_index]  = (tempHot* Pamb - Phot* tempAmb) / (Phot - Pamb)
                TsysSpec[ant_index, pol_index] = (Psky* tempAmb) / (Pamb - Psky)
                print '%s SPW=%d pol=%s: Trx=%5.1f Tsys=%5.1f' % (antList[ant_index], TsysSPW[spw_index], pol[pol_index], np.median(TrxSpec[ant_index, pol_index]), np.median(TsysSpec[ant_index, pol_index]))
            #
        #
        TrxList.append(TrxSpec)
        TsysList.append(TsysSpec)
    #
    return TrxList, TsysList
#
#-------- GridPoint
def GridPoint( value, samp_x, samp_y, point_x, point_y, kernel ):
    #---- Check NaN and replace with 0
    nan_index = np.where( value != value )[0]
    value[nan_index] = 0.0
    #---- Distance from gridding points
    dist_sq = (samp_x - point_x)**2 + (samp_y - point_y)**2
    dist_thresh = 9.0 * kernel**2
    index = np.where( dist_sq < dist_thresh)[0]
    wt = exp( -0.5* dist_sq[index] / kernel**2 )
    nan_index = np.where( value[index] != value[index] )[0]
    wt[nan_index] = 0.0
    sumWt = np.sum(wt)
    if sumWt < 1.0e-3:
        return 0.0
    #
    return np.sum(value[index]* wt) / sumWt
#
#-------- GridData
def GridData( value, samp_x, samp_y, grid_x, grid_y, kernel ):
    gridNum = len(grid_x)
    results = np.zeros(gridNum)
    for index in range(gridNum):
        results[index] = GridPoint( value, samp_x, samp_y, grid_x[index], grid_y[index], kernel)
    #
    return results
#
msfile = prefix + '.ms'
antList = GetAntName(msfile)
antNum  = len(antList)
polNum  = len(pol)
scanNum  = len(scan)
spwNum  = len(spw_ACA)
logfile = open(prefix + '_BBLOG.log', 'w')
#-------- Tsys spectrum for specified antennas
chNum, chWid, Freq = GetChNum(msfile, spw_ACA[0])
wavelength = constants.c / np.median(Freq)
FWHM = 1.13* 180.0* 3600.0* wavelength / (12.0* pi) # Gaussian beam for 12-m antenna, in unit of arcsec 
Freq = Freq* 1.0e-9  # [GHz]
TrxList, TsysList = tsysSpec( prefix, TsysScan, TsysSPW )   # List of Trx[antNum, polNum, chNum], Tsys[antNum, polNum, chNum]
antListInACA = antIndex(prefix, antList)
chRange = range( int(chNum* 0.04), int(chNum* 0.98))
#-------- Az, El scan pattern
scanTime, AntID, az, el = GetAzEl(msfile)
interval, timeStamp = GetTimerecord(msfile, 0, 0, 0, spw_ACA[0], scan[0])
timeNum = len(timeStamp)
ScanAz = np.zeros([timeNum]); ScanEl = np.zeros([timeNum])
for time_index in range(timeNum):
    ScanAz[time_index], ScanEl[time_index] = AzElMatch( timeStamp[time_index], scanTime, np.median(interval), az, el)
#
LST = mjd2gmst( timeStamp/86400.0, 30.0) + ALMA_long
ScanRA, ScanDEC = azel2radec( ScanAz, ScanEl, LST, ALMA_lat)
#-------- SubScan Pattern
ST_index, ED_index = scanPattern(timeStamp, 1.5* np.median(interval))
SSnum = len(ST_index)
OffIndex = []; OnIndex = []; knots = []
for ss_index in range(SSnum):
    knots.append( 0.5*( timeStamp[ST_index[ss_index]] + timeStamp[ST_index[ss_index] + 1]) )
    knots.append( 0.5*( timeStamp[ED_index[ss_index]] + timeStamp[ED_index[ss_index] - 1]) )
#
for ss_index in range(0, SSnum, 2):
    OffIndex = OffIndex + range(ST_index[ss_index], ED_index[ss_index])
#
for ss_index in range(1, SSnum, 2):
    OnIndex = OnIndex + range(ST_index[ss_index], ED_index[ss_index])
#
refRA, refDEC = mean(ScanRA[OnIndex]), mean(ScanDEC[OnIndex])
ScanRA  = 202164.8*(ScanRA  - refRA)*cos(refDEC); ScanDEC = 202164.8*(ScanDEC - refDEC)
OffIndex = np.where( ScanRA**2 + ScanDEC**2 > (1.8*FWHM)**2)[0]
interp1d( timeStamp[OffIndex], ScanAz[OffIndex], kind='cubic')
#
#-------- SourceMapping
for scan_index in range(scanNum):
    print 'SCAN=%d' % scan[scan_index]
    text_sd = 'ANT POL SPW peak-peak   TaBB  Ta%s Err%%' % (corrLabel)
    print text_sd; logfile.write(text_sd + '\n')
    for spw_index in range(spwNum):
        text_sd = '%s Scan=%d SPW=%d' % (prefix, scan[scan_index], spw_ACA[spw_index])
        for ant_index in range(antNum):
            for pol_index in range(polNum):
                fig  = plt.figure(figsize = (8, 11))
                timeACA, dataACA = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_ACA[spw_index], scan[scan_index])
                #-------- BB detector
                timeBB, dataBB = GetVisibility(msfile, ant_index, ant_index, pol_index, spw_BB[spw_index], scan[scan_index])
                timeMatchedBB = np.zeros(timeNum)
                plotBB = abs(dataBB[0])/ mean(abs(dataBB[0]))
                for index in range(timeNum):
                    indexRange = timeRange(timeACA[index], timeBB, 0.144)
                    if( len(indexRange) > 0):
                        timeMatchedBB[index] = np.mean(plotBB[indexRange])
                    #
                #
                indexAvail = np.where( timeMatchedBB > 0.0)[0]
                plotACA = np.mean( dataACA[chRange].real, axis=0) / mean( dataACA[chRange].real)
                skyBB  = LSQUnivariateSpline(timeStamp[OffIndex], timeMatchedBB[OffIndex], t=knots)
                skyACA = LSQUnivariateSpline(timeStamp[OffIndex], plotACA[OffIndex], t=knots)
                ACA_BB_ratio = plotACA[indexAvail]/timeMatchedBB[indexAvail]
                #
                #-------- Plot BB and ACA ratio
                plt.subplot(2, 2, 1)
                plt.plot(timeACA[indexAvail], plotACA[indexAvail], ls='steps-mid', label=corrLabel)
                plt.plot(timeACA[indexAvail], timeMatchedBB[indexAvail], '.', label="BB")
                plt.plot(timeACA[indexAvail], ACA_BB_ratio, ls='steps-mid', label="Ratio")
                plt.text(0.8*np.max(timeACA)+0.2*np.min(timeACA), 0.9*np.max(plotBB) + 0.1*np.min(plotBB), 'Pol=' + pol[pol_index], size='x-small')
                plt.title('Time-Power')
                text_sd = '%s %s %d %5.3f %5.3f ' % (antList[ant_index], pol[pol_index], spw_ACA[spw_index], np.max(ACA_BB_ratio-1.0)*100.0, np.min(ACA_BB_ratio-1.0)*100.0)
                print text_sd,
                logfile.write(text_sd)
                plt.xlabel('Time', fontsize=9)
                if pol_index == 0 :
                    plt.ylabel('Scaled Power [a.u.]', fontsize=9)
                #
                plt.legend(loc = 'lower left', prop={'size':9})
                #-------- Plot BB and ACA Comparison
                plt.subplot(2, 2, 2, aspect=1)
                plt.plot(timeMatchedBB[indexAvail], plotACA[indexAvail], '.')
                plt.plot( np.array([min(timeMatchedBB[indexAvail]), max(timeMatchedBB[indexAvail])]), np.array([min(timeMatchedBB[indexAvail]), max(timeMatchedBB[indexAvail])]) )
                plt.xlabel('BB Power'); plt.ylabel('ACA power')

                #-------- Plot Uranus BB Map
                plt.subplot(2, 2, 3, aspect=1)
                GridWidth = max(ScanRA[OnIndex])
                xi, yi = np.mgrid[ -GridWidth:GridWidth:128j, -GridWidth:GridWidth:128j]
                TABB  = mean(TsysList[spw_index][ant_index, pol_index])* GridData( timeMatchedBB[OnIndex] - skyBB(timeStamp[OnIndex]), ScanRA[OnIndex], ScanDEC[OnIndex], xi.reshape(xi.size), yi.reshape(xi.size), 3 ).reshape(len(xi), len(xi))
                TAACA = mean(TsysList[spw_index][ant_index, pol_index])* GridData( plotACA[OnIndex] - skyACA(timeStamp[OnIndex]), ScanRA[OnIndex], ScanDEC[OnIndex], xi.reshape(xi.size), yi.reshape(xi.size), 3 ).reshape(len(xi), len(xi))
                GaussBB = simple2DGaussFit((timeMatchedBB[OnIndex] - skyBB(timeStamp[OnIndex])), ScanRA[OnIndex], ScanDEC[OnIndex] )
                text_sd = 'BB: Ta* = %5.3f K' % (GaussBB[0]*  mean(TsysList[spw_index][ant_index, pol_index]))
                #plt.contourf(xi, yi, TABB, np.linspace(-3, 57, 16, endpoint=True)); plt.colorbar()
                plt.contourf(xi, yi, TABB, np.linspace(-0.1*imageMax, imageMax, 12, endpoint=True)); plt.colorbar()
                plt.text(0, 0.8*GridWidth, text_sd, size='x-small', color='yellow')
                #plt.title('Saturn BB Pol=' + pol[pol_index])
                plt.title(srcName +' BB Pol=' + pol[pol_index])
                text_sd = ' %5.3f ' % (GaussBB[0]*  mean(TsysList[spw_index][ant_index, pol_index])); print text_sd,
                logfile.write(text_sd)

                #-------- Plot Uranus ACA Map
                plt.subplot(2, 2, 4, aspect=1)
                #plt.contourf(xi, yi, TAACA, np.linspace(-2, 46, 17, endpoint=True)); plt.colorbar()
                plt.contourf(xi, yi, TAACA, np.linspace(-0.1*imageMax, imageMax, 12, endpoint=True)); plt.colorbar()
                GaussACA = simple2DGaussFit((plotACA[OnIndex] - skyACA(timeStamp[OnIndex])), ScanRA[OnIndex], ScanDEC[OnIndex] )
                text_sd = '%s: Ta* = %5.3f K' % (corrLabel, GaussACA[0]*  mean(TsysList[spw_index][ant_index, pol_index]))
                plt.text(0, 0.8*GridWidth, text_sd, size='x-small', color='yellow')
                #plt.title('Saturn ' + corrLabel + ' Pol=' + pol[pol_index])
                plt.title(srcName + ' ' + corrLabel + ' Pol=' + pol[pol_index])
                text_sd = ' %5.3f ' % (GaussACA[0]*  mean(TsysList[spw_index][ant_index, pol_index])); print text_sd,
                logfile.write(text_sd)
                text_sd = ' %6.3f ' % (100.0 * (GaussACA[0]/GaussBB[0] - 1.0)); print text_sd
                logfile.write(text_sd + '\n')
                plt.suptitle(prefix + ' ' + antList[ant_index] + ' Pol=' + pol[pol_index] + ' Spw=' + `spw_ACA[spw_index]`)
                plt.savefig( 'TP_' +prefix + '_' + antList[ant_index] + '_Pol' + pol[pol_index] + '_SPW' + `spw_ACA[spw_index]` + '.pdf', form='pdf')
                plt.close(fig)
            #
        #
    #
#
logfile.close()
