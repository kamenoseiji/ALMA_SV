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
msBLC = prefixBLC + '.ms'
msACA = prefixACA + '.ms'
antNum  = len(antList)
polNum  = len(pol)
scanNum  = len(scan)
spwNum  = len(spw_ACA)
logfile = open(prefixBLC + '-' + `spw_ACA[0]` + '-BBLOG.log', 'w')
#-------- Tsys spectrum for specified antennas
chNum, chWid, Freq = GetChNum(msBLC, BLCTsysSPW[0])
BLCchWid = np.median(chWid)* 1.0e-6
wavelength = constants.c / np.median(Freq)
FWHM = 1.13* 180.0* 3600.0* wavelength / (GetAntD(antList[0])* pi) # Gaussian beam for 12-m antenna, in unit of arcsec 
BLCFreq = Freq* 1.0e-9  # [GHz]
BLCTrxList, BLCTsysList = tsysSpec( prefixBLC, TsysScan, BLCTsysSPW )   # List of Trx[antNum, polNum, chNum], Tsys[antNum, polNum, chNum]
antListInBLC = antIndex(prefixBLC, antList)
antListInACA = antIndex(prefixACA, antList)
chNum, chWid, Freq = GetChNum(msACA, ACATsysSPW[0])
ACAchWid = np.median(chWid)* 1.0e-6
ACAFreq = Freq* 1.0e-9  # [GHz]
ACATrxList, ACATsysList = tsysSpec( prefixACA, TsysScan, ACATsysSPW )   # List of Trx[antNum, polNum, chNum], Tsys[antNum, polNum, chNum]
#-------- Az, El scan pattern
scanTime, AntID, az, el = GetAzEl(msBLC)
interval, timeStamp = GetTimerecord(msBLC, 0, 0, 0, spw_BLC[0], scan[0])
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
OffIndex = []; OnIndex = []
for ss_index in range(0, SSnum, 2):
    OffIndex = OffIndex + range(ST_index[ss_index], ED_index[ss_index])
#
for ss_index in range(1, SSnum, 2):
    OnIndex = OnIndex + range(ST_index[ss_index], ED_index[ss_index])
#
refRA, refDEC = mean(ScanRA[OnIndex]), mean(ScanDEC[OnIndex])
ScanRA  = 202164.8*(ScanRA  - refRA)*cos(refDEC); ScanDEC = 202164.8*(ScanDEC - refDEC)
BLClineCH = np.where((BLCFreq >= BLClineRange[0]) & (BLCFreq <= BLClineRange[1]))[0]
BLCfreeCH = np.where(((BLCFreq >= BLClineFree[0]) & (BLCFreq <= BLClineFree[1])) | ((BLCFreq >= BLClineFree[2]) & (BLCFreq <= BLClineFree[3])))[0]
ACAlineCH = np.where((ACAFreq >= ACAlineRange[0]) & (ACAFreq <= ACAlineRange[1]))[0]
ACAfreeCH = np.where(((ACAFreq >= ACAlineFree[0]) & (ACAFreq <= ACAlineFree[1])) | ((ACAFreq >= ACAlineFree[2]) & (ACAFreq <= ACAlineFree[3])))[0]
#
GridWidth = max(ScanRA[OnIndex])
xi, yi = np.mgrid[ -GridWidth:GridWidth:128j, -GridWidth:GridWidth:128j]
#-------- SourceMapping
for scan_index in range(scanNum):
    print 'SCAN=%d' % scan[scan_index]
    text_sd = 'ANT POL SPW  TaBLC  TaACA Error%'
    print text_sd; logfile.write(text_sd + '\n')
    for spw_index in range(spwNum):
        text_sd = '%s Scan=%d SPW=%d' % (prefixBLC, scan[scan_index], spw_BLC[spw_index])
        for ant_index in range(antNum):
            fig  = plt.figure(figsize = (8, 11))
            for pol_index in range(polNum):
                timeBLC, dataBLC = GetVisibility(msBLC, antListInBLC[ant_index], antListInBLC[ant_index], pol_index, spw_BLC[spw_index], scan[scan_index])
                timeACA, dataACA = GetVisibility(msACA, antListInACA[ant_index], antListInACA[ant_index], pol_index, spw_ACA[spw_index], scan[scan_index])
                #-------- BP calib
                BPBLC = np.mean(dataBLC.real[:, OffIndex], axis=1); BPCaledBLC = (dataBLC.real.T / BPBLC).T
                BPACA = np.mean(dataACA.real[:, OffIndex], axis=1); BPCaledACA = (dataACA.real.T / BPACA).T
                #--------
                skyBLC = np.mean( BPCaledBLC[BLCfreeCH,:].real, axis=0 )
                skyACA = np.mean( BPCaledACA[ACAfreeCH,:].real, axis=0 )
                plotBLC = BLCchWid* mean(BLCTsysList[spw_index][antListInBLC[ant_index], pol_index])* np.sum(BPCaledBLC[BLClineCH,:].real - skyBLC, axis=0) / skyBLC
                plotACA = ACAchWid* mean(ACATsysList[spw_index][antListInACA[ant_index], pol_index])* np.sum(BPCaledACA[ACAlineCH,:].real - skyACA, axis=0) / skyACA
                TABLC = GridData( plotBLC[OnIndex], ScanRA[OnIndex], ScanDEC[OnIndex], xi.reshape(xi.size), yi.reshape(xi.size), FWHM/8 ).reshape(len(xi), len(xi))
                TAACA = GridData( plotACA[OnIndex], ScanRA[OnIndex], ScanDEC[OnIndex], xi.reshape(xi.size), yi.reshape(xi.size), FWHM/8 ).reshape(len(xi), len(xi))
                #
                #-------- Plot BLC Map
                plt.subplot(3, 2, pol_index + 1, aspect=1)
                plt.contourf(xi, yi, TABLC, cntrRange); plt.colorbar()
                plt.axis([GridWidth, -GridWidth, -GridWidth, GridWidth])
                GaussBLC = simple2DGaussFit(plotBLC[OnIndex], ScanRA[OnIndex], ScanDEC[OnIndex] )
                text_sd = '%s: Ta* = %5.1f K MHz' % (corrLabelBLC, GaussBLC[0])
                plt.text(0, 0.8*GridWidth, text_sd, size='x-small', color='black')
                plt.title(srcName + ' ' + corrLabelBLC + ' Pol=' + pol[pol_index])
                #-------- Plot ACA Map
                plt.subplot(3, 2, pol_index + 3, aspect=1)
                plt.contourf(xi, yi, TAACA, cntrRange); plt.colorbar()
                plt.axis([GridWidth, -GridWidth, -GridWidth, GridWidth])
                GaussACA = simple2DGaussFit(plotACA[OnIndex], ScanRA[OnIndex], ScanDEC[OnIndex] )
                text_sd = '%s: Ta* = %5.1f K MHz' % (corrLabelACA, GaussACA[0])
                plt.text(0, 0.8*GridWidth, text_sd, size='x-small', color='black')
                plt.title(srcName + ' ' + corrLabelACA + ' Pol=' + pol[pol_index])
                #
                #-------- Plot Ratiop
                plt.subplot(3, 2, pol_index + 5, aspect=1)
                plt.plot(plotBLC, plotACA, '.')
                plt.plot( np.array([min(plotBLC), max(plotBLC)]), np.array([min(plotBLC), max(plotBLC)]) )
                slope = np.dot(plotBLC, plotACA) / np.dot(plotBLC, plotBLC)
                text_sd = 'Slope = %5.3f' % ( slope )
                plt.text(0, 0.8*max(plotBLC), text_sd, size='small', color='black')
                #
                #text_sd = '%s %s %d  %5.1f  %5.1f  %5.2f' % (antList[ant_index], pol[pol_index], spw_BLC[spw_index], GaussBLC[0], GaussACA[0], 100.0*(GaussACA[0]/GaussBLC[0] - 1.0) ); print text_sd
                text_sd = '%s %s %d  %5.1f  %5.1f  %5.2f' % (antList[ant_index], pol[pol_index], spw_BLC[spw_index], GaussBLC[0], GaussACA[0], 100.0*(slope - 1.0) ); print text_sd
                logfile.write(text_sd + '\n')
            #
            plt.suptitle(prefixBLC + ' ' + antList[ant_index] + ' Spw=' + `spw_BLC[spw_index]`)
            plt.savefig(prefixBLC + '.' + antList[ant_index] + '.SPW' + `spw_BLC[spw_index]` + '.pdf', form='pdf')
            plt.close(fig)
        #
    #
#
logfile.close()
