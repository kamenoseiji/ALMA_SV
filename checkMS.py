execfile(SCR_DIR + 'interferometry.py')
ALMA_lat = -0.4019318734    # Latitude = -23.029 deg
RADDEG   = 180.0 / math.pi
#
def antRefScan( msfile ):
    antList = GetAntName(msfile)
    antNum = len(antList)
    scanRange = np.zeros(antNum)
    Time, AntID, Offset = GetAzOffset(msfile)
    for ant_index in range(antNum):
        time_index = np.where( AntID == ant_index )[0]
        scanRange[ant_index] = max( Offset[0, time_index] ) - min( Offset[0, time_index] )
    #
    refAntIndex  = np.where( scanRange == 0.0 )[0]
    scanAntIndex = np.where( scanRange >  0.0 )[0]
    return refAntIndex.tolist(), scanAntIndex.tolist()
#
msNum = len(prefix)
#-------- Listobs
for file_index in range(msNum):
    vis = wd + prefix[file_index] + '.ms'
    listfile = prefix[file_index] + '.listobs'
    verbose = True
    listobs()
#
#-------- ScanPattern
for file_index in range(msNum):
    msfile = wd + prefix[file_index] + '.ms'
    #-------- Antenna List
    antList = GetAntName(msfile)
    antNum = len(antList)
    #
    #-------- Reference and Scan Antennas
    refAnt, scanAnt = antRefScan(msfile)
    print '-------- Reference Antennas (' + `len(refAnt)` + ')----'
    for ant_index in refAnt:
        sd_text = ' %02d %s' % (ant_index, antList[ant_index]); print sd_text
    #
    print '-------- Scanning Antennas (' + `len(scanAnt)` + ')----'
    for ant_index in scanAnt:
        sd_text = ' %02d %s' % (ant_index, antList[ant_index]); print sd_text
    #
    #-------- AZEL
    Time, AntID, Az, El = GetAzEl(msfile)
    text_sd = 'Az %5.1f - %5.1f  El %5.1f - %5.1f' % (np.min(Az)*RADDEG, np.max(Az)*RADDEG, np.min(El)*RADDEG, np.max(El)*RADDEG)
    print text_sd
    AzElfig = plt.figure(figsize = (8,11))
    AzElfig.text(0.45, 0.95, prefix[file_index])
    for ant_index in range(antNum):
        time_index = np.where( AntID == ant_index )[0].tolist()
        AzPL = AzElfig.add_subplot( antNum, 2,  2* ant_index + 1 )
        AzPL.plot(Time[time_index], Az[time_index]* RADDEG, '.', ms = 0.2, label = antList[ant_index])
        AzPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        ElPL = AzElfig.add_subplot( antNum, 2,  2* ant_index + 2 )
        ElPL.plot(Time[time_index], El[time_index]* RADDEG, '.', ms = 0.2, label = antList[ant_index])
        ElPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    #
    AzElfig.text(0.05, 0.45, 'Az / El Offset [arcsec]', rotation=90)
    AzElfig.text(0.45, 0.05, 'MJD [sec]')
    AzElfig.savefig(prefix[file_index] + '_AzEl.pdf')
    #-------- Scan
    Time, AntID, Offset = GetAzOffset(vis)
    Scanfig = plt.figure(figsize = (8,11))
    Scanfig.text(0.45, 0.95, prefix[file_index])
    for ant_index in range(antNum):
        time_index = np.where( AntID == ant_index )[0].tolist()
        num_Xpanel = int( sqrt(antNum) ); num_Ypanel = int(antNum / num_Xpanel)
        if num_Xpanel * num_Ypanel < antNum:
            num_Ypanel = num_Ypanel + 1
        #
        ScanPL = Scanfig.add_subplot( num_Ypanel, num_Xpanel,  ant_index, aspect=1 )
        ScanPL.plot( Offset[0, time_index], Offset[1, time_index], '.', ms = 0.2, label = antList[ant_index])
        ScanPL.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    #
    ScanPL.set_xlabel('Az Offset [arcsec]')
    Scanfig.text(0.05, 0.45, 'El Offset [arcsec]', rotation=90)
    Scanfig.savefig(prefix[file_index] + '_scanPattern.pdf')
    plt.close('all')
#
