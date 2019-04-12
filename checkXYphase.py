execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Grid.py')
from matplotlib.backends.backend_pdf import PdfPages
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
#----------------------------------------- Check SPW intents
msfile = wd + prefix + '.ms'
execfile(SCR_DIR + 'TsysCal.py')
antList = GetAntName(msfile)
antNum = len(antList)
blNum = antNum* (antNum - 1) / 2
msmd.open(msfile)
atmSPWs = list(set(msmd.tdmspws()) & set(msmd.spwsforintent("CALIBRATE_ATMOSPHERE*"))); atmSPWs.sort()
bpSPWs  = list(set( msmd.tdmspws()) & set(msmd.spwsforintent("CALIBRATE_DELAY*"))); bpSPWs.sort()
if len(bpSPWs) == 0: bpSPWs = list(set( msmd.tdmspws()) & set(msmd.spwsforintent("CALIBRATE_PHASE*"))); bpSPWs.sort()
atmspwNames, bpspwNames = msmd.namesforspws(atmSPWs), msmd.namesforspws(bpSPWs)
atmBandNames, atmPattern = [], r'RB_..'
for spwName in atmspwNames : atmBandNames = atmBandNames + re.findall(atmPattern, spwName)
UniqBands = unique(atmBandNames).tolist(); NumBands = len(UniqBands)
atmspwLists, bpspwLists, atmscanLists, bpscanLists, BandPA = [], [], [], [], []
for band_index in range(NumBands):
    BandPA = BandPA + [(BANDPA[int(UniqBands[band_index][3:5])] + 90.0)*pi/180.0]
    atmspwLists = atmspwLists + [np.array(atmSPWs)[indexList(np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()]
    bpspwLists  = bpspwLists  + [np.array(bpSPWs)[indexList( np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()]
    atmscanLists= atmscanLists+ [list(set(msmd.scansforintent("CALIBRATE_ATMOSPHERE*")) & set(msmd.scansforspw(atmspwLists[band_index][0])))]
    bpscanLists = bpscanLists + [msmd.scansforspw(bpspwLists[band_index][0]).tolist()]
    print ' ',
    print UniqBands[band_index] + ': atmSPW=' + `atmspwLists[band_index]` + ' bpSPW=' + `bpspwLists[band_index]`
#
#-------- Check source list
print '---Checking source list'
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList); numSource = len(sourceList)
for source in sourceList: print source,
print ''
#-------- Check Scans of BandPass, EQualization, and FluxScaling
polNum = msmd.ncorrforpol(msmd.polidfordatadesc(bpspwLists[0][0]))
PolScans = msmd.scansforintent("CALIBRATE_POLARIZATION#ON_SOURCE")
print '---SPWs and Scans for each receiver band'
msmd.done()
for band_index in range(NumBands):
    bandName = UniqBands[band_index]; bandID = int(UniqBands[band_index][3:5])
    msmd.open(msfile)
    #-------- Check Calibrators
    PolScan = np.array(bpscanLists[band_index])[indexList( PolScans, np.array(bpscanLists[band_index]))]
    onsourceScans = PolScan.tolist()
    scanNum = len(onsourceScans)
    fp = open('CalQU.data')
    lines = fp.readlines()
    fp.close()
    for eachLine in lines:
        catalogStokesI[eachLine.split()[0]] = float(eachLine.split()[1])
        catalogStokesQ[eachLine.split()[0]] = float(eachLine.split()[2])
        catalogStokesU[eachLine.split()[0]] = float(eachLine.split()[3])
    #
    #-------- Check AZEL
    azelTime, AntID, AZ, EL = GetAzEl(msfile)
    azelTime_index = np.where( AntID == 0 )[0].tolist()
    azel = np.r_[AZ[azelTime_index], EL[azelTime_index]].reshape(2, len(azelTime_index))
    refTime, OnAZ, OnEL, OnPA, BPquality, EQquality, sourceIDscan, FLscore = [], [], [], [], [], [], [], np.zeros(scanNum)
    for scan_index in range(scanNum):
        sourceIDscan.append( msmd.sourceidforfield(msmd.fieldsforscan(onsourceScans[scan_index])[0]))
        interval, timeStamp = GetTimerecord(msfile, 0, 0, bpspwLists[band_index][0], onsourceScans[scan_index])
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, 0, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; dPA = np.std(np.sin(PA)) #dPA = abs(np.sin(max(PA) - min(PA)))
        OnAZ.append(np.median(AzScan)); OnEL.append(np.median(ElScan)); OnPA.append(np.median(PA))
        refTime = refTime + [np.median(timeStamp)]
    #
    EQScan = BPScan
    BPEL = OnEL[onsourceScans.index(BPScan)]; EQEL = BPEL
    BPcal = sourceList[sourceIDscan[onsourceScans.index(BPScan)]]; EQcal = BPcal
    BPcalText = 'Use %s [EL = %4.1f deg] as Bandpass Calibrator' % (BPcal, 180.0* OnEL[onsourceScans.index(BPScan)]/np.pi); print BPcalText
    EQcalText = 'Use %s [EL = %4.1f deg] as Gain Equalizer' % (EQcal, 180.0* OnEL[onsourceScans.index(EQScan)]/np.pi); print EQcalText
    atmspw = atmspwLists[band_index]; spwNum = len(atmspw)
    scnspw = bpspwLists[band_index]; scnspwNum = len(scnspw)
    pPol, cPol = [0,3], [1,2]
    execfile(SCR_DIR + 'aprioriStokes.py')
    del(BPcal)
    del(EQcal)
#
del(spwList)
