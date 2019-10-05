#---- Script for Band-3 Astroholograpy Data
from scipy import stats
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
#-------- Procedures
msmd.open(msfile)
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile)
antNum = len(antList)
#-------- Check Scans
if 'scanList' not in locals():
    print '---Checking Scans with CALIBRATE_PHASE'
    scanList = msmd.scansforintent("CALIBRATE_PHASE#ON_SOURCE").tolist()
#
scanNum = len(scanList)
#-------- Check SPWs
if 'SPWList' not in locals():
    print '---Checking spectral windows'
    spwList = list( (set(msmd.tdmspws()) | set(msmd.fdmspws())) & set(msmd.spwsforscan(scanList[0])) )
#
spwNum = len(spwList)
msmd.close()
msmd.done()
#-------- Time Records
timeList = []
for scan_index in range(scanNum):
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], scanList[scan_index])
    timeList = timeList + timeStamp.tolist()
#
timeNum = len(timeList)
#-------- SPW records
AC = []
chNumList, freqList = [], []
for spw_index in range(spwNum):
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index])
    chNumList = chNumList + [chNum]
    freqList = freqList + [Freq* 1.0e-9]
    AC = AC + [np.zeros([antNum, timeNum, 2, chNum])]
#
#-------- Load autocorrelation power spectra
for ant_index in range(antNum):
    for spw_index in range(spwNum):
        timePointer = 0
        for scan_index in range(scanNum):
            timeStamp, Pspec = GetPSpecScan(msfile, ant_index, spwList[spw_index], scanList[scan_index])    # Pspec[pol, ch, time]
            recNum = len(timeStamp)
            AC[spw_index][ant_index, timePointer:(timePointer + recNum)] = Pspec.transpose(2, 0, 1) # AC[spw][ant, time, pol, ch] 
            timePointer += recNum
        #
    #
#
#-------- Plot BP
plotAC(prefix, antList, spwList, freqList, AC)
