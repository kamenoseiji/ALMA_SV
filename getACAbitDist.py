import numpy as np
from numpy import *
import urllib2
import re
from StringIO import StringIO
import gzip
import datetime

entryKey = 'calc3bitAndDeltaRequantCorrectionValues'
CAIentryKey = 'Antenna-CAI'
def webFileList( URL ):
    fileList = []
    request =  urllib2.Request(URL)
    response = urllib2.urlopen(request)
    http = response.read()
    httpList = http.split('\"')
    for fileName in httpList:
        if fileName.find("acsStartContainer") == 0:
            fileList.append(fileName)
        #
    #
    return fileList
#
def webLogGZ( URI ):
    request =  urllib2.Request(URI)
    request.add_header('Accept-encoding', 'gzip')
    response = urllib2.urlopen(request)
    gzBuf = StringIO( response.read())
    buf = gzip.GzipFile( fileobj = gzBuf )
    textData = buf.read()
    return(textData.split('\n'))
#

def bitDistEntries( URI ):
    Xc = []
    Yc = []
    BBID = []
    CAI = []
    mjdSec = []
    #
    logList = webLogGZ( URI )
    for logEntry in logList:
        if logEntry.find(entryKey) != -1:
            bbid, cai, mjdsec, xc, yc = parseBitDist( logEntry )
            if len(xc) == 8:
                Xc.append(xc)
                Yc.append(yc)
                BBID.append(bbid)
                CAI.append(cai)
                mjdSec.append(mjdsec)
            #
        #
    #
    return CAI, BBID, np.array(mjdSec), np.array(Xc), np.array(Yc)
#
def parseBitDist( logLine ):
    if logLine.find(entryKey) == -1:
        return()
    #
    logItems = filter(lambda w: len(w) > 0, re.split(r'\s|,', logLine))
    try:
        X_index = np.where( np.array(logItems) == 'X' )[0][0]
        Y_index = np.where( np.array(logItems) == 'Y' )[0][0]
    except:
        return [], [], [], [], []
    #
    BBID = int(logItems[X_index - 4])
    CAI  = int(logItems[X_index - 3])
    mjdSec = qa.convert(logItems[X_index - 1], 's')['value']
    if mjdSec < 4.0e9:
        mjdSec = qa.convert(logItems[0], 's')['value']
    #
    Xcounts = map(int, np.array(logItems)[range((X_index+9), (X_index+13)) + range((X_index+5), (X_index+9))])
    Ycounts = map(int, np.array(logItems)[range((Y_index+9), (Y_index+13)) + range((Y_index+5), (Y_index+9))])
    return BBID, CAI, mjdSec, Xcounts, Ycounts
#
def parseCAImap( logLine ):
    if logLine.find(CAIentryKey) == -1:
        return()
    #
    logItems = filter(lambda w: len(w) > 0, re.split(r'\s|,|\/|\(|\)', logLine))
    Ref_index = np.where(np.array(logItems) == '=')[0][0]
    if logItems[Ref_index + 2] == 'AntennaName':
        Ref_index += 6
    #
    return qa.convert(logItems[0], 's')['value'], logItems[Ref_index + 1], int(logItems[Ref_index + 2])
#
def CAImapEntries( URI ):
    mjdSec = []
    CAI = []
    antName = []
    logList = webLogGZ( URI )
    for logEntry in logList:
        if logEntry.find(CAIentryKey) != -1:
            mjdsec, antname, cai = parseCAImap( logEntry )
            mjdSec.append(mjdsec)
            CAI.append(cai)
            antName.append(antname)
        #
    #
    return np.array(mjdSec), CAI, antName 
#


Server = 'http://computing-logs.aiv.alma.cl/'
DIR1   = 'AOS/CONTAINER/'
DIR2   = '/alma/logs/coj-cpm-1/ACACORR/CDPMIF/'
DIR2a  = '/alma/logs/COJ-CPM-1/ACACORR/CDPMIF/'
DIR3   = '/alma/logs/coj-cc-1/ACACORR/OBSERVATION_CONTROL/'
DIR3a  = '/alma/logs/COJ-CC-1/ACACORR/OBSERVATION_CONTROL/'

initDate = datetime.datetime.strptime(DATE, '%Y-%m-%d')
dateText = []
for dateIndex in range(Days):
    currentDate = initDate + datetime.timedelta(days = dateIndex)
    dateText.append(currentDate.strftime('%Y-%m-%d'))
#
#-------- CAI Mapping
CAIMAP  = []
ANTNAME = []
mjdSec  = np.zeros(0)
for dateIndex in range(Days):
    try:
        URL = Server + DIR1 + dateText[dateIndex] + DIR3
        FileList = webFileList(URL)
    except:
        try:
            URL = Server + DIR1 + dateText[dateIndex] + DIR3a
            FileList = webFileList(URL)
        except:
            continue
    #
    print 'CAI Mapping on ' + dateText[dateIndex]
    for fileName in FileList:
        mjdsec, cai, antname = CAImapEntries(URL + fileName)
        if len(cai) > 0:
            mjdSec = np.r_[mjdSec, mjdsec]
            CAIMAP.extend(cai)
            ANTNAME.extend(antname)
        #
    #
#
np.save('CIT.' + DATE + '.npy', mjdSec)
np.save('CCM.' + DATE + '.npy', CAIMAP)
np.save('CAN.' + DATE + '.npy', ANTNAME)
#-------- Level Histogram
CAI = []
BBID = []
mjdSec = np.zeros(0)
XC = np.zeros([0,8])
YC = np.zeros([0,8])
for dateIndex in range(Days):
    try:
        URL = Server + DIR1 + dateText[dateIndex] + DIR2
        FileList = webFileList(URL)
    except:
        try:
            URL = Server + DIR1 + dateText[dateIndex] + DIR2a
            FileList = webFileList(URL)
        except:
            continue
    #
    print 'BitDist on ' + dateText[dateIndex]
    for fileName in FileList:
        cai, bbid, mjdsec, xcount, ycount = bitDistEntries(URL + fileName)
        if len(cai) > 0:
            CAI.extend(cai)
            BBID.extend(bbid)
            mjdSec = np.r_[mjdSec, mjdsec]
            XC = np.vstack([XC, xcount])
            YC = np.vstack([YC, xcount])
        #
    #
#
np.save('XLC.' + DATE + '.npy', XC)
np.save('YLC.' + DATE + '.npy', YC)
np.save('TIM.' + DATE + '.npy', mjdSec)
np.save('CAI.' + DATE + '.npy', CAI)
np.save('BBI.' + DATE + '.npy', BBID)
