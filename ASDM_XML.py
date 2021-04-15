# Module to read XML file in ASDM
#
import xml.etree.ElementTree as ET

def BBLOfreq( ASDM ):
    #-------- BB-SPWID connection
    SPW_XML = ASDM + '/' + 'SpectralWindow.xml'
    tree = ET.parse(SPW_XML)
    root = tree.getroot()
    spwList, BBList, spwIDList = [], [], []
    for row in root.findall('row'):
        #---- Check by BB name
        for BBID in row.findall('basebandName'): BBname = BBID.text
        if BBname == 'NOBB': continue
        BBindex = int(BBname.split('_')[1]) - 1
        #---- Check by SPW name
        for spwName in row.findall('name'): SPWname = spwName.text
        if 'FULL_RES' not in SPWname: continue
        #---- Check SPWID
        for spwID in row.findall('spectralWindowId'): spw = int(spwID.text.split('_')[1])
        BBList = BBList + [BBindex]
        spwList = spwList + [spw]
    #
    UniqBBList = sorted(set(BBList), key=BBList.index) 
    for BB in UniqBBList: spwIDList = spwIDList + [spwList[max(indexList(np.array([UniqBBList[BB]]), np.array(BBList)))]]
    #
    #-------- Check LO frequencies 
    spwNum = len(spwIDList)
    RB_XML = ASDM + '/' + 'Receiver.xml'
    tree = ET.parse(RB_XML)
    root = tree.getroot()
    LO2  = np.zeros(spwNum)
    for row in root.findall('row'):
        #-------- Avoid WVR and SQLD
        for sideBand in row.findall('receiverSideband'): SBname = sideBand.text
        if 'NOSB' not in SBname: continue
        #-------- Identify BB from SPWID
        for spwID in row.findall('spectralWindowId'):
            spwIDnumber = int(spwID.text.split('_')[1])    
            if spwIDnumber in spwIDList:
                BB_index = spwIDList.index(spwIDnumber)
                for freq in row.findall('freqLO'):
                    freqList = freq.text.split()
                    LO1 = float(freqList[2])
                    LO2[BB_index] = float(freqList[3])
                #
            #
        #
    #
    return LO1, LO2.tolist()
#
