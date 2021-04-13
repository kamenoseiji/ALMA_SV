# Module to read XML file in ASDM
#
import xml.etree.ElementTree as ET

def spwLOfreq( ASDM ):
    RB_XML = ASDM + '/' + 'Receiver.xml'
    tree = ET.parse(RB_XML)
    root = tree.getroot()
    spwList, LO1List, LO2List  = [], [], []
    for row in root.findall('row'):
        for spwID in row.findall('spectralWindowId'):
            spwList = spwList + [int(spwID.text.split('_')[1])]
        #
        for freq in row.findall('freqLO'):
            freqList = freq.text.split()
            LO1List = LO1List + [float(freqList[2])]
            LO2List = LO2List + [float(freqList[3])]
        #
    #
    return spwList, LO1List, LO2List
#
