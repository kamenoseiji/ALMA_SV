import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
execfile(SCR_DIR + 'Plotters.py')
#
msmd.open(msfile)
#-------- Flux models for solar system objects
SSONum = len(BandSSOList)
timeLabel = qa.time('%fs' % (timeXY[0]), form='ymd')[0]
SSOflux0, SSOshape, centerFreqList = [], [], []
primaryBeam = np.ones([UseBlNum])
for bl_index in range(UseBlNum):
    beam0, beam1 = 1.0/antDia[ant0[bl_index]], 1.0/antDia[ant1[bl_index]] 
    primaryBeam[bl_index] = np.sqrt(2.0/ ((beam0)**2 + (beam1)**2 )) * beam0* beam1
#
for spw_index in range(spwNum): 
    chNum, chWid, Freq = GetChNum(msfile, spw[spw_index])
    centerFreqList.append( np.median(Freq)*1.0e-9 )
#
for ssoIndex in range(SSONum):
    for spw_index in range(spwNum): 
        text_Freq = '%6.2fGHz' % (centerFreqList[spw_index])
        SSOmodel = predictcomp(objname=sourceList[BandSSOList[ssoIndex]], standard="Butler-JPL-Horizons 2012", minfreq=text_Freq, maxfreq=text_Freq, nfreqs=1, prefix="", antennalist="aca.cycle3.cfg", epoch=timeLabel, showplot=T)
        SSOflux0.append(SSOmodel['spectrum']['bl0flux']['value'])
    #
    MajMinPA = np.array([SSOmodel['shape']['majoraxis']['value']* pi / 21600.0, SSOmodel['shape']['minoraxis']['value']* pi / 21600.0, SSOmodel['shape']['positionangle']['value']* pi / 180.0])
    SSOshape.append(MajMinPA)   # arcmin -> rad, diameter -> radius
#
plt.close('all')
SSOflux0= np.array(SSOflux0).reshape(SSONum, spwNum)     # [SSO, spw]
uvFlag = np.ones([SSONum, spwNum, UseBlNum])
SSOmodelVis, SSOscanID = [], []
for ssoIndex in range(SSONum):
    UVlimit = 0.32 / SSOshape[ssoIndex][0]  # Maximum uv distane(lambda) available for the SSO size
    try:
        scanID = list(set( msmd.scansforfield(BandSSOList[ssoIndex]).tolist()) & set(onsourceScans))[0]; SSOscanID.append(scanID)
    except:
        continue
    if( scanID == FCScan):
        FCS_ID = ssoIndex
        text_sd = 'Flux Calibrator is %s at %s' % (sourceList[BandSSOList[ssoIndex]], timeLabel); logfile.write(text_sd + '\n'); print text_sd
    #
    timeStamp, UVW = GetUVW(msfile, spw[spw_index], scanID)
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    for spw_index in range(spwNum):
        uvWave = uvw[0:2,:]* centerFreqList[spw_index] / 0.299792458    # UV distance in wavelength
        uvFlag[ssoIndex, spw_index, np.where( uvDist* centerFreqList[spw_index] / 0.299792458 > UVlimit )[0].tolist()] = 0.0
        SSOmodelVis = SSOmodelVis + [diskVisBeam(SSOshape[ssoIndex], uvWave[0], uvWave[1], 1.13* 0.299792458* primaryBeam/centerFreqList[spw_index])]
        #-------- for debug
        #for bl_index in np.where(uvFlag[ssoIndex, spw_index] == 0.0)[0].tolist():
        #    print '%s SPW=%d : Flagged BL(%s - %s)' % (sourceList[BandSSOList[ssoIndex]], spw[spw_index], antList[antMap[ant0[bl_index]]],antList[antMap[ant1[bl_index]]]) 
        #
    #
#
SSOmodelVis = np.array(SSOmodelVis).reshape(SSONum, spwNum, UseBlNum)
FCSmodelVis = SSOmodelVis[FCS_ID]
FCSFlag     = uvFlag[FCS_ID]
#
msmd.done()
