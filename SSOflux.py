import sys
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
execfile(SCR_DIR + 'interferometry.py')
#
msmd.open(msfile)
#-------- Flux models for solar system objects
SSONum = len(BandSSOList)
timeLabel = qa.time('%fs' % (timeXY[0]), form='ymd')[0]
SSOflux0, SSOshape, centerFreqList = [], [], []
#-------- Primary beam for each baseline
primaryBeam = np.ones([UseBlNum])
for bl_index in range(UseBlNum):
    beam0, beam1 = 1.0/antDia[antMap[ant0[bl_index]]], 1.0/antDia[antMap[ant1[bl_index]]] 
    primaryBeam[bl_index] = np.sqrt(2.0/ ((beam0)**2 + (beam1)**2 )) * beam0* beam1
#
#-------- Center frequency of each SPW
for spw_index in range(spwNum): 
    chNum, chWid, Freq = GetChNum(msfile, spw[spw_index])
    centerFreqList.append( np.median(Freq)*1.0e-9 )
#
#-------- SSO Model
for ssoIndex in range(SSONum):
    for spw_index in range(spwNum): 
        text_Freq = '%6.2fGHz' % (centerFreqList[spw_index])
        SSOmodel = predictcomp(objname=sourceList[BandSSOList[ssoIndex]], standard="Butler-JPL-Horizons 2012", minfreq=text_Freq, maxfreq=text_Freq, nfreqs=1, prefix="", antennalist="aca.cycle3.cfg", epoch=timeLabel, showplot=True)
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
        # scanID = list(set( msmd.scansforfield(BandSSOList[ssoIndex]).tolist()) & set(onsourceScans))[0]; SSOscanID.append(scanID)
        scanID = list(set( msmd.scansforfield(sourceList[BandSSOList[ssoIndex]]).tolist()) & set(onsourceScans))[0]; SSOscanID.append(scanID)
    except:
        continue
    #
    FLScaleText = 'Flux Calibrator is %s at %s' % (sourceList[BandSSOList[ssoIndex]], timeLabel);  print FLScaleText
    timeStamp, UVW = GetUVW(msfile, spw[spw_index], scanID)
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    for spw_index in range(spwNum):
        uvWave = uvw[0:2,:]* centerFreqList[spw_index] / 0.299792458    # UV distance in wavelength
        uvFlag[ssoIndex, spw_index, np.where( uvDist* centerFreqList[spw_index] / 0.299792458 > UVlimit )[0].tolist()] = 0.0
        SSOmodelVis = SSOmodelVis + [diskVisBeam(SSOshape[ssoIndex], uvWave[0], uvWave[1], 1.13* 0.299792458* primaryBeam/centerFreqList[spw_index])]
        #-------- for debug
        text_sd = 'SPW=%d uv limit = %5.0f klambda' % (spw[spw_index], UVlimit*1.0e-3); print text_sd
        for ant0_index in range(1, UseAntNum):
            text_sd = antList[antMap[ant0_index]] + ' : '; print text_sd,
            blList = np.where(np.array(ant0) == ant0_index)[0].tolist()
            for bl_index in blList:
                uvLambda = uvDist[bl_index]* centerFreqList[spw_index] / 0.299792458
                if uvFlag[ssoIndex, spw_index, bl_index] < 1.0: text_sd = '\033[91m%4.0f\033[0m' % (uvLambda*1.0e-3)
                else: text_sd = '%4.0f' % (uvLambda*1.0e-3)
                print text_sd,
            #
            print ''
        #
        print '       ',
        for ant0_index in range(0, UseAntNum - 1):
            print antList[antMap[ant0_index]],
        #
        print ''
    #
#
SSOmodelVis = np.array(SSOmodelVis).reshape(SSONum, spwNum, UseBlNum)
#FCSmodelVis = SSOmodelVis[FCS_ID]
#FCSFlag     = uvFlag[FCS_ID]
#
msmd.done()
