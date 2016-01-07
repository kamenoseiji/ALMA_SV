#---- Script for Band-3 Astroholograpy Data
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
execfile(SCR_DIR + 'interferometry.py')
RADDEG = 180.0/pi
RADSEC = 3600* 180.0/pi
IMSIZE = 256
HALFSIZE = IMSIZE/2
CELLSIZE = 0.2
#-------- Load BP and Delay
def BP_load( prefix, spw ):
    try: 
        BP_ant = np.load( prefix + '-SPW' + `spw` + '-BPant.npy' )
        return BP_ant
    except:
        return False
    #
#
#BP_ant = BP_load( wd + BPprefix, 10)
#-------- Definitions
antNum = len(refant)
blNum = antNum* (antNum - 1) / 2 
spwNum  = len(PointSPW)
polNum  = len(pol)
ant0, ant1 = ANT0[0:blNum], ANT1[0:blNum]
#
#-------- Procedures for point source
PointMSfile = wd + PointPrefix + '.ms'
PlanetMSfile = wd + PlanetPrefix + '.ms'
antList = GetAntName(PointMSfile)[refant]
blMap = range(blNum)
blInv = [False]* blNum		# True -> inverted baseline
for bl_index in range(blNum):
	ants = Bl2Ant(bl_index)
	blMap[bl_index], blInv[bl_index]  = Ant2BlD(refant[ants[0]], refant[ants[1]])
#
#-------- Procedures
#interval, timeStamp = GetTimerecord(PointMSfile, 0, 0, pol[0], PointSPW[0], PointScan)
#
#-------- Prepare Plots
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index, figsize = (8, 11))
    figAnt.suptitle(PointPrefix + ' ' + antList[ant_index] + ' Scan = ' + `PointScan`)
    figAnt.text(0.45, 0.05, 'MJD [sec]')
    figAnt.text(0.03, 0.45, 'Gain Amplitude and Phase', rotation=90)
#
for spw_index in range(spwNum):
    figMap = plt.figure(antNum + spw_index, figsize = (11, 8))
    figVisD= plt.figure(antNum + spwNum + spw_index, figsize = (11, 8))
    figVisR= plt.figure(antNum + 2* spwNum + spw_index, figsize = (11, 8))
    figVisC= plt.figure(antNum + 3* spwNum + spw_index, figsize = (11, 8))
    figMap.suptitle(PlanetPrefix + ' ' + ' SPW = ' + `PlanetSPW[spw_index]`)
    figVisD.suptitle(PlanetPrefix + ' ' + ' SPW = ' + `PlanetSPW[spw_index]`)
    figVisR.suptitle(PointPrefix + ' ' + ' SPW = ' + `PointSPW[spw_index]`)
    figVisC.suptitle(PointPrefix + ' ' + ' SPW = ' + `PointSPW[spw_index]`)
#
#-------- Loop for SPW
for spw_index in range(spwNum):
    #-------- Baseline-based bandpass
    chNum, chWid, Freq = GetChNum(PointMSfile, PointSPW[spw_index]); Freq = 1.0e-9* Freq  # GHz
    chRange = range(int(round(chNum/chBunch * 0.05)), int(round(chNum/chBunch * 0.95)))
    BP_bl = np.ones([blNum, polNum, chNum], dtype=complex)
    if BPCAL:
        BP_ant = BP_load( wd + BPprefix, PointSPW[spw_index])
        BP_bl *= BP_ant[ant1].conjugate()
        BP_bl *= BP_ant[ant0]
    #
    #-------- Load Visibilities
    print ' Loading SPW = ' + `PointSPW[spw_index]`
    timeStamp, UVW = GetUVW(PointMSfile, PointSPW[spw_index], PointScan)    # UVW[3, blNum, timeNum]
    timeNum = len(timeStamp)
    UVD = np.sqrt(UVW[0]**2 + UVW[1]**2).reshape(blNum*timeNum)
    timeStamp, Pspec, Xspec = GetVisAllBL(PointMSfile, PointSPW[spw_index], PointScan)
    figVisR= plt.figure(antNum + 2* spwNum + spw_index, figsize = (11, 8))
    figVisC= plt.figure(antNum + 3* spwNum + spw_index, figsize = (11, 8))
    Xspec = Xspec[pol]; Xspec = Xspec[:,:,blMap]
    #
    #-------- Baseline-based cross power spectra
    Ximag = Xspec.transpose(0,1,3,2).imag * (-2.0* np.array(blInv) + 1.0)
    Xspec.imag = Ximag.transpose(0,1,3,2)
    Xspec = (Xspec.transpose(3,2,0,1) / BP_bl).transpose(2, 3, 1, 0)
    VisAmp = abs(np.mean(Xspec[:,chRange], axis=1)).reshape(2, blNum*timeNum)
    scaleFact = 1.0 / np.median( VisAmp, axis=1 )
    VisAmpPL = figVisR.add_subplot( 1, 1, 1 )
    for pol_index in range(polNum):
        VisAmpPL.plot( UVD, scaleFact[pol_index]* VisAmp[pol_index], '.', label = polName[pol_index])
    #
    VisAmpPL.set_ylim(0.0, 1.25)
    VisAmpPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
    VisAmpPL.set_xlabel("Baseline Length [m]"); VisAmpPL.set_ylabel("Visibility Amplitude [scaled]")
    figVisR.savefig('VisUncal_' + PointPrefix + '_SPW' + `PointSPW[spw_index]` + '.pdf')
    #-------- Tsys calibration
    TsysBL = np.ones([blNum, polNum, chNum])    # Baseline-based SEFD with Ae=1.0
    if TSYSCAL:
        TrxSP, TsysSP = TsysSpec(PointMSfile, pol, PointTsysScan, PointSPW[spw_index], F)
        TsysBL *= sqrt( TsysSP[ant0] * TsysSP[ant1] )
    #
    Xspec = (Xspec.transpose(3,2,0,1) * TsysBL).transpose(2, 3, 1, 0)
    #-------- Gain Cal
    VisAmpPL = figVisC.add_subplot( 1, 1, 1 )
    Gain_ant = np.ones([antNum, polNum, timeNum], dtype=complex)
    for pol_index in range(polNum):
        vis_bl  = np.mean(Xspec[pol_index, chRange], axis=0)
        Gain_ant[:, pol_index] = np.apply_along_axis(gainComplex, 0, vis_bl )
        CaledVis = gainCalVis(vis_bl, Gain_ant[:, pol_index], Gain_ant[:, pol_index])
        VisAmp = CaledVis.reshape(blNum*timeNum)
        VisAmpPL.plot( UVD, VisAmp, '.', label = polName[pol_index])
    #
    VisAmpPL.set_ylim(0.0, 1.25)
    VisAmpPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
    VisAmpPL.set_xlabel("Baseline Length [m]"); VisAmpPL.set_ylabel("Visibility Amplitude [scaled]")
    figVisC.savefig('VisCaled_' + PointPrefix + '_SPW' + `PointSPW[spw_index]` + '.pdf')
    plotMax = 2.0* np.median(abs(Gain_ant))
    #-------- Plot Gain
    for ant_index in range(antNum):
        figAnt = plt.figure(ant_index)
        for pol_index in range(polNum):
            GAampPL = figAnt.add_subplot( 4, 1, 2* pol_index + 1 )
            GAphsPL = figAnt.add_subplot( 4, 1, 2* pol_index + 2 )
            plotGA = Gain_ant[ant_index, pol_index]
            GAampPL.plot( timeStamp, abs(plotGA), ls='steps-mid', label = 'SPW=' + `PointSPW[spw_index]`)
            GAphsPL.plot( timeStamp, np.angle(plotGA), '.', label = 'SPW=' + `PointSPW[spw_index]`)
            text_sd = '%s %s %02d : %5.2f' % (antList[ant_index], polName[pol_index], PointSPW[spw_index], 100.0* np.max(abs(Gain_ant[ant_index, pol_index])**2))
            print text_sd
            if spw_index == spwNum - 1:
                GAampPL.set_ylim(0.0, 1.25* plotMax )
                GAphsPL.set_ylim(-math.pi, math.pi)
                GAampPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
                GAphsPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
                GAampPL.text( np.min(timeStamp), 1.1* plotMax, polName[pol[pol_index]] + ' Amp')
                GAphsPL.text( np.min(timeStamp), 2.5, polName[pol[pol_index]] + ' Phase')
                figAnt.savefig('GA_' + PointPrefix + '_' + antList[ant_index] + '_Scan' + `PointScan` + '.pdf')
            #
        #
    #
#
#-------- Visibilities for Planet
    print ' Loading SPW = ' + `PlanetSPW[spw_index]`
    timeStamp, Pspec, Xspec = GetVisAllBL(PlanetMSfile, PlanetSPW[spw_index], PlanetScan)
    timeNum = len(timeStamp)
    Xspec = Xspec[pol]; Xspec = Xspec[:,:,blMap]
    #-------- Baseline-based cross power spectra
    Ximag = Xspec.transpose(0,1,3,2).imag * (-2.0* np.array(blInv) + 1.0)
    Xspec.imag = Ximag.transpose(0,1,3,2)
    Xspec = (Xspec.transpose(3,2,0,1) / BP_bl).transpose(2, 3, 1, 0)
    #-------- Tsys calibration
    TsysBL = np.ones([blNum, polNum, chNum])    # Baseline-based SEFD with Ae=1.0
    if TSYSCAL:
        TrxSP, TsysSP = TsysSpec(PlanetMSfile, pol, PlanetTsysScan, PlanetSPW[spw_index], F)
        TsysBL *= sqrt( TsysSP[ant0] * TsysSP[ant1] )
    #
    Xspec = (Xspec.transpose(3,2,0,1) * TsysBL).transpose(2, 3, 1, 0)
    #-------- Gain Cal
    chNum, chWid, Freq = GetChNum(PlanetMSfile, PlanetSPW[spw_index]); Freq = 1.0e-9* Freq  # GHz
    lambdaInv = np.mean(Freq)*1.0e9/constants.c     # 1/wavelength [m^-1]
    for pol_index in range(polNum):
        GainAmp = np.outer(np.mean( abs(Gain_ant[:, pol_index]), axis=1), np.ones([timeNum], dtype=complex))
        vis_bl  = np.mean(Xspec[pol_index, chRange], axis=0)
        for time_index in range(timeNum):
            antPhase, PhsError = clphase_solve( np.angle(vis_bl[:, time_index]), abs(vis_bl[:, time_index]))
            GainAmp[:,time_index] *= exp(1.0j* antPhase)
        #
        CaledVis = gainCalVis(vis_bl, GainAmp, GainAmp)/ scaleFact[pol_index]**2 * np.mean(scaleFact)**2
        #
        #-------- Set visibilities
        timeStamp, UVW = GetUVW(PlanetMSfile, PlanetSPW[spw_index], PlanetScan)    # UVW[3, blNum, timeNum]
        UVW *= lambdaInv
        UVmax = sqrt(np.max( UVW[0]**2 + UVW[1]**2 ))
        UVmin = sqrt(np.min( UVW[0]**2 + UVW[1]**2 ))
        GridMax = RADSEC/CELLSIZE
        xi, yi = np.mgrid[ -GridMax:GridMax:((IMSIZE + 1)* 1j), -GridMax:GridMax:((IMSIZE + 1)* 1j)]       # (u, v) sampling points
        xi, yi = xi[0:IMSIZE, 0:IMSIZE], yi[0:IMSIZE, 0:IMSIZE]
        U = np.r_[UVW[0], -UVW[0]].reshape(2*blNum*timeNum)
        V = np.r_[UVW[1], -UVW[1]].reshape(2*blNum*timeNum)
        reCaledVis = np.r_[ CaledVis.real,  CaledVis.real ].reshape(2*blNum*timeNum)
        imCaledVis = np.r_[ CaledVis.imag, -CaledVis.imag ].reshape(2*blNum*timeNum)
        #
        #-------- Imaging
        GridVis =  GridData( reCaledVis, U, V, xi.reshape(xi.size), yi.reshape(xi.size), UVmin/2 ).reshape(len(xi), len(xi)) + 1.0j* GridData( imCaledVis, U, V, xi.reshape(xi.size), yi.reshape(xi.size), UVmin/2 ).reshape(len(xi), len(xi))
        tempVis = np.zeros([IMSIZE, IMSIZE], dtype=complex)
        GridMap = np.zeros([IMSIZE, IMSIZE])
        tempVis[0:HALFSIZE, 0:HALFSIZE]   = GridVis[HALFSIZE:IMSIZE, HALFSIZE:IMSIZE]
        tempVis[0:HALFSIZE, HALFSIZE:IMSIZE] = GridVis[HALFSIZE:IMSIZE, 0:HALFSIZE]
        tempVis[HALFSIZE:IMSIZE, 0:HALFSIZE] = GridVis[0:HALFSIZE,HALFSIZE:IMSIZE]
        tempVis[HALFSIZE:IMSIZE, HALFSIZE:IMSIZE] = GridVis[0:HALFSIZE,0:HALFSIZE]
        tempMap = np.fft.fft2(tempVis).real / len(np.where(tempVis.real > 0.5)[0])
        GridMap[0:HALFSIZE, 0:HALFSIZE]   = tempMap[HALFSIZE:IMSIZE, HALFSIZE:IMSIZE]
        GridMap[0:HALFSIZE, HALFSIZE:IMSIZE] = tempMap[HALFSIZE:IMSIZE, 0:HALFSIZE]
        GridMap[HALFSIZE:IMSIZE, 0:HALFSIZE] = tempMap[0:HALFSIZE,HALFSIZE:IMSIZE]
        GridMap[HALFSIZE:IMSIZE, HALFSIZE:IMSIZE] = tempMap[0:HALFSIZE,0:HALFSIZE]
        xi, yi = np.mgrid[ -(HALFSIZE* CELLSIZE):(HALFSIZE* CELLSIZE):((IMSIZE + 1)* 1j), -(HALFSIZE* CELLSIZE):(HALFSIZE* CELLSIZE):((IMSIZE + 1)* 1j)] 
        xi, yi = xi[0:IMSIZE, 0:IMSIZE], yi[0:IMSIZE, 0:IMSIZE]
        figMap = plt.figure(antNum + spw_index)
        MapPL = figMap.add_subplot( 1, 2, pol_index + 1, aspect=1 )
        mapImage = MapPL.contourf(xi, yi, GridMap.real, np.linspace(-0.4, 4.8, 14))
        figMap.colorbar(mapImage, ax=MapPL)
        # MapPL.axis('equal')
        #-------- Visibility map
        #figVis = plt.figure(spw_index)
        #ax = Axes3D(figVis)
        #ax.set_xlabel("U [lambda]"); ax.set_ylabel("V [lambda]"); ax.set_zlabel("Visibility Amplitude [scaled by the calibrator]")
        #ax.scatter3D( U, V, reCaledVis)
        #ax.plot( U, V, reCaledVis, '.')
        #figVis.savefig('VIS_' + PlanetPrefix + '_SPW_' + `PlanetSPW[spw_index]` + 'Pol' + polName[pol_index] + '.pdf')
        #
        #-------- Visibility Amp
        UVD = np.sqrt(U**2 + V**2)
        figVisD = plt.figure(antNum + spwNum + spw_index)
        VisDPL = figVisD.add_subplot( 1, 1, 1 )
        VisDPL.plot( UVD, reCaledVis, '.')
        VisDPL.set_xlim(0.0, 1.2* max(UVD) ); VisDPL.set_ylim(0.0, 1.2* max(reCaledVis) )
        VisDPL.set_xlabel("UV distance [lambda]"); VisDPL.set_ylabel("Visibility Amplitude [scaled by the calibrator]")
    #
    figMap.savefig('MAP_' + PlanetPrefix + '_SPW_' + `PlanetSPW[spw_index]` + '.pdf')
    figVisD.savefig('VISD_' + PlanetPrefix + '_SPW_' + `PlanetSPW[spw_index]` + '.pdf')
#
plt.close('all')
#
#np.save(PointPrefix + '.Ant.npy', antList) 
#np.save(PointPrefix + 'Scn' + `PointScan` + '.GA.npy', Gain_ant) 
#plt.close('all')
"""
    #-------- Source Model
    if not PointSource:
        lambdaInv = np.mean(Freq)*1.0e9/constants.c     # 1/wavelength [m^-1]
        timeStamp, UVW = GetUVW(msfile, spw[spw_index], TGscan)
        UVD = np.sqrt(UVW[0]**2 + UVW[1]**2)*lambdaInv  # [BL, Time], in unit of wavelength
        for bl_index in range(blNum):
            ants = Bl2Ant(bl_index)
            antD = sqrt(GetAntD(antList[ants[0]])* GetAntD(antList[ants[1]]))
            FWHM = GetFWHM(msfile, spw[spw_index], antD)
            zeroSpFlux = beamF( diskR / FWHM )* Tb2Flux(SourceTemp, np.mean(Freq), diskR)
            PiUD = 1.5230870989335429e-05 * UVD[bl_index]* diskR       # pi* u * D,
            VisJy = PiUD / (2.0* scipy.special.jv(1,PiUD)) 
            tempVis[:, bl_index] *= VisJy
        #
    #
"""
