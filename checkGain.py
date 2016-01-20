#---- Script for Band-3 Astroholograpy Data
from scipy import stats
import matplotlib.pyplot as plt
execfile(SCR_DIR + 'interferometry.py')
<<<<<<< HEAD
#-------- Load BP and Delay
if BPCAL:
    try: 
        BP_ant = np.load( wd + BPprefix + '-BPant.npy' )
    except:
        BPCAL = False
    #
#
=======
>>>>>>> ef7dadde44a7768568b5e50915a426e8dd40756a
#-------- Definitions
antNum = len(refant)
blNum = antNum* (antNum - 1) / 2 
spwNum  = len(spw)
polNum  = len(pol)
#
#-------- Load BP and Delay
"""
if BPCAL:
    for spw_index in range(spwNum):
    BPdata = np.load(wd.kk
    BP_ant = np.load( wd + BPprefix + '.BPant.npy' )
    except:
        BPCAL = False
    #
#
"""
#-------- Procedures
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile)[refant]
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
blMap = range(blNum)
blInv = [False]* blNum		# True -> inverted baseline
for bl_index in range(blNum):
	ants = Bl2Ant(bl_index)
	blMap[bl_index], blInv[bl_index]  = Ant2BlD(refant[ants[0]], refant[ants[1]])
#
#-------- Procedures
interval, timeStamp = GetTimerecord(msfile, 0, 0, pol[0], spw[0], TGscan)
timeNum = len(timeStamp)
#
#-------- Prepare Plots
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index, figsize = (8, 11))
    figAnt.suptitle(prefix + ' ' + antList[ant_index] + ' Scan = ' + `TGscan`)
    figAnt.text(0.45, 0.05, 'MJD [sec]')
    figAnt.text(0.03, 0.45, 'Gain Amplitude and Phase', rotation=90)
#
Gain_ant = np.ones([antNum, spwNum, polNum, timeNum], dtype=complex)
#-------- Baseline-based bandpass
"""
if BPCAL:
    BP_bl = np.ones([spwNum, polNum, chNum, blNum], dtype=complex)
    for bl_index in range(blNum):
        ants = Bl2Ant(bl_index)
        BP_bl[:,:,:,bl_index] = BP_ant[ants[0]]* BP_ant[ants[1]].conjugate()
    #
#
"""
#-------- Loop for SPW
CaledVis = np.zeros([polNum, spwNum, blNum, timeNum], dtype=complex)
for spw_index in range(spwNum):
    print ' Loading SPW = ' + `spw[spw_index]`
    chNum, chWid, Freq = GetChNum(msfile, spw[spw_index]); Freq = 1.0e-9* Freq  # GHz
    chRange = range(int(round(chNum/chBunch * 0.05)), int(round(chNum/chBunch * 0.95)))
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw[spw_index], TGscan)
    Xspec = Xspec[pol]; Xspec = Xspec[:,:,blMap]
    #-------- Baseline-based cross power spectra
    Ximag = Xspec.transpose(0,1,3,2).imag * (-2.0* np.array(blInv) + 1.0)
    Xspec.imag = Ximag.transpose(0,1,3,2)
    #-------- Tsys calibration
    if TSYSCAL:
        TrxSP, TsysSP = TsysSpec(msfile, pol, TsysScan, spw[spw_index], F)
        TsysBL = np.ones([polNum, chNum, blNum])
        for bl_index in range(blNum):
            ants = Bl2Ant(bl_index)
            TsysBL[:,:,bl_index] = sqrt( TsysSP[ants[0]]* TsysSP[ants[1]] )
        #
        for time_index in range(timeNum):
            Xspec[:,:,:,time_index] *= TsysBL
        #
    #
    #-------- Delay Determination and calibration
    if BPCAL:
        BP_ant = np.load( BPfile[spw_index] )
        tempVis = np.mean( Xspec.transpose(3,0,1,2) / (BP_ant[ant1].conjugate()* BP_ant[ant0]).transpose(1,2,0) , axis=2 ).transpose(1, 2, 0)
        # tempVis = np.mean( Xspec.transpose(3,0,1,2) / BP_bl[spw_index], axis=2).transpose(1,2,0)
    elif chNum == 1:
        tempVis = Xspec[:,0,:,:]
    else:
        tempVis = np.mean(Xspec, axis=1) 
    #
    #
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
    for pol_index in range(polNum):
        vis_bl  = tempVis[pol_index]
        Gain_ant[:, spw_index, pol_index] = np.apply_along_axis(gainComplex, 0, vis_bl )
        CaledVis[pol_index, spw_index] = gainCalVis(vis_bl, Gain_ant[:, spw_index, pol_index], Gain_ant[:, spw_index, pol_index])
    #
#
#plotMax = 2.0* np.median(abs(Gain_ant))
plotMax = 1.2* np.median(abs(Gain_ant))
#-------- Plot Gain
for ant_index in range(antNum):
    figAnt = plt.figure(ant_index)
    for pol_index in range(polNum):
        GAampPL = figAnt.add_subplot( 4, 1, 2* pol_index + 1 )
        GAphsPL = figAnt.add_subplot( 4, 1, 2* pol_index + 2 )
        for spw_index in range(spwNum):
            plotGA = Gain_ant[ant_index, spw_index, pol_index]
            GAampPL.plot( timeStamp, abs(plotGA), ls='steps-mid', label = 'SPW=' + `spw[spw_index]`)
            GAphsPL.plot( timeStamp, np.angle(plotGA), '.', label = 'SPW=' + `spw[spw_index]`)
            text_sd = '%s %s %02d : %5.2f' % (antList[ant_index], polName[pol_index], spw[spw_index], 100.0* np.max(abs(Gain_ant[ant_index, spw_index, pol_index])**2))
            print text_sd
        #
        GAampPL.set_ylim(0.0, 1.25* plotMax )
        GAphsPL.set_ylim(-math.pi, math.pi)
        GAampPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
        GAphsPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
        GAampPL.text( np.min(timeStamp), 1.1* plotMax, polName[pol[pol_index]] + ' Amp')
        GAphsPL.text( np.min(timeStamp), 2.5, polName[pol[pol_index]] + ' Phase')
    #
    figAnt.savefig('GA_' + prefix + '_' + antList[ant_index] + '_Scan' + `TGscan` + '.pdf')
#
np.save(prefix + '.Ant.npy', antList) 
np.save(prefix + 'Scn' + `TGscan` + '.GA.npy', Gain_ant) 
plt.close('all')
