execfile(SCR_DIR + 'interferometry.py')
from matplotlib.backends.backend_pdf import PdfPages
import sys
ALMA_lat = -23.029/180.0*pi     # ALMA AOS Latitude
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
#
#----------------------------------------- Antenna Mapping
# antList : (eg.) array(['DA41', 'DA42', 'DV01', ... ])     : antenna name list ordered in MS
# antID   : [0,1,2,3,4,5,6,...]                             : one-by-one number on antList 
# fragAnt : (eg.) [2,15]                                    : subset of antID to flag out
# antMap  : (eg.) [32,0,1,4, ... ]                          : subset of antID by canonical order, fragged antennas are not included
# trkAntMap : (eg.) [32,0,1,4, ...]                         : subset of antID for tracking antennas by canonical order
# refantID :  (eg.) 32                                      : antID for the reference antenna
# scnAnt  : (eg.) [33,3,38,...]                             : subset of antID for scanning antennas
#----------------------------------------- Procedures
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile); antNum = len(antList); blNum = antNum* (antNum - 1)/2
DantList, noDlist = [], []
flagAnt = indexList(antFlag, antList)
UseAnt = sort(list(set(range(antNum)) - set(flagAnt))).tolist()
UseAntNum = len(UseAnt); UseBlNum  = UseAntNum* (UseAntNum - 1) / 2
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
#----------------------------------------- Find BP table
BPpath = BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-BPant.npy'
XYpath = BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYspec.npy'
QUpath = QUprefix + '-SPW' + `spw` + '-' + refantName +'.QUXY.npy'
if not os.path.exists(BPpath): sys.exit('No BP table [%s]' % (BPpath))
if not os.path.exists(XYpath): sys.exit('No XY table [%s]' % (XYpath))
if not os.path.exists(QUpath): sys.exit('No QU table [%s]' % (QUpath))
refAntID = np.where(antList == refantName)[0][0]
#----------------------------------------- Find BP and D-term references
for ant_index in range(UseAntNum):
    Dxpath = Dprefix + '-SPW' + `spw` + '-' + antList[UseAnt[ant_index]] + '.DxSpec.npy'
    if os.path.exists(Dxpath): DantList += [UseAnt[ant_index]]
    else: noDlist += [UseAnt[ant_index]]
    #
#
if not refAntID in DantList:  sys.exit('No D-term file for the refant [%s]' % (refantName) )
antMap = [refAntID] + list(set(DantList) - set([refAntID])) + noDlist
DantNum, noDantnum = len(DantList), len(noDlist)
DblNum = DantNum* (DantNum - 1) / 2
DxList, DyList = [], []
if (DantNum * noDantnum) !=  0:
    print 'Antennas with D-term (%d):' % DantNum,
    for ant_index in range(DantNum):
        print '%s' % (antList[antMap[ant_index]]),
        Dxpath = Dprefix + '-SPW' + `spw` + '-' + antList[antMap[ant_index]] + '.DxSpec.npy'
        Dypath = Dprefix + '-SPW' + `spw` + '-' + antList[antMap[ant_index]] + '.DySpec.npy'
        DxList = DxList+[np.load(Dxpath)[1]]
        DyList = DyList+[np.load(Dypath)[1]]
    #
    print ''
    print 'Antennas without D-term (%d):' % noDantnum,
    for ant_index in range(DantNum, UseAntNum):
        print '%s' % (antList[antMap[ant_index]]),
    #
    print ''
else: sys.exit('Numbers of [Dant] and [noDant] = %d and %d' % (DantNum, noDantnum))
#----------------------------------------- QU table
QUsol = np.load(QUpath)     # [Q, U]
#QUsol = np.array([0.0, 0.0])     # [Q, U]
#----------------------------------------- BP table
blMap, blInv= range(UseBlNum), [False]* UseBlNum
for bl_index in range(UseBlNum): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print '  ' + `len(np.where( blInv )[0])` + ' baselines are inverted.'
BPantList, BP_ant, XYspec = np.load(BPprefix + '-REF' + refantName + '.Ant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-BPant.npy'), np.load(BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYspec.npy')
BP_ant = BP_ant[indexList(antList[antMap], BPantList)]      # BP antenna mapping
BP_ant[:,1] *= XYspec                                       # XY phase cal
BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()    # Baseline-based bandpass table
#----------------------------------------- Access to MS
msfile = wd + prefix + '.ms'
chNum, chWid, Freq = GetChNum(msfile, spw); chRange = range(int(0.05*chNum), int(0.95*chNum)); Freq *= 1.0e-9   # [Hz]-> [GHz]
#----------------------------------------- D-term table to store
Dx, Dy = np.zeros([UseAntNum, chNum], dtype=complex), np.zeros([UseAntNum, chNum], dtype=complex)
Dx[0:DantNum], Dy[0:DantNum] = np.array(DxList), np.array(DyList)
#----------------------------------------- Load visibilities
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
#-------- AZ, EL, PA
azelTime, AntID, az, el = GetAzEl(msfile)
AZ, EL = AzElMatch(timeStamp, azelTime, AntID, refAntID, az, el); PA = AzEl2PA(AZ, EL, ALMA_lat) + BANDPA; PAnum = len(PA)
#----------------------------------------- Visibility Calibration
tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3, 2, 0, 1)
BPCaledXspec = (tempSpec / BP_bl).transpose(2,3,1,0)            # BP Cal
chAvgVis = np.mean(BPCaledXspec[:,chRange], axis=1)
GainX, GainY = polariGain( chAvgVis[0], chAvgVis[3], PA, QUsol[0], QUsol[1])
Gain = np.array([GainX, GainY])
caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
Gain[1] *= np.sign(np.cos( XY2Phase(PA, QUsol[0], QUsol[1], np.mean(caledVis, axis=1)[[1,2]]) ))
VisSpec = BPCaledXspec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
#----------------------------------------- D-term-corrected Stokes visibilities
ant0, ant1 = ANT0[0:DblNum], ANT1[0:DblNum]
PS = InvPAVector(PA, np.ones(PAnum))
StokesVis = np.zeros([4, chNum], dtype=complex)
for ch_index in range(chNum):
    Minv = InvMullerVector(Dx[ant1, ch_index], Dy[ant1, ch_index], Dx[ant0, ch_index], Dy[ant0, ch_index], np.ones(DblNum))
    StokesVis[:,ch_index] = PS.reshape(4, 4*PAnum).dot(Minv.reshape(4, 4*DblNum).dot(VisSpec[ch_index][:,0:DblNum].reshape(4*DblNum, PAnum)).reshape(4*PAnum)) / (PAnum* DblNum)
QUsol = np.mean(StokesVis[[1,2]][:,chRange], axis=1).real
print 'Stokes Measurement: (Q, U) = (%6.3f, %6.3f)' % (QUsol[0], QUsol[1])
#----------------------------------------- Improved Gain using improved (Q,U)
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
GainX, GainY = polariGain( chAvgVis[0], chAvgVis[3], PA, QUsol[0], QUsol[1])
Gain = np.array([GainX, GainY])
caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
Gain[1] *= np.sign(np.cos( XY2Phase(PA, QUsol[0], QUsol[1], np.mean(caledVis, axis=1)[[1,2]]) ))
VisSpec = BPCaledXspec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
#-------- get D-term spectra
PS = PSvector(PA, np.array([1.0, QUsol[0], QUsol[1], 0.0])).real
DxTemp, DyTemp = np.zeros([chNum, PAnum], dtype=complex), np.zeros([chNum, PAnum], dtype=complex)
for ant_index in range(DantNum, UseAntNum):
    print 'Determining D-term of ' + antList[antMap[ant_index]]
    DtransBL = range(ant_index* (ant_index - 1) / 2, ant_index* (ant_index - 1) / 2 + DantNum)
    for ch_index in range(chNum):
        Dx[ant_index, ch_index], Dy[ant_index, ch_index] = TransferD(VisSpec[ch_index][:, DtransBL], Dx[0:DantNum, ch_index], Dy[0:DantNum, ch_index], PS) 
    #
#
spwNum = 1
#-------- Plot D-term spectra
pp = PdfPages('D_' + prefix + '-Dspec.pdf')
#for ant_index in range(UseAntNum):
for ant_index in range(DantNum, UseAntNum):
    figAnt = plt.figure(ant_index, figsize = (11, 8))
    figAnt.suptitle(prefix + ' ' + antList[antMap[ant_index]])
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'D-term Spectra (Real and Imaginary)', rotation=90)
#
#-------- Plot D-spec
#for ant_index in range(UseAntNum):
for ant_index in range(DantNum, UseAntNum):
    figAnt = plt.figure(ant_index)
    DxPL = figAnt.add_subplot( 2, 1, 1 )
    DyPL = figAnt.add_subplot( 2, 1, 2 )
    #
    plotDx, plotDy = Dx[ant_index], Dy[ant_index]
    DxPL.plot( Freq, plotDx.real, ls='steps-mid', label = 'reDx')
    DxPL.plot( Freq, plotDx.imag, ls='steps-mid', label = 'imDx')
    DxPL.axis([np.min(Freq), np.max(Freq), -0.12, 0.12])
    #
    DyPL.plot( Freq, plotDy.real, ls='steps-mid', label = 'reDy')
    DyPL.plot( Freq, plotDy.imag, ls='steps-mid', label = 'imDy')
    DyPL.axis([np.min(Freq), np.max(Freq), -0.12, 0.12])
    #
    DxPL.text( np.min(Freq), 0.09, 'SPW=' + `spw`)
    #
    DxPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
    DyPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
    figAnt.savefig(pp, format='pdf')
#
#
plt.close('all')
pp.close()
