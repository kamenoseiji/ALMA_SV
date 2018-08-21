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
#----------------------------------------- SPW and Band name
msmd.open(msfile)
spwName = msmd.namesforspws([spw])
spwPattern = r'RB_..'
BandName = re.findall(spwPattern, spwName[0])[0]
BandPA = (BANDPA[int(BandName[3:5])] + 90.0)*pi/180.0
msmd.close()
msmd.done()
#----------------------------------------- Pol-Cal query
interval, timeStamp = GetTimerecord(msfile, 0, 0, spw, scan)
os.system('rm -rf CalQU.data')
text_sd = R_DIR + 'Rscript %spolQuery.R -D%s -F%f' % (SCR_DIR, qa.time('%fs' % (timeStamp[0]), form='ymd')[0], BANDFQ[int(BandName[3:5])])
text_sd = text_sd + ' ' + PolCal
os.system(text_sd)
fp = open('CalQU.data')
lines = fp.readlines()
fp.close()
IQU = range(3)
for index in range(3): IQU[index] = float( lines[0].split()[index + 1] )
QUsol = np.array([IQU[1]/IQU[0], IQU[2]/IQU[0]])
#----------------------------------------- Find BP table
BPpath = BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-BPant.npy'
XYpath = BPprefix + '-REF' + refantName + '-SPW' + `spw` + '-XYspec.npy'
#QUpath = QUprefix + '-SPW' + `spw` + '-' + refantName +'.QUXY.npy'
if not os.path.exists(BPpath): sys.exit('No BP table [%s]' % (BPpath))
if not os.path.exists(XYpath): sys.exit('No XY table [%s]' % (XYpath))
refAntID = np.where(antList == refantName)[0][0]
#----------------------------------------- Find BP and D-term references
for ant_index in range(UseAntNum):
    Dpath = Dprefix + '-SPW' + `spw` + '-' + antList[UseAnt[ant_index]] + '.DSpec.npy'
    if os.path.exists(Dpath): DantList += [UseAnt[ant_index]]
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
        Dpath = Dprefix + '-SPW' + `spw` + '-' + antList[antMap[ant_index]] + '.DSpec.npy'
        Dspec = np.load(Dpath)
        DxList = DxList + [Dspec[1] + (0.0+1.0j)* Dspec[2]]
        DyList = DyList + [Dspec[3] + (0.0+1.0j)* Dspec[4]]
    #
    print ''
    print 'Antennas without D-term (%d):' % noDantnum,
    for ant_index in range(DantNum, UseAntNum):
        print '%s' % (antList[antMap[ant_index]]),
    #
    print ''
else: sys.exit('Numbers of [Dant] and [noDant] = %d and %d' % (DantNum, noDantnum))
#----------------------------------------- QU table
#if os.path.exists(QUpath): QUsol = np.load(QUpath)     # [Q, U]
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
chOut = sort(list(set(range(chNum)) - set(chRange)))
#----------------------------------------- D-term table to store
Dx, Dy = np.zeros([UseAntNum, chNum], dtype=complex), np.zeros([UseAntNum, chNum], dtype=complex)
Dx[0:DantNum], Dy[0:DantNum] = np.array(DxList), np.array(DyList)
#----------------------------------------- Load visibilities
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
#-------- AZ, EL, PA
azelTime, AntID, az, el = GetAzEl(msfile)
AZ, EL = AzElMatch(timeStamp, azelTime, AntID, refAntID, az, el); PA = AzEl2PA(AZ, EL, ALMA_lat) + BandPA; PAnum = len(PA)
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
#for ch_index in range(chNum):
for ch_index in chRange:
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
for ant_index in range(DantNum, UseAntNum):
    print 'Determining D-term of ' + antList[antMap[ant_index]]
    DtransBL = range(ant_index* (ant_index - 1) / 2, ant_index* (ant_index - 1) / 2 + DantNum)
    #for ch_index in range(chNum):
    for ch_index in chRange:
        Dx[ant_index, ch_index], Dy[ant_index, ch_index] = TransferD(VisSpec[ch_index][:, DtransBL], Dx[0:DantNum, ch_index], Dy[0:DantNum, ch_index], PS) 
    #
    Dx[ant_index, chOut] = np.median( Dx[ant_index, chRange] )
    Dy[ant_index, chOut] = np.median( Dy[ant_index, chRange] )
#
#-------- Save D-term spectra
np.save(prefix + '-SPW' + `spw` + '-' + refantName + '.Ant.npy', antList[antMap])
np.save(prefix + '-SPW' + `spw` + '-' + refantName + '.QUXY.npy', QUsol )
for ant_index in range(UseAntNum):
    DtermFile = np.array([Freq, Dx[ant_index].real, Dx[ant_index].imag, Dy[ant_index].real, Dy[ant_index].imag])
    np.save(prefix + '-SPW' + `spw` + '-' + antList[antMap[ant_index]] + '.DSpec.npy', DtermFile)
#
#-------- Plot D-term spectra
pp = PdfPages('D_' + prefix + '-SPW' + `spw` + '-Dspec.pdf')
for ant_index in range(DantNum, UseAntNum):
    figAnt = plt.figure(ant_index, figsize = (11, 8))
    figAnt.suptitle(prefix + ' ' + antList[antMap[ant_index]])
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'D-term Spectra (Real and Imaginary)', rotation=90)
#
#-------- Plot D-spec
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
