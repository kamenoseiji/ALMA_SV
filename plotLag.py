execfile('/users/skameno/ALMA_SV/interferometry.py')
#prefix = 'uid___A002_Xdf15dc_X307'
prefix = 'uid___A002_Xe3a5fd_X15e3d'
msfile = prefix + '.ms'
spw = 31
Scan = 3
timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, Scan)
antList = GetAntName(msfile)
DA48 = np.where(antList == 'DA48')[0][0]
bl0DA48 = np.where(np.array(ANT0) == DA48)[0].tolist()
bl1DA48 = np.where(np.array(ANT1) == DA48)[0].tolist()
Lag = arange(-1920,1920)
XX = np.mean(Xspec[0], axis=2)  # Time Average
#XY = np.mean(Xspec[1], axis=2)  # Time Average
#YX = np.mean(Xspec[2], axis=2)  # Time Average
YY = np.mean(Xspec[1], axis=2)  # Time Average
CorrXX = np.apply_along_axis(spec2corr, 0, XX)
#CorrXY = np.apply_along_axis(spec2corr, 0, XY)
#CorrYX = np.apply_along_axis(spec2corr, 0, YX)
CorrYY = np.apply_along_axis(spec2corr, 0, YY)


plt.plot(Lag, CorrXX[:, bl1DA48[2]].real, ls='steps-mid', label='Real')
plt.plot(Lag, CorrXX[:, bl1DA48[2]].imag, ls='steps-mid', label='Imag')
#plt.plot(Lag, CorrXY[:, bl1DA48[2]].real, ls='steps-mid', label='XY-Real')
#plt.plot(Lag, CorrXY[:, bl1DA48[2]].imag, ls='steps-mid', label='XY-Imag')
#plt.plot(Lag, CorrYX[:, bl1DA48[2]].real, ls='steps-mid', label='YX-Real')
#plt.plot(Lag, CorrYX[:, bl1DA48[2]].imag, ls='steps-mid', label='YX-Imag')
#plt.title(prefix + ' DA48-DA48 SPW43 Scan4 Correlation Function')
plt.title(prefix + ' DA48 Correlation Function')
plt.xlabel('Lag')
plt.ylabel('Corrlation Coeff.')
plt.legend(loc='lower left')

'''
plt.plot(Lag, CorrXX[:, 0].real, ls='steps-mid', label='Real')
plt.plot(Lag, CorrXX[:, 0].imag, ls='steps-mid', label='Imag')
plt.title(prefix + ' DA48-DA42 SPW11 Scan3 Correlation Function')
plt.xlabel('Lag')
plt.ylabel('Corrlation Coeff.')
plt.legend(loc='lower left')
'''
