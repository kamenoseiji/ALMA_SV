fig = plt.figure(figsize = (8,11))
#text_sd = '%s %s - %s %s%d %s' % (prefix, antList[1], antList[2], 'Pol=', pol_index, 'SPW=9')
text_sd = '%s %s - %s %s%d %s' % (prefix, antList[0], antList[1], 'Pol=', pol_index, 'SPW=9')
fig.text(0.25, 0.95, text_sd)
for time_index in range(10):
	plt.subplot( 10, 2, time_index* 2 + 1)
	xlim, ylim = [min(freq), max(freq)], [0.0, 0.004]
	plt.plot(freq, abs(Xspec[0, :, 0, time_index]), ls='steps-mid')
	text_sd = 'Amplitude Time Index = %d / %d' % (time_index, 10)
	plt.text(xlim[1]*0.1, ylim[1]*0.1, text_sd, size='x-small')
	plt.xticks(fontsize=6); plt.yticks(fontsize=6)
	plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]], fontsize=6)

	plt.subplot( 10, 2, time_index* 2 + 2)
	xlim, ylim = [min(freq), max(freq)], [-pi, pi]
	plt.plot(freq, np.angle(Xspec[0, :, 0, time_index]), '.')
	text_sd = 'Phase Time Index = %d / %d' % (time_index, 10)
	plt.text(xlim[1]*0.1, ylim[1]*0.7, text_sd, size='x-small')
	plt.xticks(fontsize=6); plt.yticks(fontsize=6)
	plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]], fontsize=6)

