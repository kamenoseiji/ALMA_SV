execfile('/users/skameno/Scripts/gaussNbit.py')
def readBitDistBLC(fname):
	file = open(fname)
	readLines = file.readlines()
	file.close()
	#
	bitCounts = np.zeros([4, 2, 8])  # 4BB x 2pol x 8 levels
	Fil_index = -1
	for index in range(len(readLines)):
		words = re.split(r' +', readLines[index])
		if( re.search("Fil", words[0])):
			Fil_index = Fil_index + 1
			BB_index = int(Fil_index/2)
			pol_index = Fil_index % 2
		#
		#-------- parse bit distribution
		if(re.search("S\d\d", words[0])): 
			for level_index in range(8):
				bitCounts[BB_index, pol_index, level_index] = bitCounts[BB_index, pol_index, level_index] + float(words[8 - level_index])
			#
		#
	#
	return 125e3* bitCounts
#

"""
def readBitDistBLC(fname):
	file = open(fname)
	readLines = file.readlines()
	file.close()
	#
	bitCounts = np.zeros([4, 2, 8])  # 4BB x 2pol x 8 levels
	BB_index = -1
	for index in range(len(readLines)):
		#
		words = re.split(r' +', readLines[index])
		#-------- parse BB
		if( re.search("BBpr", words[0])):
			BB_index = int(words[0][4])
		#
		#-------- parse polarization
		#if( words[0] == 'Pol'):
		#    pol_index = np.where( np.array(['X', 'Y']) == words[1][0])[0]
		#
		#-------- parse bit distribution
		if(re.search("S\d\d", words[0])): 
			for level_index in range(8):
				bitCounts[BB_index, pol_index, level_index] = bitCounts[BB_index, pol_index, level_index] + float(words[8 - level_index])
			#
		#
	#
	return 125e3* bitCounts
"""
#
def readBitDistACA(fname):
    file = open(fname)
    readLines = file.readlines()
    file.close()
    #
    bitCounts = np.zeros([4, 2, 8])  # 4BB x 2pol x 8 levels
    for index in range(len(readLines)):
        #
        words = re.split(r' +', readLines[index])
        #-------- parse BB
        if( words[0] == '#' ):
            BB_index = -int(words[2])
        #
        #-------- parse polarization
        if( words[0] == 'pol'):
            pol_index = np.where( np.array(['X', 'Y']) == words[2][0])[0]
        #
        #-------- parse bit distribution
        if( words[0] == 'lvl'): 
            level_index = int( words[2] )
            bitCounts[BB_index, pol_index, level_index] = int(words[5])
        #
    #
    return bitCounts
#
