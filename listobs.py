# wd = '/Volumes/SSD/ALMA_SV/POLBEAM/'
# prefix = ['uid___A002_X9b92b0_X408', 'uid___A002_X9b92b0_X27e', 'uid___A002_X9b92b0_X457', 'uid___A002_X9b92b0_X22f', 'uid___A002_X9b92b0_Xa5', 'uid___A002_X9b92b0_X56']
for index in range(len(prefix)):
	vis = wd + prefix[index] + '.ms'
	listfile = prefix[index] + '.listobs'
	verbose = True
	listobs()
#


