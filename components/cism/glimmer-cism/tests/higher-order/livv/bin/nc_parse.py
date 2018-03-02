import sys
import subprocess


def emptycheck(checkpath):
        noplot = 0
	file = 'temp.txt'
	comline = "ncdump -c " + checkpath + "> temp.txt"
	subprocess.call(comline, shell = True)

	input = open("temp.txt", 'r')
	line = ''

	for line in input:
		if line.find("// (0 currently)") != -1:
			noplot = 1
        return noplot


checkpath = '/tmp/work/ab3/higher-order/reg_test/gis_10km/data/gis_10km.seacism.nc'

noplot = emptycheck(checkpath)



if noplot == 1:
        print "no data"
else:
        print "data exists"
