#!/usr/bin/env

import sys
import os
from optparse import OptionParser
import subprocess




file_path = '/tmp/work/ab3/higher-order/reg_test/gis_10km/data/gis_10km.JFNK.trilinos.10.gnu.out'

flag = 0


try:
        file = open(file_path, "r")

except:
        print "error reading in file"
        sys.exit(1)
        raise

for line in file:
    if ('FATAL ERROR' in line):
        flag = 1


if flag == 1:
    print 'Fatal Error in output file'


file.close()
