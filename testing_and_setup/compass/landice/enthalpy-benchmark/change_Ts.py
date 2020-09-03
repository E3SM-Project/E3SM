#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser

parser = OptionParser(description='change surface air temperature')
parser.add_option("-f", "--file", dest="input_file", help="the input file")
parser.add_option("-v", "--value", dest="change_value", help="the temperature value in Kelvin")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

data = Dataset(options.input_file, 'r+')

T_new = options.change_value

data.variables['surfaceAirTemperature'][0,:] = float(T_new)

data.close()
