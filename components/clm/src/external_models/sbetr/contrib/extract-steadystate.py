#!/usr/bin/env python
"""FIXME: A nice python program to do something useful.

Author: Ben Andre <andre@ucar.edu>

"""

from __future__ import print_function

import sys

if sys.hexversion < 0x02070000:
    print(70 * "*")
    print("ERROR: {0} requires python >= 2.7.x. ".format(sys.argv[0]))
    print("It appears that you are running python {0}".format(
        ".".join(str(x) for x in sys.version_info[0:3])))
    print(70 * "*")
    sys.exit(1)

#
# built-in modules
#
import argparse
import os
import traceback

if sys.version_info[0] == 2:
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser

#
# installed dependencies
#
import numpy as np
from scipy.io import netcdf

#
# other modules in this package
#

# -------------------------------------------------------------------------------
#
# User input
#
# -------------------------------------------------------------------------------

def commandline_options():
    """Process the command line arguments.

    """
    parser = argparse.ArgumentParser(
        description='FIXME: python program template.')

    parser.add_argument('--backtrace', action='store_true',
                        help='show exception backtraces as extra debugging '
                        'output')

    parser.add_argument('--debug', action='store_true',
                        help='extra debugging output')

    parser.add_argument('--data-dir', nargs=1, default=['.'],
                        help='path to directory containing the data')

    parser.add_argument('--config', nargs=1, required=True,
                        help='path to config file')

    options = parser.parse_args()
    return options


def read_config_file(filename):
    """Read the configuration file and process

    """
    print("Reading configuration file : {0}".format(filename))

    cfg_file = os.path.abspath(filename)
    if not os.path.isfile(cfg_file):
        raise RuntimeError("Could not find config file: {0}".format(cfg_file))

    config = config_parser()
    config.read(cfg_file)

    return config

# -------------------------------------------------------------------------------
#
# FIXME: work functions
#
# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
#
# main
#
# -------------------------------------------------------------------------------

def main(options):
    config = read_config_file(options.config[0])
    data_dir = options.data_dir[0]
    input_filename = config.get('netcdf', 'input')
    input_filename = os.path.join(os.path.abspath(data_dir), input_filename)
    print('Opening input file {0}'.format(input_filename))
    infile = netcdf.netcdf_file(input_filename, 'r', mmap=False)
    
    output_filename = config.get('netcdf', 'output')
    output_filename = os.path.join(os.path.abspath(data_dir), output_filename)
    print('Opening output file {0}'.format(output_filename))
    outfile = netcdf.netcdf_file(output_filename, 'w', mmap=False)
    
    size_list = config.items('sizes')
    sizes = {}
    for name, value in size_list:
        sizes[name] = int(value)
    print('sizes = {0}'.format(sizes))

    for dim in sizes:
        outfile.createDimension(dim, sizes[dim])

    time_in = infile.variables['time']
    time_out = outfile.createVariable('time', time_in.typecode(), ('time', ))
    time_out[:] = 0.0

    for d in ('2d', '3d'):
        dimensions = config.get(d, 'dimensions').split()
        variables = config.get(d, 'variables').split()
        print('{0} : dimensions = {1}'.format(d, dimensions))
        print('{0} : variables = {1}'.format(d, variables))
        for var in variables:
            variable_in = infile.variables[var]
            variable_out = outfile.createVariable(var, variable_in.typecode(), dimensions)
            variable_out._attributes.update(variable_in._attributes)
            if '3d' == d:
                variable_out.data[0, :, 0] = variable_in[0, :, 0]
            elif '2d' in d:
                variable_out.data[0, 0] = variable_in.data[0, 0]
                print('{0} in : {1}'.format(var, variable_in.__dict__))
                print('{0} out : {1}'.format(var, variable_out.__dict__))
            
    infile.close()
    outfile.close()
    return 0


if __name__ == "__main__":
    options = commandline_options()
    try:
        status = main(options)
        sys.exit(status)
    except Exception as error:
        print(str(error))
        if options.backtrace:
            traceback.print_exc()
        sys.exit(1)
