#!/usr/bin/env python
"""Generate plots comparing the betr reactive transport solutions with
the compariable analytical solution.

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
import matplotlib.pyplot as plt

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

    parser.add_argument('--datafile', required=True,
                        help='path to the input data file.')

    # parser.add_argument('--config', nargs=1, required=True,
    #                    help='path to config file')

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
# FIXME: do something
#
# -------------------------------------------------------------------------------
def extract_from_netcdf(filename):
    """
    """
    with netcdf.netcdf_file(filename, 'r', mmap=False) as f:
        for v in f.variables:
            print(v)
            print('    {0} = {1}'.format(
                f.variables[v].dimensions, f.variables[v].shape))
            if len(f.variables[v].shape) == 3:
                ts = 1
                lev = ':'
                col = 0
                value = f.variables[v].data[ts, :, col]
                print('    ({ts}, {lev}, {col}) = {value}'.format(
                    ts=ts, lev=lev, col=col, value=value))

                with plt.xkcd():
                    ts = np.arange(1, 51, 10)
                    value = f.variables[v].data[ts, :, col]
                    print(value)
                    fig = plt.figure()
                    plt.xlabel('z soil * 10 [m]')
                    plt.ylabel('{0}'.format(v))

                    plt.plot(np.transpose(value))
                    plt.pause(10)
                    


# -------------------------------------------------------------------------------
#
# main
#
# -------------------------------------------------------------------------------
def main(options):
    # config = read_config_file(options.config[0])
    filename = options.datafile
    print('Reading {0}'.format(filename))
    extract_from_netcdf(filename)
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
