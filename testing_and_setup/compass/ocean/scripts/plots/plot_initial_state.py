#!/usr/bin/env python
"""
This script creates plots of the initial condition.
"""
# import modules
# {{{
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
# }}}


def main():
    # parser
    # {{{
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-i', '--input_file_name', dest='input_file_name',
        default='initial_state.nc',
        help='MPAS file name for input of initial state.')
    parser.add_argument(
        '-o', '--output_file_name', dest='output_file_name',
        default='initial_state.png',
        help='File name for output image.')
    args = parser.parse_args()
    # }}}

    # load mesh variables
    ncfile = Dataset(args.input_file_name, 'r')
    nCells = ncfile.dimensions['nCells'].size
    nEdges = ncfile.dimensions['nEdges'].size
    nVertLevels = ncfile.dimensions['nVertLevels'].size

    fig = plt.figure()
    fig.set_size_inches(16.0, 8.0)
    plt.clf()

    print('plotting histograms of the initial condition')
    txt = \
        'number cells: {}\n'.format(nCells) + \
        'number cells, millions: {:6.3f}\n'.format(nCells / 1.e6) + \
        'number layers: {}\n'.format(nVertLevels)
    print(txt)
    plt.subplot(2, 2, 1)
    plt.text(0, 0, txt, fontsize=12)
    plt.axis('off')

    maxLevelCell = ncfile.variables['maxLevelCell']
    plt.subplot(2, 2, 2)
    plt.hist(maxLevelCell, bins=nVertLevels - 4)
    plt.xlabel('maxLevelCell')
    plt.ylabel('frequency')

    bottomDepth = ncfile.variables['bottomDepth']
    plt.subplot(2, 2, 4)
    plt.hist(bottomDepth, bins=nVertLevels - 4)
    plt.xlabel('bottomDepth [m]')
    plt.ylabel('frequency')

    rx1Edge = ncfile.variables['rx1Edge']
    plt.subplot(2, 2, 3)
    plt.hist(np.squeeze(rx1Edge[0, :]))
    plt.xlabel('Haney Number')
    plt.ylabel('frequency')
    plt.title('Haney Number, max={:4.2f}'.format(
        np.max(np.squeeze(rx1Edge[0, :, :]))))

    plt.savefig(args.output_file_name)


if __name__ == '__main__':
    # If called as a primary module, run main
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
