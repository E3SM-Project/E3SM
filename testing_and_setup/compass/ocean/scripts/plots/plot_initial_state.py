#!/usr/bin/env python
"""
This script creates histogram plots of the initial condition.
"""
# import modules
from netCDF4 import Dataset
import numpy as np
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


def main():
    # parser
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

    # load mesh variables
    ncfile = Dataset(args.input_file_name, 'r')
    nCells = ncfile.dimensions['nCells'].size
    nEdges = ncfile.dimensions['nEdges'].size
    nVertLevels = ncfile.dimensions['nVertLevels'].size

    fig = plt.figure()
    fig.set_size_inches(16.0, 12.0)
    plt.clf()

    print('plotting histograms of the initial condition')
    print('see: init/initial_state/initial_state.png')
    d = datetime.datetime.today()
    txt = \
        'MPAS-Ocean initial state\n' + \
        'date: {}\n'.format(d.strftime('%m/%d/%Y')) + \
        'number cells: {}\n'.format(nCells) + \
        'number cells, millions: {:6.3f}\n'.format(nCells / 1.e6) + \
        'number layers: {}\n\n'.format(nVertLevels) + \
        '  min val   max val  variable name\n'

    plt.subplot(3, 3, 2)
    varName = 'maxLevelCell'
    var = ncfile.variables[varName]
    maxLevelCell = var[:]
    plt.hist(var, bins=nVertLevels - 4)
    plt.ylabel('frequency')
    plt.xlabel(varName)
    txt = '{}{:9.2e} {:9.2e} {}\n'.format(txt, np.amin(var), np.amax(var), varName)

    plt.subplot(3, 3, 3)
    varName = 'bottomDepth'
    var = ncfile.variables[varName]
    plt.hist(var, bins=nVertLevels - 4)
    plt.xlabel(varName)
    txt = '{}{:9.2e} {:9.2e} {}\n'.format(txt, np.amin(var), np.amax(var), varName)

    cellsOnEdge = ncfile.variables['cellsOnEdge']
    cellMask = np.zeros((nCells, nVertLevels), bool)
    edgeMask = np.zeros((nEdges, nVertLevels), bool)
    for k in range(nVertLevels):
        cellMask[:, k] = k < maxLevelCell
        cell0 = cellsOnEdge[:, 0]-1
        cell1 = cellsOnEdge[:, 1]-1
        edgeMask[:, k] = np.logical_and(np.logical_and(cellMask[cell0, k],
                                                       cellMask[cell1, k]),
                                        np.logical_and(cell0 >= 0,
                                                       cell1 >= 0))

    plt.subplot(3, 3, 4)
    varName = 'temperature'
    var = ncfile.variables[varName][0, :, :][cellMask]
    plt.hist(var, bins=100, log=True)
    plt.ylabel('frequency')
    plt.xlabel(varName)
    txt = '{}{:9.2e} {:9.2e} {}\n'.format(txt, np.amin(var), np.amax(var), varName)

    plt.subplot(3, 3, 5)
    varName = 'salinity'
    var = ncfile.variables[varName][0, :, :][cellMask]
    plt.hist(var, bins=100, log=True)
    plt.xlabel(varName)
    txt = '{}{:9.2e} {:9.2e} {}\n'.format(txt, np.amin(var), np.amax(var), varName)

    plt.subplot(3, 3, 6)
    varName = 'layerThickness'
    var = ncfile.variables[varName][0, :, :][cellMask]
    plt.hist(var, bins=100, log=True)
    plt.xlabel(varName)
    txt = '{}{:9.2e} {:9.2e} {}\n'.format(txt, np.amin(var), np.amax(var), varName)

    rx1Edge = ncfile.variables['rx1Edge']
    plt.subplot(3, 3, 7)
    varName = 'rx1Edge'
    var = ncfile.variables[varName][0, :, :][edgeMask]
    plt.hist(var, bins=100, log=True)
    plt.ylabel('frequency')
    plt.xlabel('Haney Number, max={:4.2f}'.format(
        np.max(rx1Edge[:].ravel())))
    txt = '{}{:9.2e} {:9.2e} {}\n'.format(txt, np.amin(var), np.amax(var), varName)

    font = FontProperties()
    font.set_family('monospace')
    font.set_size(12)
    print(txt)
    plt.subplot(3, 3, 1)
    plt.text(0, 1, txt, verticalalignment='top', fontproperties=font)
    plt.axis('off')

    plt.savefig(args.output_file_name)


if __name__ == '__main__':
    # If called as a primary module, run main
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
