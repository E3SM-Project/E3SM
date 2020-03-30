#!/usr/bin/env python
'''
This script creates the vertical grid for MPAS-Ocean and writes it to a netcdf file.
'''
# import modules
# {{{
from netCDF4 import Dataset
import numpy as np
import argparse
# }}}


def main():
    # parser
    # {{{
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-o', '--output_file_name', dest='output_filename_name',
        default='MPAS-Ocean_vertical_grid.nc',
        help='MPAS file name for output of vertical grid.',
        metavar='NAME')
    parser.add_argument(
        '-p', '--plot_vertical_grid', dest='plot_vertical_grid',
        action='store_true')
    parser.add_argument(
        '-bd', '--bottom_depth', dest='bottom_depth',
        default=5000,
        help='bottom depth for the chosen vertical coordinate [m]',
        type=float)
    parser.add_argument(
        '-nz', '--num_vert_levels', dest='nz',
        default=64,
        help='Number of vertical levels for the grid',
        type=int)
    parser.add_argument(
        '-dz1', '--layer1_thickness', dest='dz1',
        default=2.0,
        help='Target thickness of the first layer [m]',
        type=float)
    parser.add_argument(
        '-dz2', '--maxLayer_thickness', dest='dz2',
        default=250.0,
        help='Target maximum thickness in column [m]',
        type=float)
    parser.add_argument(
        '-eps', '--error_tolerance', dest='epsilon',
        default=1e-2,
        help='Threshold for iterations', type=float)
    parser.add_argument(
        '-maxit',
        '--max_iterations',
        dest='maxit',
        default=1000,
        help='maximum number of iterations for grid convergences',
        type=int)
    args = parser.parse_args()

    create_vertical_grid(args.bottom_depth, args.nz, args.dz1,
                         args.dz2, args.plot_vertical_grid,
                         maxit=args.maxit, epsilon=args.epsilon,
                         outFile=args.output_filename_name)
# }}}


def create_vertical_grid(
        bottom_depth,
        nz,
        dz1_in,
        dz2_in,
        plot_vertical_grid,
        maxit,
        epsilon,
        outFile):
    # {{{
    print('Creating mesh with ', nz, ' layers...')
    dz1 = dz1_in
    dz2 = dz2_in
    # open a new netCDF file for writing.
    ncfile = Dataset(outFile, 'w')
    # create the depth_t dimension.
    ncfile.createDimension('nVertLevels', nz)

    refBottomDepth = ncfile.createVariable(
        'refBottomDepth', np.dtype('float64').char, ('nVertLevels'))
    refMidDepth = ncfile.createVariable(
        'refMidDepth', np.dtype('float64').char, ('nVertLevels'))
    refLayerThickness = ncfile.createVariable(
        'refLayerThickness', np.dtype('float64').char, ('nVertLevels'))

    Hmax = bottom_depth
    nLayers = 0

    layerThickness = np.zeros(nz)
    dz = [epsilon]
    z = [0]
    count = 0

    while nLayers != nz and count < maxit:
        zval = -epsilon
        dz = [epsilon]
        z = [0]
        nLayers = 0
        while zval > -Hmax:
            difference = dz_z(zval, Hmax, epsilon, dz2) - zval
            while abs(difference) > 0.3:
                zval -= epsilon
                difference = dz_z(zval, Hmax, epsilon, dz2) - z[nLayers] + zval
            z.append(zval)
            dz.append(dz_z(zval, Hmax, epsilon, dz2))
            nLayers += 1
            zval -= epsilon

        dz_arr = np.asarray(dz)
        ind = abs(dz_arr - dz1).argmin()

        dztemp = dz_arr[ind:]
        nLayers = len(dztemp)
        change = nz - nLayers
        dz2 -= float(change)
        count += 1

    layerThickness[:nLayers] = dztemp
    nVertLevels = nz
    botDepth = np.zeros(nVertLevels)
    midDepth = np.zeros(nVertLevels)
    botDepth[0] = layerThickness[0]
    midDepth[0] = 0.5 * layerThickness[0]

    for i in range(1, nVertLevels):
        botDepth[i] = botDepth[i - 1] + layerThickness[i]
        midDepth[i] = midDepth[i - 1] + 0.5 * \
            (layerThickness[i] + layerThickness[i - 1])

    if count >= maxit:
        print('Error: grid did not converge, adjust parameters')
    else:
        refBottomDepth[:] = botDepth
        refMidDepth[:] = midDepth
        refLayerThickness[:] = layerThickness[:nVertLevels]
        ncfile.close()

    if plot_vertical_grid:
        import matplotlib
        import matplotlib.pyplot as plt
        matplotlib.use('Agg')
        fig = plt.figure()
        fig.set_size_inches(16.0, 8.0)
        zInd = np.arange(1, nVertLevels + 1)
        plt.clf()

        plt.subplot(2, 2, 1)
        plt.plot(zInd, midDepth, '.')
        plt.gca().invert_yaxis()
        plt.xlabel('vertical index (one-based)')
        plt.ylabel('layer mid-depth [m]')
        plt.grid()

        plt.subplot(2, 2, 2)
        plt.plot(layerThickness, midDepth, '.')
        plt.gca().invert_yaxis()
        plt.xlabel('layer thickness [m]')
        plt.ylabel('layer mid-depth [m]')
        plt.grid()

        plt.subplot(2, 2, 3)
        plt.plot(zInd, layerThickness, '.')
        plt.xlabel('vertical index (one-based)')
        plt.ylabel('layer thickness [m]')
        plt.grid()

        txt = \
            'number layers: ' + str(nz) + '\n' + \
            'bottom depth requested:  ' + '{:8.2f}'.format(bottom_depth) + '\n' +  \
            'bottom depth actual:     ' + '{:8.2f}'.format(np.amax(botDepth)) + '\n' +  \
            'min thickness reqeusted: ' + '{:8.2f}'.format(dz1_in) + '\n' + \
            'min thickness actual:    ' + '{:8.2f}'.format(np.amin(layerThickness[:])) + '\n' + \
            'max thickness reqeusted: ' + '{:8.2f}'.format(dz2_in) + '\n' + \
            'max thickness actual:    ' + '{:8.2f}'.format(np.amax(layerThickness[:]))
        print(txt)
        plt.subplot(2, 2, 4)
        plt.text(0, 0, txt, fontsize=12)
        plt.axis('off')
        plt.savefig('vertical_grid.png')

# }}}


def dz_z(z, Hmax, epsilon, dzmax):
    return dzmax * np.tanh(-z * np.pi / Hmax) + epsilon


if __name__ == '__main__':
    # If called as a primary module, run main
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
