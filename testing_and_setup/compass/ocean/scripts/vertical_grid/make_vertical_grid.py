#!/usr/bin/env python
"""
This script creates the vertical grid for MPAS-Ocean and writes it to a netcdf
file.
"""
# import modules
# {{{
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import argparse
from scipy.optimize import root_scalar
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
        '-nz', '--num_vert_levels', dest='nz',
        default=64,
        help='Number of vertical levels for the grid',
        type=int)
    parser.add_argument(
        '-bd', '--bottom_depth', dest='bottom_depth',
        default=5000.0,
        help='bottom depth for the chosen vertical coordinate [m]',
        type=float)
    parser.add_argument(
        '-dz1', '--min_layer_thickness', dest='dz1',
        default=2.0,
        help='Target thickness of the first layer [m]',
        type=float)
    parser.add_argument(
        '-dz2', '--maxLayer_thickness', dest='dz2',
        default=250.0,
        help='Target maximum thickness in column [m]',
        type=float)
    parser.add_argument(
        '-p', '--plot_vertical_grid', dest='plot_vertical_grid',
        action='store_true')
    parser.add_argument(
        '-o', '--output_file_name', dest='output_file_name',
        default='MPAS-Ocean_vertical_grid.nc',
        help='MPAS file name for output of vertical grid.',
        metavar='NAME')
    args = parser.parse_args()

    create_vertical_grid(args.bottom_depth, args.nz, args.dz1,
                         args.dz2, args.plot_vertical_grid,
                         outfile=args.output_file_name)
# }}}


def create_vertical_grid(
        num_vert_levels=64,
        bottom_depth=5000.0,
        min_layer_thickness=2.0,
        max_layer_thickness=250.0,
        plot_vertical_grid=False,
        outfile='MPAS-Ocean_vertical_grid.nc'):
    # {{{
    """
    This function creates the vertical grid for MPAS-Ocean and writes it to a
    NetCDF file.

    Parameters
    ----------

    bottom_depth : float, optional
        bottom depth for the chosen vertical coordinate [m]

    num_vert_levels : int, optional
        Number of vertical levels for the grid

    min_layer_thickness : float, optional
        Target thickness of the first layer [m]

    max_layer_thickness : float, optional
        Target maximum thickness in column [m]

    plot_vertical_grid : bool, optional
        Whether to plot the vertical grid

    outfile : str, optional
        MPAS file name for output of vertical grid
    """

    print('Creating mesh with ', num_vert_levels, ' layers...')
    nz = num_vert_levels
    dz1 = min_layer_thickness
    dz2 = max_layer_thickness
    # open a new netCDF file for writing.
    ncfile = Dataset(outfile, 'w')
    # create the depth_t dimension.
    ncfile.createDimension('nVertLevels', nz)

    refBottomDepth = ncfile.createVariable(
        'refBottomDepth', np.dtype('float64').char, ('nVertLevels',))
    refMidDepth = ncfile.createVariable(
        'refMidDepth', np.dtype('float64').char, ('nVertLevels',))
    refLayerThickness = ncfile.createVariable(
        'refLayerThickness', np.dtype('float64').char, ('nVertLevels',))

    # the bracket here is large enough that it should hopefully encompass any
    # reasonable value of delta, the characteristic length scale over which
    # dz varies.  The args are passed on to the match_bottom function below,
    # and the root finder will determine a value of delta (sol.root) such that
    # match_bottom is within a tolerance of zero, meaning the bottom of the
    # coordinate computed by cumsum_z hits bottom_depth almost exactly
    sol = root_scalar(match_bottom, method='brentq',
                      bracket=[dz1, 10 * bottom_depth],
                      args=(nz, dz1, dz2, bottom_depth))

    delta = sol.root
    layerThickness, z = cumsum_z(delta, nz, dz1, dz2)
    nVertLevels = nz
    botDepth = -z[1:]
    midDepth = -0.5 * (z[0:-1] + z[1:])

    refBottomDepth[:] = botDepth
    refMidDepth[:] = midDepth
    refLayerThickness[:] = layerThickness[:nVertLevels]
    ncfile.close()

    if plot_vertical_grid:
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
            'number layers: {}\n'.format(nz) + \
            'bottom depth requested:  {:8.2f}\n'.format(bottom_depth) +  \
            'bottom depth actual:     {:8.2f}\n'.format(np.amax(botDepth)) +  \
            'min thickness reqeusted: {:8.2f}\n'.format(min_layer_thickness) + \
            'min thickness actual:    {:8.2f}\n'.format(np.amin(layerThickness[:])) + \
            'max thickness reqeusted: {:8.2f}\n'.format(max_layer_thickness) + \
            'max thickness actual:    {:8.2f}'.format(np.amax(layerThickness[:]))
        print(txt)
        plt.subplot(2, 2, 4)
        plt.text(0, 0, txt, fontsize=12)
        plt.axis('off')
        plt.savefig('vertical_grid.png')

# }}}


def match_bottom(delta, nz, dz1, dz2, bottom_depth):
    """
    Compute the difference between the bottom depth computed with the given
    parameters and the target ``bottom_depth``, used in the root finding
    algorithm to determine which value of ``delta`` to use.

    Parameters
    ----------
    delta : float
        The characteristic length scale over which dz varies (this parameter
        will be optimized to hit a target depth in a target number of layers)

    nz : int
        The number of layers

    dz1 : float
        The layer thickness at the top of the ocean (z = 0)

    dz2 : float
        The layer thickness at z --> -infinity

    bottom_depth: float
        depth of the bottom of the ocean that should match the bottom layer
        interface.  Note: the bottom_depth is positive, whereas the layer
        interfaces are negative.

    Returns
    -------
    diff : float
        The computed bottom depth minus the target ``bottom_depth``.  ``diff``
        should be zero when we have found the desired ``delta``.
    """
    _, z = cumsum_z(delta, nz, dz1, dz2)
    diff = -bottom_depth - z[-1]
    return diff


def cumsum_z(delta, nz, dz1, dz2):
    """
    Compute layer interface depths and layer thicknesses over ``nz`` layers

    Parameters
    ----------
    delta : float
        The characteristic length scale over which dz varies (this parameter
        will be optimized to hit a target depth in a target number of layers)

    nz : int
        The number of layers

    dz1 : float
        The layer thickness at the top of the ocean (z = 0)

    dz2 : float
        The layer thickness at z --> -infinity

    Returns
    -------
    dz : numpy.ndarray
        The layer thicknesses for each layer

    z : numpy.ndarray
        The depth (positive up) of each layer interface (``nz + 1`` total
        elements)
    """
    dz = np.zeros(nz)
    z = np.zeros(nz + 1)
    for zindex in range(nz):
        dz[zindex] = dz_z(z[zindex], dz1, dz2, delta)
        z[zindex + 1] = z[zindex] - dz[zindex]
    return dz, z


def dz_z(z, dz1, dz2, delta):
    """
    layer thickness as a funciton of depth

    Parameters
    ----------
    z : float
        Depth coordinate (positive up) at which to find the layer thickness

    dz1 : float
        The layer thickness at the top of the ocean (z = 0)

    dz2 : float
        The layer thickness at z --> -infinity

    delta : float
        The characteristic length scale over which dz varies (this parameter
        will be optimized to hit a target depth in a target numer of layers)

    Returns
    -------
    dz : float
        The layer thickness
    """
    return (dz2 - dz1) * np.tanh(-z * np.pi / delta) + dz1


if __name__ == '__main__':
    # If called as a primary module, run main
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
