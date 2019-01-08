#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import netCDF4 as nc4
import argparse


def inject_preserve_floodplain(mesh_file, floodplain_elevation):

    nc_mesh = nc4.Dataset(mesh_file, 'r+')
    nc_vars = nc_mesh.variables.keys()

    if 'cellSeedMask' not in nc_vars:
        nc_mesh.createVariable('cellSeedMask', 'i', ('nCells'))
    nc_mesh.variables['cellSeedMask'][:] = \
        nc_mesh.variables['bottomDepthObserved'][:] < floodplain_elevation

    nc_mesh.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('mesh_file', action='store', type=str)
    parser.add_argument('floodplain_elevation', action='store', type=float)
    cl_args = parser.parse_args()

    inject_preserve_floodplain(cl_args.mesh_file, cl_args.floodplain_elevation)
