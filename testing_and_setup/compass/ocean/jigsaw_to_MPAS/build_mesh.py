#!/usr/bin/env python
"""
This script performs the first step of initializing the global ocean.  This
includes:
Step 1. Build cellWidth array as function of latitude and longitude
Step 2. Build mesh using JIGSAW
Step 3. Convert triangles from jigsaw format to netcdf
Step 4. Convert from triangles to MPAS mesh
Step 5. Create vtk file for visualization
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import scipy.io as sio
import subprocess
import os
import xarray
import argparse

from mpas_tools.conversion import convert
from mpas_tools.io import write_netcdf


from jigsaw_to_MPAS.jigsaw_driver import jigsaw_driver
from jigsaw_to_MPAS.triangle_jigsaw_to_netcdf import jigsaw_to_netcdf
from jigsaw_to_MPAS.inject_bathymetry import inject_bathymetry
from jigsaw_to_MPAS.inject_meshDensity import inject_meshDensity
from jigsaw_to_MPAS.inject_preserve_floodplain import \
    inject_preserve_floodplain

from define_base_mesh import define_base_mesh


def build_mesh(preserve_floodplain=False, floodplain_elevation=20.0,
               do_inject_bathymetry=False):

    print('Step 1. Build cellWidth array as function of latitude and '
          'longitude')
    cellWidth, lon, lat = define_base_mesh.cellWidthVsLatLon()
    sio.savemat(
        'cellWidthVsLatLon.mat', {
            'cellWidth': cellWidth, 'lon': lon, 'lat': lat})

    print('Step 2. Generate mesh with JIGSAW')
    jigsaw_driver(cellWidth, lon, lat)

    print('Step 3. Convert triangles from jigsaw format to netcdf')
    jigsaw_to_netcdf(msh_filename='mesh-MESH.msh',
                     output_name='mesh_triangles.nc', on_sphere=True)

    print('Step 4. Convert from triangles to MPAS mesh')
    write_netcdf(convert(xarray.open_dataset('mesh_triangles.nc')),
                 'base_mesh.nc')

    print('Step 5. Inject correct meshDensity variable into base mesh file')
    inject_meshDensity(mat_filename='cellWidthVsLatLon.mat',
                       mesh_filename='base_mesh.nc')

    if do_inject_bathymetry:
        print('Step 6. Injecting bathymetry')
        inject_bathymetry(mesh_file='base_mesh.nc')

    if preserve_floodplain:
        print('Step 7. Injecting flag to preserve floodplain')
        inject_preserve_floodplain(mesh_file='base_mesh.nc',
                                   floodplain_elevation=floodplain_elevation)

    print('Step 8. Create vtk file for visualization')
    args = ['paraview_vtk_field_extractor.py',
            '--ignore_time',
            '-l',
            '-d', 'maxEdges=0',
            '-v', 'allOnCells',
            '-f', 'base_mesh.nc',
            '-o', 'base_mesh_vtk']
    print("running", ' '.join(args))
    subprocess.check_call(args, env=os.environ.copy())

    print("***********************************************")
    print("**    The global mesh file is base_mesh.nc   **")
    print("***********************************************")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--preserve_floodplain', action='store_true')
    parser.add_argument('--floodplain_elevation', action='store',
                        type=float, default=20.0)
    parser.add_argument('--inject_bathymetry', action='store_true')
    cl_args = parser.parse_args()
    build_mesh(cl_args.preserve_floodplain, cl_args.floodplain_elevation,
               cl_args.inject_bathymetry)
