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

import xarray
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy

from mpas_tools.conversion import convert
from mpas_tools.io import write_netcdf
from mpas_tools.viz.paraview_extractor import extract_vtk

from jigsaw_to_MPAS.jigsaw_driver import jigsaw_driver
from jigsaw_to_MPAS.triangle_jigsaw_to_netcdf import jigsaw_to_netcdf
from jigsaw_to_MPAS.inject_bathymetry import inject_bathymetry
from jigsaw_to_MPAS.inject_meshDensity import inject_meshDensity
from jigsaw_to_MPAS.inject_preserve_floodplain import \
    inject_preserve_floodplain
from jigsaw_to_MPAS.mesh_definition_tools import register_sci_viz_colormaps

import define_base_mesh

def build_mesh(
        preserve_floodplain=False,
        floodplain_elevation=20.0,
        do_inject_bathymetry=False,
        geometry='sphere',
        plot_cellWidth=True):

    if geometry == 'sphere':
        on_sphere = True
    else:
        on_sphere = False

    print('Step 1. Build cellWidth array as function of horizontal coordinates')
    if on_sphere:
        cellWidth, lon, lat = define_base_mesh.cellWidthVsLatLon()
        da = xarray.DataArray(cellWidth,
                              dims=['lat', 'lon'],
                              coords={'lat': lat, 'lon': lon},
                              name='cellWidth')
        cw_filename = 'cellWidthVsLatLon.nc'
        da.to_netcdf(cw_filename)
        plot_cellWidth = True
        if plot_cellWidth:
            register_sci_viz_colormaps()
            fig = plt.figure(figsize=[16.0, 8.0])
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.set_global()
            im = ax.imshow(cellWidth, origin='lower',
                           transform=ccrs.PlateCarree(),
                           extent=[-180, 180, -90, 90], cmap='3Wbgy5',
                           zorder=0)
            ax.add_feature(cartopy.feature.LAND, edgecolor='black', zorder=1)
            gl = ax.gridlines(
                crs=ccrs.PlateCarree(),
                draw_labels=True,
                linewidth=1,
                color='gray',
                alpha=0.5,
                linestyle='-', zorder=2)
            gl.xlabels_top = False
            gl.ylabels_right = False
            plt.title(
                'Grid cell size, km, min: {:.1f} max: {:.1f}'.format(
                cellWidth.min(),cellWidth.max()))
            plt.colorbar(im, shrink=.60)
            fig.canvas.draw()
            plt.tight_layout()
            plt.savefig('cellWidthGlobal.png', bbox_inches='tight', dpi=300)
            plt.close()

    else:
        cellWidth, x, y, geom_points, geom_edges = define_base_mesh.cellWidthVsXY()
        da = xarray.DataArray(cellWidth,
                              dims=['y', 'x'],
                              coords={'y': y, 'x': x},
                              name='cellWidth')
        cw_filename = 'cellWidthVsXY.nc'
        da.to_netcdf(cw_filename)

    print('Step 2. Generate mesh with JIGSAW')
    if on_sphere:
        jigsaw_driver(cellWidth, lon, lat)
    else:
        jigsaw_driver(
            cellWidth,
            x,
            y,
            on_sphere=False,
            geom_points=geom_points,
            geom_edges=geom_edges)

    print('Step 3. Convert triangles from jigsaw format to netcdf')
    jigsaw_to_netcdf(msh_filename='mesh-MESH.msh',
                     output_name='mesh_triangles.nc', on_sphere=on_sphere)

    print('Step 4. Convert from triangles to MPAS mesh')
    write_netcdf(convert(xarray.open_dataset('mesh_triangles.nc')),
                 'base_mesh.nc')

    print('Step 5. Inject correct meshDensity variable into base mesh file')
    inject_meshDensity(cw_filename=cw_filename,
                       mesh_filename='base_mesh.nc', on_sphere=on_sphere)

    if do_inject_bathymetry:
        print('Step 6. Injecting bathymetry')
        inject_bathymetry(mesh_file='base_mesh.nc')

    if preserve_floodplain:
        print('Step 7. Injecting flag to preserve floodplain')
        inject_preserve_floodplain(mesh_file='base_mesh.nc',
                                   floodplain_elevation=floodplain_elevation)

    print('Step 8. Create vtk file for visualization')
    extract_vtk(ignore_time=True, lonlat=True, dimension_list=['maxEdges='],
                variable_list=['allOnCells'], filename_pattern='base_mesh.nc',
                out_dir='base_mesh_vtk')

    print("***********************************************")
    print("**    The global mesh file is base_mesh.nc   **")
    print("***********************************************")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--preserve_floodplain', action='store_true')
    parser.add_argument('--floodplain_elevation', action='store',
                        type=float, default=20.0)
    parser.add_argument('--inject_bathymetry', action='store_true')
    parser.add_argument('--geometry', default='sphere')
    parser.add_argument('--plot_cellWidth', action='store_true')
    cl_args = parser.parse_args()
    build_mesh(cl_args.preserve_floodplain, cl_args.floodplain_elevation,
               cl_args.inject_bathymetry, cl_args.geometry,
               cl_args.plot_cellWidth)
