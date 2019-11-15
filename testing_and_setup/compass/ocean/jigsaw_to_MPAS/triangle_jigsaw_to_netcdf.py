#!/usr/bin/env python
"""
Name: triangle_jigsaw_to_netcdf.py
Authors: Phillip J. Wolfram and Matthew Hoffman

Last Modified: 01/12/2018

This script converts Triangle and JIGSAW output into a netCDF file for use with
the MPASMeshConverter.x to produce a valid MPAS mesh.

# Mesh Conversion Steps

## JIGSAW mesh
1. Produce a JIGSAW mesh, e.g., example.msh, from
   https://github.com/dengwirda/jigsaw-geo-matlab
2. `./triangle_jigsaw_to_netcdf.py -m example.msh -s`
3. `MpasMeshConverter.x grid.nc`
4. Final mesh mesh.nc can then be used to create our initial condition files.

## TRIANGLE mesh
1. Produce a TRIANGLE mesh, e.g., produced from
   http://www.netlib.org/voronoi/triangle.zip
2. `./triangle -p example.poly`
3. `./triangle_jigsaw_to_netcdf.py -n example.node -e example.ele`
4. `MpasMeshConverter.x grid.nc`
5. Final mesh mesh.nc can then be used to create our initial condition files.

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np

from netCDF4 import Dataset as NetCDFFile
from jigsaw_to_MPAS.open_msh import readmsh

import argparse

import collections
point = collections.namedtuple('Point', ['x', 'y', 'z'])


def circumcenter(on_sphere, x1, y1, z1, x2, y2, z2, x3, y3, z3):  # {{{
    p1 = point(x1, y1, z1)
    p2 = point(x2, y2, z2)
    p3 = point(x3, y3, z3)
    if on_sphere:
        a = (p2.x - p3.x)**2 + (p2.y - p3.y)**2 + (p2.z - p3.z)**2
        b = (p3.x - p1.x)**2 + (p3.y - p1.y)**2 + (p3.z - p1.z)**2
        c = (p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2

        pbc = a * (-a + b + c)
        apc = b * (a - b + c)
        abp = c * (a + b - c)

        xv = (pbc * p1.x + apc * p2.x + abp * p3.x) / (pbc + apc + abp)
        yv = (pbc * p1.y + apc * p2.y + abp * p3.y) / (pbc + apc + abp)
        zv = (pbc * p1.z + apc * p2.z + abp * p3.z) / (pbc + apc + abp)
    else:
        d = 2 * (p1.x * (p2.y - p3.y) + p2.x *
                 (p3.y - p1.y) + p3.x * (p1.y - p2.y))

        xv = ((p1.x**2 + p1.y**2) * (p2.y - p3.y) + (p2.x**2 + p2.y**2)
              * (p3.y - p1.y) + (p3.x**2 + p3.y**2) * (p1.y - p2.y)) / d
        yv = ((p1.x**2 + p1.y**2) * (p3.x - p2.x) + (p2.x**2 + p2.y**2)
              * (p1.x - p3.x) + (p3.x**2 + p3.y**2) * (p2.x - p1.x)) / d
        zv = 0.0

        # Optional method to use barycenter instead.
        # xv = p1.x + p2.x + p3.x
        # xv = xv / 3.0
        # yv = p1.y + p2.y + p3.y
        # yv = yv / 3.0
    return point(xv, yv, zv)

# }}}


def triangle_to_netcdf(node, ele, output_name):
    on_sphere = False
    grid = NetCDFFile(output_name, 'w', format='NETCDF3_CLASSIC')

    # Get dimensions
    # Get nCells
    cell_info = open(node, 'r')
    nCells = -1  # There is one header line
    for block in iter(lambda: cell_info.readline(), ""):
        if block.startswith("#"):
            continue  # skip comment lines
        nCells = nCells + 1
    cell_info.close()

    # Get vertexDegree and nVertices
    cov_info = open(ele, 'r')
    vertexDegree = 3  # always triangles with Triangle!
    nVertices = -1  # There is one header line
    for block in iter(lambda: cov_info.readline(), ""):
        if block.startswith("#"):
            continue  # skip comment lines
        nVertices = nVertices + 1
    cov_info.close()

    if vertexDegree != 3:
        ValueError("This script can only compute vertices with triangular "
                   "dual meshes currently.")

    grid.createDimension('nCells', nCells)
    grid.createDimension('nVertices', nVertices)
    grid.createDimension('vertexDegree', vertexDegree)

    # Create cell variables and sphere_radius
    xCell_full = np.zeros((nCells,))
    yCell_full = np.zeros((nCells,))
    zCell_full = np.zeros((nCells,))

    cell_info = open(node, 'r')
    cell_info.readline()  # read header
    i = 0
    for block in iter(lambda: cell_info.readline(), ""):
        block_arr = block.split()
        if block_arr[0] == "#":
            continue  # skip comment lines
        xCell_full[i] = float(block_arr[1])
        yCell_full[i] = float(block_arr[2])
        zCell_full[i] = 0.0  # z-position is always 0.0 in a planar mesh
        i = i + 1
    cell_info.close()

    grid.on_a_sphere = "NO"
    grid.sphere_radius = 0.0

    cellsOnVertex_full = np.zeros(
        (nVertices, vertexDegree), dtype=np.int32)

    cov_info = open(ele, 'r')
    cov_info.readline()  # read header
    iVertex = 0
    for block in iter(lambda: cov_info.readline(), ""):
        block_arr = block.split()
        if block_arr[0] == "#":
            continue  # skip comment lines
        cellsOnVertex_full[iVertex, :] = int(-1)
        # skip the first column, which is the triangle number, and then
        # only get the next 3 columns
        for j in np.arange(0, 3):
            cellsOnVertex_full[iVertex, j] = int(block_arr[j + 1])

        iVertex = iVertex + 1

    cov_info.close()

    # Create vertex variables
    xVertex_full = np.zeros((nVertices,))
    yVertex_full = np.zeros((nVertices,))
    zVertex_full = np.zeros((nVertices,))

    for iVertex in np.arange(0, nVertices):
        cell1 = cellsOnVertex_full[iVertex, 0]
        cell2 = cellsOnVertex_full[iVertex, 1]
        cell3 = cellsOnVertex_full[iVertex, 2]

        x1 = xCell_full[cell1 - 1]
        y1 = yCell_full[cell1 - 1]
        z1 = zCell_full[cell1 - 1]
        x2 = xCell_full[cell2 - 1]
        y2 = yCell_full[cell2 - 1]
        z2 = zCell_full[cell2 - 1]
        x3 = xCell_full[cell3 - 1]
        y3 = yCell_full[cell3 - 1]
        z3 = zCell_full[cell3 - 1]

        pv = circumcenter(on_sphere, x1, y1, z1, x2, y2, z2, x3, y3, z3)
        xVertex_full[iVertex] = pv.x
        yVertex_full[iVertex] = pv.y
        zVertex_full[iVertex] = pv.z

    meshDensity_full = grid.createVariable(
        'meshDensity', 'f8', ('nCells',))

    meshDensity_full[0:nCells] = 1.0

    var = grid.createVariable('xCell', 'f8', ('nCells',))
    var[:] = xCell_full
    var = grid.createVariable('yCell', 'f8', ('nCells',))
    var[:] = yCell_full
    var = grid.createVariable('zCell', 'f8', ('nCells',))
    var[:] = zCell_full
    var = grid.createVariable('xVertex', 'f8', ('nVertices',))
    var[:] = xVertex_full
    var = grid.createVariable('yVertex', 'f8', ('nVertices',))
    var[:] = yVertex_full
    var = grid.createVariable('zVertex', 'f8', ('nVertices',))
    var[:] = zVertex_full
    var = grid.createVariable(
        'cellsOnVertex', 'i4', ('nVertices', 'vertexDegree',))
    var[:] = cellsOnVertex_full

    grid.sync()
    grid.close()


def jigsaw_to_netcdf(msh_filename, output_name, on_sphere):
    grid = NetCDFFile(output_name, 'w', format='NETCDF3_CLASSIC')

    # Get dimensions
    # Get nCells
    msh = readmsh(msh_filename)
    nCells = msh['POINT'].shape[0]

    # Get vertexDegree and nVertices
    vertexDegree = 3  # always triangles with JIGSAW output
    nVertices = msh['TRIA3'].shape[0]

    if vertexDegree != 3:
        ValueError("This script can only compute vertices with triangular "
                   "dual meshes currently.")

    grid.createDimension('nCells', nCells)
    grid.createDimension('nVertices', nVertices)
    grid.createDimension('vertexDegree', vertexDegree)

    # Create cell variables and sphere_radius
    sphere_radius = 6371000
    xCell_full = msh['POINT'][:, 0]
    yCell_full = msh['POINT'][:, 1]
    zCell_full = msh['POINT'][:, 2]
    for cells in [xCell_full, yCell_full, zCell_full]:
        assert cells.shape[0] == nCells, 'Number of anticipated nodes is' \
                                         ' not correct!'
    if on_sphere:
        grid.on_a_sphere = "YES"
        grid.sphere_radius = sphere_radius
    else:
        grid.on_a_sphere = "NO"
        grid.sphere_radius = 0.0

    # Create cellsOnVertex
    cellsOnVertex_full = msh['TRIA3'][:, :3] + 1
    assert cellsOnVertex_full.shape == (nVertices, vertexDegree), \
        'cellsOnVertex_full is not the right shape!'

    # Create vertex variables
    xVertex_full = np.zeros((nVertices,))
    yVertex_full = np.zeros((nVertices,))
    zVertex_full = np.zeros((nVertices,))

    for iVertex in np.arange(0, nVertices):
        cell1 = cellsOnVertex_full[iVertex, 0]
        cell2 = cellsOnVertex_full[iVertex, 1]
        cell3 = cellsOnVertex_full[iVertex, 2]

        x1 = xCell_full[cell1 - 1]
        y1 = yCell_full[cell1 - 1]
        z1 = zCell_full[cell1 - 1]
        x2 = xCell_full[cell2 - 1]
        y2 = yCell_full[cell2 - 1]
        z2 = zCell_full[cell2 - 1]
        x3 = xCell_full[cell3 - 1]
        y3 = yCell_full[cell3 - 1]
        z3 = zCell_full[cell3 - 1]

        pv = circumcenter(on_sphere, x1, y1, z1, x2, y2, z2, x3, y3, z3)
        xVertex_full[iVertex] = pv.x
        yVertex_full[iVertex] = pv.y
        zVertex_full[iVertex] = pv.z

    meshDensity_full = grid.createVariable(
        'meshDensity', 'f8', ('nCells',))

    for iCell in np.arange(0, nCells):
        meshDensity_full[iCell] = 1.0

    del meshDensity_full

    var = grid.createVariable('xCell', 'f8', ('nCells',))
    var[:] = xCell_full
    var = grid.createVariable('yCell', 'f8', ('nCells',))
    var[:] = yCell_full
    var = grid.createVariable('zCell', 'f8', ('nCells',))
    var[:] = zCell_full
    var = grid.createVariable('xVertex', 'f8', ('nVertices',))
    var[:] = xVertex_full
    var = grid.createVariable('yVertex', 'f8', ('nVertices',))
    var[:] = yVertex_full
    var = grid.createVariable('zVertex', 'f8', ('nVertices',))
    var[:] = zVertex_full
    var = grid.createVariable(
        'cellsOnVertex', 'i4', ('nVertices', 'vertexDegree',))
    var[:] = cellsOnVertex_full

    grid.sync()
    grid.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-n",
        "--node",
        dest="node",
        help="input .node file generated by Triangle.",
        metavar="FILE")
    parser.add_argument(
        "-e",
        "--ele",
        dest="ele",
        help="input .ele file generated by Triangle.",
        metavar="FILE")
    parser.add_argument(
        "-m",
        "--msh",
        dest="msh",
        help="input .msh file generated by JIGSAW.",
        metavar="FILE")
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="output file name.",
        metavar="FILE")
    parser.add_argument(
        "-s",
        "--spherical",
        dest="spherical",
        action="store_true",
        default=False,
        help="Determines if the input/output should be spherical or not.")

    options = parser.parse_args()

    if not options.msh:
        if not options.node:
            parser.error("A .node file is required.")

        if not options.ele:
            parser.error("A .ele file is required.")

    options.density = False  # I'm not sure if we can get this or not...
    if not options.density:
        const_dens = True
    else:
        const_dens = False

    if not options.output:
        output_name = "grid.nc"
    else:
        output_name = options.output

    if options.msh:
        on_sphere = options.spherical
    else:
        # These will always be planar meshes for non-JIGSAW inputs
        on_sphere = False

    if options.msh:
        jigsaw_to_netcdf(options.msh, output_name, on_sphere)
    else:
        triangle_to_netcdf(options.node, options.ele, output_name)

# vim: ai ts=4 sts=4 et sw=4 ft=python
