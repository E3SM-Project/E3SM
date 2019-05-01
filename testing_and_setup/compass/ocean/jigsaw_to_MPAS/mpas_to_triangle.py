#!/usr/bin/env python
'''
Script to convert from MPAS netCDF format to the Triangle format:
https://www.cs.cmu.edu/~quake/triangle.node.html
https://www.cs.cmu.edu/~quake/triangle.ele.html

Only works for planar meshes.
'''
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import sys
from netCDF4 import Dataset as NetCDFFile
from optparse import OptionParser


def mpas_to_triangle(mpasfile, trifile):

    fin = NetCDFFile(options.mpasfile, 'r')
    if fin.on_a_sphere == "YES":
        sys.abort("ERROR: This script only works for planar meshes!")

    if len(fin.dimensions['vertexDegree']) != 3:
        sys.abort("ERROR: This script only works for vertexDegree of 3!")

    nCells = len(fin.dimensions['nCells'])
    nVertices = len(fin.dimensions['nVertices'])

    xCell = fin.variables['xCell'][:]
    yCell = fin.variables['yCell'][:]
    ConC = fin.variables['cellsOnCell'][:]
    nConC = fin.variables['nEdgesOnCell'][:]
    ConV = fin.variables['cellsOnVertex'][:]

    # create node file
    fnode = open(options.trifile + ".node", 'w')
    # write node file header:   First line: <# of vertices> <dimension (must
    # be 2)> <# of attributes> <# of boundary markers (0 or 1)>
    fnode.write("{:d} 2 0 1\n".format(nCells))
    # Remaining lines: <vertex #> <x> <y> [attributes] [boundary marker]
    for i in range(nCells):
        if ConC[i, 0:nConC[i]].min() == 0:
            isBdy = 1
        else:
            isBdy = 0
        fnode.write(
            "{:d} {:f} {:f} {:d}\n".format(
                i + 1,
                xCell[i],
                yCell[i],
                isBdy))
    fnode.write("# Generated from MPAS file: {}\n".format(options.mpasfile))
    fnode.close()

    # create ele file
    fele = open(options.trifile + ".ele", "w")

    # calculate number of non-degenerate triangles
    numtri = 0
    for i in range(nVertices):
        if ConV[i, :].min() > 0:
            numtri += 1

    # write ele file header:  First line: <# of triangles> <nodes per
    # triangle> <# of attributes>
    fele.write("{:d} 3 0\n".format(numtri))
    # Remaining lines: <triangle #> <node> <node> <node> ... [attributes]
    cnt = 0
    for i in range(nVertices):
        # write non-generate triangles only
        if ConV[i, :].min() > 0:
            cnt += 1
            fele.write("{:d} {:d} {:d} {:d}\n".format(
                cnt, ConV[i, 0], ConV[i, 1], ConV[i, 2]))
    fele.write("# Generated from MPAS file: {}\n".format(options.mpasfile))
    fele.close()

    fin.close()
    print("Conversion complete.")


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option(
        "-m",
        "--mpas",
        dest="mpasfile",
        help="input MPAS netCDF file.",
        metavar="FILE")
    parser.add_option(
        "-t",
        "--triangle",
        dest="trifile",
        help="output file name template to be in triangle format (FILE.1.node,"
             " FILE.1.ele).",
        metavar="FILE")

    options, args = parser.parse_args()

    if not options.mpasfile:
        parser.error("An input MPAS file is required.")

    if not options.trifile:
        parser.error("A output Triangle format file name is required.")

    mpas_to_triangle(mpasfile=options.mpasfile, trifile=options.trifile)
