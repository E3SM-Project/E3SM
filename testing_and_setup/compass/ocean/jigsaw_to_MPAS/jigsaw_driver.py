from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import jigsawpy


def jigsaw_driver(cellWidth, x, y, on_sphere=True, geom_points=None, geom_edges=None):
    '''
    A function for building a jigsaw mesh

    Parameters
    ----------
    cellWidth : ndarray
        The size of each cell in the resulting mesh as a function of space

    x, y : ndarray
        The x and y coordinates of each point in the cellWidth array (lon and lat for spherical mesh)

    on_sphere : logical
        Whether this mesh is spherical or planar

    geom_points : list of point coordinates for bounding polygon for planar mesh

    geom_edges : list of edges between points in geom_points that define the bounding polygon

    '''
    # Authors
    # -------
    # Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    # setup files for JIGSAW
    opts = jigsawpy.jigsaw_jig_t()
    opts.geom_file = 'mesh.msh'
    opts.jcfg_file = 'mesh.jig'
    opts.mesh_file = 'mesh-MESH.msh'
    opts.hfun_file = 'mesh-HFUN.msh'

    # save HFUN data to file
    hmat = jigsawpy.jigsaw_msh_t()
    if on_sphere:
       hmat.mshID = 'ELLIPSOID-GRID'
       hmat.xgrid = numpy.radians(x)
       hmat.ygrid = numpy.radians(y)
    else:
       hmat.mshID = 'EUCLIDEAN-GRID'
       hmat.xgrid = x
       hmat.ygrid = y
    hmat.value = cellWidth
    jigsawpy.savemsh(opts.hfun_file, hmat)

    # define JIGSAW geometry
    geom = jigsawpy.jigsaw_msh_t()
    if on_sphere:
       geom.mshID = 'ELLIPSOID-MESH'
       geom.radii = 6371.*numpy.ones(3, float)
    else:
       geom.mshID = 'EUCLIDEAN-MESH'
       geom.vert2 = geom_points
       geom.edge2 = geom_edges
       #geom.edge2.index = geom_edges
       print (geom_points)
    jigsawpy.savemsh(opts.geom_file, geom)

    # build mesh via JIGSAW!
    mesh = jigsawpy.jigsaw_msh_t()
    opts.hfun_scal = 'absolute'
    opts.hfun_hmax = float("inf")
    opts.hfun_hmin = 0.0
    opts.mesh_dims = +2  # 2-dim. simplexes
    opts.optm_qlim = 0.9375
    opts.verbosity = +1

    jigsawpy.cmd.jigsaw(opts)
