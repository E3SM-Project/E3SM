from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import jigsawpy


def jigsaw_driver(cellWidth, lon, lat):
    '''
    A function for building a jigsaw mesh

    Parameters
    ----------
    cellWidth : ndarray
        The size of each cell in the resulting mesh as a function of space

    lon, lat : ndarray
        The lon. and lat. of each point in the cellWidth array
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
    hmat.mshID = 'ELLIPSOID-GRID'
    hmat.xgrid = numpy.radians(lon)
    hmat.ygrid = numpy.radians(lat)
    hmat.value = cellWidth
    jigsawpy.savemsh(opts.hfun_file, hmat)
   
    # define JIGSAW geometry
    geom = jigsawpy.jigsaw_msh_t()
    geom.mshID = 'ELLIPSOID-MESH'
    geom.radii = 6371.*numpy.ones(3, float)
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
