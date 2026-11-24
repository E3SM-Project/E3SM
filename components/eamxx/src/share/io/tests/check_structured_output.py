#!/usr/bin/env python3

"""
Checks that output on a structured grid (lat-lon or curvilinear) matches expected values.

The tgt grid in the map file is an (ny,nx) array, with nx the fastest varying index
(dst_grid_dims=[nx,ny]). The remapped fields are compared against a mat-vec computed
independently, directly from the map file. The output file is expected to have
  - lat-lon grids: dims (lat,lon), with 1d lat(lat) and lon(lon) coordinate variables
  - curvilinear grids: dims (lat,lon) too (like ncremap does), but of length ny and nx, and 2d
    lat(lat,lon) and lon(lat,lon) variables
"""

import argparse, sys, pathlib
from netCDF4 import Dataset
import numpy as np

###############################################################################
def parse_command_line(args, description):
###############################################################################
    parser = argparse.ArgumentParser(
        usage="""\n{0} <param>=<val> [<param>=<val>] ...
OR
{0} --help
""".format(pathlib.Path(args[0]).name),
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--map-file",required=True,type=str,
        help="The path to the map file used to generate the output",
    )
    parser.add_argument("--out-file",required=True,type=str,
        help="The name of the output file to check",
    )
    parser.add_argument("--grid-kind",required=True,type=str,choices=["latlon","curvilinear"],
        help="The kind of structured grid the map file tgt grid is expected to be",
    )

    return parser.parse_args(args[1:])

###############################################################################
def check_structured_output(map_file,out_file,grid_kind):
###############################################################################
    with Dataset(map_file, 'r') as map_ds:
        xc_a = map_ds.variables['xc_a'][:]
        yc_a = map_ds.variables['yc_a'][:]
        xc_b = map_ds.variables['xc_b'][:]
        yc_b = map_ds.variables['yc_b'][:]
        row  = map_ds.variables['row'][:] - 1   # map files are 1-based
        col  = map_ds.variables['col'][:] - 1
        S    = map_ds.variables['S'][:]
        nx, ny = [int(n) for n in map_ds.variables['dst_grid_dims'][:]]  # x is the fastest varying

    ncols = nx*ny
    assert xc_b.size==ncols, "nx*ny does not match the tgt grid size in the map file"

    # Independent reference mat-vec
    def matvec(src):
        tgt = np.zeros(ncols)
        np.add.at(tgt,row,S*src[col])
        return tgt.reshape(ny,nx)

    lat2d = yc_b.reshape(ny,nx)
    lon2d = xc_b.reshape(ny,nx)

    # Sanity check: make sure the map file is of the kind we are expecting, otherwise the test would be moot
    is_latlon = len(np.unique(np.round(lat2d,4)[:,0]))==ny and \
                np.allclose(lat2d,lat2d[:,:1],atol=1e-3) and np.allclose(lon2d,lon2d[:1,:],atol=1e-3)
    assert is_latlon==(grid_kind=="latlon"), \
        f"The tgt grid in the map file is {'' if is_latlon else 'not '}a lat-lon grid, but --grid-kind={grid_kind}"
    if grid_kind=="curvilinear":
        # nx should not be the number of lons, and ny should not be the number of lats
        assert len(np.unique(np.round(xc_b,4)))>nx and len(np.unique(np.round(yc_b,4)))>ny

    # Names of the dims in the file (slow dim first), which are the same for both kinds of grids
    ydim, xdim = "lat", "lon"
    # Dims that must NOT be in the file
    bad_dims = ["ncol","ny","nx","x","y"]

    ok = True
    with Dataset(out_file, 'r') as out_ds:
        for d,n in [(ydim,ny),(xdim,nx)]:
            if d not in out_ds.dimensions or out_ds.dimensions[d].size!=n:
                print(f"Dimension '{d}' is missing, or has the wrong size (expected {n})")
                ok = False
        for d in bad_dims:
            if d in out_ds.dimensions:
                print(f"Dimension '{d}' should NOT be in the output file")
                ok = False
        if not ok:
            return False

        for name,src in [('f_lat',yc_a),('f_lon',xc_a)]:
            var = out_ds.variables[name]
            if var.dimensions!=('time',ydim,xdim):
                print(f"Variable '{name}' has wrong dims: {var.dimensions}")
                ok = False
                continue
            expected = matvec(src)
            if not np.allclose(var[0,:,:], expected, rtol=1e-5, atol=1e-4):
                print(f"'{name}' does not match the expected remapped values")
                print(var[0,:,:]); print(expected)
                ok = False

        # The geo data lat/lon should have the proper layout, CF units, and match the map file coordinates
        if grid_kind=="latlon":
            geo = [('lat',('lat',),lat2d[:,0]),('lon',('lon',),lon2d[0,:])]
        else:
            geo = [('lat',('lat','lon'),lat2d),('lon',('lat','lon'),lon2d)]
        for name,dims,expected in geo:
            if name not in out_ds.variables:
                print(f"Geo variable '{name}' not found in the output file")
                ok = False
                continue
            var = out_ds.variables[name]
            units = {'lat':'degrees_north','lon':'degrees_east'}[name]
            if var.dimensions!=dims:
                print(f"Geo variable '{name}' has wrong dims: {var.dimensions}")
                ok = False
            elif getattr(var,'units','') != units:
                print(f"Geo variable '{name}' has non-CF units: {getattr(var,'units','')}")
                ok = False
            elif not np.allclose(var[:], expected, rtol=1e-6, atol=1e-4):
                print(f"Geo variable '{name}' does not match the map file tgt coordinates")
                ok = False

    return ok

###############################################################################
def _main_func(description):
###############################################################################
    success = check_structured_output(**vars(parse_command_line(sys.argv, description)))
    sys.exit(0 if success else 1)

###############################################################################

if (__name__ == "__main__"):
    _main_func(__doc__)
