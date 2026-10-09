#!/usr/bin/env python3

"""
Generate a map file from the ne4pg2 grid to a small, tilted, 2d *rectilinear*
(but not lat-lon) grid over North America. The tgt grid has nx*ny points,
stored with x as the fastest-varying index (gid = iy*nx + ix), and
(lat,lon) of the points are NOT separable, i.e., every tgt point has a
distinct lat and a distinct lon. The grid is meant to be used in tests.

The src grid (ne4pg2) is extracted from an existing ncremap map file
(such as map_ne4pg2_to_10x20_20260112.nc), so no src scrip file is needed.

The weights are linear (barycentric) interpolation weights on the Delaunay
triangulation of the src cell centers on the sphere. We do not use ncremap/TempestRemap,
since the latter cannot generate overlap meshes when the tgt cells are much
smaller than the src cells (as is the case here).

Requires: python3 with netCDF4, numpy and scipy
"""

import argparse
import numpy as np
from netCDF4 import Dataset
from scipy.spatial import ConvexHull

###############################################################################
def to_xyz(lat, lon):
###############################################################################
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    return np.stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)], axis=-1)

###############################################################################
def tgt_geometry(nx, ny, lat0, lon0, dx, dy, rot_deg):
###############################################################################
    """
    Map (i,j) "index space" (nodes are at integers, cell centers at half-integers)
    to (lat,lon) via a rotated, locally-cartesian frame. Returns arrays with
    x (i) as fastest-varying index.
    """
    c, s = np.cos(np.deg2rad(rot_deg)), np.sin(np.deg2rad(rot_deg))
    def ll(i, j):
        u = (i - nx/2)*dx
        v = (j - ny/2)*dy
        lat = lat0 + c*v + s*u
        lon = lon0 + (c*u - s*v) / np.cos(np.deg2rad(lat))
        return lat, lon

    ix, iy = np.meshgrid(np.arange(nx), np.arange(ny))   # shape (ny,nx)
    ix, iy = ix.ravel(), iy.ravel()
    clat, clon = ll(ix+0.5, iy+0.5)
    # counter-clockwise corners
    corners = [(0,0), (1,0), (1,1), (0,1)]
    vlat = np.stack([ll(ix+a, iy+b)[0] for a, b in corners], axis=1)
    vlon = np.stack([ll(ix+a, iy+b)[1] for a, b in corners], axis=1)
    return clat, clon, vlat, vlon

###############################################################################
def tri_area(a, b, c):
###############################################################################
    # Area of spherical triangles with unit-vector vertices (Van Oosterom-Strackee)
    num = np.abs(np.einsum('...i,...i', a, np.cross(b, c)))
    den = 1 + np.einsum('...i,...i', a, b) + np.einsum('...i,...i', b, c) \
            + np.einsum('...i,...i', c, a)
    return 2*np.arctan2(num, den)

###############################################################################
def compute_weights(src_lat, src_lon, tgt_lat, tgt_lon):
###############################################################################
    """
    Return (row,col,S), 0-based, with one barycentric triplet per tgt point.
    """
    pts = to_xyz(src_lat, src_lon)
    tris = ConvexHull(pts).simplices          # Delaunay triangulation on the sphere
    tpts = to_xyz(tgt_lat, tgt_lon)

    A = pts[tris]                             # (ntri,3,3)
    rows, cols, wgts = [], [], []
    for i, p in enumerate(tpts):
        # Solve p ~ w0*A0 + w1*A1 + w2*A2 (central projection onto each triangle plane)
        w = np.linalg.solve(np.transpose(A, (0, 2, 1)), p)   # (ntri,3)
        # Keep triangles containing p, on the same side of the sphere center (w sums > 0)
        inside = np.where(np.all(w > -1e-12, axis=1) & (w.sum(axis=1) > 0))[0]
        assert len(inside)>0, f"tgt point {i} not found in any src triangle"
        k = inside[0]
        wk = np.clip(w[k], 0, None)
        wk /= wk.sum()
        rows += [i]*3; cols += list(tris[k]); wgts += list(wk)
    return np.array(rows), np.array(cols), np.array(wgts)

###############################################################################
def main():
###############################################################################
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--src-map", required=True,
                   help="Existing map file, from which the src (ne4pg2) grid is extracted")
    p.add_argument("--out", required=True, help="Output map file")
    p.add_argument("--nx", type=int, default=10)
    p.add_argument("--ny", type=int, default=20)
    p.add_argument("--lat0", type=float, default=40.0, help="Center latitude")
    p.add_argument("--lon0", type=float, default=257.0, help="Center longitude (deg east)")
    p.add_argument("--dx", type=float, default=2.0, help="Grid spacing along x (deg)")
    p.add_argument("--dy", type=float, default=1.5, help="Grid spacing along y (deg)")
    p.add_argument("--rot", type=float, default=10.0, help="Grid tilt (deg)")
    args = p.parse_args()

    n_b = args.nx*args.ny
    with Dataset(args.src_map) as src:
        g = lambda n: np.array(src[n][:])
        n_a = src.dimensions["n_a"].size
        src_data = {n: g(n) for n in ["xc_a","yc_a","xv_a","yv_a","area_a","mask_a","src_grid_dims"]}

    clat, clon, vlat, vlon = tgt_geometry(args.nx, args.ny, args.lat0, args.lon0,
                                          args.dx, args.dy, args.rot)
    # Area of each (spherical) quad cell as the sum of two spherical triangles
    v = to_xyz(vlat, vlon)
    area_b = tri_area(v[:,0],v[:,1],v[:,2]) + tri_area(v[:,0],v[:,2],v[:,3])

    row, col, S = compute_weights(src_data["yc_a"], src_data["xc_a"], clat, clon)
    n_s = len(S)

    with Dataset(args.out, "w", format="NETCDF4") as ds:
        ds.Title = "Barycentric ne4pg2 to rectilinear (non lat-lon) test map"
        ds.history = "Generated by components/eamxx/scripts/gen-rectilinear-test-map.py " + \
                     f"(nx={args.nx}, ny={args.ny}, lat0={args.lat0}, lon0={args.lon0}, " + \
                     f"dx={args.dx}, dy={args.dy}, rot={args.rot})"
        for n, l in [("src_grid_rank",1),("dst_grid_rank",2),("n_a",n_a),("n_b",n_b),
                     ("nv_a",4),("nv_b",4),("n_s",n_s)]:
            ds.createDimension(n, l)
        def var(name, dtype, dims, data, **atts):
            x = ds.createVariable(name, dtype, dims)
            x[:] = data
            x.setncatts(atts)
        var("src_grid_dims","i4",("src_grid_rank",),src_data["src_grid_dims"])
        var("dst_grid_dims","i4",("dst_grid_rank",),[args.nx,args.ny])
        var("yc_a","f8",("n_a",),src_data["yc_a"],units="degrees")
        var("yc_b","f8",("n_b",),clat,units="degrees")
        var("xc_a","f8",("n_a",),src_data["xc_a"],units="degrees")
        var("xc_b","f8",("n_b",),clon,units="degrees")
        var("yv_a","f8",("n_a","nv_a"),src_data["yv_a"],units="degrees")
        var("yv_b","f8",("n_b","nv_b"),vlat,units="degrees")
        var("xv_a","f8",("n_a","nv_a"),src_data["xv_a"],units="degrees")
        var("xv_b","f8",("n_b","nv_b"),vlon,units="degrees")
        var("area_a","f8",("n_a",),src_data["area_a"],units="steradians")
        var("area_b","f8",("n_b",),area_b,units="steradians")
        var("mask_a","i4",("n_a",),src_data["mask_a"],units="unitless")
        var("mask_b","i4",("n_b",),np.ones(n_b,dtype=int),units="unitless")
        var("frac_a","f8",("n_a",),np.zeros(n_a),units="unitless")
        var("frac_b","f8",("n_b",),np.ones(n_b),units="unitless")
        var("row","i4",("n_s",),row+1,first_index=1)
        var("col","i4",("n_s",),col+1,first_index=1)
        var("S","f8",("n_s",),S)

if __name__ == "__main__":
    main()
