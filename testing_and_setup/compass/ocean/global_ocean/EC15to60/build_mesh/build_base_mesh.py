#!/usr/bin/env python
import numpy as np
import mpas_tools.mesh.creation.mesh_definition_tools as mdt
from mpas_tools.ocean import build_spherical_mesh
from mpas_tools.mesh.creation.signed_distance import \
    signed_distance_from_geojson
from mpas_tools.cime.constants import constants
from geometric_features import read_feature_collection

import xarray


def cellWidthVsLatLon():
    """
    Create cell width array for this mesh on a regular latitude-longitude grid.
    Returns
    -------
       cellWidth : ndarray
            m x n array, entries are desired cell width in km

       lat : ndarray
            latitude, vector of length m, with entries between -90 and 90,
            degrees

       lon : ndarray
            longitude, vector of length n, with entries between -180 and 180,
            degrees
    """

    dlon = 0.1
    dlat = dlon
    nlon = int(360./dlon) + 1
    nlat = int(180./dlat) + 1
    lon = np.linspace(-180., 180., nlon)
    lat = np.linspace(-90., 90., nlat)

    cellWidthVsLat = mdt.EC_CellWidthVsLat(lat)

    _, cellWidth = np.meshgrid(lon, cellWidthVsLat)

    fc = read_feature_collection('north_mid_res_region.geojson')

    earth_radius = constants['SHR_CONST_REARTH']

    mr_signed_distance = signed_distance_from_geojson(fc, lon, lat,
                                                      earth_radius,
                                                      max_length=0.25)

    fc = read_feature_collection('arctic_high_res_region.geojson')

    hr_signed_distance = signed_distance_from_geojson(fc, lon, lat,
                                                      earth_radius,
                                                      max_length=0.25)

    frac = (-mr_signed_distance / (-mr_signed_distance + hr_signed_distance))

    frac = np.maximum(0., np.minimum(1., frac))

    dx_min = 15.
    dx_max = 30.

    arctic_widths = dx_max + (dx_min - dx_max) * frac

    trans_width = 1000e3
    trans_start = 0.

    weights = 0.5 * (1 + np.tanh((mr_signed_distance - trans_start) /
                                 trans_width))

    cellWidth = arctic_widths * (1 - weights) + cellWidth * weights

    return cellWidth, lon, lat


def main():
    cellWidth, lon, lat = cellWidthVsLatLon()
    build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc')


if __name__ == '__main__':
    main()

