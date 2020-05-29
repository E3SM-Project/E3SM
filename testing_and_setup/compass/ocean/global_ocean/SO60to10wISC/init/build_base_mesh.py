#!/usr/bin/env python
import numpy as np
import mpas_tools.mesh.creation.mesh_definition_tools as mdt
from mpas_tools.mesh.creation.signed_distance import \
    signed_distance_from_geojson
from geometric_features import read_feature_collection
from mpas_tools.ocean import build_spherical_mesh
from mpas_tools.cime.constants import constants


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
    #earth_radius = constants['SHR_CONST_REARTH']
    earth_radius = 6371.0e3
    nlon = int(360./dlon) + 1
    nlat = int(180./dlat) + 1
    lon = np.linspace(-180., 180., nlon)
    lat = np.linspace(-90., 90., nlat)

    cellWidthSouth = mdt.EC_CellWidthVsLat(lat, cellWidthEq=30.,
                                           cellWidthMidLat=45.,
                                           cellWidthPole=45.,
                                           latPosEq=7.5, latWidthEq=3.0)

    # Transition at Equator
    cellWidthNorth = mdt.EC_CellWidthVsLat(lat, cellWidthEq=30.,
                                           cellWidthMidLat=60.,
                                           cellWidthPole=60.,
                                           latPosEq=7.5, latWidthEq=3.0)
    latTransition = 0.0
    latWidthTransition = 2.5
    cellWidthVsLat = mdt.mergeCellWidthVsLat(
        lat,
        cellWidthSouth,
        cellWidthNorth,
        latTransition,
        latWidthTransition)

    _, cellWidth = np.meshgrid(lon, cellWidthVsLat)

    # now, add the high-res region
    fc = read_feature_collection('high_res_region.geojson')

    so_signed_distance = signed_distance_from_geojson(fc, lon, lat,
                                                      earth_radius,
                                                      max_length=0.25)

    # Equivalent to 20 degrees latitude
    trans_width = 1600e3
    trans_start = -500e3
    dx_min = 12.

    weights = 0.5 * (1 + np.tanh((so_signed_distance - trans_start) /
                                 trans_width))

    cellWidth = dx_min * (1 - weights) + cellWidth * weights

    fc = read_feature_collection('north_mid_res_region.geojson')

    ar_signed_distance = signed_distance_from_geojson(fc, lon, lat,
                                                      earth_radius,
                                                      max_length=0.25)

    fc = read_feature_collection('greenland.geojson')

    gr_signed_distance = signed_distance_from_geojson(fc, lon, lat,
                                                      earth_radius,
                                                      max_length=0.25)

    frac = (-ar_signed_distance/(-ar_signed_distance + gr_signed_distance))

    frac = np.maximum(0., np.minimum(1., frac))

    dx_min = 15.
    dx_max = 35.

    arctic_widths = dx_max + (dx_min - dx_max)*frac

    trans_width = 1000e3
    trans_start = 0.

    weights = 0.5 * (1 + np.tanh((ar_signed_distance - trans_start) /
                                 trans_width))

    cellWidth = arctic_widths * (1 - weights) + cellWidth * weights

    return cellWidth, lon, lat


def main():
    cellWidth, lon, lat = cellWidthVsLatLon()
    build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc')


if __name__ == '__main__':
    main()
