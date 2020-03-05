#!/usr/bin/env python
"""
% Create cell width array for this mesh on a regular latitude-longitude grid.
% Outputs:
%    cellWidth - m x n array, entries are desired cell width in km
%    lat - latitude, vector of length m, with entries between -90 and 90, degrees
%    lon - longitude, vector of length n, with entries between -180 and 180, degrees
"""
import numpy as np
import jigsaw_to_MPAS.mesh_definition_tools as mdt


def cellWidthVsLatLon():

    ddeg = 0.1

    lat = np.arange(-90, 90.01, ddeg)
    lon = np.arange(-180, 180.01, ddeg)

    cellWidthSouth = 15 * np.ones(lat.size)
    cellWidthNorth = 60 * np.ones(lat.size)
    latTransition = -30
    latWidthTransition = 5

    cellWidthVsLat = mdt.mergeCellWidthVsLat(
        lat,
        cellWidthSouth,
        cellWidthNorth,
        latTransition,
        latWidthTransition)

    cellWidth = np.ones((lat.size, lon.size))
    for i in range(lon.size):
        cellWidth[:, i] = cellWidthVsLat

    return cellWidth, lon, lat
