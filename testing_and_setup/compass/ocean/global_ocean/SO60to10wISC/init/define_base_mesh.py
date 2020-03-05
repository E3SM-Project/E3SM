# /usr/bin/env python
"""
% Create cell width array for this mesh on a regular latitude-longitude grid.
% Outputs:
%    cellWidth - m x n array, entries are desired cell width in km
%    lat - latitude, vector of length m, with entries between -90 and 90, degrees
%    lon - longitude, vector of length n, with entries between -180 and 180, degrees
"""
import numpy as np
import jigsaw_to_MPAS.mesh_definition_tools as mdt

# Uncomment to plot the cell size distribution.
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt


def cellWidthVsLatLon():
    lat = np.arange(-90, 90.01, 0.1)
    # Note that longitude step is 10 degrees, but should only be used if mesh does not 
    # vary with longitude. Otherwise, set to 0.1 degrees.
    lon = np.arange(-180, 180.01, 10.0)

    # define uniform distributions
    cellWidth10 = 10.0 * np.ones(lat.size)
    cellWidth30 = 30.0 * np.ones(lat.size)

    # Southern transition
    latTransition = -48.0
    latWidthTransition = 10.0
    cellWidthSouth = mdt.mergeCellWidthVsLat(
        lat,
        cellWidth10,
        cellWidth30,
        latTransition,
        latWidthTransition)

    # Transition at Equator
    cellWidthNorth = mdt.EC_CellWidthVsLat(lat)
    latTransition = 0.0
    latWidthTransition = 5.0
    cellWidthVsLat = mdt.mergeCellWidthVsLat(
        lat,
        cellWidthSouth,
        cellWidthNorth,
        latTransition,
        latWidthTransition)

    # Uncomment to plot the cell size distribution.
    #plt.plot(lat,cellWidthVsLat)
    #plt.grid(True)
    #plt.xlabel('latitude')
    #plt.ylabel('grid cell size')
    #plt.title('SO60to10wISC, transition at 55S')
    #plt.savefig('cellWidthVsLat.pdf')
    #plt.savefig('cellWidthVsLat.png')

    cellWidth = np.ones((lat.size, lon.size))
    for i in range(lon.size):
        cellWidth[:, i] = cellWidthVsLat

    return cellWidth, lon, lat
