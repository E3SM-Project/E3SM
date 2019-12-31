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


def cellWidthVsLatLon():
    lat = np.arange(-90, 90.01, 0.1)
    lon = np.arange(-180, 180.01, 0.1)

    QU1 = np.ones(lat.size)
    EC60to30 = mdt.EC_CellWidthVsLat(lat)
    RRS30to10 = mdt.RRS_CellWidthVsLat(lat, 30, 10)

    AtlNH = RRS30to10
    AtlGrid = mdt.mergeCellWidthVsLat(lat, EC60to30, AtlNH, 0, 6)

    PacNH = mdt.mergeCellWidthVsLat(lat, 30 * QU1, RRS30to10, 50, 10)
    PacGrid = mdt.mergeCellWidthVsLat(lat, EC60to30, PacNH, 0, 6)

    cellWidth = mdt.AtlanticPacificGrid(lat, lon, AtlGrid, PacGrid)

    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    plt.clf()
    plt.plot(lat, AtlGrid, label='Atlantic')
    plt.plot(lat, PacGrid, label='Pacific')
    plt.grid(True)
    plt.xlabel('latitude')
    plt.title('Grid cell size, km')
    plt.legend()
    plt.savefig('cellWidthVsLat.png')

    return cellWidth, lon, lat
