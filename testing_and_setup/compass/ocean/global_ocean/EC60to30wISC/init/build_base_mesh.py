#!/usr/bin/env python
import numpy as np
import mpas_tools.mesh.creation.mesh_definition_tools as mdt
from mpas_tools.ocean import build_spherical_mesh


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
    lat = np.arange(-90, 90.01, 0.1)
    lon = np.arange(-180, 180.01, 10.0)

    cellWidthVsLat = mdt.EC_CellWidthVsLat(lat)
    cellWidth = np.outer(cellWidthVsLat, np.ones([1, lon.size]))

    return cellWidth, lon, lat


def main():
    cellWidth, lon, lat = cellWidthVsLatLon()
    build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc')


if __name__ == '__main__':
    main()
