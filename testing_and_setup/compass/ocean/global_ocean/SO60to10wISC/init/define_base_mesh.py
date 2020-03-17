import numpy as np
import jigsaw_to_MPAS.mesh_definition_tools as mdt
from jigsaw_to_MPAS.coastal_tools import signed_distance_from_geojson, \
    compute_cell_width
from geometric_features import read_feature_collection
import xarray

# Uncomment to plot the cell size distribution.
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt


def cellWidthVsLatLon():
    """
    Create cell width array for this mesh on a regular latitude-longitude grid.
    Returns
    -------
       cellWidth : numpy.ndarray
            m x n array, entries are desired cell width in km
       lat : numpy.ndarray
            latitude, vector of length m, with entries between -90 and 90,
            degrees
       lon : numpy.ndarray
            longitude, vector of length n, with entries between -180 and 180,
            degrees
    """
    dlon = 0.1
    dlat = dlon
    nlon = int(360./dlon) + 1
    nlat = int(180./dlat) + 1
    lon = np.linspace(-180., 180., nlon)
    lat = np.linspace(-90., 90., nlat)

    cellWidthSouth = 30. * np.ones((len(lat)))

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

    _, cellWidth = np.meshgrid(lon, cellWidthVsLat)

    # now, add the high-res region
    fc = read_feature_collection('high_res_region.geojson')

    signed_distance = signed_distance_from_geojson(fc, lon, lat,
                                                   max_length=0.25)

    da = xarray.DataArray(signed_distance,
                          dims=['y', 'x'],
                          coords={'y': lat, 'x': lon},
                          name='signed_distance')
    cw_filename = 'signed_distance.nc'
    da.to_netcdf(cw_filename)

    # multiply by 5 because transition_width gets multiplied by 0.2 in
    # compute_cell_width
    # Equivalent to 10 degrees latitude
    trans_width = 5*1100e3
    # The last term compensates for the offset in compute_cell_width.
    # The middle of the transition is ~2.5 degrees (300 km) south of the
    # region boundary to best match previous transition at 48 S. (The mean lat
    # of the boundary is 45.5 S.)
    trans_start = -300e3 - 0.5 * trans_width
    dx_min = 10.

    cellWidth = compute_cell_width(signed_distance, cellWidth, lon,
                                   lat, dx_min, trans_start, trans_width,
                                   restrict_box={'include': [], 'exclude': []})

    # Uncomment to plot the cell size distribution.
    # Lon, Lat = np.meshgrid(lon, lat)
    # ax = plt.subplot(111)
    # plt.pcolormesh(Lon, Lat, cellWidth)
    # plt.colorbar()
    # ax.set_aspect('equal')
    # ax.autoscale(tight=True)
    # plt.tight_layout()
    # plt.savefig('cellWidthVsLat.png', dpi=200)

    return cellWidth, lon, lat
