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
    dlon = 1.0
    dlat = dlon
    nlon = int(360./dlon) + 1
    nlat = int(180./dlat) + 1
    lon = np.linspace(-180., 180., nlon)
    lat = np.linspace(-90., 90., nlat)

    QU1 = np.ones(lat.size)

# change back later:
    #EC60to30 = mdt.EC_CellWidthVsLat(lat)
    #RRS30to10 = mdt.RRS_CellWidthVsLat(lat, 30, 10)
    #AtlNH = RRS30to10
    #AtlGrid = mdt.mergeCellWidthVsLat(lat, EC60to30, AtlNH, 0, 6)
    #PacNH = mdt.mergeCellWidthVsLat(lat, 30 * QU1, RRS30to10, 50, 10)
    #PacGrid = mdt.mergeCellWidthVsLat(lat, EC60to30, PacNH, 0, 6)

    AtlVsLat= QU1
    PacVsLat= 2*QU1

    _, AtlGrid = np.meshgrid(lon, AtlVsLat)
    _, PacGrid = np.meshgrid(lon, PacVsLat)

    fc = read_feature_collection('Atlantic_region.geojson')

    signed_distance = signed_distance_from_geojson(fc, lon, lat,
                                                   max_length=0.25)

    mask = 0.5*(1+np.tanh(signed_distance/1000.0e3)) 
    
    cellWidth = AtlGrid*mask + PacGrid*(1-mask)

    # save signed distance to a file
    #da = xarray.DataArray(signed_distance,
    #                      dims=['y', 'x'],
    #                      coords={'y': lat, 'x': lon},
    #                      name='signed_distance')
    #cw_filename = 'signed_distance.nc'
    #da.to_netcdf(cw_filename)

# replace this:
    #cellWidth = mdt.AtlanticPacificGrid(lat, lon, AtlGrid, PacGrid)

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
