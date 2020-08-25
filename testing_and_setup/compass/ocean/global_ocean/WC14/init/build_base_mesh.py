#!/usr/bin/env python
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from mpas_tools.ocean import build_spherical_mesh
import mpas_tools.mesh.creation.mesh_definition_tools as mdt
from mpas_tools.mesh.creation.signed_distance import \
    signed_distance_from_geojson, mask_from_geojson
from mpas_tools.viz.colormaps import register_sci_viz_colormaps
from mpas_tools.cime.constants import constants
from geometric_features import read_feature_collection

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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
    # To speed up for testing, set following line to 1.0 degrees
    dlon = 0.1
    dlat = dlon
    #earth_radius = constants['SHR_CONST_REARTH']
    earth_radius = 6371.0e3
    print('\nCreating cellWidth on a lat-lon grid of: {0:.2f} x {0:.2f} '
          'degrees'.format(dlon, dlat))
    print('This can be set higher for faster test generation\n')
    nlon = int(360. / dlon) + 1
    nlat = int(180. / dlat) + 1
    lon = np.linspace(-180., 180., nlon)
    lat = np.linspace(-90., 90., nlat)
    km = 1.0e3

    print('plotting ...')
    fig = plt.figure()
    plt.clf()
    fig.set_size_inches(10.0, 14.0)
    register_sci_viz_colormaps()

    # Create cell width vs latitude for Atlantic and Pacific basins
    EC60to30 = mdt.EC_CellWidthVsLat(lat)
    EC60to30Narrow = mdt.EC_CellWidthVsLat(lat, latPosEq=8.0, latWidthEq=3.0)

    # Expand from 1D to 2D
    _, cellWidth = np.meshgrid(lon, EC60to30Narrow)
    plot_cartopy(2, 'narrow EC60to30', cellWidth, '3Wbgy5')
    plotFrame = 3

    # global settings for regionally refines mesh
    highRes = 14.0  # [km]

    fileName = 'region_Central_America'
    transitionWidth = 800.0*km
    transitionOffset = 0.0
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, earth_radius,
                                                  max_length=0.25)
    mask = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) /
                              (transitionWidth/2.)))
    cellWidth = 30.0 * mask + cellWidth * (1 - mask)

    fileName = 'coastline_CUSP'
    distanceToTransition = 600.0*km
    # transitionWidth is distance from 0.07 to 0.03 of transition within tanh
    transitionWidth = 600.0*km
    transitionOffset = distanceToTransition + transitionWidth/2.0
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, earth_radius,
                                                  max_length=0.25)
    mask = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) /
                              (transitionWidth/2.)))
    cellWidth = highRes * mask + cellWidth * (1 - mask)
    plot_cartopy(plotFrame, fileName + ' mask', mask, 'Blues')
    plot_cartopy(plotFrame+1, 'cellWidth ', cellWidth, '3Wbgy5')
    plotFrame += 2

    fileName = 'region_Gulf_of_Mexico'
    transitionOffset = 600.0*km
    transitionWidth = 600.0*km
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, earth_radius,
                                                  max_length=0.25)
    maskSmooth = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) /
                                    (transitionWidth/2.)))
    maskSharp = 0.5 * (1 + np.sign(-signedDistance))
    fc = read_feature_collection('land_mask_Mexico.geojson')
    signedDistance = signed_distance_from_geojson(fc, lon, lat, earth_radius,
                                                  max_length=0.25)
    landMask = 0.5 * (1 + np.sign(-signedDistance))
    mask = maskSharp * landMask + maskSmooth * (1-landMask)
    cellWidth = highRes * mask + cellWidth * (1 - mask)
    plot_cartopy(plotFrame, fileName + ' mask', mask, 'Blues')
    plot_cartopy(plotFrame+1, 'cellWidth ', cellWidth, '3Wbgy5')
    plotFrame += 2

    fileName = 'region_Bering_Sea'
    transitionOffset = 0.0*km
    transitionWidth = 600.0*km
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, earth_radius,
                                                  max_length=0.25)
    maskSmoothEast = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) /
                                        (transitionWidth/2.)))

    fc = read_feature_collection('region_Bering_Sea_reduced.geojson')
    signedDistance = signed_distance_from_geojson(fc, lon, lat, earth_radius,
                                                  max_length=0.25)
    maskSmoothWest = 0.5 * (1 + np.tanh((transitionOffset - signedDistance) /
                                        (transitionWidth / 2.)))

    fc = read_feature_collection('land_mask_Kamchatka.geojson')
    maskWest = mask_from_geojson(fc, lon, lat)
    mask = maskSmoothWest * maskWest + maskSmoothEast * (1 - maskWest)
    cellWidth = highRes * mask + cellWidth * (1 - mask)
    plot_cartopy(plotFrame, fileName + ' mask', mask, 'Blues')
    plot_cartopy(plotFrame+1, 'cellWidth ', cellWidth, '3Wbgy5')
    plotFrame += 2

    fileName = 'region_Arctic_Ocean'
    transitionOffset = 0.0*km
    transitionWidth = 600.0*km
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, earth_radius,
                                                  max_length=0.25)
    mask = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) /
                              (transitionWidth/2.)))
    cellWidth = highRes * mask + cellWidth * (1 - mask)
    plot_cartopy(plotFrame, fileName + ' mask', mask, 'Blues')
    plot_cartopy(plotFrame+1, 'cellWidth ', cellWidth, '3Wbgy5')
    plotFrame += 2

    fileName = 'region_Gulf_Stream_extension'
    transitionOffset = 0.0*km
    transitionWidth = 600.0*km
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, earth_radius,
                                                  max_length=0.25)
    mask = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) /
                              (transitionWidth/2.)))
    cellWidth = highRes * mask + cellWidth * (1 - mask)
    plot_cartopy(plotFrame, fileName + ' mask', mask, 'Blues')
    plot_cartopy(plotFrame+1, 'cellWidth ', cellWidth, '3Wbgy5')
    plotFrame += 2

    # save signed distance to a file
    # da = xarray.DataArray(signedDistance,
    #                      dims=['y', 'x'],
    #                      coords={'y': lat, 'x': lon},
    #                      name='signedDistance')
    # cw_filename = 'signedDistance.nc'
    # da.to_netcdf(cw_filename)

    ax = plt.subplot(6, 2, 1)
    ax.plot(lat, EC60to30, label='original EC60to30')
    ax.plot(lat, EC60to30Narrow, label='narrow EC60to30')
    ax.grid(True)
    plt.title('Grid cell size [km] versus latitude')
    plt.legend(loc="upper left")

    plt.savefig('mesh_construction.png', dpi=300)

    return cellWidth, lon, lat


def plot_cartopy(nPlot, varName, var, map_name):
    ax = plt.subplot(6, 2, nPlot, projection=ccrs.PlateCarree())
    ax.set_global()

    im = ax.imshow(var,
                   origin='lower',
                   transform=ccrs.PlateCarree(),
                   extent=[-180, 180, -90, 90], cmap=map_name,
                   zorder=0)
    ax.add_feature(cfeature.LAND, edgecolor='black', zorder=1)
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=1,
        color='gray',
        alpha=0.5,
        linestyle='-', zorder=2)
    ax.coastlines()
    gl.top_labels = False
    gl.bottom_labels = False
    gl.right_labels = False
    gl.left_labels = False
    plt.colorbar(im, shrink=.9)
    plt.title(varName)


def main():
    cellWidth, lon, lat = cellWidthVsLatLon()
    build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc')


if __name__ == '__main__':
    main()
