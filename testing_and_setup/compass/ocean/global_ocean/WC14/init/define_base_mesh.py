import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import jigsaw_to_MPAS.mesh_definition_tools as mdt
from jigsaw_to_MPAS.coastal_tools import signed_distance_from_geojson, \
    mask_from_geojson, distance_from_geojson
from geometric_features import read_feature_collection
import xarray
import matplotlib
import matplotlib.colors as colors
matplotlib.use('Agg')


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
    # To speed up for testing, set following line to 1.0 degrees
    dlon = 0.1
    dlat = dlon
    print('\nCreating cellWidth on a lat-lon grid of: {0:.2f} x {0:.2f} degrees'.format(dlon,dlat))
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
    mdt.register_sci_viz_colormaps()
    #cmapBluesHalf = truncate_colormap(cmapIn='Blues', minval=0.0, maxval=0.7)

    # Create cell width vs latitude for Atlantic and Pacific basins
    QU1 = np.ones(lat.size)
    EC60to30 = mdt.EC_CellWidthVsLat(lat)
    EC60to30Narrow = mdt.EC_CellWidthVsLat(lat, latPosEq = 8.0, latWidthEq = 3.0)

    # Expand from 1D to 2D
    _, cellWidth = np.meshgrid(lon, EC60to30Narrow)
    plot_cartopy(2, 'narrow EC60to30', cellWidth, '3Wbgy5')
    plotFrame = 3

    # global settings for regionally refines mesh
    highRes = 14.0 #[km]

    fileName = 'region_Central_America'
    transitionWidth = 800.0*km
    transitionOffset = 0.0
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, max_length=0.25)
    mask = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) / (transitionWidth/2.)))
    cellWidth = 30.0 * mask + cellWidth * (1 - mask)

    fileName = 'coastline_CUSP'
    distanceToTransition = 600.0*km
    # transitionWidth is distance from 0.07 to 0.03 of transition within tanh
    transitionWidth = 600.0*km
    transitionOffset = distanceToTransition + transitionWidth/2.0
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, max_length=0.25)
    mask = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) / (transitionWidth/2.)))
    cellWidth = highRes * mask + cellWidth * (1 - mask)
    plot_cartopy(plotFrame, fileName + ' mask', mask, 'Blues')
    plot_cartopy(plotFrame+1, 'cellWidth ', cellWidth, '3Wbgy5')
    plotFrame += 2

    fileName = 'region_Gulf_of_Mexico'
    transitionOffset = 600.0*km
    transitionWidth = 600.0*km
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, max_length=0.25)
    maskSmooth = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) / (transitionWidth/2.)))
    maskSharp = 0.5 * (1 + np.sign(-signedDistance))
    fc = read_feature_collection('land_mask_Mexico.geojson')
    signedDistance = signed_distance_from_geojson(fc, lon, lat, max_length=0.25)
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
    signedDistance = signed_distance_from_geojson(fc, lon, lat, max_length=0.25)
    maskSmoothEast = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) /
                                        (transitionWidth/2.)))

    fc = read_feature_collection('region_Bering_Sea_reduced.geojson')
    signedDistance = signed_distance_from_geojson(fc, lon, lat, max_length=0.25)
    maskSmoothWest = 0.5 * (1 + np.tanh((transitionOffset - signedDistance) /
                                        (transitionWidth / 2.)))

    fc = read_feature_collection('land_mask_Kamchatka.geojson')
    maskWest = mask_from_geojson(fc, lon, lat)
    mask = maskSmoothWest * maskWest + maskSmoothEast * (1 - maskWest)
    cellWidth = highRes * mask + cellWidth * (1 - mask)
    ds = xarray.Dataset()
    ds['maskSmoothEast'] = xarray.DataArray(
        maskSmoothEast, dims=['y', 'x'], coords={'y': lat, 'x': lon})
    ds['maskSmoothWest'] = xarray.DataArray(
        maskSmoothWest, dims=['y', 'x'], coords={'y': lat, 'x': lon})
    ds['maskWest'] = xarray.DataArray(
        maskWest, dims=['y', 'x'], coords={'y': lat, 'x': lon})
    ds['mask'] = xarray.DataArray(
        mask, dims=['y', 'x'], coords={'y': lat, 'x': lon})
    ds['cellWidth'] = xarray.DataArray(
        cellWidth, dims=['y', 'x'], coords={'y': lat, 'x': lon})
    ds.to_netcdf('bering.nc')
    plot_cartopy(plotFrame, fileName + ' mask', mask, 'Blues')
    plot_cartopy(plotFrame+1, 'cellWidth ', cellWidth, '3Wbgy5')
    plotFrame += 2

    fileName = 'region_Arctic_Ocean'
    transitionOffset = 0.0*km
    transitionWidth = 600.0*km
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, max_length=0.25)
    mask = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) / (transitionWidth/2.)))
    cellWidth = highRes * mask + cellWidth * (1 - mask)
    plot_cartopy(plotFrame, fileName + ' mask', mask, 'Blues')
    plot_cartopy(plotFrame+1, 'cellWidth ', cellWidth, '3Wbgy5')
    plotFrame += 2

    fileName = 'region_Gulf_Stream_extension'
    transitionOffset = 0.0*km
    transitionWidth = 600.0*km
    fc = read_feature_collection('{}.geojson'.format(fileName))
    signedDistance = signed_distance_from_geojson(fc, lon, lat, max_length=0.25)
    mask = 0.5 * (1 + np.tanh((transitionOffset-signedDistance) / (transitionWidth/2.)))
    cellWidth = highRes * mask + cellWidth * (1 - mask)
    plot_cartopy(plotFrame, fileName + ' mask', mask, 'Blues')
    plot_cartopy(plotFrame+1, 'cellWidth ', cellWidth, '3Wbgy5')
    plotFrame += 2

    # save signed distance to a file
    # da = xarray.DataArray(signedDistance,
    #                      dims=['y', 'x'],
    #                      coords={'y': lat, 'x': lon},
    #                      name='signedDistance')
    #cw_filename = 'signedDistance.nc'
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
    gl.xlabels_top = False
    gl.xlabels_bottom = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    plt.colorbar(im, shrink=.9)
    plt.title(varName)

def truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100):
    '''truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100)'''
    cmapIn = plt.get_cmap(cmapIn)

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmapIn.name, a=minval, b=maxval),
        cmapIn(np.linspace(minval, maxval, n)))

    return new_cmap

if __name__ == '__main__':
    cellWidthVsLatLon()
