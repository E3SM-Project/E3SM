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
    dlon = 1.0
    dlat = dlon
    #earth_radius = constants['SHR_CONST_REARTH']
    earth_radius = 6378206.4
    print('\nCreating cellWidth on a lat-lon grid of: {0:.2f} x {0:.2f} '
          'degrees'.format(dlon,dlat))
    print('This can be set higher for faster test generation\n')
    nlon = int(360. / dlon) + 1
    nlat = int(180. / dlat) + 1
    lon = np.linspace(-180., 180., nlon)
    lat = np.linspace(-90., 90., nlat)
    km = 1.0e3

    print('plotting ...')
    fig = plt.figure()
    plt.clf()
    fig.set_size_inches(10.0, 10.0)
    register_sci_viz_colormaps()

    # Create cell width vs latitude for Atlantic and Pacific basins
    QU1 = np.ones(lat.size)
    EC60to30 = mdt.EC_CellWidthVsLat(lat)
    RRS30to10 = mdt.RRS_CellWidthVsLat(lat, 30, 10)
    AtlNH = RRS30to10
    AtlVsLat = mdt.mergeCellWidthVsLat(lat, EC60to30, AtlNH, 0, 6)
    PacNH = mdt.mergeCellWidthVsLat(lat, 30 * QU1, RRS30to10, 50, 10)
    PacVsLat = mdt.mergeCellWidthVsLat(lat, EC60to30, PacNH, 0, 6)

    # Expand from 1D to 2D
    _, AtlGrid = np.meshgrid(lon, AtlVsLat)
    _, PacGrid = np.meshgrid(lon, PacVsLat)

    # Signed distance of Atlantic region
    fc = read_feature_collection('Atlantic_region.geojson')
    signedDistance = signed_distance_from_geojson(fc, lon, lat, earth_radius,
                                                  max_length=0.25)

    # Merge Atlantic and Pacific distrubutions smoothly
    transitionWidth = 500.0*km
    maskSmooth = 0.5 * (1 + np.tanh(signedDistance / transitionWidth))
    cellWidthSmooth = PacGrid * maskSmooth + AtlGrid * (1 - maskSmooth)

    # Merge Atlantic and Pacific distrubutions with step function
    maskSharp = 0.5 * (1 + np.sign(signedDistance))
    cellWidthSharp = PacGrid * maskSharp + AtlGrid * (1 - maskSharp)

    # Create a land mask that is 1 over land
    fc = read_feature_collection('Americas_land_mask.geojson')
    Americas_land_mask = mask_from_geojson(fc, lon, lat)
    fc = read_feature_collection('Europe_Africa_land_mask.geojson')
    Europe_Africa_land_mask = mask_from_geojson(fc, lon, lat)
    landMask = np.fmax(Americas_land_mask, Europe_Africa_land_mask)

    # Merge: step transition over land, smooth transition over water
    cellWidth = cellWidthSharp * landMask + cellWidthSmooth * (1 - landMask)

    ax = plt.subplot(4, 2, 1)
    ax.plot(lat, AtlVsLat, label='Atlantic')
    ax.plot(lat, PacVsLat, label='Pacific')
    ax.grid(True)
    plt.title('Grid cell size [km] versus latitude')
    plt.legend()

    varNames = [
        'signedDistance',
        'maskSmooth',
        'cellWidthSmooth',
        'maskSharp',
        'cellWidthSharp',
        'landMask',
        'cellWidth']
    j = 2
    for varName in varNames:
        plot_cartopy(j, varName, vars()[varName], '3Wbgy5')
        j += 1
    fig.canvas.draw()
    plt.tight_layout()

    plt.savefig('mesh_construction.png')

    return cellWidth, lon, lat


def plot_cartopy(nPlot, varName, var, map_name):
    ax = plt.subplot(4, 2, nPlot, projection=ccrs.PlateCarree())
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
