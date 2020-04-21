#!/usr/bin/env python
'''
name: mesh_definition_tools

These functions are tools used to define the cellWidth variable on
regular lat/lon grids.  The cellWidth variable is a jigsaw input that
defines the mesh.
'''
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
import xml.etree.ElementTree as ET
import pkg_resources
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt


##########################################################################
# Functions
##########################################################################


def register_sci_viz_colormaps():
    """Register all SciVisColor colormaps with matplotlib"""

    for mapName in ['3wave-yellow-grey-blue', '3Wbgy5',
                    '4wave-grey-red-green-mgreen', '5wave-yellow-brown-blue',
                    'blue-1', 'blue-3', 'blue-6', 'blue-8', 'blue-orange-div',
                    'brown-2', 'brown-5', 'brown-8', 'green-1', 'green-4',
                    'green-7', 'green-8', 'orange-5', 'orange-6',
                    'orange-green-blue-gray', 'purple-7', 'purple-8', 'red-1',
                    'red-3', 'red-4', 'yellow-1', 'yellow-7']:

        xmlFile = pkg_resources.resource_filename(
            __name__, 'SciVisColorColormaps/{}.xml'.format(mapName))
        _read_xml_colormap(xmlFile, mapName)


def mergeCellWidthVsLat(
        lat,
        cellWidthInSouth,
        cellWidthInNorth,
        latTransition,
        latWidthTransition):
    '''
    mergeCellWidthVsLat: combine two cell width distributions using a tanh function.
    This is inted as part of the workflow to make an MPAS global mesh.

    Syntax: cellWidthOut = mergeCellWidthVsLat(lat, cellWidthInSouth, cellWidthInNorth, latTransition, latWidthTransition)

    Inputs:
       lat - vector of length n, with entries between -90 and 90, degrees
       cellWidthInSouth - vector of length n, first distribution
       cellWidthInNorth - vector of length n, second distribution

    Optional inputs:
       latTransition = 0 # lat to change from cellWidthInSouth to cellWidthInNorth, degrees
       latWidthTransition = 0 # width of lat transition, degrees

    Outputs:
       cellWidthOut - vector of length n, entries are cell width as a function of lat
    '''
    # Assign defaults
    # latTransition = 0 # lat to change from cellWidthInSouth to cellWidthInNorth, degrees
    # latWidthTransition = 0 # width of lat transition, degrees

    cellWidthOut = np.ones(lat.size)
    if latWidthTransition == 0:
        for j in range(lat.size):
            if lat[j] < latTransition:
                cellWidthOut[j] = cellWidthInSouth[j]
            else:
                cellWidthOut[j] = cellWidthInNorth[j]
    else:
        for j in range(lat.size):
            weightNorth = 0.5 * \
                (np.tanh((lat[j] - latTransition) / latWidthTransition) + 1.0)
            weightSouth = 1.0 - weightNorth
            cellWidthOut[j] = weightSouth * cellWidthInSouth[j] + \
                weightNorth * cellWidthInNorth[j]

    return cellWidthOut


def EC_CellWidthVsLat(lat, cellWidthEq=30.0, cellWidthMidLat=60.0,
                      cellWidthPole=35.0, latPosEq=15.0, latPosPole=73.0,
                      latTransition=40.0, latWidthEq=6.0, latWidthPole=9.0):
    """
    Create Eddy Closure spacing as a function of lat. This is intended as part
    of the workflow to make an MPAS global mesh.

    Parameters
    ----------
    lat : numpy.ndarray
       vector of length n, with entries between -90 and 90, degrees

    cellWidthEq : float, optional
       Cell width in km at the equator

    cellWidthMidLat : float, optional
       Cell width in km at mid latitudes

    cellWidthPole : float, optional
       Cell width in km at the poles

    latPosEq : float, optional
       Latitude in degrees of center of the equatorial transition region

    latPosPole : float, optional
       Latitude in degrees of center of the polar transition region

    latTransition : float, optional
       Latitude in degrees of the change from equatorial to polar function

    latWidthEq : float, optional
       Width in degrees latitude of the equatorial transition region

    latWidthPole : float, optional
       Width in degrees latitude of the polar transition region

    Returns
    -------

    cellWidthOut : numpy.ndarray
         1D array of same length as ``lat`` with entries that are cell width as
         a function of lat

    Examples
    --------
    Default

    >>> EC60to30 = EC_CellWidthVsLat(lat)

    Half the default resolution:

    >>> EC120to60 = EC_CellWidthVsLat(lat, cellWidthEq=60., cellWidthMidLat=120., cellWidthPole=70.)
    """

    minCellWidth = min(cellWidthEq, min(cellWidthMidLat, cellWidthPole))
    densityEq = (minCellWidth / cellWidthEq)**4
    densityMidLat = (minCellWidth / cellWidthMidLat)**4
    densityPole = (minCellWidth / cellWidthPole)**4
    densityEqToMid = ((densityEq - densityMidLat) * (1.0 + np.tanh(
        (latPosEq - np.abs(lat)) / latWidthEq)) / 2.0) + densityMidLat
    densityMidToPole = ((densityMidLat - densityPole) * (1.0 + np.tanh(
        (latPosPole - np.abs(lat)) / latWidthPole)) / 2.0) + densityPole
    mask = np.abs(lat) < latTransition
    densityEC = np.array(densityMidToPole)
    densityEC[mask] = densityEqToMid[mask]
    cellWidthOut = minCellWidth / densityEC**0.25

    return cellWidthOut


def RRS_CellWidthVsLat(lat, cellWidthEq, cellWidthPole):
    '''
    RRS_CellWidthVsLat - Create Rossby Radius Scaling as a function of lat.
    This is inted as part of the workflow to make an MPAS global mesh.

    Syntax: cellWidthOut = RRS_CellWidthVsLat(lat, cellWidthEq, cellWidthPole)

    Inputs:
       lat - vector of length n, with entries between -90 and 90, degrees
       cellWidthEq - Cell width at the equator, km
       cellWidthPole - Cell width at the poles, km

    Outputs:
       RRS_CellWidth - vector of length n, entries are cell width as a function of lat

    Example:
       RRS18to6 = RRS_CellWidthVsLat(lat,18,6)
    '''

    # ratio between high and low resolution
    gamma = (cellWidthPole / cellWidthEq)**4.0

    densityRRS = (1.0 - gamma) * \
        np.power(np.sin(np.deg2rad(np.absolute(lat))), 4.0) + gamma
    cellWidthOut = cellWidthPole / np.power(densityRRS, 0.25)
    return cellWidthOut


def AtlanticPacificGrid(lat, lon, cellWidthInAtlantic, cellWidthInPacific):
    '''
    AtlanticPacificGrid: combine two cell width distributions using a tanh function.

    Inputs:
      lon - vector of length m, with entries between -180, 180, degrees
      lat - vector of length n, with entries between -90, 90, degrees
      cellWidthInAtlantic - vector of length n, cell width in Atlantic as a function of lon, km
      cellWidthInPacific - vector of length n, cell width in Pacific as a function of lon, km

    Optional inputs:

    Outputs:
      cellWidthOut - m by n array, grid cell width on globe, km
    '''
    cellWidthOut = np.zeros((lat.size, lon.size))
    for i in range(lon.size):
        for j in range(lat.size):
            # set to Pacific mask as default
            cellWidthOut[j, i] = cellWidthInPacific[j]
            # test if in Atlantic Basin:
            if lat[j] > 65.0:
                if (lon[i] > -150.0) & (lon[i] < 170.0):
                    cellWidthOut[j, i] = cellWidthInAtlantic[j]
            elif lat[j] > 20.0:
                if (lon[i] > -100.0) & (lon[i] < 35.0):
                    cellWidthOut[j, i] = cellWidthInAtlantic[j]
            elif lat[j] > 0.0:
                if (lon[i] > -2.0 * lat[j] - 60.0) & (lon[i] < 35.0):
                    cellWidthOut[j, i] = cellWidthInAtlantic[j]
            else:
                if (lon[i] > -60.0) & (lon[i] < 20.0):
                    cellWidthOut[j, i] = cellWidthInAtlantic[j]
    return cellWidthOut


def _read_xml_colormap(xmlFile, mapName):
    """Read in an XML colormap"""

    xml = ET.parse(xmlFile)

    root = xml.getroot()
    colormap = root.findall('ColorMap')
    if len(colormap) > 0:
        colormap = colormap[0]
        colorDict = {'red': [], 'green': [], 'blue': []}
        for point in colormap.findall('Point'):
            x = float(point.get('x'))
            color = [float(point.get('r')), float(point.get('g')),
                     float(point.get('b'))]
            colorDict['red'].append((x, color[0], color[0]))
            colorDict['green'].append((x, color[1], color[1]))
            colorDict['blue'].append((x, color[2], color[2]))
        cmap = LinearSegmentedColormap(mapName, colorDict, 256)

        _register_colormap_and_reverse(mapName, cmap)


def _register_colormap_and_reverse(mapName, cmap):
    plt.register_cmap(mapName, cmap)
    plt.register_cmap('{}_r'.format(mapName), cmap.reversed())

##############################################################
