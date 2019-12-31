#!/usr/bin/env python
'''
name: mesh_definition_tools
authors: Mark R. Petersen

last modified: 07/09/2018

These groups of functions are used to define the cellWidth variable on
regular lat/lon grids.  The cellWidth variable is a jigsaw input that
defines the mesh.
'''
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')


##########################################################################
# Functions
##########################################################################


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

    Author: Mark Petersen
    Los Alamos National Laboratory
    March 2018 # Last revision: 4/20/2018
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


def EC_CellWidthVsLat(lat):
    '''
    EC_CellWidthVsLat - Create Eddy Closure spacing as a function of lat.
    This is inted as part of the workflow to make an MPAS global mesh.

    Syntax: cellWidthOut = EC_CellWidthVsLat(lat, cellWidthEq, cellWidthMidLat, cellWidthPole,
                                             latPosEq, latPosPole, latTransition,
                                             latWidthEq, latWidthPole)
    Inputs:
       lat - vector of length n, with entries between -90 and 90, degrees

    Optional inputs:
       Default values for Cell width, km
       cellWidthEq = 30.0 # Eq is equatorial lat
       cellWidthMidLat = 60.0 # MidLat is mid lat
       cellWidthPole = 35.0 # Pole is polar lat

       Default values for lat positions in degrees
       latPosEq = 15.0 # position of center of transition region
       latPosPole = 73.0 # position of center of transition region
       latTransition = 40 # lat to change from Eq to Pole function
       latWidthEq = 6.0 # width of transition region
       latWidthPole = 9.0 # width of transition region

    Outputs:
       cellWidthOut - vector of length n, entrie are cell width as a function of lat

    Example:
       EC60to30 = EC_CellWidthVsLat(lat)
       EC120to60 = EC_CellWidthVsLat(lat,60,120,70)

    Author: Mark Petersen
    Los Alamos National Laboratory
    March 2018 # Last revision: 4/20/2018
    '''

    # Default values for Cell width, km
    cellWidthEq = 30.0  # Eq is equatorial lat
    cellWidthMidLat = 60.0  # MidLat is mid lat
    cellWidthPole = 35.0  # Pole is polar lat

    # Default values for lat positions in degrees
    latPosEq = 15.0  # position of center of transition region
    latPosPole = 73.0  # position of center of transition region
    latTransition = 40  # lat to change from Eq to Pole function
    latWidthEq = 6.0  # width of transition region
    latWidthPole = 9.0  # width of transition region

    # try
    #  cellWidthEq =     varargin{1} #
    #  cellWidthMidLat = varargin{2} #
    #  cellWidthPole =   varargin{3} #
    #  latPosEq =        varargin{4} #
    #  latPosPole =      varargin{5} #
    #  latTransition =   varargin{6} #
    #  latWidthEq =      varargin{7} #
    #  latWidthPole =    varargin{8} #

    minCellWidth = min(cellWidthEq, min(cellWidthMidLat, cellWidthPole))
    densityEq = (minCellWidth / cellWidthEq)**4
    densityMidLat = (minCellWidth / cellWidthMidLat)**4
    densityPole = (minCellWidth / cellWidthPole)**4
    densityEC = np.zeros(lat.shape)
    cellWidthOut = np.zeros(lat.shape)
    for j in range(lat.size):
        if np.abs(lat[j]) < latTransition:
            densityEC[j] = ((densityEq - densityMidLat) * (1.0 + np.tanh(
                (latPosEq - np.abs(lat[j])) / latWidthEq)) / 2.0) + densityMidLat
        else:
            densityEC[j] = ((densityMidLat - densityPole) * (1.0 + np.tanh(
                (latPosPole - np.abs(lat[j])) / latWidthPole)) / 2.0) + densityPole
        cellWidthOut[j] = minCellWidth / densityEC[j]**0.25

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

    Author: Mark Petersen
    Los Alamos National Laboratory
    March 2018 # Last revision: 4/20/2018
    '''

    # ratio between high and low resolution
    gamma = (cellWidthPole / cellWidthEq)**4.0

    densityRRS = (1.0 - gamma) * \
        np.power(np.deg2rad(np.sin(np.absolute(lat))), 4.0) + gamma
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
            cellWidthOut[j,i] = cellWidthInPacific[j]
            # test if in Atlantic Basin:
            if lat[j]>65.0:
                if (lon[i]>-150.0) & (lon[i]<170.0):
                    cellWidthOut[j,i] = cellWidthInAtlantic[j]
            elif lat[j]>20.0:
                if (lon[i]>-100.0) & (lon[i]<35.0):
                    cellWidthOut[j,i] = cellWidthInAtlantic[j]
            elif lat[j]>0.0:
                if (lon[i]>-2.0*lat[j]-60.0) & (lon[i]<35.0):
                    cellWidthOut[j,i] = cellWidthInAtlantic[j]
            else:
                if (lon[i]>-60.0) & (lon[i]<20.0):
                    cellWidthOut[j,i] = cellWidthInAtlantic[j]
    return cellWidthOut

#
# def  circleOnGrid(lon, lat, centerLon, centerLat, radius, tanhWidth)
# '''
# circleOnGrid: combine two cell width distributions using a tanh function.
# This is inted as part of the workflow to make an MPAS global mesh.
#
# Syntax: cellWidthOut = circleOnGrid(lat, lon, cellWidthInAtlantic, cellWidthInPacific)
#
# Inputs:
#   lon - vector of length m, with entries between -180, 180, degrees
#   lat - vector of length n, with entries between -90, 90, degrees
#
# Optional inputs:
#
# Outputs:
#   cellWidthOut - m by n array, grid cell width on globe, km
#
# Example:
#   RRS18to6 = circleOnGrid(lat,18,6)
#
# See also:
#
# Author: Mark Petersen
# Los Alamos National Laboratory
# March 2018 # Last revision: 3/27/2018
# '''
#
#cellWidthOut = zeros(lon.size,lat.size))
# for i in range(lon.size)
#  for j in range(lat.size)
#    [dist d2km]=lldistkm([centerLat, centerLon], [lat[j], lon[i]])
#    cellWidthOut(i,j) = 0.5*(-tanh((dist - radius)/tanhWidth) + 1.0)
#
#
#
# def lldistkm(latlon1,latlon2)
# '''
# format: [d1km d2km]=lldistkm(latlon1,latlon2)
# Distance:
# d1km: distance in km based on Haversine formula
# (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
# d2km: distance in km based on Pythagoras theorem
# (see: http://en.wikipedia.org/wiki/Pythagorean_theorem)
# After:
# http://www.movable-type.co.uk/scripts/latlong.html
#
# --Inputs:
#  latlon1: latlon of origin point [lat lon]
#  latlon2: latlon of destination point [lat lon]
#
# --Outputs:
#  d1km: distance calculated by Haversine formula
#  d2km: distance calculated based on Pythagoran theorem
#
# --Example 1, short distance:
#  latlon1=[-43 172]
#  latlon2=[-44  171]
#  [d1km d2km]=distance(latlon1,latlon2)
#  d1km =
#          137.365669065197 (km)
#  d2km =
#          137.368179013869 (km)
#  d1km approximately equal to d2km
#
# --Example 2, longer distance:
#  latlon1=[-43 172]
#  latlon2=[20  -108]
#  [d1km d2km]=distance(latlon1,latlon2)
#  d1km =
#          10734.8931427602 (km)
#  d2km =
#          31303.4535270825 (km)
#  d1km is significantly different from d2km (d2km is not able to work
#  for longer distances).
#
# First version: 15 Jan 2012
# Updated: 17 June 2012
# --------------------------------------------------------------------------
# '''
#
# radius=6371
# lat1=latlon1(1)*pi/180
# lat2=latlon2(1)*pi/180
# lon1=latlon1(2)*pi/180
# lon2=latlon2(2)*pi/180
# deltaLat=lat2-lat1
# deltaLon=lon2-lon1
#a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2
# c=2*atan2(sqrt(a),sqrt(1-a))
# d1km=radius*c #    Haversine distance
#
# x=deltaLon*cos((lat1+lat2)/2)
# y=deltaLat
# d2km=radius*sqrt(x*x + y*y) # Pythagoran distance
#
#

##############################################################
