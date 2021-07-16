from __future__ import print_function
from netCDF4 import Dataset
import netCDF4
import numpy as np
import os
import sys
import ConfigParser
import math
from scipy.interpolate import griddata
from create_forcing import create_scrip_grid_file, get_mpas_grid_info, create_scrip_file_MPAS, write_scrip_in_file, create_output_times, get_remapping_data

#-------------------------------------------------------------------------------

def latlon_to_xyz(lat, lon):

    x = math.cos(lat) * math.cos(lon)
    y = math.cos(lat) * math.sin(lon)
    z = math.sin(lat)

    return x, y, z

#-------------------------------------------------------------------------------

def xyz_to_latlon(x, y, z):

    lon = 0.0
    if (x != 0.0 or y != 0.0): lon = math.atan2(y,x)

    lat = math.asin(z / math.sqrt(x*x + y*y + z*z))

    return lat, lon

#-------------------------------------------------------------------------------

def create_scrip_file_gx1(filenameScrip, filenameGx1Grid):

    filenameGx1Grid = "/Users/akt/Work/Forcing/gx1/grid_info/global_gx1.nc"
    fileIn = Dataset(filenameGx1Grid,"r")

    nx = len(fileIn.dimensions["nx"])
    ny = len(fileIn.dimensions["ny"])

    ULONin = fileIn.variables["ULON"][:]
    ULATin = fileIn.variables["ULAT"][:]

    KMT = fileIn.variables["KMT"][:]

    fileIn.close()

    nCells = nx * ny
    gridDims = [nx, ny]

    gridImask = np.ones(nCells,dtype="i")

    ULAT = np.zeros((ny+1,nx+1))
    ULON = np.zeros((ny+1,nx+1))

    ULAT[1:,1:] = ULATin[:,:]
    ULON[1:,1:] = ULONin[:,:]

    ULAT[:,0] = ULAT[:,-1]
    ULON[:,0] = ULON[:,-1]

    ULON[0,:] = ULON[1,:]
    ULAT[0,:] = ULAT[1,:] - math.pi / 180.0

    cornerLat = np.zeros((4,nCells))
    cornerLon = np.zeros((4,nCells))

    for i in range(0,nx):
        for j in range(0,ny):

            ii = i + 1
            jj = j + 1

            iCell = ii + nx * (jj-1) - 1

            i1 = ii-1 ; j1 = jj-1
            i2 = ii   ; j2 = jj-1
            i3 = ii   ; j3 = jj
            i4 = ii-1 ; j4 = jj

            cornerLat[0,iCell] = ULAT[j1,i1]
            cornerLat[1,iCell] = ULAT[j2,i2]
            cornerLat[2,iCell] = ULAT[j3,i3]
            cornerLat[3,iCell] = ULAT[j4,i4]

            cornerLon[0,iCell] = ULON[j1,i1]
            cornerLon[1,iCell] = ULON[j2,i2]
            cornerLon[2,iCell] = ULON[j3,i3]
            cornerLon[3,iCell] = ULON[j4,i4]

    centerLat = np.zeros(nCells)
    centerLon = np.zeros(nCells)

    for i in range(0,nx):
        for j in range(0,ny):

            ii = i + 1
            jj = j + 1

            iCell = ii + nx * (jj-1) - 1

            x1,y1,z1 = latlon_to_xyz(cornerLat[0,iCell],cornerLon[0,iCell])
            x2,y2,z2 = latlon_to_xyz(cornerLat[1,iCell],cornerLon[1,iCell])
            x3,y3,z3 = latlon_to_xyz(cornerLat[2,iCell],cornerLon[2,iCell])
            x4,y4,z4 = latlon_to_xyz(cornerLat[3,iCell],cornerLon[3,iCell])

            x0 = 0.25 * (x1 + x2 + x3 + x4)
            y0 = 0.25 * (y1 + y2 + y3 + y4)
            z0 = 0.25 * (z1 + z2 + z3 + z4)

            centerLat[iCell], centerLon[iCell] = xyz_to_latlon(x0, y0, z0)

    create_scrip_grid_file(filenameScrip, nCells, 4, 2, gridDims, centerLat, centerLon, gridImask, cornerLat, cornerLon, "gx1")

#-------------------------------------------------------------------------------

def fill_array(arrayIn):

    nTimes = arrayIn.shape[0]
    nx = arrayIn.shape[1]
    ny = arrayIn.shape[2]

    arrayOut = np.zeros((nTimes,nx,ny))
    arrayOut[:] = arrayIn[:]

    grid_x, grid_y = np.mgrid[0:nx, 0:ny]

    for iTime in range(0,nTimes):

        array = np.zeros((nx,3*ny))

        array[:,   0:  ny] = arrayIn[iTime,:,:]
        array[:,  ny:2*ny] = arrayIn[iTime,:,:]
        array[:,2*ny:3*ny] = arrayIn[iTime,:,:]

        pointsGood = []
        valuesGood = []

        pointsBad = []

        for i in range(0,nx):
            for j in range(0,ny):
                if (array[i,j] > -900.0):
                    pointsGood.append((i,j))
                    valuesGood.append(array[i,j])
                else:
                    pointsBad.append((i,j))

        pointsGood = np.array(pointsGood)
        valuesGood = np.array(valuesGood)
        pointsBad = np.array(pointsBad)

        valuesBad = griddata(pointsGood, valuesGood, (grid_x, grid_y), method='nearest')

        for iBad in range(0,pointsBad.shape[0]):
            i = pointsBad[iBad,0]
            j = pointsBad[iBad,1]
            arrayOut[iTime,i,j] = valuesBad[i,j]

    return arrayOut

#-------------------------------------------------------------------------------

def interpolate_array(nCells, remapMatrix, arrayIn):

    arrayOut = np.zeros((12,nCells))

    for iTime in range(0,12):

        arrayInTime = arrayIn[iTime,:,:].flatten()

        arrayOut[iTime,:] = remapMatrix.dot(arrayInTime)

    return arrayOut

#-------------------------------------------------------------------------------

def create_forcing(\
        filenameIn, \
        filenameOut, \
        nCells, \
        remapMatrix):

    fileIn = Dataset(filenameIn,"r")

    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells",nCells)
    fileOut.createDimension("StrLen",64)
    fileOut.createDimension("Time",None)

    # time
    xtimes = create_output_times(12, 0)
    varXtime = fileOut.createVariable("xtime","c",dimensions=["Time","StrLen"])
    for iTime in range(0,12):
        varXtime[iTime,0:19] = netCDF4.stringtochar(np.array(xtimes[iTime], 'S19'))
        varXtime[iTime,19:] = " "*45

    varSST  = fileOut.createVariable("seaSurfaceTemperature",   "d",dimensions=["Time","nCells"])
    varSSS  = fileOut.createVariable("seaSurfaceSalinity",      "d",dimensions=["Time","nCells"])
    varU    = fileOut.createVariable("uOceanVelocity",          "d",dimensions=["Time","nCells"])
    varV    = fileOut.createVariable("vOceanVelocity",          "d",dimensions=["Time","nCells"])
    varDhdx = fileOut.createVariable("seaSurfaceTiltU",         "d",dimensions=["Time","nCells"])
    varDhdy = fileOut.createVariable("seaSurfaceTiltV",         "d",dimensions=["Time","nCells"])
    varHblt = fileOut.createVariable("oceanMixedLayerDepth",    "d",dimensions=["Time","nCells"])
    varQdp  = fileOut.createVariable("oceanHeatFluxConvergence","d",dimensions=["Time","nCells"])

    print("Interpolate seaSurfaceTemperature")
    arrayIn = fileIn.variables["T"][:]
    arrayIn = fill_array(arrayIn)
    arrayOut = interpolate_array(nCells, remapMatrix, arrayIn)
    varSST[:] = arrayOut[:]

    print("Interpolate seaSurfaceSalinity")
    arrayIn = fileIn.variables["S"][:]
    arrayIn = fill_array(arrayIn)
    arrayOut = interpolate_array(nCells, remapMatrix, arrayIn)
    varSSS[:] = arrayOut[:]

    print("Interpolate uOceanVelocity")
    arrayIn = fileIn.variables["U"][:,0,:,:]
    arrayIn = fill_array(arrayIn)
    arrayOut = interpolate_array(nCells, remapMatrix, arrayIn)
    varU[:] = arrayOut[:]

    print("Interpolate vOceanVelocity")
    arrayIn = fileIn.variables["V"][:,0,:,:]
    arrayIn = fill_array(arrayIn)
    arrayOut = interpolate_array(nCells, remapMatrix, arrayIn)
    varV[:] = arrayOut[:]

    print("Interpolate seaSurfaceTiltU")
    arrayIn = fileIn.variables["dhdx"][:]
    arrayIn = fill_array(arrayIn)
    arrayOut = interpolate_array(nCells, remapMatrix, arrayIn)
    varDhdx[:] = arrayOut[:]

    print("Interpolate seaSurfaceTiltV")
    arrayIn = fileIn.variables["dhdy"][:]
    arrayIn = fill_array(arrayIn)
    arrayOut = interpolate_array(nCells, remapMatrix, arrayIn)
    varDhdy[:] = arrayOut[:]

    print("Interpolate oceanMixedLayerDepth")
    arrayIn = fileIn.variables["hblt"][:]
    arrayIn = fill_array(arrayIn)
    arrayOut = interpolate_array(nCells, remapMatrix, arrayIn)
    varHblt[:] = arrayOut[:]

    print("Interpolate oceanHeatFluxConvergence")
    arrayIn = fileIn.variables["qdp"][:]
    arrayIn = fill_array(arrayIn)
    arrayOut = interpolate_array(nCells, remapMatrix, arrayIn)
    varQdp[:] = arrayOut[:]

    fileIn.close()
    fileOut.close()

#-------------------------------------------------------------------------------

def perform_remapping(\
        filenameMPASGrid, \
        filenameGx1Grid, \
        filenameGx1OceanMixed, \
        filenameMPASOceanMixed, \
        scripDir):

    # create MPAS scrip grid file
    print("create_scrip_file_MPAS")
    scripGridFilename = "remap_grid_MPAS_tmp.nc"
    create_scrip_file_MPAS(filenameMPASGrid, scripGridFilename)

    # create gx1 scrip grid file
    print("create_scrip_file_gx1")
    scripGx1Filename = "remap_grid_gx1_tmp.nc"
    create_scrip_file_gx1(scripGx1Filename, filenameGx1Grid)

    # create input scrip file
    print("write_scrip_in_file")
    write_scrip_in_file("gx1")

    # run scrip to generate weights
    print("SCRIP")
    cmd = scripDir + "/scrip"
    os.system(cmd)

    # get remapping weights
    print("get_remapping_data")
    filenameRemapping = "remap_gx1_to_MPAS_tmp.nc"
    remapMatrix, dstGridSize = get_remapping_data(filenameRemapping)

    print("create_forcing ocean climatology")
    # combined ocean climatology
    create_forcing(\
            filenameGx1OceanMixed, \
            filenameMPASOceanMixed, \
            dstGridSize, \
            remapMatrix)

#-------------------------------------------------------------------------------

'''
create_ocean_forcing.py
=======================

Usage
-----

This script creates ocean forcing using CESM output.

Usage: python create_ocean_forcing.py configFilename

where configFilename is a python config file with the following example format:

[forcing_generation]
filenameMPASGrid = /location/of/MPAS/grid
filenameGx1Grid = /location/of/gx1/grid
filenameGx1OceanMixed = /location/of/gx1/ocean_mixed_file
filenameMPASOceanMixed = /location/of/output/ocean_mixed_file
scripDir = /location/of/SCRIP/executable

SCRIP
-----

This script requires the SCRIP package to be installed.
SCRIP is a software package which computes addresses and weights for remapping
and interpolating fields between grids in spherical coordinates. It can be
obtained from https://github.com/SCRIP-Project/SCRIP

gx1 input data
--------------

This script requires a gx1 grid file and ocean mixed file as input. These can be
obtained from:
https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/forcing/
MPAS-Seaice_clim_data.tar.gz
'''

if (len(sys.argv) != 2):
    print("Usage: python create_ocean_forcing.py configFilename")
    sys.exit()

config = ConfigParser.ConfigParser()
config.read(sys.argv[1])

filenameMPASGrid       = config.get('forcing_generation','filenameMPASGrid')
filenameGx1Grid        = config.get('forcing_generation','filenameGx1Grid')
filenameGx1OceanMixed  = config.get('forcing_generation','filenameGx1OceanMixed')
filenameMPASOceanMixed = config.get('forcing_generation','filenameMPASOceanMixed')
scripDir               = config.get('forcing_generation','scripDir')

perform_remapping(\
        filenameMPASGrid, \
        filenameGx1Grid, \
        filenameGx1OceanMixed, \
        filenameMPASOceanMixed, \
        scripDir)
