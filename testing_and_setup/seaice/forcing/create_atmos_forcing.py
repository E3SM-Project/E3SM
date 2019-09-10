from __future__ import print_function
from netCDF4 import Dataset
import netCDF4
import numpy as np
import os
import sys
import glob
import ConfigParser
from create_forcing import create_scrip_grid_file, get_mpas_grid_info, create_scrip_file_MPAS, write_scrip_in_file, create_output_times, get_remapping_data

#-------------------------------------------------------------------------------

def create_T62_remap_file(filenameScrip, title, dataDirSixHourly):

    nLat = 94
    nLon = 192

    filenames = sorted(glob.glob(dataDirSixHourly+"/t_10/t_10.*.nc"))
    filenameT62 = filenames[0]
    fileT62 = Dataset(filenameT62,"r")

    LATin = fileT62.variables["LAT"][:]
    LONin = fileT62.variables["LON"][:]

    fileT62.close()

    LAT = np.zeros(nLat+2)
    LON = np.zeros(nLon+2)

    LAT[1:-1] = LATin[:]
    LON[1:-1] = LONin[:]

    # center
    latCenter = np.zeros((nLat,nLon))
    lonCenter = np.zeros((nLat,nLon))

    for iLon in range(0,nLon):
       latCenter[:,iLon] = LATin[:]

    for iLat in range(0,nLat):
       lonCenter[iLat,:] = LONin[:]

    latCenter = np.radians(latCenter)
    lonCenter = np.radians(lonCenter)

    # corners
    latCorner = np.zeros((nLat,nLon,4))
    lonCorner = np.zeros((nLat,nLon,4))

    for iLon in range(0,nLon):

        iLon2 = iLon + 1

        lonCorner[:,iLon,0] = 0.5 * (LON[iLon2] + LON[iLon2-1])
        lonCorner[:,iLon,1] = 0.5 * (LON[iLon2] + LON[iLon2+1])
        lonCorner[:,iLon,2] = 0.5 * (LON[iLon2] + LON[iLon2+1])
        lonCorner[:,iLon,3] = 0.5 * (LON[iLon2] + LON[iLon2-1])

    lonCorner = (np.where(lonCorner < 0.0, lonCorner + 360.0, lonCorner))
    lonCorner = np.radians(lonCorner)

    for iLat in range(0,nLat):

        iLat2 = iLat + 1

        latCorner[iLat,:,0] = 0.5 * (LAT[iLat2] + LAT[iLat2-1])
        latCorner[iLat,:,1] = 0.5 * (LAT[iLat2] + LAT[iLat2-1])
        latCorner[iLat,:,2] = 0.5 * (LAT[iLat2] + LAT[iLat2+1])
        latCorner[iLat,:,3] = 0.5 * (LAT[iLat2] + LAT[iLat2+1])

    latCorner = np.radians(latCorner)

    # create file
    nGridSize = nLat * nLon
    nGridCorners = 4
    gridRank = 2
    gridDims = np.array([nLon,nLat])
    gridImask = np.ones(nGridSize,dtype="i")

    create_scrip_grid_file(filenameScrip, nGridSize, nGridCorners, gridRank, gridDims, latCenter.flatten(), lonCenter.flatten(), gridImask, latCorner.flatten(), lonCorner.flatten(), title)

#-------------------------------------------------------------------------------

def create_forcing(\
        filenameOutTemplate, \
        varnameInput, \
        yearStart, \
        yearStop, \
        filenameInputTemplate, \
        varnameOutput, \
        inputTimesPerYear, \
        yearStartData, \
        remapMatrix, \
        dstGridSize):

    # loop over years
    for year in range(yearStart,yearStop+1):

        print("  Year: %i of %i to %i" %(year, yearStart, yearStop))

        # create output file
        filenameOut = filenameOutTemplate.replace("$Y",str(year))
        fileForcing = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

        # dimensions
        nCells = fileForcing.createDimension("nCells",dstGridSize)
        StrLen = fileForcing.createDimension("StrLen",64)
        Time   = fileForcing.createDimension("Time",)

        # time
        xtimes = create_output_times(inputTimesPerYear, year)
        varXtime = fileForcing.createVariable("xtime","c",dimensions=["Time","StrLen"])
        for iTime in range(0,inputTimesPerYear):
            varXtime[iTime,0:19] = netCDF4.stringtochar(np.array(xtimes[iTime], 'S19'))
            varXtime[iTime,19:] = " "*45

        # loop over variables
        for iVariable in range(0,len(varnameInput)):

            print("    Variable: %s to %s" %(varnameInput[iVariable], varnameOutput[iVariable]))

            # open input file
            filenameInput = sorted(glob.glob(filenameInputTemplate[iVariable].replace("$Y",str(year))))[0]
            fileInput = Dataset(filenameInput,"r")
            arrayIn = fileInput.variables[varnameInput[iVariable]][:]
            fileInput.close()

            arrayOut = np.zeros((inputTimesPerYear,dstGridSize))

            # loop over times
            for iTime in range(0,inputTimesPerYear):

                arrayInTime = arrayIn[iTime,:,:].flatten()
                arrayOut[iTime,:] = remapMatrix.dot(arrayInTime)

            # output variable to netcdf file
            var = fileForcing.createVariable(varnameOutput[iVariable],"d",dimensions=["Time","nCells"])
            var[:] = arrayOut[:]

        # close forcing file
        fileForcing.close()

#-------------------------------------------------------------------------------

def write_scrip_in_file(srcTitle):

    scripFile = open("scrip_in","w")

    scripFile.write("&remap_inputs\n")
    scripFile.write("    num_maps = 1\n")
    scripFile.write("    grid1_file = 'remap_grid_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    grid2_file = 'remap_grid_MPAS_tmp.nc'\n")
    scripFile.write("    interp_file1 = 'remap_%s_to_MPAS_tmp.nc'\n" %(srcTitle))
    scripFile.write("    interp_file2 = 'remap_MPAS_to_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    map1_name = '%s to MPAS bilinear mapping'\n" %(srcTitle))
    scripFile.write("    map2_name = 'MPAS to %s bilinear mapping'\n" %(srcTitle))
    scripFile.write("    map_method = 'bilinear'\n")
    scripFile.write("    normalize_opt = 'frac'\n")
    scripFile.write("    output_opt = 'scrip'\n")
    scripFile.write("    restrict_type = 'latitude'\n")
    scripFile.write("    num_srch_bins = 90 \n")
    scripFile.write("    luse_grid1_area = .false.\n")
    scripFile.write("    luse_grid2_area = .false.\n")
    scripFile.write("/\n")

    scripFile.close()

#-------------------------------------------------------------------------------

def perform_remapping(\
        filenameMPASGrid, \
        outputDir, \
        startYear, \
        endYear, \
        dataDirSixHourly, \
        dataDirMonthly, \
        scripDir):

    # create MPAS scrip grid file
    print("create_scrip_file_MPAS")
    scripGridFilename  = "remap_grid_MPAS_tmp.nc"
    create_scrip_file_MPAS(filenameMPASGrid, scripGridFilename)

    # create T62 remapping file
    print("create_T62_remap_file")
    scripT62Filename = "remap_grid_T62_tmp.nc"
    create_T62_remap_file(scripT62Filename, "T62", dataDirSixHourly)

    # create input scrip file
    print("write_scrip_in_file")
    write_scrip_in_file("T62")

    # run scrip to generate weights
    print("SCRIP")
    cmd = scripDir + "/scrip"
    os.system(cmd)

    # get remapping weights
    print("get_remapping_data")
    filenameRemapping = "remap_T62_to_MPAS_tmp.nc"
    remapMatrix, dstGridSize = get_remapping_data(filenameRemapping)

    print("create_forcing six hourly")
    # combined six hourly file
    create_forcing(\
       outputDir+"/LYq_six_hourly.$Y.nc", \
       ["T_10_MOD","Q_10_MOD","U_10_MOD","V_10_MOD"], \
       startYear, \
       endYear, \
       [dataDirSixHourly+"/t_10/t_10.$Y.*.nc", \
        dataDirSixHourly+"/q_10/q_10.$Y.*.nc", \
        dataDirSixHourly+"/u_10/u_10.$Y.*.nc", \
        dataDirSixHourly+"/v_10/v_10.$Y.*.nc"], \
       ["airTemperature", \
        "airSpecificHumidity", \
        "uAirVelocity", \
        "vAirVelocity"], \
       1460, \
       1948, \
       remapMatrix, \
       dstGridSize)

    print("create_forcing monthly")
    # combined monthly file
    create_forcing(\
       outputDir+"/LYq_monthly.nc", \
       ["cldf","prec"], \
       0, \
       0, \
       [dataDirMonthly+"/cldf.omip.nc", \
        dataDirMonthly+"/prec.nmyr.nc"], \
       ["cloudFraction", \
        "rainfallRate"], \
       12, \
       0, \
       remapMatrix, \
       dstGridSize)

#-------------------------------------------------------------------------------

'''
create_atmos_forcing.py
=======================

Usage
-----

This script creates atmospheric forcing using six hourly CORE-II data and
monthly AOMIP climatologies.

Usage: python create_atmos_forcing.py configFilename

where configFilename is a python config file with the following example format:

[forcing_generation]
filenameMPASGrid = /location/of/MPAS/grid
outputDir = /location/to/put/output/forcing
startYear = 1948
endYear = 2007
dataDirSixHourly = /location/of/CORE-II/data
dataDirMonthly = /location/of/AOMIP/climatologies
scripDir = /location/of/SCRIP/executable

SCRIP
-----

This script requires the SCRIP package to be installed.
SCRIP is a software package which computes addresses and weights for remapping
and interpolating fields between grids in spherical coordinates. It can be
obtained from https://github.com/SCRIP-Project/SCRIP

CORE-II data
------------

Six-hourly air temperature, velocity and specific humidity comes from CORE-II.
Data files can be obtained from
https://data1.gfdl.noaa.gov/nomads/forms/core/COREv2/CIAF_v2.html
To generate forcing for a given year YYYY, the following files are required:
${dataDirSixHourly}/t_10/t_10.YYYY.*.nc
${dataDirSixHourly}/q_10/q_10.YYYY.*.nc
${dataDirSixHourly}/u_10/u_10.YYYY.*.nc
${dataDirSixHourly}/v_10/v_10.YYYY.*.nc
where ${dataDirSixHourly} is the local location of the six hourly data.

AOMIP climatologies
-------------------

Monthly climatologies of cloudiness and precipitation comes from AOMIP.
The following data files are required:
${dataDirMonthly}/cldf.omip.nc
${dataDirMonthly}/prec.nmyr.nc
where ${dataDirMonthly} is the local location of the monthly data.
These files can be obtained from:
https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/forcing/
MPAS-Seaice_clim_data.tar.gz
'''

if (len(sys.argv) != 2):
    print("Usage: python create_atmos_forcing.py configFilename")
    sys.exit()

config = ConfigParser.ConfigParser()
config.read(sys.argv[1])

filenameMPASGrid = config.get   ('forcing_generation','filenameMPASGrid')
outputDir        = config.get   ('forcing_generation','outputDir')
startYear        = config.getint('forcing_generation','startYear')
endYear          = config.getint('forcing_generation','endYear')
dataDirSixHourly = config.get   ('forcing_generation','dataDirSixHourly')
dataDirMonthly   = config.get   ('forcing_generation','dataDirMonthly')
scripDir         = config.get   ('forcing_generation','scripDir')

perform_remapping(\
        filenameMPASGrid, \
        outputDir, \
        startYear, \
        endYear, \
        dataDirSixHourly, \
        dataDirMonthly, \
        scripDir)
