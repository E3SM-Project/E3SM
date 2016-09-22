#!/usr/bin/env python

"""
Read MPAS Landice files from MISMIP+ experiments, and convert to netCDF grounding-line
files following the protocol in Asay-Davis et al. (2016, GMD).

The following input fields are needed to create the grounding-line file
(assuming GL points are tracked at cell edges)::

     Dimensions: 
         nCells
         nEdges
         nVertLevels
         time

     Mesh quantities:
         xEdge, yEdge
         cellsOnEdge
         areaCell
         layerThicknessFractions
       
     Prognostic quantities:
         edgeMask
         cellMask
         thickness
         bedTopography
         uReconstructZonal
         uReconstructMeridional

Authors William Lipscomb and Matthew Hoffman

"""

import os, sys
import numpy as np
import netCDF4
import datetime
from optparse import OptionParser
import matplotlib.pyplot as plt

model = "_MPASLI"
GLbit = 256
Icebit = 32
Floatbit = 4
secInYr = 3600.0 * 24.0 * 365.0  # Note: this may be slightly wrong for some calendar types!
rhoi = 918.
rhow = 1028.

parser = OptionParser()

parser.add_option("-f", "--file", dest="filename", type='string', default='output_00000.nc', help="output file to analyze", metavar="FILE")
parser.add_option("-x", "--expt", dest="experiment", type='string', default = 'all', help="MISMIP+ experiment(s) to set up", metavar="EXPT")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"

options, args = parser.parse_args()

if options.experiment:
    if options.experiment == 'all':
        # Read output data for all experiments
        experiments = ['Ice0', 'Ice1r', 'Ice1ra', 'Ice1rr', 'Ice1rax', 'Ice1rrx', 'Ice2r', 'Ice2ra', 'Ice2rr', 'Ice2rax', 'Ice2rrx']
    else:
        experiments = [options.experiment]
else:
    sys.exit('Error: No experiment specified.  Please specify experiment(s) with the -x option')

################### DEFINE FUNCTIONS ######################

def xtime2numtime(xtime):
    """Define a function to convert xtime character array to numeric time values using datetime objects"""
    # First parse the xtime character array into a string 
    xtimestr = netCDF4.chartostring(xtime) # convert from the character array to an array of strings using the netCDF4 module's function

    dt = []
    for stritem in xtimestr:
        itemarray = stritem.strip().replace('_', '-').replace(':', '-').split('-')  # Get an array of strings that are Y,M,D,h,m,s
	results = [int(i) for i in itemarray]
	if (results[0] < 1900):  # datetime has a bug where years less than 1900 are invalid on some systems
            results[0] += 1900
	    dt.append( datetime.datetime(*results) ) # * notation passes in the array as arguments

	    # use the netCDF4 module's function for converting a datetime to a time number
            numtime = netCDF4.date2num(dt, units='seconds since '+str(dt[0]))   
	    numtime /= (3600.0 * 24.0 * 365.0)
	    numtime -= numtime[0]  # return years from start
    return numtime

def xtimeGetYear(xtime):
    """Get an array of years from an xtime array, ignoring any partial year information"""
    # First parse the xtime character array into a string 
    xtimestr = netCDF4.chartostring(xtime) # convert from the character array to an array of strings using the netCDF4 module's function
    years = np.zeros( (len(xtimestr),) )
    for i in range(len(xtimestr)):
        years[i] = ( int(xtimestr[i].split('-')[0]) ) # Get the year part and make it an integer
    return years

def glplot(ncfile, times, colora, label):
    """
    add a plot of grounding line points to current axes.
    uses the numpy.ma.MaskedArray when reading xGL,yGL
    Adapted from Steph Cornford's script in Asay-Davis et al. (2016)
    """
    ncid = netCDF4.Dataset(ncfile, 'r')
    ltime = ncid.variables["time"][:]

    print 'ltime:', ltime[:]

    lxmax = 0.0
    lxmin = 800.0
    for i in range(0, len(times)):
        seq = (ltime == times[i])
        xGL = ncid.variables["xGL"][:, seq]*1e-3
#        print 'i, xGL, max:', i, xGL[:], np.max(xGL)
        lxmax = max(np.max(xGL), lxmax)
        lxmin = min(np.min(xGL), lxmin)
        yGL = ncid.variables["yGL"][:, seq]*1e-3
        plt.plot(xGL, yGL, 's', ms=3, mfc=colora[i],
                 mec=colora[i], label=label + ', t = ' + format(times[i]))
    ncid.close()
    return lxmin, lxmax

outputFile = options.filename

for expt in experiments:

    print '\n Looking for output file', outputFile, 'in directory', expt

    try:
        os.chdir(expt)
    except:
        sys.exit('Could not find a directory for this experiment')

    try:
        ncfile = netCDF4.Dataset(outputFile, 'r')
    except:
        sys.exit('Could not find the output file in this directory')

    # Get dimensions
    nTime = len(ncfile.dimensions['Time'])
    nCells = len(ncfile.dimensions['nCells'])
    nEdges = len(ncfile.dimensions['nEdges'])
    nVertLevels = len(ncfile.dimensions['nVertLevels'])
    nVertInterfaces = nVertLevels + 1

    xtime = ncfile.variables['xtime'][:]
    years = xtimeGetYear(xtime)

    # Get mesh fields
    layerThicknessFractions = ncfile.variables['layerThicknessFractions'][:]
    areaCell = ncfile.variables['areaCell'][:]
    xEdge = ncfile.variables['xEdge'][:]
    yEdge = ncfile.variables['yEdge'][:]
    cellsOnEdge = ncfile.variables['cellsOnEdge'][:,:]

    # Get prognostic fields on cells
    cellMask = ncfile.variables['cellMask'][:,:]
    thickness = ncfile.variables['thickness'][:,:]
    bedTopography = ncfile.variables['bedTopography'][:,:]
    uReconstructZonal = ncfile.variables['uReconstructZonal'][:,:,:]
    uReconstructMeridional = ncfile.variables['uReconstructMeridional'][:,:,:]
    
    # Get prognostic fields on edges
    edgeMask = ncfile.variables['edgeMask'][:,:]

    # Compute vertical mean velocity
    uMeanZonal = np.zeros((nTime,nCells))
    uMeanMeridional = np.zeros((nTime,nCells))

    for k in range(nVertLevels):
        uMeanZonal[:,:] += layerThicknessFractions[k] * 0.5 * (uReconstructZonal[:,:,k] + uReconstructZonal[:,:,k+1])
        uMeanMeridional[:,:] += layerThicknessFractions[k] * 0.5 * (uReconstructMeridional[:,:,k] + uReconstructMeridional[:,:,k+1])

    # Create GL file 
    GLfile = expt + model + '.nc'
    ncfile = netCDF4.Dataset(GLfile, 'w')

    print 'Creating output GL file', GLfile

    # Set dimensions                                                                                  
    glptdim = ncfile.createDimension('nPointGL', size = None)
    timedim = ncfile.createDimension('nTime', size = nTime)

    # Create variables
    xGL = ncfile.createVariable('xGL', 'f4', ('nPointGL', 'nTime'))
    yGL = ncfile.createVariable('yGL', 'f4', ('nPointGL', 'nTime'))
    time = ncfile.createVariable('time', 'f4', ('nTime'))

    iceThicknessGL = ncfile.createVariable('iceThicknessGL', 'f4', ('nPointGL', 'nTime'))
    uSurfaceGL = ncfile.createVariable('uSurfaceGL', 'f4', ('nPointGL', 'nTime'))
    vSurfaceGL = ncfile.createVariable('vSurfaceGL', 'f4', ('nPointGL', 'nTime'))
    uBaseGL = ncfile.createVariable('uBaseGL', 'f4', ('nPointGL', 'nTime'))
    vBaseGL = ncfile.createVariable('vBaseGL', 'f4', ('nPointGL', 'nTime'))
    uMeanGL = ncfile.createVariable('uMeanGL', 'f4', ('nPointGL', 'nTime'))
    vMeanGL = ncfile.createVariable('vMeanGL', 'f4', ('nPointGL', 'nTime'))

    iceVolume = ncfile.createVariable('iceVolume', 'f4', ('nTime'))
    iceVAF = ncfile.createVariable('iceVAF', 'f4', ('nTime'))
    groundedArea = ncfile.createVariable('groundedArea', 'f4', ('nTime'))

    iceVolume = np.zeros((nTime,)) 
    iceVAF = np.zeros((nTime,)) 
    groundedArea = np.zeros((nTime,)) 

    print 'Created output variables'

    # Loop over time slices
    for iTime in range(nTime):

        time[iTime] = years[iTime]
        print 'iTime, time =', iTime, time[iTime] 

        # Loop over edges to gather GL info
        nGL = 0

        for iEdge in range(nEdges):

            # Identify grounding-line edges
            #WHL - The following commented line is from the MISMIP3d script, but seems too inclusive.
#            if np.nonzero( np.logical_and( ( (edgeMask[iTime,iEdge] & GLbit) / GLbit == 1), (xEdge[iEdge] > 0.0) ) ):
            #WHL - This logic seems to work.
            if np.logical_and( ( (edgeMask[iTime,iEdge] & GLbit) == GLbit), (xEdge[iEdge] > 0.0) ):

                nGL += 1
                m = nGL - 1   # so indexing starts at 0

                # find indices of adjacent cells
                # subtract 1 from cellsOnEdge because indexing of cell-centered fields starts at 0
                iCell1 = cellsOnEdge[iEdge,0] - 1  
                iCell2 = cellsOnEdge[iEdge,1] - 1

                xGL[m,iTime] = xEdge[iEdge]
                yGL[m,iTime] = yEdge[iEdge]

                # average quantities from neighboring cell centers to the edge
                # convert velocity units to m/yr

                iceThicknessGL[m,iTime] = 0.5 * (thickness[iTime,iCell1] + thickness[iTime,iCell2])

                iLev = 0    # subtract 1 from level 1 for 0-based indexing
                uSurfaceGL[m,iTime] = 0.5 * (uReconstructZonal[iTime,iCell1,iLev] + uReconstructZonal[iTime,iCell2,iLev]) * secInYr
                vSurfaceGL[m,iTime] = 0.5 * (uReconstructMeridional[iTime,iCell1,1] + uReconstructMeridional[iTime,iCell2,1]) * secInYr

                iLev = nVertInterfaces - 1   # subtract 1 for 0-based indexing
                uBaseGL[m,iTime] = 0.5 * (uReconstructZonal[iTime,iCell1,iLev] + uReconstructZonal[iTime,iCell2,iLev]) * secInYr
                vBaseGL[m,iTime] = 0.5 * (uReconstructMeridional[iTime,iCell1,iLev] + uReconstructMeridional[iTime,iCell2,iLev]) * secInYr

                uMeanGL[m,iTime] = 0.5 * (uMeanZonal[iTime,iCell1] + uMeanZonal[iTime,iCell2]) * secInYr
                vMeanGL[m,iTime] = 0.5 * (uMeanMeridional[iTime,iCell1] + uMeanMeridional[iTime,iCell2]) * secInYr

        # compute ice volume, ice volume-above-flotation, and grounded area
        # Note: thickness above flotation = H - Hf, where Hf = max( (-rhow/rhoi)*topg, 0.)

        for iCell in range(nCells):

            iceVolume[iTime] += areaCell[iCell] * thickness[iTime,iCell]

            if (cellMask[iTime,iCell] & Floatbit) != Floatbit:    # not floating
                if bedTopography[iTime,iCell] < 0.:
                    iceVAF[iTime] += areaCell[iCell] * (thickness[iTime,iCell] + (rhow/rhoi)*bedTopography[iTime,iCell])
                else:
                    iceVAF[iTime] += areaCell[iCell] * thickness[iTime,iCell]

            if np.logical_and( (cellMask[iTime,iCell] & Icebit), 
                               ( (cellMask[iTime,iCell] & Floatbit) != Floatbit) ):   # ice is present and not floating
                groundedArea[iTime] += areaCell[iCell]

        print 'ice volume (m^3) =', iceVolume[iTime]
        print 'ice VAF (m^3) =', iceVAF[iTime]
        print 'grounded area (m^2) =', groundedArea[iTime]

        ncfile.variables['iceVolume'][:] = iceVolume[:]
        ncfile.variables['iceVAF'][:] = iceVAF[:]
        ncfile.variables['groundedArea'][:] = groundedArea[:]

    # Close file for this experiment
    ncfile.close()

    # Make sure the grounded area is reasonable
    ncfile = netCDF4.Dataset(GLfile, 'r')
    groundedArea = ncfile.variables['groundedArea'][:]
    print 'groundedArea =', groundedArea[:]
    ncfile.close()

    # Create a test plot from the data in this file
    print 'Making plot from file', GLfile

    if expt == 'Ice0' or expt == 'Ice1r' or expt == 'Ice2r':
        timeList = [0, 20, 40, 60, 80, 100]
    elif expt == 'Ice1ra' or expt == 'Ice1rr' or expt == 'Ice2ra' or expt == 'Ice2rr':
        timeList = [100, 120, 140, 160, 180, 200]
    elif expt == 'Ice1rax' or expt == 'Ice1rrx' or expt == 'Ice2rax' or expt == 'Ice2rrx':
        timeList = [200, 300, 400, 600, 800, 1000]  # hardwired to 6 points for now

    xmin, xmax = glplot(GLfile, timeList,
                        ['black', 'red', 'orange', 'green', 'blue', 'purple'], 'Test')

    plt.legend(loc='center right', ncol=1, frameon=True, borderaxespad=0)

    plotname = 'plot_' + expt + model + '.pdf'
    plt.savefig(plotname)
    plt.clf()

    print 'Created test plot', plotname

    # Change to parent directory to process the next experiment
    os.chdir('..')
