#! /usr/bin/env python
"""
    Take an input netcdf database and extract particle information, converting to
    vtu format for input into ParaView.

    Usage exampe:
    ./animated_particles.py -f free_floating.nc

    produces

    particles_free_floating*.vtu
    particles_free_floating.pvd

    Phillip Wolfram
    LANL
    07/15/2014

"""
from pyevtk.vtk import VtkGroup
from pyevtk.hl import pointsToVTK
import numpy as np

import netCDF4
import os
import pdb

def cleanup(data):
    if type(data) == np.ma.core.MaskedArray:
        data.data[data.mask] = data.fill_value
        data = data.data
    return data

def build_particle_file(fname_in, fname_outbase):
    """
    Take existing netCDF4 input f_in and produce
    vtu files for f_outbase with all the particle data

    Phillip Wolfram
    LANL
    07/15/2014
    """
    ## initialize the starting points
    fullpath = fname_in.split(os.sep)

    # don't want to place this in starting directory
    #g = VtkGroup(os.sep.join(fullpath[:-1]) + '/' + fname_outbase + '_' + os.path.splitext(fullpath[-1])[0] )
    g = VtkGroup(fname_outbase + '_' + os.path.splitext(fullpath[-1])[0] )

    # load the netCDF database

    ## get time
    #time = 0.0


    # open the database
    f_in = netCDF4.Dataset(fname_in, 'r')

    # num particles and data frames
    Npoints = len(f_in.dimensions['nParticles'])
    Ntime = len(f_in.dimensions['Time'])

    # initialize empty dictionary
    particleData = dict()

    for t in np.arange(Ntime):

        # get the points locations
        #x = f_in.variables['xParticle'][t] #+ f_in.variables['zLevelParticle'][t]
        #y = f_in.variables['yParticle'][t] #+ f_in.variables['zLevelParticle'][t]
        #z = f_in.variables['zParticle'][t] #+ f_in.variables['zLevelParticle'][t]
        x = f_in.variables['xParticle'][t]
        y = f_in.variables['yParticle'][t]
        z = f_in.variables['zParticle'][t]
        print t
        print np.vstack((x,y,z)).T

        particleType = f_in.variables['xParticle'].dimensions
        particleTypeStatic = f_in.variables['indexToParticleID'].dimensions

        # get particle data
        particleData.clear()
        coordinateVars = ['xParticle', 'yParticle', 'zParticle']
        nonCoords = [i for i in f_in.variables if i not in coordinateVars]
        for v in nonCoords:
            dim = f_in.variables[v].dimensions
            if dim == particleType:
                particleData[str(v)] = cleanup(f_in.variables[v][t])
            if dim == particleTypeStatic:
                particleData[str(v)] = cleanup(f_in.variables[v][:])

        #print particleData

        # set the new time (just use index for simplicity for now)
        time = 2*t
        #time = ''.join(f_in.variables['xtime'][t,0,:]).replace(' ','')

        # file name
        #f_out = fname_outbase + '_' + os.path.splitext(fname_in)[0] + '_' + str(t)
        f_out = fname_outbase + '_' + os.path.splitext(fullpath[-1])[0] + '_' + str(t)
        # output file
        pointsToVTK(f_out, x, y, z, data = particleData)
        # make sure the file is added to the time vector
        g.addFile(filepath = f_out + '.vtu', sim_time = time)


    # save the animation file
    g.save()


if __name__ == "__main__":
    from optparse import OptionParser

    # Get command line parameters
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="inputfilename",
            help="file to open for appending \
                    particle data 'particle_' extension",
                    metavar="FILE")

    parser.add_option("-o", "--output", dest="outputfilename",
            help="output file name base for export to *.vtu \
                    file which stores particle data",
                    metavar="FILE")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename is a required input.")

    if not options.outputfilename:
        parser.error("Output filename base for vtu files is required input.")

    build_particle_file(options.inputfilename, options.outputfilename)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
