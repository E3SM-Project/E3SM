#! /usr/bin/env python
"""
    Take an input netcdf database and extract particle information,
    converting to vtu format for input into ParaView.

    Usage exampe:
    ./animated_particles.py -f free_floating.nc

    produces

    particles_free_floating*.vtu
    particles_free_floating.pvd

    Phillip Wolfram
    LANL
    07/15/2014

    Riley Brady
    CU Boulder/LANL
    05/15/2019
    * Compatible with python3, code cleanup.
    * VTK time now decimal format and compatible with `annotate_date`.
    * --ignore_xtime option if `xtime` not in particle output.
"""
from pyevtk.vtk import VtkGroup
from pyevtk.hl import pointsToVTK
import numpy as np
import netCDF4
import os
from datetime import datetime
from netCDF4 import date2num
import argparse


def cleanup(data):
    if type(data) == np.ma.core.MaskedArray:
        data.data[data.mask] = data.fill_value
        data = data.data
    return data


def build_particle_file(fname_in, fname_outbase, ignore_xtime):
    """
    Take existing netCDF4 input f_in and produce
    vtu files for f_outbase with all the particle data

    Phillip Wolfram
    LANL
    07/15/2014
    """
    # initialize the starting points
    fullpath = fname_in.split(os.sep)

    # don't want to place this in starting directory
    g = VtkGroup(fname_outbase + '_' + os.path.splitext(fullpath[-1])[0])

    # open the database
    f_in = netCDF4.Dataset(fname_in, 'r')

    # num particles and data frames
    Ntime = len(f_in.dimensions['Time'])

    # initialize empty dictionary
    particleData = dict()

    for t in np.arange(Ntime):
        # get the points locations
        x = f_in.variables['xParticle'][t]
        y = f_in.variables['yParticle'][t]
        z = f_in.variables['zParticle'][t]
        print(f'{t+1}/{Ntime}...')

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

        if (ignore_xtime) or ('xtime' not in f_in.variables):
            # Set time to a simple index.
            time = 2*t
        else:
            # Set decimal time similar to that in `paraview_vtk_extractor` so
            # that the `annotate_date` macro from MPAS-Tools can be used.
            xtime = f_in.variables['xtime'][t].tostring().decode('utf-8') \
                        .strip()
            date = datetime(int(xtime[0:4]), int(xtime[5:7]),
                            int(xtime[8:10]), int(xtime[11:13]),
                            int(xtime[14:16]), int(xtime[17:19]))
            time = date2num(date, units='days since 0000-01-01',
                            calendar='noleap')/365

        # file name
        f_out = f'{fname_outbase}_{os.path.splitext(fullpath[-1])[0]}_{str(t)}'

        # output file
        pointsToVTK(f_out, x, y, z, data=particleData)
        # make sure the file is added to the time vector
        g.addFile(filepath=f_out + '.vtu', sim_time=time)

    # save the animation file
    g.save()


if __name__ == "__main__":
    # Get command line parameters
    ap = argparse.ArgumentParser(
            description="Converts raw LIGHT output to VTK files for analysis \
            in ParaView.")
    ap.add_argument("-f", "--file", dest="inputfilename", required=True,
                    type=str, help="file to open for appending particle data \
                            'particle_' extension")
    ap.add_argument("-o", "--output", dest="outputfilename", required=True,
                    type=str, help="output file name base for export to *.vtu \
                            file which stores particle data")
    ap.add_argument("--ignore_xtime", dest="ignore_xtime", required=False,
                    action="store_true",
                    help="Ignore the `xtime` variable and just pass integer "
                         "time to ParaView. (Default is to convert `xtime` "
                         "to decimal time for `annotate_date` to work "
                         "from MPAS-Tools.")
    args = vars(ap.parse_args())
    build_particle_file(args['inputfilename'], args['outputfilename'],
                        args['ignore_xtime'])

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
