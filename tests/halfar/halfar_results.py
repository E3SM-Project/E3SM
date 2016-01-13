#!/usr/bin/env python
# A script to compare CISM model output to the Halfar analytic solution of SIA evolution of a dome.

# Matt Hoffman, LANL, October 2013
# Reconfigured by Joseph H Kennedy at ORNL on August 7, 2015 to work with the regression testing
#     NOTE: Did not adjust inner workings except where needed.


import os
import sys
import glob
import subprocess

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

from math import tan, pi, sin, cos
from netCDF import *
from ConfigParser import ConfigParser

from runHalfar import halfarDome   # located in current directory

import argparse
parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# The command line options
# ------------------------
parser.add_argument('-o', '--output-dir', default='./output',
        help="Directory containing the tests output files. Warning: if there is a" \
        +"path passed via the -f/--out-file option, this argument will be" \
        +"ignored.")

parser.add_argument('-f', '--output-file', 
        help="The tests output file you would like to plot. If a path is" \
        +"passed via this option, the -o/--output-dir option will be ignored.")

parser.add_argument('-t', '--time-level', type=int, default=-1,
        help="Which time level to use.")


# ===========================================================
# Define some variables and functions used in the main script
# ===========================================================

def get_in_file():
    if args.output_file:
        out_d, out_f = os.path.split(args.output_file)
        if out_d:
            args.output_dir = out_d
            args.output_file = out_f
    
        print("\nUsing "+os.path.join(args.output_dir, args.output_file)+"\n")
        
    else:
        outpath = os.path.join(args.output_dir, '*.out.nc')
        matching = glob.glob(outpath)
        if len(matching) == 1:
            newest = matching[0]
            print("\nUsing "+newest+"\n")
       
        elif len(matching) > 1:
            newest = max(matching, key=os.path.getmtime)
            print("\nWARNING: MULTIPLE *.out.nc FILES DETECTED!")
            print(  "==========================================")
            print(  "Ploting the most recently modified file in the output directory:")
            print(  "    "+newest)
            print(  "To plot another file, specify it with the -f/--outfile option.\n")
            
        else:
            print("\nERROR: NO *.out.nc FILES DETECTED!")
            print(  "==================================")
            print(  "Either specify a location to look for the test output")
            print(  "files with the -o/--output-dir option, or the test output")
            print(  "file with the -f/--output-file option.\n")
            sys.exit(1)

        args.output_file = os.path.basename(newest)

    filein = NetCDFFile(os.path.join(args.output_dir, args.output_file),'r')
     
    return filein


# =========================
# Actual script starts here
# =========================
def main():
    """
    Plot the slab test results.
    """

    print("WARNING: THIS TEST CASE IS IN DEVELOPMENT. USE AT YOUR OWN RISK!")


    filein = get_in_file()    

    # Open config file for reading
    config_parser = ConfigParser()
    config_file = os.path.join(args.output_dir, args.output_file.replace('out.nc','config'))
    config_parser.read(config_file)

    # Get the value of flwa specified in the default_flwa parameter in the config file.
    # This is the only way this test case supports specifying flwa.
    try: 
        flwa = float(config_parser.get('parameters','default_flwa'))
        print 'Parameter used: ' + config_file + ' has specified a flwa value of ' + str(flwa)
        flow_law = int(config_parser.get('options','flow_law'))
        if flow_law != 0:
            sys.exit('Error: The option "flow_law" must be set to 0 for the test case to work properly.')
    except:
        raise
        sys.exit('Error: problem getting default_flwa parameter value from the config file')

    # Try to get ice density used by the model
    try:
        rhoi = float( subprocess.check_output( 'grep "real(dp),parameter :: rhoi =" ../../libglimmer/glimmer_physcon.F90 | cut -d " " -f 7 | cut -d "." -f 1', shell='/bin/bash' ) )
        print 'Parameter used: ../../libglimmer/glimmer_physcon.F90 has specified a rhoi value of ' + str(rhoi)
    except:
        print 'Warning: problem getting ice density value from ../../../libglimmer/glimmer_physcon.F90  Assuming 910.0 kg/m^3 as a default value.'
        rhoi = 910.0



    # open supplied file and get thickness slice needed
    x1 = filein.variables['x1'][:]
    dx = x1[1]-x1[0]
    y1 = filein.variables['y1'][:]
    ny = y1.size
    time = filein.variables['time'][:]
    thk = filein.variables['thk'][:]
    if netCDF_module == 'Scientific.IO.NetCDF':
        thk = thk * filein.variables['thk'].scale_factor

    # Call the halfar function
    thkHalfar = halfarDome(time[args.time_level]-time[0], x1, y1, flwa, rhoi)

    thkDiff = thk[args.time_level, :, :] - thkHalfar
    thkDiffIce = thkDiff[ np.where( thk[args.time_level,:,:] > 0.0) ]  # Restrict to cells modeled to have ice
    RMS = ( (thkDiffIce**2).sum() / float(len(thkDiffIce)) )**0.5

    # Print some stats about the error
    print '\nError statistics for cells modeled to have ice (in m):'
    print '* RMS error = ' + str( RMS )
    print '* Maximum error is ' + str( thkDiffIce.max() )
    print '* Minimum error is ' + str( thkDiffIce.min() )
    print '* Mean error is ' + str( thkDiffIce.mean() )
    print '* Median error is ' + str( np.median(thkDiffIce) )
    print '* Mean absolute error = ' + str( np.absolute(thkDiffIce).mean() )
    print '* Median absolute error = ' + str( np.median(np.absolute(thkDiffIce)) )
    print ''


    # ================
    # Plot the results
    # ================
    fig = plt.figure(1, facecolor='w', figsize=(10, 4), dpi=100)
    gray = np.ones(3)*0.8

    fig.add_subplot(1,3,1)
    plt.pcolor(x1/1000.0,y1/1000.0,  ma.masked_values(thk[args.time_level,:,:], 0.0) )
    plt.colorbar()
    plt.axis('equal')
    plt.title('Modeled thickness (m) \n at time ' + str(time[args.time_level]) )
    plt.xlabel('x (km)'); plt.ylabel('y (km)')

    fig.add_subplot(1,3,2)
    plt.pcolor(x1/1000.0,y1/1000.0, ma.masked_values(thkHalfar, 0.0) )
    plt.colorbar()
    plt.axis('equal')
    plt.title('Analytic thickness (m) \n at time ' + str(time[args.time_level]) )
    plt.xlabel('x (km)'); plt.ylabel('y (km)')

    fig.add_subplot(1,3,3)
    plt.pcolor(x1/1000.0,y1/1000.0, ma.masked_values(thkDiff, 0.0))
    plt.colorbar()
    plt.axis('equal')
    plt.title('Modeled thickness - Analytic thickness \n at time ' + str(time[args.time_level]) ) 
    plt.xlabel('x (km)'); plt.ylabel('y (km)')


    # optional second figure - cross section through center of dome
    # -------------------------------------------------------------

    #fig = plt.figure(2, facecolor='w', figsize=(10, 4), dpi=100)
    #yind = ny//2
    #print yind, y1[yind]

    #x1dense = np.linspace(x1[0], x1[-1], 1000)
    #thkHalfarDense = halfarDome(time[args.time_level]-time[0], x1dense, y1[[0, yind, -1]], flwa, rhoi)

    #plt.step(x1/1000.0 + 0.5*dx/1000.0, thk[args.time_level,yind,:], '-r', label='model')
    #plt.plot(x1/1000.0, thk[args.time_level,yind,:], '.r')

    #plt.plot(x1dense/1000.0, thkHalfarDense[1,:], '-k', label='analytic')

    #plt.stem(x1/1000.0, (thk[args.time_level,yind,:] - thkHalfar[yind,:]) * 10.0, '-.b', label='error x10.0')

    #plt.xlabel('x (km)'); plt.ylabel('Elevation (m)')
    #plt.legend()



    plt.draw()
    plt.show()

    filein.close()

# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())

