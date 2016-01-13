#!/usr/bin/env python

"""
This script plots the results of an experiment with an ice stream.
"""

# Reconfigured by Joseph H Kennedy at ORNL on August 7, 2015 to work with the regression testing
#     NOTE: Did not adjust inner workings except where needed.

import os
import sys
import glob
import numpy

import matplotlib.pyplot as plt

from netCDF import *
from math import tan, pi, sin, cos
from runStream import *  # Get all the parameter values and functions defined in the run script

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


    y0 = filein.variables['y0'][:]
    ny = y0.shape[0] + 1  # the actual ny dimension, not the size of y0
    dy = y0[1]-y0[0]

    uvel = filein.variables['uvel'][0,:,:,:]  # get all levels at time 0
    vvel = filein.variables['vvel'][0,:,:,:]  # get all levels at time 0
    mintauf = filein.variables['tauf'][0,:,:]
    #btractx = filein.variables['btractx'][0,:,:]
    #btracty = filein.variables['btracty'][0,:,:]
    if netCDF_module == 'Scientific.IO.NetCDF':
        uvel = uvel * filein.variables['uvel'].scale_factor
        vvel = vvel * filein.variables['vvel'].scale_factor
        mintauf = mintauf * filein.variables['tauf'].scale_factor
    #    btractx = btractx * filein.variables['btractx'].scale_factor
    #    btracty = btracty * filein.variables['btracty'].scale_factor
    #btract = (btractx**2 + btracty**2)**0.5

    x0 = filein.variables['x0'][:]
    xpos = x0.shape[0]/2   # integer division on x-length to get the middle column of the domain

    ypos = y0.shape[0]/2   # integer division on y-length to get the middle row of the domain

    # Calculate analytic velocity profile - the analytic functions are in runStream.py
    if analytic_solution == 'raymond':
        uvel_analytic_profile = raymond_uvel(y0)
        analytic_name = 'raymond analytic solution'
    elif analytic_solution == 'schoof':
        uvel_analytic_profile = schoof_uvel(y0)
        analytic_name = 'schoof analytic solution'
    else:
        sys.exit("Error: Invalid value for 'analytic_solution'.")

    # ===================
    # Plot other diagnostics 
    # SP: moved this to Fig. 1 so that it lies on bottom)
    fig = plt.figure(1, facecolor='w', figsize=(16, 10), dpi=80)

    fig.add_subplot(2,2,1)
    colors = plt.cm.get_cmap('jet',len(y0))
    for j in range(len(y0)):
        plt.plot(x0/1000.0, uvel[ 0,j,:] - uvel[ 0,j,:].mean(), '.-', color=colors(j), label=str(j))
    plt.xlabel('distance along flow (km)')
    plt.ylabel('surface ALONG flow velocity, demeaned (m/a)')
    plt.title('Longit. x-sect. of model vel. for all rows \n(should be constant-valued)')
    # kludgy code to get a colorbar for the lines
    sm = plt.cm.ScalarMappable(cmap=colors, norm=plt.Normalize(vmin=0, vmax=len(y0)))
    sm._A = []; 
    cb=plt.colorbar(sm)
    cb.set_label('y-index')

    fig.add_subplot(2,2,2)
    for j in range(len(y0)):
        plt.plot(x0/1000.0, uvel[ -1,j,:] - uvel[ -1,j,:].mean(), '.-', color=colors(j), label=str(j))
    plt.xlabel('distance along flow (km)')
    plt.ylabel('basal ALONG flow velocity, demeaned (m/a)')
    #plt.title('Longit. x-sect. of model vel. for all rows (should be constant-valued)')
    # kludgy code to get a colorbar for the lines
    sm = plt.cm.ScalarMappable(cmap=colors, norm=plt.Normalize(vmin=0, vmax=len(y0)))
    sm._A = []; 
    cb=plt.colorbar(sm)
    cb.set_label('y-index')

    fig.add_subplot(2,2,3)
    for j in range(len(y0)):
        plt.plot(x0/1000.0, vvel[ 0,j,:] - vvel[ 0,j,:].mean(), '.-', color=colors(j), label=str(j))
    plt.xlabel('distance along flow (km)')
    plt.ylabel('surface ACROSS flow velocity, demeaned (m/a)')
    #plt.title('Longit. x-sect. of model vel. for all rows (should be constant-valued)')
    # kludgy code to get a colorbar for the lines
    sm = plt.cm.ScalarMappable(cmap=colors, norm=plt.Normalize(vmin=0, vmax=len(y0)))
    sm._A = []
    cb=plt.colorbar(sm)
    cb.set_label('y-index')

    fig.add_subplot(2,2,4)
    for j in range(len(y0)):
        plt.plot(x0/1000.0, vvel[-1,j,:] - vvel[-1,j,:].mean(), '.-', color=colors(j), label=str(j))
    plt.xlabel('distance along flow (km)')
    plt.ylabel('basal ACROSS flow velocity, demeaned (m/a)')
    #plt.title('Longit. x-sect. of model vel. for all rows (should be constant-valued)')
    # kludgy code to get a colorbar for the lines
    sm = plt.cm.ScalarMappable(cmap=colors, norm=plt.Normalize(vmin=0, vmax=len(y0)))
    sm._A = []
    cb=plt.colorbar(sm)
    cb.set_label('y-index')

    # ===================
    # Setup plot of uvel cross-section 
    # SP: moved this to Fig. 2 so that it lies on top)
    fig = plt.figure(2, facecolor='w', figsize=(12, 10), dpi=80)

    fig.add_subplot(2,1,1)
    plt.plot(y0/1000.0, uvel_analytic_profile, '-or', label=analytic_name)
    plt.plot(y0/1000.0, uvel[ 0,:,xpos], '-xk', label='CISM surface')
    plt.plot(y0/1000.0, uvel[-1,:,xpos], '-^k', label='CISM basal')

    plt.xlabel('distance across flow (km)')
    plt.ylabel('along flow velocity (m/a)')
    plt.title(analytic_name +  ' at x=%.1f km'%(x0[xpos]/1000.0))
    plt.legend()

    fig.add_subplot(2,1,2)
    plt.plot(y0/1000.0, mintauf[:, xpos], '-bo', label='Yield stress')
    plt.plot(y0/1000.0, numpy.ones(y0.size) * taud, '-k', label='Driving stress')
    #plt.plot(y0/1000.0, btract[:, xpos], '-go', label='basal traction')
    plt.xlabel('distance across flow (km)')
    plt.ylabel('stress (Pa)')
    plt.legend()

    plt.draw()
    plt.show()

    filein.close()


# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())






