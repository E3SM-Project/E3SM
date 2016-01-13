#!/usr/bin/env python

"""
This script runs an experiment for an ice sheet with a "dome" shape on a flat
base.  The time evolution of this dome shape using the shallow-ice
approximation has an analytic solution.  For details, see: 
    Halfar, P. 1983. On the Dynamics of the Ice Sheets 2.  Journal of Geophysical
    Research, 88, 6043-6051.
"""

# Authors
# -------
# Modified from dome.py script written by Glen Granzow at the University of Montana on April 13, 2010
# Modified for Halfar test case by Matt Hoffman, October 2013.
# Reconfigured by Joseph H Kennedy at ORNL on August 12, 2015 to work with the regression testing

import os
import sys
import errno
import subprocess
import ConfigParser 

import numpy
import netCDF
from math import sqrt


# Parse the command line options
# ------------------------------
import argparse
parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# small helper function so argparse will understand unsigned integers
def unsigned_int(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("This argument is an unsigned int type! Should be an integer greater than zero.")
    return x

parser.add_argument('-c','--config', default='./halfar.config', 
        help="The configure file used to setup the test case and run CISM")
parser.add_argument('-e','--executable', default='./cism_driver', 
        help="The CISM driver")
parser.add_argument('--hpc', action='store_true',
        help="Shortcuts parallel run command lookup for High Performance Computing Systems. Will set run command to `time apirun -n N`.")
parser.add_argument('-m', '--modifier', metavar='MOD', default='',
        help="Add a modifier to file names. FILE.EX will become FILE.MOD.EX")
parser.add_argument('-n','--parallel', metavar='N', type=unsigned_int, default=0, 
        help="Run in parallel using N processors.")
parser.add_argument('-o', '--output-dir', default='./output',
        help="Write all created files here.")
parser.add_argument('-q', '--quiet', action='store_true',
        help="Run the CISM process quietly.")
parser.add_argument('-s','--setup-only', action='store_true',
        help="Set up the test, but don't actually run it.")

# Additional test specific options:


# Some useful functions
# ---------------------

# function to make a directory, and not worry if it exists.
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


# prep the command line functions
def prep_commands(args, config_name):
    driver = os.path.abspath(args.executable)
   
    quiet_mod = ''
    if args.quiet:
        quiet_mod = ' > '+config_name+'.oe'

    commands = []
    mkdir_p(args.output_dir)
    commands.append("cd "+os.path.abspath(args.output_dir))
    
    if args.hpc and (args.parallel > 0):
        mpiexec = 'aprun -n ' + str(args.parallel)+" "
    elif (args.parallel > 0):
        # These calls to os.system will return the exit status: 0 for success (the command exists), some other integer for failure
        if os.system('which openmpirun > /dev/null') == 0:
            mpiexec = 'openmpirun -np ' + str(args.parallel)+" "
        elif os.system('which mpirun > /dev/null') == 0:
            mpiexec = 'mpirun -np ' + str(args.parallel)+" "
        elif os.system('which aprun > /dev/null') == 0:
            mpiexec = 'aprun -n ' + str(args.parallel)+" "
        elif os.system('which mpirun.lsf > /dev/null') == 0:
            # mpirun.lsf does NOT need the number of processors (options.parallel)
            mpiexec = 'mpirun.lsf '
        else:
            print("Unable to execute parallel run!")
            print("   Please edit the script to use your MPI run command, or run the model manually with")
            print("   something like: mpirun -np 4 ./cism_driver halfar.config")
            sys.exit(1)
    else:
        mpiexec = ''

    commands.append(mpiexec+driver+" "+config_name+quiet_mod)

    return commands


# Define the function to calculate the Halfar thickness
# Halfar, P. 1983. On the Dynamics of the Ice Sheets 2.  Journal of Geophysical Research, 88, 6043-6051.
def halfarDome(t,x,y,flwa,rhoi):
  # Input: t - time in years
  # Input: x - 1-d array of cell center x-positions
  # Input: y - 1-d array of cell center y-positions
  # Input: flwa - flow law parameter A in units of Pa^-3 yr^-1
  # Input: rhoi - ice density in kg/m3

  # Initial radius and central thickness of dome
  R0 = 60000.0 * numpy.sqrt(0.125)
  H0 = 2000.0 * numpy.sqrt(0.125)

  n = 3.0
  grav = 9.8101
  alpha = 1.0/9.0
  beta = 1.0/18.0
  secpera = 31556926.0
  Gamma = 2.0/(n+2.0) * flwa * (rhoi * grav)**n

  xcenter = max(x)/2.0
  ycenter = max(y)/2.0

  t0 = (beta/Gamma) * (7.0/4.0)**3 * (R0**4/H0**7)  # Note: this line assumes n=3!
  tr=(t+t0)/t0 

  H=numpy.zeros((len(y), len(x)))
  for i in range(len(x)):
    for j in range(len(y)):
      r = numpy.sqrt( (x[i]-xcenter)**2 + (y[j]-ycenter)**2)
      r=r/R0
      inside = max(0.0, 1.0 - (r / tr**beta)**((n+1.0) / n))

      H[j,i] = H0 * inside**(n / (2.0*n+1.0)) / tr**alpha
  return H.astype(numpy.float32)


# the main script function
# ------------------------
def main():
    """
    Run the test.
    """
    
    # check that file name modifier, if it exists, starts with a '-'
    if not (args.modifier == '') and not args.modifier.startswith('-') :
        args.modifier = '-'+args.modifier
         
    # get the configuration
    # ---------------------
    try:
        config_parser = ConfigParser.SafeConfigParser()
        config_parser.read( args.config )
        
        nz = int(config_parser.get('grid','upn'))
        
        nx = int(config_parser.get('grid','ewn'))
        ny = int(config_parser.get('grid','nsn'))
        dx = float(config_parser.get('grid','dew'))
        dy = float(config_parser.get('grid','dns'))
        
        flwa = float(config_parser.get('parameters','default_flwa'))
        flow_law = int(config_parser.get('options','flow_law'))
        if flow_law != 0:
            print('Error: The option "flow_law" must be set to 0 for the test case to work properly.')
            sys.exit(1)
       
        file_name = config_parser.get('CF input', 'name')
        root, ext = os.path.splitext(file_name)

    except ConfigParser.Error as error:
        print("Error parsing " + args.config )
        print("   "), 
        print(error)
        sys.exit(1)
    
    res = str(nx).zfill(4)
    if args.parallel > 0:
        mod = args.modifier+'.'+res+'.p'+str(args.parallel).zfill(3)
    else:
        mod = args.modifier+'.'+res
   
    file_name = root+mod+ext
    config_name = root+mod+'.config'
    out_name = root+mod+'.out'+ext

    
    # create the new config file
    # --------------------------
    if not args.quiet: 
        print("\nCreating config file: "+config_name)
    
    config_parser.set('grid', 'udn', str(nz))
    
    config_parser.set('grid', 'ewn', str(nx))
    config_parser.set('grid', 'nsn', str(ny))
    config_parser.set('grid', 'dew', str(dx))
    config_parser.set('grid', 'dns', str(dy))

    config_parser.set('parameters', 'default_flwa', str(flwa))
    config_parser.set('options', 'flow_law', str(flow_law))
    
    config_parser.set('CF input', 'name', file_name)
    config_parser.set('CF output', 'name', out_name)
    config_parser.set('CF output', 'xtype', 'double')
   
    with open(config_name, 'wb') as config_file:
        config_parser.write(config_file)


    # create the input netCDF file
    # ----------------------------
    if not args.quiet: 
        print("\nCreating halfar netCDF file: "+file_name)
    try:
        nc_file = netCDF.NetCDFFile(file_name,'w',format='NETCDF3_CLASSIC')
    except TypeError:
        nc_file = netCDF.NetCDFFile(file_name,'w')

    nc_file.createDimension('time',1)
    nc_file.createDimension('x1',nx)
    nc_file.createDimension('y1',ny)
    nc_file.createDimension('level',nz)
    nc_file.createDimension('x0',nx-1) # staggered grid 
    nc_file.createDimension('y0',ny-1)

    x = dx*numpy.arange(nx,dtype='float32')
    y = dx*numpy.arange(ny,dtype='float32')

    nc_file.createVariable('time','f',('time',))[:] = [0]
    nc_file.createVariable('x1','f',('x1',))[:] = x
    nc_file.createVariable('y1','f',('y1',))[:] = y
    nc_file.createVariable('x0','f',('x0',))[:] = dx/2 + x[:-1] # staggered grid
    nc_file.createVariable('y0','f',('y0',))[:] = dy/2 + y[:-1]

    # Calculate values for the required variables.
    thk  = numpy.zeros([1,ny,nx],dtype='float32')
    topg = numpy.zeros([1,ny,nx],dtype='float32')


    # Try to get ice density used by the model
    # NOTE: might be better to just hard-code the same value here, and put a
    # note in glimmer_physcon.F90 to also change it here when changed there. --JHK 
    try:
       rhoi = float(subprocess.check_output( 
           'grep "real(dp),parameter :: rhoi =" ../../libglimmer/glimmer_physcon.F90 | cut -d " " -f 7 | cut -d "." -f 1', 
           shell='/bin/bash' 
           ))

       print('Parameter used: ../../libglimmer/glimmer_physcon.F90 has specified a rhoi value of '+str(rhoi))

    except:
       print('Warning: problem getting ice density value from ../../../libglimmer/glimmer_physcon.F90  Assuming 910.0 kg/m^3 as a default value.')
       rhoi = 910.0


    # Calculate the thickness of the halfar dome of ice
    thk = halfarDome(0.0, x, y, flwa, rhoi)  # Get the initial time shape from the halfar function
    # Note: The halfar solution will assume flwa = 1.0e-16, 
    #   so don't modify the default temperature settings.

    # Create the required variables in the netCDF file.
    nc_file.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
    nc_file.createVariable('topg','f',('time','y1','x1'))[:] = topg

    nc_file.close()
    
    mkdir_p(args.output_dir)
    subprocess.check_call("cp *rilinosOptions.xml "+args.output_dir, shell=True)
    subprocess.check_call("mv "+file_name+" "+args.output_dir, shell=True)
    subprocess.check_call("mv "+config_name+" "+args.output_dir, shell=True)

    # Run CISM
    # --------
    command_list =  prep_commands(args, config_name) 
    commands_all = ["# HALFAR"+mod+" test"]
    commands_all.extend(command_list)
   
    result_mv = "mv results "+root+mod+".results 2>/dev/null"
    timing_mv = "for file in cism_timing*; do mv $file "+root+mod+".$file 2>/dev/null; done"
    commands_all.append(result_mv)
    commands_all.append(timing_mv)
    commands_all.append(" ")

    if not args.setup_only:
        if not args.quiet: 
            print("\nRunning CISM halfar test")
            print(  "========================\n")

        process = subprocess.check_call(str.join("; ",command_list), shell=True)
   
        try:
            subprocess.check_call("cd "+args.output_dir+"; "+result_mv, shell=True)
        except subprocess.CalledProcessError:
            pass 

        try:
            subprocess.check_call("cd "+args.output_dir+"; "+timing_mv, shell=True)
        except subprocess.CalledProcessError:
            pass

        if not args.quiet: 
            print("\nFinished running the CISM halfar test")
            print(  "=====================================\n")
    else:
        run_script = args.output_dir+os.sep+root+mod+".run" 
        
        run_file = open(run_script,'w') 
        
        run_file.write('#!/bin/bash \n')
        for command in commands_all:
            run_file.write(command+" \n")

        run_file.close()
        os.chmod(run_script, 0o755)   # uses an octal number!

        if not args.quiet:
            print("\nFinished setting up the CISM halfar test")
            print(  "========================================")
            print(  "   To run the test, use: "+run_script)


# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())

