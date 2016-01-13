#!/usr/bin/env python

#FIXME: More detailed description of this test case!!!
"""
Run an experiment with an ice "slab". 
"""

# Authors
# -------
# Modified from dome.py by Matt Hoffman, Dec. 16, 2013
#    Test case described in sections 5.1-2 of:
#    J.K. Dukoqicz, 2012. Reformulating the full-Stokes ice sheet model for a 
#    more efficient computational solution. The Cryosphere, 6, 21-34. www.the-cryosphere.net/6/21/2012/
# Reconfigured by Joseph H Kennedy at ORNL on April 27, 2015 to work with the regression testing

import os
import sys
import errno
import subprocess
import ConfigParser 

import numpy
import netCDF
from math import sqrt, tan, pi, cos


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

parser.add_argument('-c','--config', default='./slab.config', 
        help="The configure file used to setup the test case and run CISM")
parser.add_argument('-e','--executable', default='./cism_driver', 
        help="The CISM driver")
parser.add_argument('--hpc', action='store_true',
        help="Shortcuts parallel run command lookup for High Performance Computing Systems. Will set run command to `time apirun -n N`.")
parser.add_argument('-m', '--modifier', metavar='MOD', default='',
        help="Add a modifier to file names. FILE.EX will become FILE.MOD.EX")
parser.add_argument('-n','--parallel', metavar='N', type=unsigned_int, default=0, 
        help="Run in parallel using N processors.")
parser.add_argument('-o', '--output-dir',default='./output',
        help="Write all created files here.")
parser.add_argument('-q', '--quiet', action='store_true',
        help="Run the CISM process quietly.")
parser.add_argument('-s','--setup-only', action='store_true',
        help="Set up the test, but don't actually run it.")


# Additional test specific options:
#parser.add_argument('--scale', type=unsigned_int, default=0, 
#        help="Scales the problem size by 2**SCALE. SCALE=0 creates a 31x31 grid, SCALE=1 " 
#            +"creates a 62x62 grid, and SCALE=2 creates a 124x124 grid.")


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
            # mpirun.lsf does NOT need the number of processors
            mpiexec = 'mpirun.lsf '
        else:
            print("Unable to execute parallel run!")
            print("   Please edit the script to use your MPI run command, or run the model manually with")
            print("   something like: mpirun -np 4 ./cism_driver slab.config")
            sys.exit(1)
    else:
        mpiexec = ''

    commands.append(mpiexec+driver+" "+config_name+quiet_mod)

    return commands


# Hard coded test specific parameters
# -----------------------------------
#FIXME: Some of these could just be options!

# Physical parameters 
n = 1 # flow law parameter - only the n=1 case is currently supported
# (implementing the n=3 case would probably require implementing a new efvs option in CISM)
rhoi = 910.0 # kg/m3
grav = 9.1801 # m^2/s

# Test case parameters
theta = 18  # basal inclination angle (degrees)  unpub. man. uses example with theta=18
thickness = 1000.0  # m  thickness in the rotated coordinate system, not in CISM coordinates
baseElevation = 1000.0 # arbitrary height to keep us well away from sea level

efvs = 2336041.42829      # hardcoded in CISM for constant viscosity setting (2336041.42829 Pa yr)

eta = 10.0   # unpub. man. uses example with eta=10.0
beta = eta / thickness / efvs**-n / (rhoi * grav * thickness)**(n-1)  # Pa yr m^-1
# Note: Fig. 3 in Ducowicz (2013) uses eta=18, where eta=beta*H/efvs
  
  
# the main script function
# ------------------------
def main():
    """
    Run the slab test.
    """

    print("WARNING: THIS TEST CASE IS IN DEVELOPMENT. USE AT YOUR OWN RISK!")

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

    
    # Setup the domain
    # ----------------
    offset = 1.0 * float(nx)*dx * tan(theta * pi/180.0)


    # create the new config file
    # --------------------------
    if not args.quiet: 
        print("\nCreating config file: "+config_name)
    
    config_parser.set('grid', 'upn', str(nz))
    config_parser.set('grid', 'ewn', str(nx))
    config_parser.set('grid', 'nsn', str(ny))
    config_parser.set('grid', 'dew', str(dx))
    config_parser.set('grid', 'dns', str(dy))

    config_parser.set('parameters', 'periodic_offset_ew', str(offset))
    
    config_parser.set('CF input', 'name', file_name)
    config_parser.set('CF output', 'name', out_name)
    config_parser.set('CF output', 'xtype', 'double')
    
    with open(config_name, 'wb') as config_file:
        config_parser.write(config_file)


    # create the input netCDF file
    # ----------------------------
    if not args.quiet: 
        print("\nCreating slab netCDF file: "+file_name)
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
    artm = numpy.zeros([1,ny,nx],dtype='float32')
    unstagbeta = numpy.zeros([1,ny,nx],dtype='float32')

    # Calculate the geometry of the slab of ice
    thk[:] = thickness / cos(theta * pi/180.0)
    xmax = x[:].max()
    for i in range(nx):
        topg[0,:,i] = (xmax - x[i]) * tan(theta * pi/180.0) + baseElevation
    unstagbeta[:] = beta

    # Create the required variables in the netCDF file.
    nc_file.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
    nc_file.createVariable('topg','f',('time','y1','x1'))[:] = topg
    nc_file.createVariable('unstagbeta','f',('time','y1','x1'))[:] = unstagbeta

    nc_file.close()
    mkdir_p(args.output_dir)
    subprocess.check_call("cp *rilinosOptions.xml "+args.output_dir, shell=True)
    subprocess.check_call("mv "+file_name+" "+args.output_dir, shell=True)
    subprocess.check_call("mv "+config_name+" "+args.output_dir, shell=True)


    # Run CISM
    # --------
    command_list = prep_commands(args, config_name)
    commands_all = ["# SLAB"+mod+" test"]
    commands_all.extend( command_list )
    
    result_mv = "mv results "+root+mod+".results 2>/dev/null"
    timing_mv = "for file in cism_timing*; do mv $file "+root+mod+".$file 2>/dev/null; done"
    commands_all.append(result_mv)
    commands_all.append(timing_mv)
    commands_all.append(" ")
    
    if not args.setup_only:
        if not args.quiet: 
            print("\nRunning CISM slab test")
            print(  "======================\n")

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
            print("\nFinished running the CISM slab test")
            print(  "===================================\n")
    else:
        run_script = args.output_dir+os.sep+root+mod+".run" 
        
        run_file = open(run_script,'w') 
        
        run_file.write('#!/bin/bash \n')
        for command in commands_all:
            run_file.write(command+" \n")

        run_file.close()
        os.chmod(run_script, 0o755)   # uses an octal number!

        if not args.quiet:
            print("\nFinished setting up the CISM slab test")
            print(  "======================================")
            print(  "   To run the test, use: "+run_script)

    print("WARNING: THIS TEST CASE IS IN DEVELOPMENT. USE AT YOUR OWN RISK!")

# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())

