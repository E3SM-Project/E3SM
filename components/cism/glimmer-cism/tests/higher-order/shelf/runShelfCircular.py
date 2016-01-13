#!/usr/bin/env python

#FIXME: More detailed description of this test case!!!
"""
Run an experiment with an idealized circular ice shelf. 
"""

# Authors
# -------
# Written by Glen Granzow at the University of Montana on April 9, 2010
# Slight alterations by SFP on 2-3-11 for Glimmer-CISM 2.0 release
# Reconfigured by Joseph H Kennedy at ORNL on April 27, 2015 to work with the regression testing

import os
import sys
import errno
import subprocess
import ConfigParser 

import numpy
import netCDF
from math import sqrt, exp

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

parser.add_argument('-c','--config', default='./shelf-circular.config', 
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
        help="Run the cism process quietly.")
parser.add_argument('-s','--setup-only', action='store_true',
        help="Set up the test, but don't actually run it.")


# Additional test specific options:
#parser.add_argument('--scale', type=unsigned_int, default=0, 
#        help="Scales the problem size by 2**SCALE. SCALE=0 creates a 31x31 grid, SCALE=1 " 
#            +"creates a 62x62 grid, and SCALE=2 creates a 124x124 grid.")
parser.add_argument('-b','--beta', action='store_true',
        help="Use a Guassian function for beta.")
#FIXME: more descriptive!
parser.add_argument('-a','--alpha', action='store_true',
        help="Use a conically topped ice thickness")
#FIXME: more desciptive!
parser.add_argument('-d','--dirichlet', action='store_true',
        help="Apply Dirichlet boundary condition at the center")


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
            print("   something like: mpirun -np 4 ./cism_driver shelf-circular.config")
            sys.exit(1)
    else:
        mpiexec = ''

    commands.append(mpiexec+driver+" "+config_name+quiet_mod)

    return commands


# the main script function
# ------------------------
def main():
    """
    Run the circular shelf test.
    """
    
    # check that file name modifier, if it exists, starts with a '-'
    if not (args.modifier == '') and not args.modifier.startswith('-') :
        args.modifier = '-'+args.modifier
         
    # get the configuration
    # ---------------------
    #scale_factor = 2 ** args.scale
    
    try:
        config_parser = ConfigParser.SafeConfigParser()
        config_parser.read( args.config )
        
        nz = int(config_parser.get('grid','upn'))
        
        nx = int(config_parser.get('grid','ewn'))#*scale_factor
        ny = int(config_parser.get('grid','nsn'))#*scale_factor
        dx = float(config_parser.get('grid','dew'))#/float(scale_factor)
        dy = float(config_parser.get('grid','dns'))#/float(scale_factor)
        
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
    
    config_parser.set('grid', 'ewn', str(nx))
    config_parser.set('grid', 'nsn', str(ny))
    config_parser.set('grid', 'dew', str(dx))
    config_parser.set('grid', 'dns', str(dy))

    config_parser.set('CF input', 'name', file_name)
    config_parser.set('CF output', 'name', out_name)
    config_parser.set('CF output', 'xtype', 'double')
    
    with open(config_name, 'wb') as config_file:
        config_parser.write(config_file)



    # create the input netCDF file
    # ----------------------------
    if not args.quiet: 
        print("\nCreating circular shelf netCDF file: "+file_name)
    try:
        nc_file = netCDF.NetCDFFile(file_name,'w',format='NETCDF3_CLASSIC')
    except TypeError:
        nc_file = netCDF.NetCDFFile(file_name,'w')

    nc_file.createDimension('time',1)
    nc_file.createDimension('x1',nx)
    nc_file.createDimension('y1',ny)
    nc_file.createDimension('level',nz)
    #nc_file.createDimension('staglevel',nz-1)
    #nc_file.createDimension('stagwbndlevel',nz+1)
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
    beta = numpy.zeros([1,ny-1,nx-1],dtype='float32')


    # Calculate the topography
    topg[0,:,:] = -2000
    
    
    # Calculate the thickness and beta
    # Domain size
    Lx = nx*dx
    Ly = ny*dy

    for i in range(nx):
        x = float(i)/(nx-1) - 0.5                        # -1/2 < x < 1/2 
        xx = x*Lx                                        # -L/2 < xx < L/2
        for j in range(ny):
            y = float(j)/(ny-1) - 0.5                    # -1/2 < y < 1/2
            yy = y*Ly                                    # -L/2 < yy < L/2
            
            r = sqrt(x*x+y*y)                            # radial distance from the center

            if r < 0.44:                                 # Inside a circle we have
                thk[0,j,i] = 1000                        # constant ice thickness unless
                if args.alpha:                           # command line option specifies
                    thk[0,j,i] *= (1-r)                  # conical top
            
            if args.beta and (i < nx-1) and (j < ny-1):  # use Gaussian on staggered grid if specified
                beta[0,j,i] = 1.0 + 1.0e10*exp(-(xx*xx+yy*yy)/5.0e5) # Gaussian

    if not args.beta:        # Don't use a Gaussian
        beta[0,:,:] =  0                             # beta is 0 almost everywhere
        beta[0,ny/2-1:ny/2+1,nx/2-1:nx/2+1] = 1.0e8  # but large in the center

    # Add a single bedrock spike in the domain center, to "ground" shelf for 
    # bisicles dycore
    topg[0,(ny-1)/2-1:(ny-1)/2+2,(nx-1)/2-1:(nx-1)/2+2] = -880. 

    # Create the required variables in the netCDF file.
    nc_file.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
    nc_file.createVariable('topg','f',('time','y1','x1'))[:] = topg
    nc_file.createVariable('beta','f',('time','y0','x0'))[:] = beta
    
    # Dirichlet BCs at center
    if args.dirichlet:  
        uvelbc = nc_file.createVariable('uvelbc','f',('time','level','y0','x0'))
        vvelbc = nc_file.createVariable('vvelbc','f',('time','level','y0','x0'))
        
        bc = numpy.empty(1,[ny-1,nx-1],dtype='float32')
        
        # boundary condition is NaN almost everywhere
        bc[0,:,:] = float('NaN')
        
        # boundary condition is 0 in the center
        bc[0,ny/2-1:ny/2+2,nx/2-1:nx/2+2] = 0
        for k in range(nz): # loop over levels
            uvelbc[0,k,:,:] = bc
            vvelbc[0,k,:,:] = bc

    nc_file.close()
    mkdir_p(args.output_dir)
    subprocess.check_call("cp *rilinosOptions.xml "+args.output_dir, shell=True)
    subprocess.check_call("mv "+file_name+" "+args.output_dir, shell=True)
    subprocess.check_call("mv "+config_name+" "+args.output_dir, shell=True)


    # Run CISM
    # --------
    command_list = prep_commands(args, config_name)
    commands_all = ["# SHELF-CIRCULAR"+mod+" test"]
    commands_all.extend( command_list )
    
    result_mv = "mv results "+root+mod+".results 2>/dev/null"
    timing_mv = "for file in cism_timing*; do mv $file "+root+mod+".$file 2>/dev/null; done"
    commands_all.append(result_mv)
    commands_all.append(timing_mv)
    commands_all.append(" ")
    
    if not args.setup_only:
        if not args.quiet: 
            print("\nRunning CISM circular shelf test")
            print(  "================================\n")


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
            print("\nFinished running the CISM circular shelf test")
            print(  "=============================================\n")
    else:
        run_script = args.output_dir+os.sep+root+mod+".run" 
        
        run_file = open(run_script,'w') 
        
        run_file.write('#!/bin/bash \n')
        for command in commands_all:
            run_file.write(command+" \n")

        run_file.close()
        os.chmod(run_script, 0o755)   # uses an octal number!

        if not args.quiet:
            print("\nFinished setting up the CISM circular shelf test")
            print(  "================================================")
            print(  "   To run the test, use: "+run_script)


# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())

