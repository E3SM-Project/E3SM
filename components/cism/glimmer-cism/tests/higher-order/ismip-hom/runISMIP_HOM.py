#!/usr/bin/env python

"""
Script to run the Ice Sheet Model Intercomparison Project for Higher-Order
Models (ISMIP-HOM) experiments. For more information, see the README and/or
http://homepages.ulb.ac.be/~fpattyn/ismip/. 
"""

# Authors
# -------
# Written March 2, 2010 by Glen Granzow at the University of Montana.
# Reconfigured by Joseph H Kennedy at ORNL on April 27, 2015 to work with the regression testing

import os
import sys
import errno
import subprocess
import ConfigParser 

import numpy
import netCDF
from math import tan, sin, pi, exp

defaultSizes=[20,80]
defaultExperiments=['a','c']

# Parse the command line options
# ------------------------------
import argparse
parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

def unsigned_int(x):
    """
    Allows argparse to understand unsigned integers. 
    """
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("This argument is an unsigned int type! Should be an integer greater than zero.")
    return x

parser.add_argument('-c','--config', default='./ismip-hom.config', 
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
def lower_case(mixed):
    """
    Converts a string to all lower cased. If used for a type in argparse, this 
    conversion will be applied before the possible choices are evaluated.
    """
    return mixed.lower()

parser.add_argument('--cyclic', action='store_true', 
        help="If specified, all fields (including scalars) are truly periodic across the domain. This will only be applied to experiments a and b. "
            +"Note: this model setup should NOT be compared to ISMIP-HOM results.")
parser.add_argument('-r','--experiments', nargs='*', default=defaultExperiments, choices=['a','b','c','d','e','f'], type=lower_case, metavar='{a b c d e f}',
        help="List (separated by spaces) the ISMIP-HOM experiments to run.")
parser.add_argument('--scale', type=unsigned_int, default=0, 
        help="Scales the problem size by 2**SCALE. SCALE=0 creates a 40x40 grid, SCALE=1 " 
            +"creates a 80x80 grid, and SCALE=2 creates a 160x160 grid.")
#FIXME: do we need a valid range here?
parser.add_argument('--sizes', nargs='*', default=defaultSizes, type=unsigned_int, metavar='KM', 
        help="List (separated by spaces) the domain sizes to run. Recommended sizes: 5, 10, 20, 40, 80 and 160 km. "
            +"Note: sizes will only be applied to experiments a though e. Experiment f has only one size (100 km). ")
parser.add_argument('--vertical', type=unsigned_int,
        help="Override the vertical grid size (upn) in the config file.")


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
            print("   Please edit the script to use your MPI run command, or run the model mannually with")
            print("   something like: mpirun -np 4 ./cism_driver ismip-hom.config")
            sys.exit(1)
    else:
        mpiexec = ''

    commands.append(mpiexec+driver+" "+config_name+quiet_mod)

    return commands

# the main script function
# ------------------------
def main():
    """
    Run the ISMIP-HOM tests.
    """
    
    # check that file name modifier, if it exists, starts with a '-'
    if not (args.modifier == '') and not args.modifier.startswith('-') :
        args.modifier = '-'+args.modifier
         
    # get the configuration
    # ---------------------
    scale_factor = 2 ** args.scale
    

    commands_all = []
    ran_f = False
    for experiment in args.experiments:
        for size in args.sizes:
            if experiment == 'f':
                size = 100;
                if ran_f:
                    break
                else:
                    ran_f = True
         
            try:
                config_parser = ConfigParser.SafeConfigParser()
                config_parser.read( args.config )
                
                if args.vertical:
                    nz = args.vertical
                else:
                    nz = int(config_parser.get('grid','upn'))
                
                nx = int(config_parser.get('grid','ewn'))*scale_factor
                ny = int(config_parser.get('grid','nsn'))*scale_factor
                
                file_name = config_parser.get('CF input', 'name')
                root, ext = os.path.splitext(file_name)

            except ConfigParser.Error as error:
                print("Error parsing " + args.config )
                print("   "), 
                print(error)
                sys.exit(1)
            
           
            dx = float(size)*1000./float(nx)
            dy = float(size)*1000./float(ny)

            res = str(size).zfill(4)
            if args.parallel > 0:
                mod = '-'+experiment+args.modifier+'.'+res+'.p'+str(args.parallel).zfill(3)
            else:
                mod = '-'+experiment+args.modifier+'.'+res
            
            file_name = root+mod+ext
            config_name = root+mod+'.config'
            out_name = root+mod+'.out'+ext

            # create the new config file
            # --------------------------
            if not args.quiet: 
                print("\nCreating config file: "+config_name)
            
            config_parser.set('grid','upn', str(nz))
            config_parser.set('grid', 'ewn', str(nx))
            config_parser.set('grid', 'nsn', str(ny))
            config_parser.set('grid', 'dew', str(dx))
            config_parser.set('grid', 'dns', str(dy))

            config_parser.set('CF default', 'title', 'ISMIP-HOM Experiment '+experiment.upper() )
            config_parser.set('CF input', 'name', file_name)
            config_parser.set('CF output', 'name', out_name)
            config_parser.set('CF output', 'xtype', 'double')
            
            if not args.cyclic:
                if experiment in ('a','b'):
                    offset = float(size)*1000.0 * tan(0.5 * pi/180.0)
                elif experiment in ('c','d'):
                    offset = float(size)*1000.0 * tan(0.1 * pi/180.0)
                elif experiment in ('f'):
                    offset = float(size)*1000.0 * tan(3.0 * pi/180.0)
                config_parser.set('parameters', 'periodic_offset_ew', str(offset))

            if experiment in ('c','d'):
                # These tests have beta passed in from the input file, so change option accordingly.
                config_parser.set('ho_options', 'which_ho_babc', '5')

            ##Optional: if doing experiment C, one can alternatively use the ho_babc option setup for this test case rather than passing in a beta
            #if experiment in ('c'):
            #    config_parser.set('ho_options', 'which_ho_babc', '8')

            # For test case F we need to make a few additional adjustments to the config
            if experiment in ('f'):
                # Set to efvs to be the constant value of 2336041.42829 hardcoded for this option - this corresponds to the value needed for this test case
                config_parser.set('ho_options', 'which_ho_efvs', '0')
                # 2.4 yr is the longest dt ok for diffusive CFL, when using dx=dy=2500.0
                config_parser.set('time', 'dt', '2.2')  
                # Need to run to steady-state...  
                # It's close to SS by 400 years, but there are some long-period oscillations that still appear out to 1000 yrs.  
                # Not sure yet how much longer than that to eliminate those.
                config_parser.set('time', 'tend', '400.0')
                # Include flwa, efvs and the CFL variables to the output file
                config_parser.set('CF output', 'variables', 'uvel vvel uvel_extend vvel_extend uvel_icegrid vvel_icegrid topg thk usurf wvel_ho velnorm efvs adv_cfl_dt diff_cfl_dt')  
                # we don't want to output a whole lot of time levels, but want to be able to see we've reached SS.
                config_parser.set('CF output', 'frequency', '25.0')  


            # write the config file
            with open(config_name, 'wb') as config_file:
                config_parser.write(config_file)


            # create the input netCDF file
            # ----------------------------
            #FIXME: This whole section could use a clean up. -JHK
            if not args.quiet: 
                print("Creating input netCDF file: "+file_name)
            try:
                nc_file = netCDF.NetCDFFile(file_name,'w',format='NETCDF3_CLASSIC')
            except TypeError:
                nc_file = netCDF.NetCDFFile(file_name,'w')

            
            nc_file.createDimension('time',1)
            nc_file.createDimension('x1',nx)   # unstaggered grid
            nc_file.createDimension('y1',ny)
            nc_file.createDimension('x0',nx-1) # staggered grid 
            nc_file.createDimension('y0',ny-1)
            time = nc_file.createVariable('time','f',('time',))
            x1   = nc_file.createVariable('x1','f',('x1',))
            y1   = nc_file.createVariable('y1','f',('y1',))
            x0   = nc_file.createVariable('x0','f',('x0',))
            y0   = nc_file.createVariable('y0','f',('y0',))
            thk  = nc_file.createVariable('thk' ,'f',('time','y1','x1'))
            topg = nc_file.createVariable('topg','f',('time','y1','x1'))
            if experiment in ('c','d'):
                unstagbeta = nc_file.createVariable('unstagbeta','f',('time','y1','x1'))
            time[0] = 0
            x1[:] = [(i+0.5)*dx for i in range(nx)] # unstaggered grid
            y1[:] = [(j+0.5)*dy for j in range(ny)]
            x0[:] = [(i+1)*dx for i in range(nx-1)] # staggered grid 
            y0[:] = [(j+1)*dy for j in range(ny-1)]

            #   Generate the ice thickness, bed topography, and (sometimes) 
            #   basal friction coefficient for the experiment

            thickness  = list()
            topography = list()
            basalFriction = list()

            xx = [(i+0.5)*dx for i in range(nx)]
            yy = [(j+0.5)*dy for j in range(ny)]

            if experiment in ('a','b'):
                if not args.cyclic:              
                    alpha = 0.5 * pi/180
                else: # optional flags to allow for truly periodic domain setup
                    alpha = 0.
                zz = [4000-x1[i]*tan(alpha) for i in range(nx)]
            elif experiment in ('c','d'):
                #      if not args.cyclic:              
                alpha = 0.1 * pi/180
            #      else:  # optional flags to allow for truly periodic domain setup
            #        alpha = 0.
                zz = [1000-x1[i]*tan(alpha) for i in range(nx)]
            elif experiment == 'f':
                alpha = 3.0 * pi/180
                zz = [7000-x1[i]*tan(alpha) for i in range(nx)]
                xc = (xx[0]+xx[-1])/2
                yc = (yy[0]+yy[-1])/2
                a0 = 100
                sigma2 = 10000**2

            omega = 2*pi / (size*1000)      
            for y in yy:
                row = list()
                for x in xx:
                    if experiment == 'a':
                        row.append(1000 - 500*sin(omega*x)*sin(omega*y))
                    elif experiment == 'b':
                        row.append(1000 - 500*sin(omega*x))
                    elif experiment == 'c':
            #           if args.cyclic:    # for test case w/ truly periodic domain, add some non-zero topog to force flow 
            #               row.append(1000 - 500*sin(omega*x)*sin(omega*y))
            #           elif not args.cyclic:
                        row.append(1000 + 1000*sin(omega*x)*sin(omega*y))
                    elif experiment == 'd':
                        row.append(1000 + 1000*sin(omega*x))
                    elif experiment == 'f':
                        row.append(1000 - a0*exp(-((x-xc)**2+(y-yc)**2)/sigma2))
                if experiment in ('a','b','f'):
                    thickness.append(row)
                    if not args.cyclic:        
                        topography.append([z-t for (z,t) in zip(zz,row)])
                    else:  # args to allow for truly periodic domain setup
                        topography.append([z for z in zz])
                else:
            #        basalFriction.append(row[:-1])
                    basalFriction.append(row[:])

            if experiment in ('a','b','f'):
                thk [:] = thickness
                topg[:] = topography
            elif experiment in ('c','d'):
                thk [:] = ny*[nx*[1000]]
                topg[:] = ny*[zz]
            #    if not args.cyclic:         # args to allow for truly periodic domain setup
                unstagbeta[:] = basalFriction[:]
            
            # close the new netCDF and config file, and move it to the output directory (if given)
            nc_file.close()
            mkdir_p(args.output_dir)
            subprocess.check_call("cp *rilinosOptions.xml "+args.output_dir, shell=True)
            subprocess.check_call("mv "+file_name+" "+args.output_dir, shell=True)
            subprocess.check_call("mv "+config_name+" "+args.output_dir, shell=True)


            # Run CISM
            # --------
            command_list = prep_commands(args, config_name)
            commands_all.append("# ISMIP-HOM"+mod.upper()+" test")
            commands_all.extend(command_list) 
            
            result_mv = "mv results "+root+mod+".results 2>/dev/null"
            timing_mv = "for file in cism_timing*; do mv $file "+root+mod+".$file 2>/dev/null; done"
            commands_all.append(result_mv)
            commands_all.append(timing_mv)
            commands_all.append(" ")

            if not args.setup_only:
                if not args.quiet: 
                    print("\nRunning CISM ISMIP-HOM"+mod.upper()+" test")
                    print(  "==================================\n")

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
                    print("\nFinished running the CISM ISMIP-HOM"+mod.upper()+" test")
                    print(  "===============================================\n")

    if args.setup_only:
        run_script = args.output_dir+os.sep+root+args.modifier+".run" 
        
        run_file = open(run_script,'w') 
        
        run_file.write('#!/bin/bash \n')
        for command in commands_all:
            run_file.write(command+" \n")

        run_file.close()
        os.chmod(run_script, 0o755)   # uses an octal number!

        if not args.quiet:
            print("\nFinished setting up the CISM ISMIP-HOM tests")
            print(  "============================================")
            print(  "   To run the test, use: "+run_script)


# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())
