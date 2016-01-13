#!/usr/bin/env python

#FIXME: More detailed description of this test case!!!
"""
Run an experiment with an ice "stream". 
"""

# Authors
# -------
# Original author unlisted.
# Reconfigured by Joseph H Kennedy at ORNL on April 27, 2015 to work with the regression testing

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

parser.add_argument('-c','--config', default='./stream.config', 
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
#parser.add_argument('--scale', type=unsigned_int, default=0, 
#        help="Scales the problem size by 2**SCALE. SCALE=0 creates a 31x31 grid, SCALE=1 " 
#            +"creates a 62x62 grid, and SCALE=2 creates a 124x124 grid.")
parser.add_argument('-z','--stream-size', type=unsigned_int, default=25,
        help=") The number of grid cells used to model the ice stream portion of the domain."
            +"Note: values <19 may not work properly for all problems.")    
#optparser.add_option('-s','--stream-size',dest='stream_grid_size',default=25,type='int',help='Number of cells to use to model the ice stream portion of the domain (values <19 may not work properly for all problems).')

parser.add_argument('--vertical', type=unsigned_int,
        help="Override the vertical grid size (upn) in the config file.")
#optparser.add_option('-v','--vert-grid-size',dest='vertical_grid_size',default=2,type='int',help='Number of vertical layers to use (upn); minimum value = 2')

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
            print("   something like: mpirun -np 4 ./cism_driver stream.config")
            sys.exit(1)
    else:
        mpiexec = ''

    commands.append(mpiexec+driver+" "+config_name+quiet_mod)

    return commands


# Hard coded test specific parameters
# -----------------------------------
#FIXME: Some of these could just be options!

analytic_solution = 'raymond'  # can be 'raymond' or 'schoof'
kinflag = 1    # 1=apply kinematic bc (analytic soln) at points in the domain (discussed further below); 0=the run will be doubly periodic (preferred)
fillInitialGuess = 0  # 1=use the analytic solution as the initial guess for the velocity solver to speed convergence; 0=use the default 0-velocity initial guess

# Domain parameters
streamHalfWidth = 25000.0   # ice stream half-width, in m - used for both raymond & schoof formulations
alongFlowLength = 30000.0   # the desired along-flow length of the domain, in m; set to -1 to get a square domain
H = 1000.0       # ice thickness
dsdx = -1.0e-3   # bed (and surface) slope in the x-direction (y-direction is flat)

# Physical parameters
rho = 910.0   # ice density kg/m3
g = -9.81     # gravity m/s2
n = 3         # flow law exponent
A = 1e-16     # flow rate factor in Pa^-3 yr^-1

# schoof solution parameters
m = 1.55  # schoof exponent
L = streamHalfWidth / (m+1.0)**(1.0/m)  # This comes from the line above eq. 4.3 in schoof (2006)

taud = rho * g * H * dsdx  # Driving stress
# Calculate a good size for the size of the domain outside of the stream (in m)
if analytic_solution == 'raymond':
    strongWidth = 5.0 * H  # 5 ice thicknesses should get us beyond the zone of lateral stress transfer.  Adjust as needed
elif analytic_solution == 'schoof':
    # schoof (2006) uses a domain size that is 3L on either side of the central axis
    strongWidth = 3.0 * L - streamHalfWidth


# Test specific functions
# -----------------------

# raymond yield stress
def raymond_tau(yy):
    tau0 = 5.2e3*numpy.ones(yy.shape)         # set the stream value everywhere
    tau0[numpy.absolute(yy)>=streamHalfWidth] = 0.7e5        # set a very large value  outside the stream
    return tau0

# raymond velocity solution
def raymond_uvel(yy):
    tau0r = raymond_tau(yy)
    tau0r[tau0r>taud] = taud
    ur = 2.0 * A / (n+1.0) * ( (taud - tau0r)/H )**n * ( streamHalfWidth**(n+1) - numpy.absolute(yy)**(n+1) )
    ur[ur<0.0] = 0.0
    return ur

# schoof yield stress distribution
def schoof_tau(yy):
    return taud * numpy.absolute( yy / L )**m

# schoof velocity solution
def schoof_uvel(yy):
    B = A**(-1.0/n)
    us = -2.0 * taud**3 * L**4 / (B**3 * H**3) * ( ((yy/L)**4 - (m+1.0)**(4.0/m))/4.0 - 3.0*( numpy.absolute(yy/L)**(m+4.0) \
    - (m+1.0)**(1.0+4.0/m) )/((m+1.0)*(m+4.0)) + 3.0*( numpy.absolute(yy/L)**(2.0*m+4.0) - (m+1.0)**(2.0+4.0/m) )/((m+1.0)**2*(2.0*m+4.0)) \
    - ( numpy.absolute(yy/L)**(3.0*m+4.0) - (m+1.0)**(3.0+4.0/m) )/ ( (m+1.0)**3*(3.0*m+4.0)) )

    # Some adjustments to the analytic profile - not entirely sure why these are needed.
    ind = numpy.nonzero( numpy.absolute(yy) >= streamHalfWidth )
    us[ind] = 0.0

    return us


# the main script function
# ------------------------
def main():
    """
    Run the stream test.
    """

    # check that file name modifier, if it exists, starts with a '-'
    if not (args.modifier == '') and not args.modifier.startswith('-') :
        args.modifier = '-'+args.modifier
         
    # get the configuration
    # ---------------------
    try:
        config_parser = ConfigParser.SafeConfigParser()
        config_parser.read( args.config )
        
        if args.vertical:
            nz = args.vertical
        else:
            nz = int(config_parser.get('grid','upn'))
        
        file_name = config_parser.get('CF input', 'name')
        root, ext = os.path.splitext(file_name)

    except ConfigParser.Error as error:
        print("Error parsing " + args.config )
        print("   "), 
        print(error)
        sys.exit(1)
    
    # Setup the domain
    # ----------------
    nStream = args.stream_size
    # Check domain sizes for usefulness
    if (nStream % 2) == 0 and analytic_solution == 'schoof':
        print("Warning: For the schoof version, you might want the number of cells in the stream to be an odd number.")
    
    dy = 2.0 * streamHalfWidth / float(nStream)
    dx = dy  # always want this
    
    # Figure out the number of cells we need to add to get as close t0 the 
    # desired width of the strong region as possible (note: may want to use 
    # ceil() instead of round() here)
    nStrongStrip = int(round(strongWidth / dy))  

    # a +1 is needed to convert from y0 to y1 but we leaving it off lets the 
    #stream boundaries fall on the y0 grid, which is needed to best match the 
    #analytic solution
    ny = nStream + 2 * nStrongStrip  
    if alongFlowLength < 0:
      nx = ny  # square domain
    else:
      nx = int(round(alongFlowLength / dx))

    offset = -dsdx * dx * nx
    
    if not args.quiet:
        print("\nDomain setup:")
        print(  "=============" )
        print(  "Number of cells for stream (N-S): "+str(nStream))
        print(  "Number of cells for entire domain (N-S): "+str(ny))
        print(  "dy=dx= "+str(dy))
        print(  "Domain N-S (across-flow) width (m): "+str(ny*dy))
        print(  "Domain E-W (along-flow) width (m): "+str(nx*dx))

    res = str(nStream).zfill(4)
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
        print("\nCreating stream netCDF file: "+file_name)
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


    x1 = dx*numpy.arange(nx,dtype='float64')
    y1 = dy*numpy.arange(ny,dtype='float64') - dy*float(ny-1)/2.0  # make the y-coordinates centered about 0

    x0 = dx/2.0 + x1[:-1] # staggered grid
    y0 = dy/2.0 + y1[:-1]

    # Make sure the edge of the stream lands on the grid cells on the y0 grid.  
    # This should always happen with the logic above, so this check should never be activated.
    if (analytic_solution == 'raymond') and (not True in (numpy.absolute(streamHalfWidth-y0) < 0.0001)):
        print("\nERROR: the stream edge does not land on the y0 grid so the stream will "
               +"not be resolved adequately for the raymond case. Adjust the domain size, "
               +"stream size, and/or horizontal resolution.")
        print(  "   Stream half width = "+str(streamHalfWidth))
        print(  "   y0 grid has values at: ")
        print(  "      "+str(y0[:]))
        sys.exit(1)

    # Make sure we have at least two non-stream rows on each side
    if (numpy.absolute(y0[:])>streamHalfWidth).sum() < 4:
        print("\nERROR: there are less than two non-stream rows on each side of the stream."
               +"Adjust the domain size, stream size, and/or horizontal resolution.")
        print(  "   Stream half width = "+str(streamHalfWidth))
        print(  "   y0 grid has values at: ")
        print(  "      "+str(y0[:]))
        sys.exit(1)

    nc_file.createVariable('time','f',('time',))[:] = [0]
    nc_file.createVariable('x1','f',('x1',))[:] = numpy.float32(x1)
    nc_file.createVariable('y1','f',('y1',))[:] = numpy.float32(y1)
    nc_file.createVariable('x0','f',('x0',))[:] = numpy.float32(x0) # staggered grid
    nc_file.createVariable('y0','f',('y0',))[:] = numpy.float32(y0)


    # Calculate values for the required variables.
    thk  = numpy.zeros([1,ny,nx],dtype='float32')
    topg = numpy.zeros([1,ny,nx],dtype='float32')
    tauf = numpy.zeros([1,ny-1,nx-1],dtype='float32')

    # Calculate input field values
    thk[:] = H  # constant thickness

    for j in range(ny):
        topg[0,j,:] = 1000.0 + dsdx * x1[:]   # sloped bed.  add 1000.0 to stay well above sea level

    if analytic_solution == 'raymond':
        tau0Profile = raymond_tau(y0)
        uvelProfile = raymond_uvel(y0)
    elif analytic_solution == 'schoof':
        tau0Profile = schoof_tau(y0)
        uvelProfile = schoof_uvel(y0)
    else:
        print("\nERROR: Invalid value for 'analytic_solution'.")
        sys.exit(1)

    for i in range(nx-1):
        tauf[0,:,i] = tau0Profile


    # =======================================
    # Save the required variables to the netCDF file.
    nc_file.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
    nc_file.createVariable('topg','f',('time','y1','x1'))[:] = topg
    nc_file.createVariable('tauf','f',('time','y0','x0'))[:] = tauf

    if kinflag == 1 or fillInitialGuess == 1:
        nc_file.createVariable('uvel','f',('time','level','y0','x0'))
        nc_file.createVariable('vvel','f',('time','level','y0','x0'))

    if kinflag == 1:
        # setup Dirichlet boundary conditions for uvel and/or vvel at points in the domain

        dudy = numpy.gradient( uvelProfile, dy )
        vvelProfile = -dudy*dy

        kinbcmask = numpy.zeros([1,ny-1,nx-1],dtype='int32')
        uvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
        vvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')


    # =================================================================
    # fill both uvel and vvel at the upstream and downstream domain ends

        # Fill first column
#        i = 0
#        uvel[0,:,:,i] = numpy.tile(uvelProfile, [nz, 1])  # uniform in the vertical
#        vvel[0,:,:,i] = -numpy.tile(vvelProfile, [nz, 1])  # uniform in the vertical
#        kinbcmask[0,:,i] = 1

        # Fill last column
#        i = nx-1 - 1
#        uvel[0,:,:,i] = numpy.tile(uvelProfile, [nz, 1])  # uniform in the vertical
#        vvel[0,:,:,i] = numpy.tile(vvelProfile, [nz, 1])  # uniform in the vertical
#        kinbcmask[0,:,i] = 1

    # =================================================================
    # fill both uvel and vvel at the upstream and downstream domain ends
        # Fill just a single across-flow profile in domain interior
        i = 2
        uvel[0,:,:,i] = numpy.tile(uvelProfile, [nz, 1])  # uniform in the vertical
#        vvel[0,:,:,i] = -numpy.tile(vvelProfile, [nz, 1])  # uniform in the vertical
        kinbcmask[0,:,i] = 1

        nc_file.variables['uvel'][:] = uvel[:]
        nc_file.variables['vvel'][:] = vvel[:]
        nc_file.createVariable('kinbcmask','i',('time','y0','x0'))[:] = kinbcmask[:]

    if fillInitialGuess == 1:
        # Fill the analytic solution into the initial guess to speed convergence
        dudy = numpy.gradient( uvelProfile, dy )
        vvelProfile = -dudy*dy

        uvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
        vvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')

        for i in range(nx-1):
            uvel[0,:,:,i] = numpy.tile(uvelProfile, [nz, 1])  # uniform in the vertical
            vvel[0,:,:,i] = numpy.tile(vvelProfile, [nz, 1])  # uniform in the vertical

        nc_file.variables['uvel'][:] = uvel[:]
        nc_file.variables['vvel'][:] = vvel[:]

    nc_file.close()
    mkdir_p(args.output_dir)
    subprocess.check_call("cp *rilinosOptions.xml "+args.output_dir, shell=True)
    subprocess.check_call("mv "+file_name+" "+args.output_dir, shell=True)
    subprocess.check_call("mv "+config_name+" "+args.output_dir, shell=True)

    # Run CISM
    # --------
    command_list = prep_commands(args, config_name)
    commands_all = ["# STREAM"+mod+" test"]
    commands_all.extend( command_list )
    
    result_mv = "mv results "+root+mod+".results 2>/dev/null"
    timing_mv = "for file in cism_timing*; do mv $file "+root+mod+".$file 2>/dev/null; done"
    commands_all.append(result_mv)
    commands_all.append(timing_mv)
    commands_all.append(" ")
    
    if not args.setup_only:
        if not args.quiet: 
            print("\nRunning CISM stream test")
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
            print("\nFinished running the CISM stream test")
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
            print("\nFinished setting up the CISM stream test")
            print(  "========================================")
            print(  "   To run the test, use: "+run_script)


# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())


