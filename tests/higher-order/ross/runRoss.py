#!/usr/bin/env python

#FIXME: More detailed description of this test case!!!
"""
Run an experiment with an ice "Ross". 
"""

# Authors
# -------
# Written March 18, 2010 by Glen Granzow at the University of Montana.
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

parser.add_argument('-c','--config', default='./ross.config', 
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

#optparser.add_option("-r", "--run", dest="doRun", default=False, action="store_true", help="Including this flag will run CISM.  Excluding it will cause the script to only setup the initial condition file")


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
            print("   something like: mpirun -np 4 ./cism_driver ross.config")
            sys.exit(1)
    else:
        mpiexec = ''

    commands.append(mpiexec+driver+" "+config_name+quiet_mod)

    return commands


# Hard coded test specific parameters
# -----------------------------------
#FIXME: Some of these could just be options!

create_files = True
use_inlets   = (False, True, 'reverse')[1]


# Test specific functions
# -----------------------

# Raymond yield stress
def createArray(nx,ny,data,dtype):
    # nx and ny dimensions are one more than the raw data, 
    # because we are choosing to use the raw data on the velocity grid.
    field = numpy.empty((ny-1,nx-1),dtype=dtype)
    for i in range(nx-1):
        for j in range(ny-1):
            field[j,i] = data[j][i]
    return field

def plot(variable): # used for debugging only
    from matplotlib import pyplot
    pyplot.figure()
    pyplot.imshow(variable,origin='lower',interpolation='nearest')
    pyplot.colorbar()
    pyplot.show()


# the main script function
# ------------------------
def main():
    """
    Run the Ross test.
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
        
        # nx and ny dimensions are one more than the raw data, 
        # because we are choosing to use the raw data on the velocity grid.
        if nx != 148:
            print("WARNING: ewn should be set to 148 in ross.config")
            raise ConfigParser.Error
        if ny != 112:
            print("WARNING: nsn should be set to 112 in ross.config")
            raise ConfigParser.Error
        if dx != 6822 or dy !=6822:
            print("WARNING: dew and dns should be set to 6822 in ross.config")
            raise ConfigParser.Error

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
    
    nz = int(config_parser.get('grid','upn'))
    config_parser.set('grid', 'ewn', str(nx))
    config_parser.set('grid', 'nsn', str(ny))
    config_parser.set('grid', 'dew', str(dx))
    config_parser.set('grid', 'dns', str(dy))

    config_parser.set('CF input', 'name', file_name)
    config_parser.set('CF output', 'name', out_name)
    config_parser.set('CF output', 'xtype', 'double')
    
    with open(config_name, 'wb') as config_file:
        config_parser.write(config_file)


    # Read the input files
    # --------------------
    if create_files:
        # Read the main data file into a dictionary mapping names to lists (of lists)
        # The dictionary keys are the headers that begin each section in the file
        inputname = os.path.join('data','111by147Grid.dat')
        inputfile = open(inputname)
        if not args.quiet:
            print("\nReading: "+inputname)
        
        currentKey = None
        currentList = list()
        data = dict()
        for line in inputfile:
            line = line.strip()
            if line.startswith('#'):
                if currentKey != None:
                    data[currentKey] = currentList
                    currentList = list()
                currentKey = line[1:].strip().lower()
            elif len(line) > 0:
                if line.find('.') > 0:
                    currentList.append([float(x) for x in line.split()])
                else:
                    currentList.append([int(x) for x in line.split()])
        data[currentKey] = currentList
        inputfile.close()

        if not args.quiet:
            print("The "+str(len(data.keys()))+" data fields read from 111by147Grid.dat are:")
            print(data.keys())


        #NOTE: optional.
        # Manually edit the existency mask so that it extends to be connected to the inlets in inlets.dat 
        data['existency table:'][109][80]=1
        data['existency table:'][97][101]=1
        data['existency table:'][98][100]=1
        data['existency table:'][96][102]=1
        data['existency table:'][96][103]=1
        data['existency table:'][78][121]=1
        data['existency table:'][77][123]=1
        data['existency table:'][79][120]=1
        data['existency table:'][80][120]=1
        data['existency table:'][80][119]=1
        data['existency table:'][51][140]=1
        data['existency table:'][53][138]=1


        # Read the kinematic boundary conditions (kbc) mask file.  
        # This is a set of (i,j) coordinates specifying where the velocity read from
        # the data file should be used as a boundary condition for the model.
        kbc_mask = numpy.zeros((111,147), dtype='i')
        inputname = os.path.join('data','kbc.dat')
        if not args.quiet:
            print("\nReading"+inputname)
        
        inputfile = open(inputname)
        for line in inputfile:
            i,j = map(int,line.split())
            kbc_mask[i,j] = 1   # 1=where we have Dirichlet; 0=otherwise
        inputfile.close()
        

        if not args.quiet:
            print(str(numpy.sum(kbc_mask))+" points were read from kbc.dat")

        # Read in the inlets file, which specifies additional Dirichlet boundary conditions
        filename = os.path.join('data','inlets.dat')
        if not args.quiet:
            print("\nReading "+filename)

        counter = [0,0]
        inputfile = open(filename)
        for line in inputfile:
            i, j, azimuth, magnitude = line.split()
            i,j = map(int,(i,j))
            if use_inlets == 'reverse':
                magnitude,azimuth = map(float,(azimuth,magnitude))
            else:
                azimuth,magnitude = map(float,(azimuth,magnitude))
            if use_inlets:
                if not args.quiet:
                    indices = '(%d,%d):' % (i,j)
                    print("Changing azimuth at "+indices+str(data["ice velocity azimuth grid"][i][j])+"->"+str(azimuth))
                    print("Changing velocity at"+indices+str(data["ice velocity magnitude"][i][j])+"->"+str(magnitude))
                data['ice velocity azimuth grid'][i][j] = azimuth
                data['ice velocity magnitude'][i][j] = magnitude
            counter[kbc_mask[i,j]] += 1
            if not args.quiet:
                print(str(i)+","+str(j)+": "+str(kbc_mask[i,j]) )
            kbc_mask[i,j] += 2
        inputfile.close()
        print("inlets.dat contains"+str(counter[0])+"points that are not in kbc.dat")
        print("inlets.dat contains"+str(counter[1])+"points that are in kbc.dat")


        # Put the data into numpy arrays with two extra rows and columns all around
        # This reproduces a previous script's output (makerossnc.py)
        mask1        = createArray(nx,ny,data['existency table:'],dtype='i')
        azimuth      = createArray(nx,ny,data['ice velocity azimuth grid'],dtype='f')
        velocity     = createArray(nx,ny,data['ice velocity magnitude'],dtype='f')
        thickness    = createArray(nx,ny,data['thickness'],dtype='f')
        #mask2        = createArray(nx,ny,data['reliable velocity obs'],dtype='i')
        seabed       = createArray(nx,ny,data['seabed depth'],dtype='f')
        mask3        = createArray(nx,ny,data['fake ice shelf region'],dtype='i')
        #accumulation = createArray(nx,ny,data['surface accumulation'],dtype='f')
        #bbar         = createArray(nx,ny,data['flowlaw'],dtype='f')
        #temperature  = createArray(nx,ny,data['surface temperature'],dtype='f')
        kbc          = createArray(nx,ny,kbc_mask,dtype='i')

        # Remove any parts of the "fake shelf" that are not in the "existency table"
        if not args.quiet:
            print("Removing"+str(numpy.sum(numpy.logical_and(mask1==0,mask3!=0)))+"points from mask3")
        mask3[mask1==0] = 0

        # Set the velocity to zero except where needed as a kinematic boundary condition.
        velocity[kbc == 0] = 0.0
        # Get the components of the velocity vector
        #NOTE: velocity1 is vvel, velocity 2 is uvel
        azimuth *= numpy.pi/180.0
        velocity1 = velocity * numpy.cos(azimuth)
        velocity2 = velocity * numpy.sin(azimuth)


        # Create a netCDF file containing the raw data
        # --------------------------------------------
        #NOTE: This is not necessary, but can be useful to debug the rest of the script, if needed.

        #filename = os.path.join('output','raw.nc')
        #print("\nWriting "+filename)
        #if netCDF_module == 'netCDF4':
        #    netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
        #else:
        #    netCDFfile = NetCDFFile(filename,'w')
        #
        #netCDFfile.createDimension('x',data['rows columns number of sub parameters'][0][1])
        #netCDFfile.createDimension('y', data['rows columns number of sub parameters'][0][0])
        #netCDFfile.createVariable('x','f',('x',))[:] = [x[0] for x in data['columns position'][:-1]]
        #netCDFfile.createVariable('y','f',('y',))[:] = [x[0] for x in data['rows position'][:-1]]
        #netCDFfile.createVariable('mask1',    'i',('y','x'))[:] = data['existency table:']
        #netCDFfile.createVariable('azimuth',  'f',('y','x'))[:] = data['ice velocity azimuth grid']
        #netCDFfile.createVariable('velocity', 'f',('y','x'))[:] = data['ice velocity magnitude']
        #netCDFfile.createVariable('thickness','f',('y','x'))[:] = data['thickness']
        #netCDFfile.createVariable('mask2',    'i',('y','x'))[:] = data['reliable velocity obs']
        #netCDFfile.createVariable('seabed',   'f',('y','x'))[:] = data['seabed depth']
        #netCDFfile.createVariable('mask3',    'i',('y','x'))[:] = data['fake ice shelf region']
        #netCDFfile.createVariable('accumulation','f',('y','x'))[:] = data['surface accumulation']
        #netCDFfile.createVariable('bbar',        'f',('y','x'))[:] = data['flowlaw']
        #netCDFfile.createVariable('temperature', 'f',('y','x'))[:] = data['surface temperature']
        #numpy.set_printoptions(threshold='nan')
        #netCDFfile.createVariable('kbc',         'i',('y','x'))[:] = kbc_mask
        #netCDFfile.close()
        #del(netCDFfile) # remove this variable from the name-space (pycdf might fail if we don't)


        ##NOTE: optional plot of kinematic bc positions
        #import matplotlib.pyplot as plt
        #plt.imshow(kbc_mask, interpolation='nearest', origin='lower')
        ##plt.imshow(mask1[:,:], interpolation='nearest', origin='lower')
        #for i in range(kbc_mask.shape[1]):
        #    for j in range(kbc_mask.shape[0]):
        #        if kbc_mask[j,i] == 1:
        #            plt.plot(i,j,'og')  # big inlets
        #        if kbc_mask[j,i] == 2:
        #            plt.plot(i,j,'oc')  # data inlets
        #        if data['ice velocity magnitude'][j][i] != 0.0:
        #            plt.plot(i,j,'xk')  # nonzero velo
        #plt.colorbar()
        #plt.axis('equal')
        #plt.show()


    # create the input netCDF file
    # ----------------------------
    if not args.quiet: 
        print("\nCreating Ross netCDF file: "+file_name)
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
    beta = numpy.zeros((ny-1,nx-1),dtype='float32') 
    uvel = numpy.array(nz*[velocity2])
    vvel = numpy.array(nz*[velocity1])
  
    mask = numpy.logical_and(velocity==0,numpy.logical_or(mask1==1,mask3==1))
    kinbcmask = numpy.int32(numpy.where(mask, 0, 1))
  
    # interpolate thk and topg from the staggered grid onto the scalar grid
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            #assert False
            if numpy.count_nonzero(mask3[j-1:j+1, i-1:i+1]) == 0:
                thk[0,j,i] = thickness[j-1:j+1, i-1:i+1].mean()
            else:
                thk[0,j,i] = 0.0
            topg[0,j,i] = -seabed[j-1:j+1, i-1:i+1].mean()
    # The above line averages the supplied topography from the staggered grid to the unstaggered grid.
    # However, it results in a grounded cell (or maybe cells) in the region specified as the ice shelf.
    # This occurs in a little region just downstream of Roosevelt Island and causes very rapid velocities
    # in that region because it creates a grounded region with a steep slope and beta of 0.
    # Because the test case is only solved for the shelf itself, and there are Dirichlet velocity b.c.
    # at all grounding lines, the topography is not technically needed.  
    # Therefore it is easier to simpy specify a really deep basal topography everywhere to avoid
    # any inadvertent groundings of the ice.
    #topg[:] = -5000.0    
    topg[0,22:42,24:47] -= 500.0

    # extrapolate the edges
    topg[0,0,:] = topg[0,1,:]
    topg[0,-1,:] = topg[0,-2,:]
    topg[0,:,0] = topg[0,:,1]
    topg[0,:,-1] = topg[0,:,-2]
    thk[0,0:2,:] = 0.0  # no ice along bottom two rows of domain to ensure ice shelf front extends to Ross Island side.
    thk[0,-1,:] = thk[0,-2,:]
    thk[0,:,0] = thk[0,:,1]
    thk[0,:,-1] = thk[0,:,-2]
  
    # Create the required variables in the netCDF file.
    nc_file.createVariable('uvel','f',('time','level','y0','x0'))[:] = uvel
    nc_file.createVariable('vvel','f',('time','level','y0','x0'))[:] = vvel
    nc_file.createVariable('kinbcmask','i',('time','y0','x0'))[:] = kinbcmask
    nc_file.createVariable('beta','f',('time','y0','x0'))[:] = beta 
    nc_file.createVariable('thk' ,'f',('time','y1','x1'))[:] = thk
    nc_file.createVariable('topg','f',('time','y1','x1'))[:] = topg
  

    ##NOTE: optional plot of kinematic bc positions
    #import matplotlib.pyplot as plt
    #fig = plt.figure(2, facecolor='w', figsize=(10, 4), dpi=100)
    #
    ##plt.imshow(kinbcmask[0,:,:], interpolation='nearest', origin='lower')
    #plt.imshow(mask1[:,:], interpolation='nearest', origin='lower')
    #for i in range(kbc.shape[1]):
    #    for j in range(kbc.shape[0]):
    #        if kbc[j,i] == 1:
    #            plt.plot(i,j,'og')  # big inlets
    #        if kbc[j,i] == 2:
    #            plt.plot(i,j,'oc')  # data inlets
    #        if velocity[j,i] != 0.0:
    #            plt.plot(i,j,'xk')  # nonzero velo
    #plt.colorbar()
    #plt.axis('equal')
    #plt.show()

    nc_file.close()
    mkdir_p(args.output_dir)
    subprocess.check_call("cp *rilinosOptions.xml "+args.output_dir, shell=True)
    subprocess.check_call("mv "+file_name+" "+args.output_dir, shell=True)
    subprocess.check_call("mv "+config_name+" "+args.output_dir, shell=True)


    # Run CISM
    # --------
    command_list = prep_commands(args, config_name)
    commands_all = ["# ROSS"+mod+" test"]
    commands_all.extend( command_list )
    
    result_mv = "mv results "+root+mod+".results 2>/dev/null"
    timing_mv = "for file in cism_timing*; do mv $file "+root+mod+".$file 2>/dev/null; done"
    commands_all.append(result_mv)
    commands_all.append(timing_mv)
    commands_all.append(" ")
    
    if not args.setup_only:
        if not args.quiet: 
            print("\nRunning CISM Ross test")
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
            print("\nFinished running the CISM Ross test")
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
            print("\nFinished setting up the CISM Ross test")
            print(  "======================================")
            print(  "   To run the test, use: "+run_script)


# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())

