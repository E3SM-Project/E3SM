#!/usr/bin/env python
# This script runs the mismip experiment on a linear downward sloping bed.
# Files are written in the "output" subdirectory.

#### TO BE REWRITTEN
# The script performs the following two steps:
# 1. Created an initial profie for the mismip experiment.
# 2. Create a netCDF input file for Glimmer.

# Parse options
from optparse import OptionParser
optparser = OptionParser()
optparser.add_option("-c", "--config", dest="configfile", type='string', default='mismipInit.config', help="Name of .config file to use to setup and run the mismip experiment", metavar="FILE")
optparser.add_option('-m','--parallel',dest='parallel',type='int', help='Number of processors to run the model with: if specified then execute run in parallel', metavar="NUMPROCS")
optparser.add_option('-e','--exec',dest='executable',default='./cism_driver',help='Set path to the CISM executable')
for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = optparser.parse_args()

import sys, os, numpy
from netCDF import *
from math import sqrt
from ConfigParser import ConfigParser



# =====================================
# Create a netCDF file according to the information in the config file.
try:
    parser = ConfigParser()
    parser.read(options.configfile)
    nx = int(parser.get('grid','ewn'))
    ny = int(parser.get('grid','nsn'))
    nz = int(parser.get('grid','upn'))
    dx = float(parser.get('grid','dew'))
    dy = float(parser.get('grid','dns'))
    A = float(parser.get('parameters','default_flwa'))
#    C = float(parser.get('parameters','coulomb_c'))
    p = float(parser.get('parameters','p_ocean_penetration'))
    filename = parser.get('CF input', 'name')
except:
    sys.exit('Error parsing ' + options.configfile)

print 'Writing', filename
try:
  netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
except TypeError:
  netCDFfile = NetCDFFile(filename,'w')

netCDFfile.createDimension('time',1)
netCDFfile.createDimension('x1',nx)
netCDFfile.createDimension('y1',ny)
netCDFfile.createDimension('level',nz)
netCDFfile.createDimension('staglevel',nz-1)
netCDFfile.createDimension('stagwbndlevel',nz+1)
netCDFfile.createDimension('x0',nx-1) # staggered grid 
netCDFfile.createDimension('y0',ny-1)

x = dx*numpy.arange(nx,dtype='float32')
y = dx*numpy.arange(ny,dtype='float32')

netCDFfile.createVariable('time','f',('time',))[:] = [0]
netCDFfile.createVariable('x1','f',('x1',))[:] = x
netCDFfile.createVariable('y1','f',('y1',))[:] = y
netCDFfile.createVariable('x0','f',('x0',))[:] = dx/2 + x[:-1] # staggered grid
netCDFfile.createVariable('y0','f',('y0',))[:] = dy/2 + y[:-1]

thk  = numpy.zeros([1,ny,nx],dtype='float32')
topg = numpy.zeros([1,ny,nx],dtype='float32')
beta = numpy.zeros([1,ny-1,nx-1],dtype='float32')
acab = numpy.zeros([1,ny,nx],dtype='float32')
uvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
vvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
kinbcmask = numpy.zeros([1,ny-1,nx-1],dtype='int32')
effecpress = numpy.zeros([1,ny,nx],dtype='float32')
tau_b = numpy.zeros([1,ny-1,nx-1],dtype='float32')

print 'p_ocean=', p

# Calculating the linear bed topography from Schoof 2007
def computeBedLinear(x):
    xBar = 1.e6
    HBar = 1.e3
    schoofx = 750000.
    slope = -778.5
    b0 = 720.

    eps_b = 1e-10
    abs_x = numpy.sqrt(x**2 + eps_b**2)
    xprime = abs_x*xBar/schoofx
    b = -(b0 + slope*xprime)/HBar

    return b


deltaX = dx/1.e6
xc = (nx-1)*deltaX
xu = numpy.linspace(0,xc,nx)
xH = xu - 0.5*deltaX
bed = computeBedLinear(xH)

    
# We are using non-dimensionalized equations at first. We will re-dimensionalized at when saving for input netcd file
# parameters for initialization (maybe to read from config file?)
sPerY = 365.25*24.*3600.
HBar = 1000.
xBar = 1000000.
aBar = .3/sPerY
uBar = aBar*xBar/HBar
tBar = xBar/uBar

rho_i = 910.
rho_w = 1028.
delta = 1.0 - rho_i/rho_w
g = 9.81
n = 3.0
a = 1.0  # (non-dimensionalized accumulation)
Ab = 3.1688e-24 # ice temperature at the bed.
A = 4.6416e-24 # first MISMIP exp value.
Cs = 7.624e6 # C-schoof 2007.

lambda_0 = 2.
m_0 = 0.5
linearSlope = 778.5

NBar = rho_i*g*HBar
taudBar = NBar*HBar/xBar
taubBar = Cs*uBar**(1./n)
Kappa = (m_0*uBar)/(lambda_0*Ab)/(NBar**n)
gamma = taubBar/taudBar

taulBar = (A**(-1./n))*(uBar/xBar)**(1./n)*HBar/xBar
epsilon = taulBar/(2*taudBar)

print "kappa=", Kappa
print "gamma=", gamma
print "Epsilon=", epsilon

toleranceInner = 1e-3

xg_init = 0.9


# Creating initial profile using shallow ice approximation (balance between tau_b and tau_d)

glIndex = numpy.argmin(numpy.abs(xH-xg_init))
xg = xH[glIndex]
Hxg = bed[glIndex]/(1.0-delta)


print "dx=", deltaX
print "xH(xg)=", xH[glIndex]
print "Hxg=",Hxg

print "xg=",xg

if (Hxg <= 0.):
  raise ValueError("Hf(xg) <= 0. Cannot produce H profile")
  print "xg = ",xg
  print "Hxg =",Hxg

uxg = a*xg/Hxg
c1 = (delta*a/(8*epsilon))**n
#operand = numpy.maximum(c1*(xH**(n+1) - xg**(n+1)),0.0) + uxg**(n+1)
operand = c1*(xH**(n+1) - xg**(n+1)) + uxg**(n+1)


HShelf = a*xH*(operand)**(-1/(n+1))
uShelf = a*xH/HShelf
    
H = HShelf.copy()
u = uShelf.copy()

#Hip1 = H[glIndex+1]
#for xIndex in range(glIndex,-1,-1):
#   Hi = Hip1
#   deltaB = (bed[glIndex+1]-bed[glIndex])
#   for iterIndex in range(100):
#     Hmid = 0.5*(Hi+Hip1)
#     umid = a*xu[xIndex]/Hmid
#     taub = -gamma*umid**(1./n)
#     HPrev = Hi
#     Hi = 0.5*(-deltaB + numpy.sqrt(deltaB**2-4*(-Hip1**2 +Hip1*deltaB + 2*deltaX*taub)))
#     deltaH = numpy.abs(Hi-HPrev)
#     if(deltaH < toleranceInner):
#        break
#   #print "deltaH:", Hi-HPrev, Hi, Hip1
#   Hip1 = Hi
#   H[xIndex] = Hi
#   u[xIndex] = umid


for xIndex in range(glIndex,-1,-1):
   Hip1 = H[xIndex+1]
   xx = xH[xIndex+1]
   H[xIndex] = -(bed[xIndex+1]-bed[xIndex])+deltaX*gamma*a*xx*numpy.abs(a*xx)**(1/n-1)/Hip1**(1/n+1)+Hip1

Hf = numpy.max(bed/(1-delta),0)
Fp = numpy.max((1-Hf/H),0)**p
N_effec = H*Fp
tauB = gamma*u**(1/n)*(N_effec**n/(Kappa*u + N_effec**n))**(1/n)


# Re-dimensionalize the variables for output file writing
H = H*HBar
u = u*uBar
xH = xH*HBar
xu = xu*HBar
bed_topo = bed*HBar
N_effec = N_effec*NBar
accumulation = 0.3
tauB = tauB*taubBar

# Calculating the ice thickness and velocity and the bed topo over the entire
for j in range(ny):
  topg[0,j,:] = -bed_topo     # takes into account the silly sign convention from Schoof 2007
#  effecpress[0,j,:] = N_effec
  thk[0,j,3:-3] = H[:-6]      # leaving 3 cells on both ends to be empty
  effecpress[0,j,3:-3] = N_effec[:-6] 
  acab[0,j,:] = accumulation  
#  tau_b[0,j,:] = tauB[:-1] 
  
acab[:,:,0:3] = 0.0   # no ice accumulates in the three first cells
acab[:,:,-3:] = 0.0   # no ice accumulates in the three last cells
kinbcmask[:,:,2] = 1


# Calculate tempstag and beta values, if desired.  See lines below to enable them being written to file.
tau_b[:] = tauB[:-1]


# Create the required variables in the netCDF file.
netCDFfile.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
netCDFfile.createVariable('topg','f',('time','y1','x1'))[:] = topg
netCDFfile.createVariable('acab', 'f',('time','y1','x1'))[:] = acab
netCDFfile.createVariable('uvel','f',('time','level','y0','x0'))[:] = uvel[:]
netCDFfile.createVariable('vvel','f',('time','level','y0','x0'))[:] = vvel[:]
netCDFfile.createVariable('kinbcmask','i',('time','y0','x0'))[:] = kinbcmask[:]
netCDFfile.createVariable('effecpress','f',('time','y1','x1'))[:] = effecpress
netCDFfile.createVariable('tau_b','f',('time','y0','x0'))[:] = tau_b


# Optional fields that could be added to the initial condition file.  
# Uncomment these lines (and modify their values above), if you want to include them
#netCDFfile.createVariable('beta','f',('time','y0','x0'))[:] = beta

netCDFfile.close()


# =====================================
# Run CISM
print 'Running CISM'
print '============\n'
if options.parallel == None:
   # Perform a serial run
   os.system(options.executable + ' ' + options.configfile)
else:
   # Perform a parallel run
   if options.parallel <= 0:
      sys.exit( 'Error: Number of processors specified for parallel run is <=0.' )
   else:
      # These calls to os.system will return the exit status: 0 for success (the command exists), some other integer for failure
      if os.system('which openmpirun > /dev/null') == 0:
         mpiexec = 'openmpirun -np ' + str(options.parallel)
      elif os.system('which mpirun > /dev/null') == 0:
         mpiexec = 'mpirun -np ' + str(options.parallel)
      elif os.system('which aprun > /dev/null') == 0:
         mpiexec = 'aprun -n ' + str(options.parallel)
      elif os.system('which mpirun.lsf > /dev/null') == 0:
         # mpirun.lsf does NOT need the number of processors (options.parallel)
         mpiexec = 'mpirun.lsf'
      else:
         sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run the model manually with something like: mpirun -np 4 ./cism_driver mismip_init.config')
      runstring = mpiexec + ' ' + options.executable + ' ' + options.configfile
      print 'Executing parallel run with:  ' + runstring + '\n\n'
      os.system(runstring)  # Here is where the parallel run is actually executed!


