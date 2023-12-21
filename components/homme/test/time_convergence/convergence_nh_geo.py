#from scipy.io import netcdf
import glob
import numpy as np
import statistics as st
import xarray as xr
import scipy.io
from scipy.io import netcdf as nc
import matplotlib.pyplot as plt

def L2diff(field1,field2,hwts,zwts):
  diff = np.subtract(field1[:,:],field2[:,:])
  rel  = field2[:,:]
  err  = np.square(diff)/np.square(rel)
  err  = np.sum(err*hwts,1)
  err  = np.sum(err*zwts,0)
  err  = np.sqrt(err/np.sum(hwts))
  return err

outputdir = '/nscratch/asteyer/test_convergence/'

refnhsol = outputdir+'refnhsol.nc'
#refhsol  = outputdir+'refhsol.nc'

fref   = scipy.io.netcdf_file(refnhsol)
time   = fref.variables['time'][:]
area   = fref.variables['area'][:]
varstr = 'geo'
var    = fref.variables[varstr][:]
vwts   = fref.variables['ilev'][:]
ht     = vwts
vwts   = vwts[1:]-vwts[0:-1]

E9nh  = np.zeros(5)
E5nh  = np.zeros(5)
E10nh = np.zeros(5)
E9h   = np.zeros(5)
E5h   = np.zeros(5)
E10h  = np.zeros(5)
for j in range(1,6):
  if (j == 1):
    fileteststr9   = '/nscratch/asteyer/test_convergence/tstep9_dt0p225.nc'
    fileteststr5   = '/nscratch/asteyer/test_convergence/tstep5_dt0p225.nc'
    fileteststr10  = '/nscratch/asteyer/test_convergence/tstep10_dt0p225.nc'
  elif (j == 2):
    fileteststr9   = '/nscratch/asteyer/test_convergence/tstep9_dt0p45.nc'
    fileteststr5   = '/nscratch/asteyer/test_convergence/tstep5_dt0p45.nc'
    fileteststr10  = '/nscratch/asteyer/test_convergence/tstep10_dt0p45.nc'
  elif (j == 3):
    fileteststr9   = '/nscratch/asteyer/test_convergence/tstep9_dt0p9.nc'
    fileteststr5   = '/nscratch/asteyer/test_convergence/tstep5_dt0p9.nc'
    fileteststr10  = '/nscratch/asteyer/test_convergence/tstep10_dt0p9.nc'
  elif (j == 4):
    fileteststr9   = '/nscratch/asteyer/test_convergence/tstep9_dt1p8.nc'
    fileteststr5   = '/nscratch/asteyer/test_convergence/tstep5_dt1p8.nc'
    fileteststr10  = '/nscratch/asteyer/test_convergence/tstep10_dt1p8.nc'
  elif (j == 5):
    fileteststr9   = '/nscratch/asteyer/test_convergence/tstep9_dt3p6.nc'
    fileteststr5   = '/nscratch/asteyer/test_convergence/tstep5_dt3p6.nc'
    fileteststr10  = '/nscratch/asteyer/test_convergence/tstep10_dt3p6.nc'

  ftest9 = scipy.io.netcdf_file(fileteststr9)
  vart9  = ftest9.variables[varstr][:]
  L2err9 =  L2diff(vart9[-1,:,:],var[-1,:,:],area,vwts)
  E9nh[j-1] = L2err9
  ftest9.close()
  ftest5 = scipy.io.netcdf_file(fileteststr5)
  vart5  = ftest5.variables[varstr][:]
  L2err5 =  L2diff(vart5[-1,:,:],var[-1,:,:],area,vwts)
  E5nh[j-1] = L2err5
  ftest5.close()
  ftest10 = scipy.io.netcdf_file(fileteststr10)
  vart10  = ftest10.variables[varstr][:]
  L2err10 = L2diff(vart10[-1,:,:],var[-1,:,:],area,vwts)
  E10nh[j-1] = L2err10
  ftest10.close()
  
plt.title('Final-time (max-norm) error vs dt for L2-spatial norm of Th, 22.8 s')
dt = [0.225, 0.45, 0.9, 1.8, 3.6]
Th9nh = plt.loglog(dt,E9nh,label='type9 Max-norm error')
Th5nh = plt.loglog(dt,E5nh,label='type5 Max-norm error')
Th10nh = plt.loglog(dt,E10nh,label='type10 Max-norm error')
slp2   = [ 1e-7/(4**4), 1e-7/(4**3), 1e-7/(4**2), 1e-7/4, 1e-7] 
slp3   = [ 1e-8/(8**4), 1e-8/(8**3), 1e-8/(8**2), 1e-8/8 , 1e-8]
slope2 = plt.loglog(dt,slp2,label='Slope 2')
slope3 = plt.loglog(dt,slp3,label='Slope 3')
plt.xticks([0.225, 0.45,0.9,1.8,3.6], ['0.225','0.45','0.9','1.8','3.6'])
plt.legend()
plt.show()
savestr = varstr+'_nherr.jpg'

                                                              
plt.savefig(savestr)      
fref.close()
