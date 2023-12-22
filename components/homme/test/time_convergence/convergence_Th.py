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

varstr = 'Th'

refnhsol = 'refnhsol.nc'
frefnh = scipy.io.netcdf_file(refnhsol)
areanh = frefnh.variables['area'][:]
varnh  = frefnh.variables[varstr][:]
vwtsnh = frefnh.variables['ilev'][:]
htnh   = vwtsnh
vwtsnh = vwtsnh[1:]-vwtsnh[0:-1]

refhsol = 'refhsol.nc'
frefh   = scipy.io.netcdf_file(refhsol)
areah   = frefh.variables['area'][:]
varh    = frefh.variables[varstr][:]
vwtsh   = frefh.variables['ilev'][:]
htnh    = vwtsh
vwtsh   = vwtsh[1:]-vwtsh[0:-1]

E9nh  = np.zeros(4)
E5nh  = np.zeros(4)
E10nh = np.zeros(4)
E9h   = np.zeros(4)
E5h   = np.zeros(4)
E10h  = np.zeros(4)
for j in range(1,5):
  if (j == 1):
    fileteststr9    = 'tstep9_dt0p225.nc'
    fileteststr5    = 'tstep5_dt0p225.nc'
    fileteststr10   = 'tstep10_dt0p225.nc'
    fileteststr9h   = 'tstep9h_dt0p225.nc'
    fileteststr5h   = 'tstep5h_dt0p225.nc'
    fileteststr10h  = 'tstep10h_dt0p225.nc'
  elif (j == 2):
    fileteststr9    = 'tstep9_dt0p45.nc'
    fileteststr5    = 'tstep5_dt0p45.nc'
    fileteststr10   = 'tstep10_dt0p45.nc'
    fileteststr9h   = 'tstep9h_dt0p45.nc'
    fileteststr5h   = 'tstep5h_dt0p45.nc'
    fileteststr10h  = 'tstep10h_dt0p45.nc'
  elif (j == 3):
    fileteststr9    = 'tstep9_dt0p9.nc'
    fileteststr5    = 'tstep5_dt0p9.nc'
    fileteststr10   = 'tstep10_dt0p9.nc'
    fileteststr9h   = 'tstep9h_dt0p9.nc'
    fileteststr5h   = 'tstep5h_dt0p9.nc'
    fileteststr10h  = 'tstep10h_dt0p9.nc'
  elif (j == 4):
    fileteststr9    = 'tstep9_dt1p8.nc'
    fileteststr5    = 'tstep5_dt1p8.nc'
    fileteststr10   = 'tstep10_dt1p8.nc'
    fileteststr9h   = 'tstep9h_dt1p8.nc'
    fileteststr5h   = 'tstep5h_dt1p8.nc'
    fileteststr10h  = 'tstep10h_dt1p8.nc'

  ftest9 = scipy.io.netcdf_file(fileteststr9)
  vart9  = ftest9.variables[varstr][:]
  L2err9 =  L2diff(vart9[-1,:,:],varnh[-1,:,:],areanh,vwtsnh)
  E9nh[j-1] = L2err9
  ftest9.close()
  ftest5 = scipy.io.netcdf_file(fileteststr5)
  vart5  = ftest5.variables[varstr][:]
  L2err5 =  L2diff(vart5[-1,:,:],varnh[-1,:,:],areanh,vwtsnh)
  E5nh[j-1] = L2err5
  ftest5.close()
  ftest10 = scipy.io.netcdf_file(fileteststr10)
  vart10  = ftest10.variables[varstr][:]
  L2err10 = L2diff(vart10[-1,:,:],varnh[-1,:,:],areanh,vwtsnh)
  E10nh[j-1] = L2err10
  ftest10.close()

  ftest9h = scipy.io.netcdf_file(fileteststr9h)
  vart9h  = ftest9h.variables[varstr][:]
  L2err9h =  L2diff(vart9h[-1,:,:],varh[-1,:,:],areah,vwtsh)
  E9h[j-1] = L2err9h
  ftest9h.close()
  ftest5h = scipy.io.netcdf_file(fileteststr5h)
  vart5h  = ftest5h.variables[varstr][:]
  L2err5h =  L2diff(vart5h[-1,:,:],varh[-1,:,:],areah,vwtsh)
  E5h[j-1] = L2err5h
  ftest5h.close()
  ftest10h = scipy.io.netcdf_file(fileteststr10h)
  vart10h  = ftest10h.variables[varstr][:]
  L2err10h = L2diff(vart10h[-1,:,:],varh[-1,:,:],areah,vwtsh)
  E10h[j-1] = L2err10h
  ftest10.close()

plt.title('Max-norm error vs dt for L2-spatial norm of Th, 22.8 s')
dt = [0.45, 0.9, 1.8, 3.6]
Th9nh = plt.loglog(dt,E9nh,label='type9 Max-norm error')
Th5nh = plt.loglog(dt,E5nh,label='type5 Max-norm error')
Th10nh = plt.loglog(dt,E10nh,label='type10 Max-norm error')
slp2   = [ .2e-7/(4**3), .2e-7/(4**2), .2e-7/4, .2e-7] 
slp3   = [ .5e-9/(8**3), .5e-9/(8**2),.5e-9/8 , .5e-9]
slope2 = plt.loglog(dt,slp2,label='Slope 2')
slope3 = plt.loglog(dt,slp3,label='Slope 3')
plt.xticks([ 0.45,0.9,1.8,3.6], ['0.45','0.9','1.8','3.6'])
plt.legend()
plt.show()
savestr = varstr+'_nherr.jpg'
plt.savefig(savestr)      
plt.close()

plt.title('Max-norm error vs dt for L2-spatial norm of Th (hydrostatic mode), 22.8 s')
dt = [0.45, 0.9, 1.8, 3.6]
Th9h = plt.loglog(dt,E9h,label='type9h Max-norm error',marker='o')
Th5h = plt.loglog(dt,E5h,label='type5h Max-norm error')
Th10h = plt.loglog(dt,E10h,label='type10h Max-norm error')
slp2   = [  1e-7/(4**3), 1e-7/(4**2), 1e-7/4, 1e-7]
slp3   = [  1e-10/(8**3), 1e-10/(8**2), 1e-10/8 , 1e-10]
slope2 = plt.loglog(dt,slp2,label='Slope 2')
slope3 = plt.loglog(dt,slp3,label='Slope 3')
plt.xticks([ 0.45,0.9,1.8,3.6], ['0.45','0.9','1.8','3.6'])
plt.legend()
plt.show()
savestr = varstr+'_herr.jpg'
plt.savefig(savestr)
plt.close()

frefnh.close()
frefh.close()
