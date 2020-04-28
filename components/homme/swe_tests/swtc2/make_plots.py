import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

DS_new = xr.open_dataset('run-new/swtc21.nc')
DS_orig = xr.open_dataset('run-orig/swtc21.nc')

geop_new = DS_new.geop
geop_orig = DS_orig.geop

u_new = DS_new.u
u_orig = DS_orig.u

v_new = DS_new.v
v_orig = DS_orig.v

def plot_3Dvar(newvar, origvar, name, level, i):
    plt.figure(figsize=(10,8))
    plt.contourf(newvar.isel(time=i,lev=level))
    plt.colorbar()
    plt.contour(newvar.isel(time=i,lev=level))
    plt.xlabel('lat')
    plt.ylabel('lon')
    plt.savefig(name+str(i)+'new.png')
    plt.close('all')

    plt.figure(figsize=(10,8))
    plt.contourf(origvar.isel(time=i,lev=level))
    plt.colorbar()
    plt.contour(origvar.isel(time=i,lev=level))
    plt.xlabel('lat')
    plt.ylabel('lon')
    plt.savefig(name+str(i)+'orig.png')
    plt.close('all')

    plt.figure(figsize=(10,8))
    plt.contourf(newvar.isel(time=i,lev=level) - origvar.isel(time=i,lev=level))
    plt.colorbar()
    plt.contour(newvar.isel(time=i,lev=level) - origvar.isel(time=i,lev=level))
    plt.xlabel('lat')
    plt.ylabel('lon')
    plt.savefig(name+str(i)+'diff.png')
    plt.close('all')

def plot_2Dvar(newvar, origvar, name, i):
    plt.figure(figsize=(10,8))
    plt.contourf(newvar.isel(time=i))
    plt.colorbar()
    plt.contour(newvar.isel(time=i))
    plt.xlabel('lat')
    plt.ylabel('lon')
    plt.savefig(name+str(i)+'new.png')
    plt.close('all')

    plt.figure(figsize=(10,8))
    plt.contourf(origvar.isel(time=i))
    plt.colorbar()
    plt.contour(origvar.isel(time=i))
    plt.xlabel('lat')
    plt.ylabel('lon')
    plt.savefig(name+str(i)+'orig.png')
    plt.close('all')

    plt.figure(figsize=(10,8))
    plt.contourf(newvar.isel(time=i) - origvar.isel(time=i))
    plt.colorbar()
    plt.contour(newvar.isel(time=i) - origvar.isel(time=i))
    plt.xlabel('lat')
    plt.ylabel('lon')
    plt.savefig(name+str(i)+'diff.png')
    plt.close('all')
#STATS PLOTTING

plot_3Dvar(geop_new, geop_orig, 'geop', 0, 12)
plot_3Dvar(u_new, u_orig, 'u', 0, 12)
plot_3Dvar(v_new, v_orig, 'v', 0, 12)
