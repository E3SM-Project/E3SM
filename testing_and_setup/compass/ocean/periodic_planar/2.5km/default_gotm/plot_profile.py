#!/usr/bin/env python

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# render statically by default
plt.switch_backend('agg')

if __name__ == "__main__":

    # constants
    kappa = 0.4
    z0b = 1.5e-3
    gssh = 1e-5
    g = 9.81
    h = 15
    # load output
    ds = xr.open_dataset('output.nc')
    # velocity
    u = ds.velocityZonal.isel(Time=-1, nCells=0).values
    # viscosity
    nu = ds.vertViscTopOfCell.isel(Time=-1, nCells=0).values
    # depth
    bottom_depth = ds.refBottomDepth.values
    z = np.zeros_like(bottom_depth)
    z[0] = -0.5*bottom_depth[0]
    z[1:] = -0.5*(bottom_depth[0:-1]+bottom_depth[1:])
    zi = np.zeros(bottom_depth.size+1)
    zi[0] = 0.0
    zi[1:] = -bottom_depth[0:]
    # analytical solution
    ustarb = np.sqrt(g*h*gssh)
    u_a = ustarb/kappa*np.log((z0b+z+h)/z0b)
    nu_a = -ustarb/h*kappa*(z0b+z+h)*z
    # infered drag coefficient
    cd = ustarb**2/u_a[-1]**2
    print('$C_d$ = {:6.4g}'.format(cd))
    # plot velocity
    plt.figure()
    plt.plot(u_a, z, 'k--', label='Analytical')
    plt.plot(u, z, 'k-', label='GOTM')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Depth (m)')
    plt.legend()
    plt.savefig('velocity_profile.png')
    # plot viscosity
    plt.figure()
    plt.plot(nu_a, z, 'k--', label='Analytical')
    plt.plot(nu, zi, 'k-', label='GOTM')
    plt.xlabel('Viscosity (m$^2$/s)')
    plt.ylabel('Depth (m)')
    plt.legend()
    plt.savefig('viscosity_profile.png')


