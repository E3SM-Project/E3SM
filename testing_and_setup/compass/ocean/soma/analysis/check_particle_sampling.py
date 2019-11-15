#!/usr/bin/env python
"""
    Verifies that temperature and salinity profiles for SOMA case from particle
    interpolations are consistent with initialized profiles.

    Phillip J. Wolfram
    04/22/2019
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

ds = xr.open_mfdataset('particleoutput.nc')

# obtain particle data
zLevelParticle = ds.zLevelParticle.values
particleTemperature = ds.particleTemperature.values
particleSalinity = ds.particleSalinity.values

# compute analytical formulations

# hard coded factors from config files
gamma = 1.0/1250.0 # psu /m
s0 = 34 # psu
alpha=0.25 # kg / m^3
beta=0.05
drho = 5 # kg/m^3
zt = 300 # m
t0 = 15 # deg C
hshelf = 100 # m
h0 = 2400 # m
z = np.linspace(0, zLevelParticle.min(), 30)

# compute profiles
salinity = s0 - gamma*z
temperature = ((1-beta)*drho*np.tanh(z/zt) + beta*drho*(z/(h0+hshelf)))/alpha + t0

# temperature comparison
plt.figure()
plt.plot(particleTemperature.ravel(), zLevelParticle.ravel(), '.',
         label='Particle value')
plt.plot(temperature, z, '-', label='Initialized scalar value')
plt.ylabel('Depth (m)')
plt.xlabel('Temperature $^\circ$C')
plt.title('Particle temperature comparison')
plt.savefig('particle_temperature_comparison.png')


# salinity comparison
plt.figure()
plt.plot(particleSalinity.ravel(), zLevelParticle.ravel(), '.',
         label='Particle value')
plt.plot(salinity, z, '-', label='Initialized scalar value')
plt.title('Particle salinity comparison')
plt.ylabel('Depth (m)')
plt.xlabel('Salinity psu')
plt.savefig('particle_salinity_comparison.png')

