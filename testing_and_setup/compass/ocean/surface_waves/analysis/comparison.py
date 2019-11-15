#!/usr/bin/env python
"""

Tidal channel comparison betewen MPAS-O and analytical forcing result.

Phillip J. Wolfram
04/12/2019

"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# render statically by default
plt.switch_backend('agg')

# analytical case
x = np.linspace(0,24,100)
y = np.sin(x*2*np.pi/24)
plt.plot(x,y, lw=3, color='black', label='analytical')

# data from MPAS-O on boundary
ds = xr.open_mfdataset('output.nc')
ds.ssh.where(ds.tidalInputMask).mean('nCells').plot(marker='o', label='MPAS-O')

plt.legend()
plt.ylabel('Tidal amplitude (m)')
plt.xlabel('Time (hrs)')

plt.savefig('tidalcomparison.png')
