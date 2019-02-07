#!/usr/bin/env python
'''
Script to plot melt perturbations for Thwaites experiment
'''

import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

perturbs = np.arange(0.0, 2.1, 0.1)
nruns = len(perturbs)

# init some vectors
meltMag = np.zeros((nruns,))
GAchange = np.zeros((nruns,))
GAchangeInst = np.zeros((nruns,))
VAFchange = np.zeros((nruns,))
VAFchangeInst = np.zeros((nruns,))
GLfchange = np.zeros((nruns,))
GLfchangeInst = np.zeros((nruns,))

for i in range(nruns):
   p = perturbs[i]
   dname = "meltfactor{}".format(p)
   print dname
   fname = dname + "/globalStats.nc"
   if os.path.isfile(fname):
      f = netCDF4.Dataset(fname, 'r')
      #meltMag[i] = f.variables['totalFloatingBasalMassBal'][1] /1.0e9 * -1.0
      meltMag[i] = f.variables['totalFloatingBasalMassBal'][1:].mean() /1.0e9 /1000.0 * -1.0
      GA = f.variables['groundedIceArea'][:] / 1000.0**2
      VAF = f.variables['volumeAboveFloatation'][:] *910.0 /1.0e9 /1000.0
      GLf = f.variables['groundingLineFlux'][:] /1.0e9 /1000.0
      dt = f.variables['deltat'][:] / (3600.0*24.0*365.0)
      GAchange[i] = GA[-1] - GA[0]
      GAchangeInst[i] = (GA[2] - GA[1]) / dt[2]
      VAFchange[i] = VAF[-1] - VAF[0]
      VAFchangeInst[i] = (VAF[2] - VAF[1]) / dt[2]
      GLfchange[i] = GLf[-1] - GLf[0]
      GLfchangeInst[i] = (GLf[2] - GLf[1]) / dt[2]
      f.close()
   else:
      print "err"

print meltMag
print GAchange

# plot ===============
ind = np.where(perturbs==1.0)[0][0]

fig = plt.figure(1, facecolor='w', figsize=(10,10))
axGA = fig.add_subplot(3,2,1)
plt.plot(meltMag, GAchange, '.')
plt.plot(meltMag[ind], GAchange[ind], 'ko')
plt.xlabel('ice shelf melt rate (Gt yr$^{-1}$)')
plt.ylabel('grounded area change\nin one year (km$^2$)')
plt.grid()

axGAinst = fig.add_subplot(3,2,2)
plt.plot(meltMag, GAchangeInst, '.')
plt.plot(meltMag[ind], GAchangeInst[ind], 'ko')
plt.xlabel('ice shelf melt rate (Gt yr$^{-1}$)')
plt.ylabel('grounded area change\nin one time step (km$^2$)')
plt.grid()

axMelt = fig.add_subplot(3,2,3)
plt.plot(meltMag, VAFchange, '.')
plt.plot(meltMag[ind], VAFchange[ind], 'ko')
plt.xlabel('ice shelf melt rate (Gt yr$^{-1}$)')
plt.ylabel('VAF change\nin one year (Gt)')
plt.grid()

axMelt = fig.add_subplot(3,2,4)
plt.plot(meltMag, VAFchangeInst, '.')
plt.plot(meltMag[ind], VAFchangeInst[ind], 'ko')
plt.xlabel('ice shelf melt rate (Gt yr$^{-1}$)')
plt.ylabel('VAF change\nin one time step (Gt)')
plt.grid()

axMelt = fig.add_subplot(3,2,5)
plt.plot(meltMag, GLfchange, '.')
plt.plot(meltMag[ind], GLfchange[ind], 'ko')
plt.xlabel('ice shelf melt rate (Gt yr$^{-1}$)')
plt.ylabel('GL flux change\nin one year (Gt yr$^{-1}$)')
plt.grid()

axMelt = fig.add_subplot(3,2,6)
plt.plot(meltMag, GLfchangeInst, '.')
plt.plot(meltMag[ind], GLfchangeInst[ind], 'ko')
plt.xlabel('ice shelf melt rate (Gt yr$^{-1}$)')
plt.ylabel('GL flux change\nin one time step (Gt yr$^{-1}$)')
plt.grid()

fig.tight_layout()


# another plot ===============

fig = plt.figure(2, facecolor='w', figsize=(7,4))
ax = fig.add_subplot(1,2,1)
rng = np.arange(0,nruns-ind)
plt.plot(perturbs[ind+rng]-perturbs[ind], ((VAFchange[ind+rng] - VAFchange[ind]) - (VAFchange[ind] - VAFchange[ind-rng])) / VAFchange[ind] * 100.0, '.')

#print ind, rng
#print ind+rng
#print ind-rng
#print perturbs[ind+rng]-perturbs[ind]
#print perturbs
#print VAFchange

plt.xlabel('fractional change in melt rate from control')
plt.ylabel('difference between + and - perturbations (%)')
plt.grid()

ax = fig.add_subplot(1,2,2)
plt.plot(perturbs[ind+rng]-perturbs[ind], (VAFchange[ind+rng] - VAFchange[ind]) - (VAFchange[ind] - VAFchange[ind-rng]) , '.')
plt.xlabel('fractional change in melt rate from control')
plt.ylabel('difference between + and - perturbations (Gt yr$^{-1}$)')
plt.grid()

fig.tight_layout()

plt.show()
