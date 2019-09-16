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

perturbs = np.arange(0.0, 2.1, 0.2)
nruns = len(perturbs)

# === set up plots ===
# another plot ===============

fig = plt.figure(2, facecolor='w', figsize=(13,5))
axVAF = fig.add_subplot(1,3,1)
plt.xlabel('ice shelf melt rate (Gt yr$^{-1}$)')
plt.ylabel('VAF change (Gt/yr)')
plt.grid()


axPctDiff = fig.add_subplot(1,3,2)
plt.xlabel('fractional change in melt rate from control')
plt.ylabel('difference between + and - perturbations (%)')
plt.grid()

axAbsDiff = fig.add_subplot(1,3,3)
plt.xlabel('fractional change in melt rate from control')
plt.ylabel('difference between + and - perturbations (Gt yr$^{-1}$)')
plt.grid()

fig.tight_layout()




yr_samples = [1.0, 2.5, 10.0]
for yr in yr_samples:
 print "Analyzing at year {}".format(yr)
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
      yrs = f.variables['daysSinceStart'][:]/365.0
      #ind = np.where(yrs == yr)[0][0]
      ind = np.argmin(np.absolute(yrs-yr))
      ind2 = np.argmin(np.absolute(yrs-(yr-1.0)))  # index to one year before 'yr'
      GAchange[i] = (GA[ind] - GA[0]) / (yrs[ind] - yrs[0])
      GAchangeInst[i] = (GA[2] - GA[1]) / dt[2]
      VAFchange[i] = (VAF[ind] - VAF[0]) / (yrs[ind] - yrs[0])
      VAFchange[i] = (VAF[ind] - VAF[ind2]) / (yrs[ind] - yrs[ind2])
      VAFchangeInst[i] = (VAF[2] - VAF[1]) / dt[2]
      GLfchange[i] = (GLf[ind] - GLf[0]) / (yrs[ind] - yrs[0])
      GLfchangeInst[i] = (GLf[2] - GLf[1]) / dt[2]
      f.close()
   else:
      print "err"
 # plot multiple years on the plots set up above
 ind = np.where(perturbs==1.0)[0][0]
 axVAF.plot(meltMag, VAFchange, '.', label="{} yr".format(yr))
 axVAF.plot(meltMag[ind], VAFchange[ind], 'k.')
 ind06 = np.argmin(np.absolute(perturbs-0.6))
 ind14 = np.argmin(np.absolute(perturbs-1.4))
 axVAF.plot(meltMag[ind], VAFchange[ [ind06, ind14] ].mean(), 'kx')
 rng = np.arange(0,nruns-ind)
 axPctDiff.plot(perturbs[ind+rng]-perturbs[ind], ((VAFchange[ind+rng] - VAFchange[ind]) - (VAFchange[ind] - VAFchange[ind-rng])) / VAFchange[ind] * 100.0, '.')
 axAbsDiff.plot(perturbs[ind+rng]-perturbs[ind], (VAFchange[ind+rng] - VAFchange[ind]) - (VAFchange[ind] - VAFchange[ind-rng]) , '.')

axVAF.legend()

# plot ===============
ind = np.where(perturbs==1.0)[0][0]

fig = plt.figure(1, facecolor='w', figsize=(10,10))
axGA = fig.add_subplot(3,2,1)
plt.plot(meltMag, GAchange, '.')
plt.plot(meltMag[ind], GAchange[ind], 'ko')
plt.xlabel('ice shelf melt rate (Gt yr$^{-1}$)')
plt.ylabel('grounded area change\n (km$^2$/yr)')
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
plt.ylabel('VAF change (Gt/yr)')
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
plt.ylabel('GL flux change\n (Gt yr$^{-1}$/yr)')
plt.grid()

axMelt = fig.add_subplot(3,2,6)
plt.plot(meltMag, GLfchangeInst, '.')
plt.plot(meltMag[ind], GLfchangeInst[ind], 'ko')
plt.xlabel('ice shelf melt rate (Gt yr$^{-1}$)')
plt.ylabel('GL flux change\nin one time step (Gt yr$^{-1}$)')
plt.grid()

fig.tight_layout()


plt.show()
