#!/usr/bin/env python
'''
Script to plot results of MISMIP3D perturbation experiments as in Pattyn et al. 2014, Fig. 5 (and possibly other figures)

'''

import numpy as np
import netCDF4
from optparse import OptionParser
import matplotlib.pyplot as plt
from collections import OrderedDict

GLbit = 256


from optparse import OptionParser
parser = OptionParser()
parser.add_option("-s", "--sliding", dest="sfilename", type='string', help="file containing output from PXXS MISMIP3D Perturbation experiment.  This should be a full width domain.", metavar="FILE")
parser.add_option("-r", "--reverse", dest="rfilename", type='string', help="file containing output from PXXR MISMIP3D Perturbation experiment.  This should be a full width domain.", metavar="FILE")
#parser.add_option("-p", "--perturb", dest="perturb", type='float', help="perturbation amount, presumably 75 or 10", metavar="P")
options, args = parser.parse_args()

# Make plot axis
fig = plt.figure(1, facecolor='w', figsize=(5.0, 7.0))

# =====================
# First plot: GL in plan view
print "=== Beginning plan view GL plot ==="
ax1 = fig.add_subplot(211)

# dictionary describing what's what for each run
results = OrderedDict()
if options.sfilename:
   results['stnd'] = {'fname':options.sfilename, 'timelev':0, 'color':'k'}
   results['S'] = {'fname':options.sfilename, 'timelev':-1, 'color':'r'}
if options.rfilename:
   results['R'] = {'fname':options.rfilename, 'timelev':-1, 'color':'c'}

# Loop from the 2 or 3 runs defined above and plot the GL
for run in results:
  t = results[run]['timelev']
  color = results[run]['color']
  f = netCDF4.Dataset(results[run]['fname'], 'r')
  xEdge = f.variables['xEdge'][:]
  yEdge = f.variables['yEdge'][:]
  xVertex = f.variables['xVertex'][:]
  yVertex = f.variables['yVertex'][:]
  vOnE = f.variables['verticesOnEdge'][:]
  edgeMask = f.variables['edgeMask'][t,:]
  xtime = f.variables['xtime'][:]

  print "For run ", run, " using file ", results[run]['fname'], " and time=", "".join(xtime[t,:])

  GLindEast = np.nonzero( np.logical_and(
                 (edgeMask[:] & GLbit) / GLbit == 1,
                 xEdge > 0.0 ) )[0]
  #plt.plot(xEdge[GLindEast] / 1000.0, yEdge[GLindEast] / 1000.0, '*', color=color)  # to just plot the edge locations
  for i in range(len(GLindEast)):
     vind = vOnE[GLindEast[i], :] - 1  # switch to 0-based
     plt.plot(xVertex[vind] / 1000.0, yVertex[vind] / 1000.0, '-', color=color)

  f.close()

plt.xlabel('X-position (km)')
plt.ylabel('Y-position (km)')

plt.ylim( (0.0, 51.0) )
x0,x1 = ax1.get_xlim()
plt.xlim( (x0, x0+60.0) )

x0,x1 = ax1.get_xlim()
y0,y1 = ax1.get_ylim()
ax1.set_aspect((x1-x0)/(y1-y0))

# =====================
# Second plot: GL position time series at two points
print "=== Beginning GL time series plot ==="
ax2 = fig.add_subplot(212)

# Find the y-position at the wall - assume it is the same for S and R files
f = netCDF4.Dataset(options.sfilename, 'r')
yEdge = f.variables['yEdge'][:]
unique_yEdge = np.array(sorted(list(set(yEdge[:]))))
wallY = unique_yEdge[-2] # get the second to last location.  Last one will be 'outside' the domain
print "For wall yEdge using value:", wallY
f.close()

# Build dictionary
# dictionary describing what's what for each run
results = OrderedDict()
if options.sfilename:
   results['Saxis'] = {'fname':options.sfilename, 'y':0.0, 'color':'r', 'reverse':False}
   results['Swall'] = {'fname':options.sfilename, 'y':wallY, 'color':'r', 'reverse':False}
if options.rfilename:
   results['Raxis'] = {'fname':options.rfilename, 'y':0.0, 'color':'c', 'reverse':True}
   results['Rwall'] = {'fname':options.rfilename, 'y':wallY, 'color':'c', 'reverse':True}

# Loop from the 2 or 3 runs defined above and plot the GL
for run in results:
  color = results[run]['color']
  yValue = results[run]['y']
  f = netCDF4.Dataset(results[run]['fname'], 'r')
  xEdge = f.variables['xEdge'][:]
  yEdge = f.variables['yEdge'][:]
  xtime = f.variables['xtime'][:]
  year = f.variables['daysSinceStart'][:] / 365.0
  if results[run]['reverse']:
     year = 100.0 - year
  nt = len(f.dimensions['Time'])

  print "For run ", run, " using file ", results[run]['fname'], " and yEdge value ", yValue

  GLx = np.zeros( (nt,) )
  for t in range(nt):
     edgeMask = f.variables['edgeMask'][t,:]
     GLindEast = np.nonzero( np.logical_and( np.logical_and(
                 (edgeMask[:] & GLbit) / GLbit == 1,
                 xEdge > 0.0 ),
                 yEdge == yValue) )[0]  # should only be 1 index
     GLx[t] = xEdge[GLindEast]

  plt.plot(year, GLx / 1000.0, '.', color=color)

  # save year, GL position
  txtfname = "MPASLI_MISMIP3D_P75_GLposition_" + run + ".txt"
  np.savetxt(txtfname, np.column_stack((year, GLx)), fmt='%10.2f', delimiter=',', newline='\n', header='', footer='', comments='# ')

  f.close()

plt.xlabel('Time (yr)')
plt.ylabel('X-position (km)')


plt.xlim( (0.0, 101.0) )
y0,y1 = ax2.get_ylim()
plt.ylim( (y0, y0+40.0) )

x0,x1 = ax2.get_xlim()
y0,y1 = ax2.get_ylim()
ax2.set_aspect((x1-x0)/(y1-y0))



plt.draw()
plt.show()

