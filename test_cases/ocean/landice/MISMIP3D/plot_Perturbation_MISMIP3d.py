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

# dictionary describing what's what for each run
results = OrderedDict()
if options.sfilename:
   results['stnd'] = {'fname':options.sfilename, 'timelev':0, 'color':'k'}
   results['S'] = {'fname':options.sfilename, 'timelev':-1, 'color':'r'}
if options.rfilename:
   results['R'] = {'fname':options.rfilename, 'timelev':-1, 'color':'c'}

# Make plot axis
fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111)

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

plt.ylim( (0.0, 50.0) )
plt.xlabel('X-position (km)')
plt.ylabel('Y-position (km)')


plt.draw()
plt.show()

