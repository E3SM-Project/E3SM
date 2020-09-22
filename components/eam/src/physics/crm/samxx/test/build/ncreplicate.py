import netCDF4
import sys
import numpy as np
import os

#######################################################################################
#######################################################################################
##
## Simple replication of unlimited dimension N times like ncrcat with the same file
## N times. unlimited dim must be slowest varying (first dim in C, Python). All
## Variable are assumed to contain the unlimited dimension
##
#######################################################################################
#######################################################################################

#Complain if there aren't two arguments
if (len(sys.argv) < 3) :
  print("Usage: python ncreplicate.py file.nc  N")
  sys.exit(1)

nrep = int(sys.argv[2])

#Open the two files
nc1 = netCDF4.Dataset(sys.argv[1])
nc2 = netCDF4.Dataset(os.path.splitext(sys.argv[1])[0]+"_"+sys.argv[2]+"x.nc","w")

for d in nc1.dimensions.keys() :
  dim = nc1.dimensions[d]
  if ( dim.isunlimited() ) :
    nc2.createDimension(dim.name,0       )
  else :
    nc2.createDimension(dim.name,dim.size)

#Loop through all variables
for v in nc1.variables.keys() :
  print("Working on variable: "+v)
  var1 = nc1.variables[v]
  var2 = nc2.createVariable(v,var1.datatype,var1.dimensions)
  if ( not nc1.dimensions[var1.dimensions[0]].isunlimited() ) :
    var2[:] = var1[:]
  else :
    ulen = var1.shape[0]
    slist = [var1.shape[i] for i in range(var1.ndim)]
    slist[0] *= nrep
    shape = tuple(slist)
    var2[:] = np.zeros(shape,var1.dtype)
    for i in range(nrep) :
      var2[i*ulen:(i+1)*ulen,...] = var1[:,...]
