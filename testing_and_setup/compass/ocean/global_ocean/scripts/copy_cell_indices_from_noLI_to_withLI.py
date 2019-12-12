#!/usr/bin/env python
'''
Script to map cell indices from MPASO noLI mesh to those of the wLI mesh in the runoff mapping file.
Start by building a runoff mapping file that has all the mesh description from wLI mapping file
but the actual mapping from the noLI mapping file:
ncks -x -v S,col,row /project/projectdirs/acme/inputdata/cpl/cpl6/map_rx1_to_oEC60to30v3wLI_smoothed.r300e600.170328.nc newfile.nc
ncks -A -v S,col,row /project/projectdirs/acme/inputdata/cpl/cpl6/map_rx1_to_oEC60to30v3_smoothed.r300e600.161222.nc newfile.nc
'''

import netCDF4
import numpy as np
file_to_modify = "map_rx1_to_oEC60to30v3wLI_smoothed.r300e600.170328.mapping_from_noLI_mesh.nc"
ltable_file = "noLI_to_wLI_lookup_table.txt"

build_table = False
if build_table:
   #noLI mesh
   fnoLI=netCDF4.Dataset('/project/projectdirs/acme/inputdata/ocn/mpas-o/oEC60to30v3/oEC60to30v3_60layer.170905.nc','r')
   noLIxCell=fnoLI.variables['xCell'][:]
   noLIyCell=fnoLI.variables['yCell'][:]
   noLInCells = len(fnoLI.dimensions['nCells'])

   #wLI mesh
   fwLI =netCDF4.Dataset('/project/projectdirs/acme/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc','r')
   wLIxCell=fwLI.variables['xCell'][:]
   wLIyCell=fwLI.variables['yCell'][:]

   # init lookup table
   lookup=np.zeros((noLInCells,), dtype=np.uint32)

   print "nCells=", noLInCells
   for i in range(noLInCells):
   #for i in range(30):
           if i%1000==0:
              print "Cell: ", i
           # find index of wLI mesh that is the same location as each cell in the noLI mesh
           lookup[i] = np.argmin( (noLIxCell[i]-wLIxCell[:])**2 + (noLIyCell[i]-wLIyCell[:])**2 )
   fnoLI.close()
   fwLI.close()
   print "Lookup table complete."
   np.savetxt(ltable_file, lookup, fmt='%d')
   print "Saved to ", ltable_file
else:
   lookup=np.loadtxt(ltable_file, dtype=np.uint32)
   print "Loaded lookup table from:", ltable_file

print "Lookup: first entries:", lookup[0:10]
print "Lookup: last entries:", lookup[-10:]

# now swap in wLI indices into the runoff mapping file
f=netCDF4.Dataset(file_to_modify, "r+")
row=f.variables['row'][:]
rownew=row*0
for i in range(len(row)):
    rownew[i] = lookup[row[i] - 1] + 1 # 1-based
f.variables['row'][:] = rownew[:]
f.close()
print "Copied over indices."