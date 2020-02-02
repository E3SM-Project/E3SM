#!/usr/bin/env python
'''
Script to map cell indices from MPASO noLI mesh to those of the wLI mesh in the runoff mapping file.
Start by building a runoff mapping file that has all the mesh description from wLI mapping file
but the actual mapping from the noLI mapping file:
ncks -x -v S,col,row /project/projectdirs/acme/inputdata/cpl/cpl6/map_rx1_to_oEC60to30v3wLI_smoothed.r300e600.170328.nc newfile.nc
ncks -A -v S,col,row /project/projectdirs/acme/inputdata/cpl/cpl6/map_rx1_to_oEC60to30v3_smoothed.r300e600.161222.nc newfile.nc
'''

# import modules # {{{ 
import netCDF4
import numpy as np
import argparse
import shutil
# }}}

# parser # {{{ 
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', dest='input_file',
                    default='map_rx1_to_oEC60to30v3wLI.nc',
                    help='Input file, original runoff mapping file'
                    )
parser.add_argument('-o', '--output_file', dest='output_file',
                    default='map_rx1_to_oEC60to30v3wLI_final.nc',
                    help='Output file, revised runoff mapping file with no runoff below ice shelf cavities'
                    )
parser.add_argument('-l', '--lookup_table_file', dest='lookup_table_file',
                    default='lookup_table.txt',
                    help='lookup table file, only used locally'
                    )
parser.add_argument('-w', '--mesh_with_ISC', dest='mesh_with_ISC',
                    default='culled_mesh.nc',
                    help='mesh file, including ice shelf cavities'
                    )
parser.add_argument('-n', '--mesh_no_ISC', dest='mesh_no_ISC',
                    default='no_ISC_culled_mesh.nc',
                    help='mesh file, but without ice shelf cavities'
                    )


input_file = parser.parse_args().input_file
output_file = parser.parse_args().output_file
lookup_table_file = parser.parse_args().lookup_table_file
shutil.copy2(input_file, output_file)
# }}}

build_table = True
if build_table:
    # noLI mesh
    mesh_no_ISC = netCDF4.Dataset(parser.parse_args().mesh_no_ISC, 'r')
    noLIxCell = mesh_no_ISC.variables['xCell'][:]
    noLIyCell = mesh_no_ISC.variables['yCell'][:]
    noLInCells = len(mesh_no_ISC.dimensions['nCells'])

    # wLI mesh
    mesh_with_ISC = netCDF4.Dataset(parser.parse_args().mesh_with_ISC, 'r')
    wLIxCell = mesh_with_ISC.variables['xCell'][:]
    wLIyCell = mesh_with_ISC.variables['yCell'][:]

    # init lookup table
    lookup = np.zeros((noLInCells,), dtype=np.uint32)

    print("nCells=", noLInCells)
    for i in range(noLInCells):
        # for i in range(30):
        if i % 1000 == 0:
            print("Cell: ", i)
        # find index of wLI mesh that is the same location as each cell in the
        # noLI mesh
        lookup[i] = np.argmin((noLIxCell[i] - wLIxCell[:])
                              ** 2 + (noLIyCell[i] - wLIyCell[:])**2)
    mesh_no_ISC.close()
    mesh_with_ISC.close()
    print( "Lookup table complete.")
    np.savetxt(lookup_table_file, lookup, fmt='%d')
    print("Saved to ", lookup_table_file)
else:
    lookup = np.loadtxt(lookup_table_file, dtype=np.uint32)
    print("Loaded lookup table from:", lookup_table_file)

print("Lookup: first entries:", lookup[0:10])
print("Lookup: last entries:", lookup[-10:])

# now swap in wLI indices into the runoff mapping file
f = netCDF4.Dataset(output_file, "r+")
row = f.variables['row'][:]
rownew = row * 0
for i in range(len(row)):
    rownew[i] = lookup[row[i] - 1] + 1  # 1-based
f.variables['row'][:] = rownew[:]
f.close()
print("Copied over indices.")

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
