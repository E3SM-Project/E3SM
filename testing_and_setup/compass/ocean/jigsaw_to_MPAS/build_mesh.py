#!/usr/bin/env python
"""
This script performs the first step of initializing the global ocean.  This
includes:
Step 1. Build cellWidth array as function of latitude and longitude
Step 2. Build mesh using JIGSAW
Step 3. Convert triangles from jigsaw format to netcdf
Step 4. Convert from triangles to MPAS mesh
Step 5. Create vtk file for visualization
"""
import sys
sys.path.append(".")

import os
import subprocess
import scipy.io as sio
import define_base_mesh


def removeFile(fileName):
    try:
        os.remove(fileName)
    except OSError:
        pass


print 'Step 1. Build cellWidth array as function of latitude and longitude'
cellWidth, lon, lat = define_base_mesh.cellWidthVsLatLon()
sio.savemat(
    'cellWidthVsLatLon.mat', {
        'cellWidth': cellWidth, 'lon': lon, 'lat': lat})

print 'Step 2. Build mesh using JIGSAW'
args = ["octave", "--silent", "--eval",
        "jigsaw_driver"]
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

print 'Step 3. Convert triangles from jigsaw format to netcdf'
args = ['./triangle_jigsaw_to_netcdf.py',
        '-s',
        '-m', 'mesh-MESH.msh',
        '-o', 'mesh_triangles.nc']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

print 'Step 4. Convert from triangles to MPAS mesh'
args = ['./MpasMeshConverter.x',
        'mesh_triangles.nc',
        'base_mesh.nc']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

print 'Step 4.1 Inject correct meshDensity variable into base mesh file'
args = ['./inject_meshDensity.py',
        'cellWidthVsLatLon.mat',
        'base_mesh.nc']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

print 'Step 5. Injecting bathymetry'
args = ['./inject_bathymetry.py',
        'base_mesh.nc']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

#print 'Step 6. Create vtk file for visualization'
# args = ['./paraview_vtk_field_extractor.py',
#        '--ignore_time',
#				'-d','maxEdges=0',
#        '-v', 'allOnCells',
#        '-f', 'base_mesh.nc',
#        '-o', 'base_mesh_vtk']
#print "running", ' '.join(args)
#subprocess.check_call(args, env=os.environ.copy())

print 'Step 7. Cull land cells'
args = ['./MpasCellCuller.x',
        'base_mesh.nc',
        'base_mesh_culled.nc']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

print 'Step 8. Injecting bathymetry'
args = ['./inject_bathymetry.py',
        'base_mesh_culled.nc']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

print "***********************************************"
print "**    The global mesh file is base_mesh.nc   **"
print "***********************************************"
