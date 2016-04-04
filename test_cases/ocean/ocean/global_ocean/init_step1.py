#!/usr/bin/env python
"""
This script performs the first step of initializing the global ocean.  This
includes:
1. combining Natural Earth land coverage north of 60S with Antarctic
   ice coverage or grounded ice coverage from Bedmap2.  
2. combining transects defining cricial passages.
3. combining points used to seed a flood fill of the global ocean.
4. create masks from land coverage.
5. create masks from transects.
6. cull cells based on land coverage but with transects present
7. create flood-fill mask based on seeds
8. cull cells based on flood-fill mask
9. create masks from transects on the final culled mesh

Optionally, the -p flag provides the path to the geometric_features 
repository, which is assumed to be the current directory by default.
Also, the optional --with_cavities flag indicates that ice-shelf cavities
are present and the grounded-ice mask from Bedmap2 should be used.
"""
import os
import os.path
import subprocess
import shutil
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--with_cavities", action="store_true", dest="with_cavities")
parser.add_option("-p", "--geom_feat_path", type="string", dest="path", default="geometric_features",
                  help="Path to the geometric_features repository.")
options, args = parser.parse_args()

path = options.path

defaultFileName = 'features.geojson'

try:
    os.remove(defaultFileName)
except OSError:
    pass

args = ['%s/difference_features.py'%path, 
        '-f', '%s/natural_earth/region/Land_Coverage/region.geojson'%path, 
        '-m', '%s/ocean/region/Global_Ocean_90S_to_60S/region.geojson'%path]
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

if options.with_cavities:
    antarcticLandCoverage = '%s/bedmap2/region/AntarcticGroundedIceCoverage/region.geojson'%path
else:
    antarcticLandCoverage = '%s/bedmap2/region/AntarcticIceCoverage/region.geojson'%path

args = ['%s/merge_features.py'%path, '-f', antarcticLandCoverage]
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

outName = 'land_coverage.geojson'
shutil.move(defaultFileName,outName)

args = ['%s/merge_features.py'%path, '-d', '%s/ocean/transect'%path, 
        '-t', 'Critical_Passage']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

outName = 'critical_passages.geojson'
shutil.move(defaultFileName,outName)

# use all points in the ocean directory, on the assumption that they are, in
# fact, in the ocean
args = ['%s/merge_features.py'%path, '-d', '%s/ocean/point'%path]
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

outName = 'seed_points.geojson'
shutil.move(defaultFileName,outName)

# Run command is:
# ./MpasMaskCreator.x  base_mesh.nc land_mask.nc -f land_coverage.geojson
subprocess.check_call(['./MpasMaskCreator.x', 'base_mesh.nc', 'land_mask.nc', 
                       '-f', 'land_coverage.geojson'], env=os.environ.copy())


# Run command is:
# ./MpasMaskCreator.x  base_mesh.nc critical_passages_mask.nc -f critical_passages.geojson
args = ['./MpasMaskCreator.x', 'base_mesh.nc', 'critical_passages_mask.nc',
        '-f', 'critical_passages.geojson']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())


# Run command is:
# ./MpasCellCuller.x  base_mesh.nc culled_mesh.nc -m land_mask.nc -p critical_passages_mask.nc
args = ['./MpasCellCuller.x', 'base_mesh.nc', 'culled_mesh.nc', 
        '-m', 'land_mask.nc', '-p', 'critical_passages_mask.nc']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())


# Run command is:
# ./MpasMaskCreator.x  culled_mesh.nc seed_mask.nc -s seed_points.geojson
args = ['./MpasMaskCreator.x', 'culled_mesh.nc', 'seed_mask.nc',
        '-s', 'seed_points.geojson']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())


# Run command is:
# ./MpasCellCuller.x  culled_mesh.nc culled_mesh_final.nc -i seed_mask.nc
args = ['./MpasCellCuller.x', 'culled_mesh.nc', 'culled_mesh_final.nc', 
        '-i', 'seed_mask.nc']

print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

# Run command is:
# ./MpasMaskCreator.x  culled_mesh_final.nc critical_passages_mask_final.nc -f critical_passages.geojson
args = ['./MpasMaskCreator.x', 'culled_mesh_final.nc', 
        'critical_passages_mask_final.nc', 
        '-f', 'critical_passages.geojson']
print "running", ' '.join(args)
subprocess.check_call(args, env=os.environ.copy())

