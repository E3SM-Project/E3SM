#!/usr/bin/env python
"""
This script performs the first step of initializing the global ocean.  This
includes:
1. combining Natural Earth land coverage north of 60S with Antarctic
   ice coverage or grounded ice coverage from Bedmap2.
2. combining transects defining critical passages.*
3. combining points used to seed a flood fill of the global ocean.
4. create masks from land coverage.
5. add land-locked cells to land coverage mask.
6. create masks from transects.*
7. cull cells based on land coverage but with transects present
8. create flood-fill mask based on seeds
9. cull cells based on flood-fill mask
10. create masks from transects on the final culled mesh*
* skipped if flag --with_critical_passages not present

Optionally, the -p flag provides the path to the geometric_data
directory from the geometric_features repository, which is assumed
to be the current directory by default. Also, the optional
--with_cavities flag indicates that ice-shelf cavities are present
and the grounded-ice mask from Bedmap2 should be used. The optional
--with_critical_passages flag indicates that critical passages are
to be opened. Otherwise, steps 2, 5 and 9 are skipped
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import os.path
import subprocess
from optparse import OptionParser

from geometric_features import GeometricFeatures


def removeFile(fileName):
    try:
        os.remove(fileName)
    except OSError:
        pass


parser = OptionParser()
parser.add_option("--with_cavities", action="store_true", dest="with_cavities")
parser.add_option("--with_critical_passages", action="store_true",
                  dest="with_critical_passages")
parser.add_option("-p", "--geom_data_path", type="string", dest="path",
                  default="geometric_data",
                  help="Path to the geometric_data from the geometric_features"
                       " repository.")
parser.add_option("--preserve_floodplain", action="store_true",
                  dest="preserve_floodplain", default=False)
options, args = parser.parse_args()

gf = GeometricFeatures(cacheLocation='{}'.format(options.path),
                       remoteBranchOrTag='0.1')

# start with the land coverage from Natural Earth
fcLandCoverage = gf.read(componentName='natural_earth', objectType='region',
                         featureNames=['Land Coverage'])

# remove the region south of 60S so we can replace it based on ice-sheet
# topography
fcSouthMask = gf.read(componentName='ocean', objectType='region',
                      featureNames=['Global Ocean 90S to 60S'])

fcLandCoverage = fcLandCoverage.difference(fcSouthMask)

# Add "land" coverage from either the full ice sheet or just the grounded
# part
if options.with_cavities:
    fcAntarcticLand = gf.read(componentName='bedmap2', objectType='region',
                              featureNames=['AntarcticGroundedIceCoverage'])
else:
    fcAntarcticLand = gf.read(componentName='bedmap2', objectType='region',
                              featureNames=['AntarcticIceCoverage'])

fcLandCoverage.merge(fcAntarcticLand)

# save the feature collection to a geojson file
fcLandCoverage.to_geojson('land_coverage.geojson')

# Create the land mask based on the land coverage, i.e. coastline data.
# Run command is:
# ./MpasMaskCreator.x  base_mesh.nc land_mask.nc -f land_coverage.geojson
args = [
    './MpasMaskCreator.x',
    'base_mesh.nc',
    'land_mask_1_from_land_coverage.nc',
    '-f', 'land_coverage.geojson']
print("running {}".format(' '.join(args)))

subprocess.check_call(args, env=os.environ.copy())

if options.with_critical_passages:
    outMaskFile = 'land_mask_2_from_land_locked_cells.nc'
else:
    outMaskFile = 'land_mask_final.nc'

# Add land-locked cells to land coverage mask.
args = ['./add_land_locked_cells_to_mask.py',
        '-f', 'land_mask_1_from_land_coverage.nc',
        '-o', outMaskFile,
        '-m', 'base_mesh.nc',
        '-l', '43.0',
        '-n', '20']
print("running {}".format(' '.join(args)))

subprocess.check_call(args, env=os.environ.copy())

# create seed points for a flood fill of the ocean
# use all points in the ocean directory, on the assumption that they are, in
# fact, in the ocean
fcSeed = gf.read(componentName='ocean', objectType='point',
                 tags=['seed_point'])
fcSeed.to_geojson('seed_points.geojson')

if options.with_critical_passages:
    # merge transects for critical passages into critical_passages.geojson
    fcCritPassages = gf.read(componentName='ocean', objectType='transect',
                             tags=['Critical_Passage'])
    fcCritPassages.to_geojson('critical_passages.geojson')

    # create masks from the transects
    # Run command is:
    # ./MpasMaskCreator.x  base_mesh.nc critical_passages_mask.nc
    # -f critical_passages.geojson
    args = ['./MpasMaskCreator.x', 'base_mesh.nc', 'critical_passages_mask.nc',
            '-f', 'critical_passages.geojson']
    print("running {}".format(' '.join(args)))

    subprocess.check_call(args, env=os.environ.copy())

    # Alter critical passages to be at least two cells wide, to avoid sea ice
    # blockage.
    args = ['./widen_transect_edge_masks.py',
            '-f', 'critical_passages_mask.nc',
            '-m', 'base_mesh.nc',
            '-l', '43.0']
    print("running {}".format(' '.join(args)))

    subprocess.check_call(args, env=os.environ.copy())

    # merge transects for critical land blockages into
    # critical_land_blockages.geojson
    fcCritBlockages = gf.read(componentName='ocean', objectType='transect',
                              tags=['Critical_Land_Blockage'])
    fcCritBlockages.to_geojson('critical_land_blockages.geojson')

    # create masks from the transects for critical land blockages
    args = ['./MpasMaskCreator.x', 'base_mesh.nc',
            'critical_land_blockages_mask.nc',
            '-f', 'critical_land_blockages.geojson']
    print("running {}".format(' '.join(args)))

    subprocess.check_call(args, env=os.environ.copy())

    # add critical land blockages to land_mask_final.nc
    args = ['./add_critical_land_blockages_to_mask.py',
            '-f', 'land_mask_2_from_land_locked_cells.nc',
            '-o', 'land_mask_final.nc',
            '-b', 'critical_land_blockages_mask.nc']
    print("running {}".format(' '.join(args)))

    subprocess.check_call(args, env=os.environ.copy())

    # Cull the mesh based on the land mask while keeping critical passages open
    # Run command is:
    # ./MpasCellCuller.x  base_mesh.nc culled_mesh_preliminary.nc -m land_mask_final.nc
    # -p critical_passages_mask.nc
    if options.preserve_floodplain:
      args = ['./MpasCellCuller.x', 'base_mesh.nc', 'culled_mesh_preliminary.nc',
              '-m', 'land_mask_final.nc', '-p', 'critical_passages_mask.nc', '-p', 'base_mesh.nc']
    else:
      args = ['./MpasCellCuller.x', 'base_mesh.nc', 'culled_mesh_preliminary.nc',
              '-m', 'land_mask_final.nc', '-p', 'critical_passages_mask.nc']
    print("running {}".format(' '.join(args)))

    subprocess.check_call(args, env=os.environ.copy())
else:

    # cull the mesh based on the land mask
    # Run command is:
    # ./MpasCellCuller.x  base_mesh.nc culled_mesh_preliminary.nc -m land_mask_final.nc
    if options.preserve_floodplain:
      args = ['./MpasCellCuller.x', 'base_mesh.nc', 'culled_mesh_preliminary.nc',
              '-m', 'land_mask_final.nc', '-p', 'base_mesh.nc']
    else:
      args = ['./MpasCellCuller.x', 'base_mesh.nc', 'culled_mesh_preliminary.nc',
              '-m', 'land_mask_final.nc']
    print("running {}".format(' '.join(args)))

    subprocess.check_call(args, env=os.environ.copy())

# create a mask for the flood fill seed points
# Run command is:
# ./MpasMaskCreator.x  culled_mesh_preliminary.nc seed_mask.nc -s seed_points.geojson
args = ['./MpasMaskCreator.x', 'culled_mesh_preliminary.nc', 'seed_mask.nc',
        '-s', 'seed_points.geojson']
print("running {}".format(' '.join(args)))

subprocess.check_call(args, env=os.environ.copy())


# cull the mesh a second time using a flood fill from the seed points
# Run command is:
# ./MpasCellCuller.x  culled_mesh_preliminary.nc culled_mesh.nc -i seed_mask.nc
args = ['./MpasCellCuller.x', 'culled_mesh_preliminary.nc', 'culled_mesh.nc',
        '-i', 'seed_mask.nc']

print("running {}".format(' '.join(args)))

subprocess.check_call(args, env=os.environ.copy())

if options.with_critical_passages:
    # make a new version of the critical passages mask on the culled mesh
    # Run command is:
    # ./MpasMaskCreator.x  culled_mesh.nc critical_passages_mask_final.nc
    # -f critical_passages.geojson
    args = ['./MpasMaskCreator.x', 'culled_mesh.nc',
            'critical_passages_mask_final.nc',
            '-f', 'critical_passages.geojson']
    print("running {}".format(' '.join(args)))

    subprocess.check_call(args, env=os.environ.copy())
