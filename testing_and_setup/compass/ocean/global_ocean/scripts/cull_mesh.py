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

The optional --with_cavities flag indicates that ice-shelf cavities are present
and the grounded-ice mask from Bedmap2 should be used. The optional
--with_critical_passages flag indicates that critical passages are
to be opened. Otherwise, steps 2, 5 and 9 are skipped
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from optparse import OptionParser
import xarray

from geometric_features import GeometricFeatures, FeatureCollection, \
    read_feature_collection
from mpas_tools.mesh import conversion
from mpas_tools.io import write_netcdf
from mpas_tools.ocean.coastline_alteration import widen_transect_edge_masks, \
    add_critical_land_blockages, add_land_locked_cells_to_mask
from mpas_tools.viz.paraview_extractor import extract_vtk


parser = OptionParser()
parser.add_option("--with_cavities", action="store_true", dest="with_cavities",
                  help="Whether the mesh should include Antarctic ice-shelf"
                       " cavities")
parser.add_option("--with_critical_passages", action="store_true",
                  dest="with_critical_passages",
                  help="Whether the mesh should open the standard critical "
                       "passages and close land blockages from "
                       "geometric_features")
parser.add_option(
    "--custom_critical_passages",
    dest="custom_critical_passages",
    help="A geojson file with critical passages to open.  This "
    "file may be supplied in addition to or instead of "
    "the default passages (--with_critical_passages)")
parser.add_option(
    "--custom_land_blockages",
    dest="custom_land_blockages",
    help="A geojson file with critical land blockages to close. "
    "This file may be supplied in addition to or instead of "
    "the default blockages (--with_critical_passages)")
parser.add_option("--preserve_floodplain", action="store_true",
                  dest="preserve_floodplain", default=False,
                  help="Whether to use the cellSeedMask field in the base "
                       "mesh to preserve a floodplain at elevations above z=0")
options, args = parser.parse_args()

# required for compatibility with MPAS
netcdfFormat = 'NETCDF3_64BIT'

critical_passages = options.with_critical_passages or \
    (options.custom_critical_passages is not None)

land_blockages = options.with_critical_passages or \
    (options.custom_land_blockages is not None)

gf = GeometricFeatures()

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
    fcAntarcticLand = gf.read(componentName='bedmachine', objectType='region',
                              featureNames=['AntarcticGroundedIceCoverage'])
else:
    fcAntarcticLand = gf.read(componentName='bedmachine', objectType='region',
                              featureNames=['AntarcticIceCoverage'])

fcLandCoverage.merge(fcAntarcticLand)

# save the feature collection to a geojson file
fcLandCoverage.to_geojson('land_coverage.geojson')

# Create the land mask based on the land coverage, i.e. coastline data.
dsBaseMesh = xarray.open_dataset('base_mesh.nc')
dsLandMask = conversion.mask(dsBaseMesh, fcMask=fcLandCoverage)

dsLandMask = add_land_locked_cells_to_mask(dsLandMask, dsBaseMesh,
                                           latitude_threshold=43.0,
                                           nSweeps=20)

# create seed points for a flood fill of the ocean
# use all points in the ocean directory, on the assumption that they are, in
# fact, in the ocean
fcSeed = gf.read(componentName='ocean', objectType='point',
                 tags=['seed_point'])

if land_blockages:
    if options.with_critical_passages:
        # merge transects for critical land blockages into
        # critical_land_blockages.geojson
        fcCritBlockages = gf.read(componentName='ocean', objectType='transect',
                                  tags=['Critical_Land_Blockage'])
    else:
        fcCritBlockages = FeatureCollection()

    if options.custom_land_blockages is not None:
        fcCritBlockages.merge(read_feature_collection(
            options.custom_land_blockages))

    # create masks from the transects
    dsCritBlockMask = conversion.mask(dsBaseMesh, fcMask=fcCritBlockages)

    dsLandMask = add_critical_land_blockages(dsLandMask, dsCritBlockMask)

fcCritPassages = FeatureCollection()
dsPreserve = []

if critical_passages:
    if options.with_critical_passages:
        # merge transects for critical passages into critical_passages.geojson
        fcCritPassages.merge(gf.read(componentName='ocean',
                                     objectType='transect',
                                     tags=['Critical_Passage']))

    if options.custom_critical_passages is not None:
        fcCritPassages.merge(read_feature_collection(
            options.custom_critical_passages))

    # create masks from the transects
    dsCritPassMask = conversion.mask(dsBaseMesh, fcMask=fcCritPassages)

    # Alter critical passages to be at least two cells wide, to avoid sea ice
    # blockage.
    dsCritPassMask = widen_transect_edge_masks(dsCritPassMask, dsBaseMesh,
                                               latitude_threshold=43.0)

    dsPreserve.append(dsCritPassMask)

if options.preserve_floodplain:
    dsPreserve.append(dsBaseMesh)


# cull the mesh based on the land mask
dsCulledMesh = conversion.cull(dsBaseMesh, dsMask=dsLandMask,
                               dsPreserve=dsPreserve)

# create a mask for the flood fill seed points
dsSeedMask = conversion.mask(dsCulledMesh, fcSeed=fcSeed)

# cull the mesh a second time using a flood fill from the seed points
dsCulledMesh = conversion.cull(dsCulledMesh, dsInverse=dsSeedMask,
                               graphInfoFileName='culled_graph.info')
write_netcdf(dsCulledMesh, 'culled_mesh.nc', format=netcdfFormat)

if critical_passages:
    # make a new version of the critical passages mask on the culled mesh
    dsCritPassMask = conversion.mask(dsCulledMesh, fcMask=fcCritPassages)
    write_netcdf(dsCritPassMask, 'critical_passages_mask_final.nc',
                 format=netcdfFormat)

if options.with_cavities:
    fcAntarcticIce = gf.read(componentName='bedmachine', objectType='region',
                             featureNames=['AntarcticIceCoverage'])
    fcAntarcticIce.to_geojson('ice_coverage.geojson')
    dsMask = conversion.mask(dsCulledMesh, fcMask=fcAntarcticIce)
    landIceMask = dsMask.regionCellMasks.isel(nRegions=0)
    dsLandIceMask = xarray.Dataset()
    dsLandIceMask['landIceMask'] = landIceMask

    write_netcdf(dsLandIceMask, 'land_ice_mask.nc', format=netcdfFormat)

    dsLandIceCulledMesh = conversion.cull(dsCulledMesh, dsMask=dsMask)
    write_netcdf(dsLandIceCulledMesh, 'no_ISC_culled_mesh.nc',
                 format=netcdfFormat)

extract_vtk(ignore_time=True, dimension_list=['maxEdges='],
            variable_list=['allOnCells'], filename_pattern='culled_mesh.nc',
            out_dir='culled_mesh_vtk')

if options.with_cavities:
    extract_vtk(ignore_time=True, dimension_list=['maxEdges='],
                variable_list=['allOnCells'],
                filename_pattern='no_ISC_culled_mesh.nc',
                out_dir='no_ISC_culled_mesh_vtk')
