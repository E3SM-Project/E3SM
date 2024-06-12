# Generate a new Grid File

[Back to Adding Support for New Grids](../adding-grid-support-step-by-step-guide.md)

In order to generate mapping files between a new atmosphere grid and the surface component grids, we need a file that describes the new grid. [TempestRemap](https://github.com/ClimateGlobalChange/tempestremap) is our preferred tool for grid file generation because it can handle the spectral element grids used by the atmosphere dycore. The initial grid file will be saved in an "exodus" file with a `.g` extension (see [Types of Grid Description Files](../adding-grid-support-grid-types.md) for more info). TempestRemap can be installed via conda.

## Generating a Standard Exodus Grid File

Once TempestRemap is in our environment we can easily generate an exodus file by calling TempestRemap directly:

```shell
GenerateCSMesh --alt --res 4 --file ne4.g
```

## Generating a Regionally Refined Grid File

For a regionally refined mesh (RRM) [SQuadGen](https://github.com/ClimateGlobalChange/squadgen) is used to define the refined area(s). [This tutorial](generate-RRM-grid-file.md) includes details and examples of using SQuadGen to generate RRM grid files.

The naming convention for RRM grid files should follow:

```shell
<refined_area_name>_<base_resolution>x<refinement_level>.g
```

For example, for a RRM with 4x refinement from ne30 to ne120 over CONUS, we should use the convention conus_ne30x4.g (note that existing meshes may have used the old naming convention `<area of refinement>x<refinement level>v<version>.g`, but future meshes should use the new naming convention).

The Exodus file contains only information about the position of the spectral element on the sphere. For SE aware utilities such as TempestRemap, they can use the polynomial order and the reference element map to fill in necessary data such as the locations of the nodal GLL points. For non-SE aware utilities, we need additional meta data, described in the next section.  

## Generating a "pg2" SCRIP Grid File

Starting in E3SMv2 the physics calculations and standard history output use a finite volume (FV) "pg2" grid. Online mapping within the component coupler between the atmosphere and surface components requires FV-to-FV type maps, and generating these maps will require a grid file for pg2 grid. These are easily generated with TempestRemap commands as follows:

```shell
NE=30
GenerateCSMesh --alt --res ${NE} --file ${GRID_FILE_PATH}/ne${NE}.g
GenerateVolumetricMesh --in ${GRID_FILE_PATH}/ne${NE}.g --out ${GRID_FILE_PATH}/ne${NE}pg2.g --np 2 --uniform
ConvertMeshToSCRIP --in ${GRID_FILE_PATH}/ne${NE}pg2.g --out ${GRID_FILE_PATH}/ne${NE}pg2_scrip.nc
```
