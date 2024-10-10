# Add New Grid Configuration to E3SM

In addition to generating input data to support a new grid, several code modifications are required before E3SM can run with the grid. However, the specific changes will depend on how the grid will be used. The intendend model configuration for the new grid will change which files need to be modified. For instance, a grid intended for aquaplanet experiments does not require as many changes as a historical AMIP-style run.

The guidelines here are meant to outline various possible changes the user should consider when adding support for a new grid. This document cannot be exhaustive, and it is important that the user understands the changes they are making. It is often useful to use a pre-existing grid configuration as a template. Note that the guidelines here are only relevant for "horizontal" grids. Similar considerations are needed to support a new vertical grid.

When setting up a new grid you will need to edit some or all of these files:

- `cime_config/config_grids.xml`
- `components/eam/bld/config_files/horiz_grid.xml`
- `components/eam/bld/namelist_files/namelist_defaults_eam.xml`
- `components/eam/bld/namelist_files/namelist_definition.xml`
- `components/elm/bld/namelist_files/namelist_definition.xml`
- `components/elm/bld/namelist_files/namelist_defaults.xml`

## Mono-Grid vs Bi-Grid vs Tri-Grid

The mono-bi-tri grid options in E3SM can be confusing, but it's important to understand what these terms mean when adding a new grid to E3SM. At the surface these terms mean that the whole model either using a single grid for all componennt models, or a combination of 2 or 3 grids shared among the component models. Note that mono-grid and bi-grid terms often ignore that the river model needs to be on i/ts own regular lat-lon grid (for now).

Historically, climate models would use a single grid for all components (i.e. mono-grid), but this is often not the case anymore. In E3SM the ocean and sea-ice components often use targeted regional refinement with special consideration of ocean mesoscale eddies, whereas the atmosphere will generally use a globally homogenous grid. In practice, the main difference between "bi" and "tri" grids often comes down to whether the land surface model shares a grid with the atmosphere or not. The component coupler is in charge of facilitating communication between component models, primarily through fluxes, and so mapping files are needed to support a combination of different grids. E3SMv3 uses a tri-grid configuration for production simulations.

## Grid Naming Conventions

The atmosphere grid name should always indicate the "ne" value and add "pg2" to indicate that the physgrid is being used. For a regionally refined mesh (RRM) the grid name should always start with `ne0` followed a descriptive string that includes the region being refined and the degree of refinement.

**Example**: `ne0np4_northamericax4v1`

Note that the example indicates a `4x` refinement, but does indicate the base resolution, which is useful to know. A better grid name would be `ne0np4_northamerica30x4v1`, because this tells us that the grid is consistent with `ne30` in the unrefined regions.

For a rectilinear lat-lon grid used by the land and/or river models the grid name should start with "r" and typically use spacing less than one degree, so they indicate the nominal grid spacing, starting with "0" and omitting the decimal.

**Examples**: `r05` is 0.5 degree spacing and `r0125` is 1/8 or 0.125 degree spacing.

For a mono-grid, which can only be used for idealized simulations such as aqua planet and RCE, the grid is written twice to indicate that both atmosphere/land and ocean/sea-ice models are on the same grid.

**Example**: `ne30pg2_ne30pg2`

Bi-grid options should indicate two different grids used for atmosphere/land and ocean/sea-ice models.

**Example**: `ne30pg2_IcoswISC30E3r5`

Tri-grid options should indicate three different grids used for atmosphere, land, and ocean/sea-ice models, with the land grid appearing in the middle.

**Example**: `ne30pg2_r05_IcoswISC30E3r5`

## Grid Definition

### Adding a New Grid Alias

Grid aliases are defined in specified in `cime_config/config_grids.xml` and are used to specify the grid for a case when calling `create_newcase` via the `--res` argument. Below is an example grid alias for the `ne30pg2_r05_IcoswISC30E3r5` grid used in E3SMv3 production simulations.

```xml
    <model_grid alias="ne30pg2_r05_IcoswISC30E3r5">
      <grid name="atm">ne30np4.pg2</grid>
      <grid name="lnd">r05</grid>
      <grid name="ocnice">IcoswISC30E3r5</grid>
      <grid name="rof">r05</grid>
      <grid name="glc">null</grid>
      <grid name="wav">null</grid>
      <mask>IcoswISC30E3r5</mask>
    </model_grid>
```

### Domain Files

Domain files are needed for each grid and are specified in the `<domains>` section of `cime_config/config_grids.xml`. The default domain files are grouped by the atmosphere grid. The section for the typical `ne30pg2` grid looks as follows:

```xml
    <domain name="ne30np4.pg2">
      <nx>21600</nx>
      <ny>1</ny>
      ...
      <file grid="atm|lnd" mask="IcoswISC30E3r5">$DIN_LOC_ROOT/share/domains/domain.lnd.ne30pg2_IcoswISC30E3r5.231121.nc</file>
      <file grid="ice|ocn" mask="IcoswISC30E3r5">$DIN_LOC_ROOT/share/domains/domain.ocn.ne30pg2_IcoswISC30E3r5.231121.nc</file>
      ...
      <desc>ne30np4.pg2 is Spectral Elem 1-deg grid w/ 2x2 FV physics grid per element:</desc>
    </domain>
```

Notice where I've used ellipses `...` to omit all entires except the lines relevant to the `ne30pg2_r05_IcoswISC30E3r5` grid. Also, note that all of these paths are relative to the input data path set as `DIN_LOC_ROOT` which has a default for each machine. See [Generating Domain Files](/generate_domain_files/) for information about creating domain files.

### Coupler Mapping Files

The mapping files used by the component coupler to communicate fluxes between the component models must be specified in the `<gridmaps>` section of `cime_config/config_grids.xml`. these are organized for specific pairs of grids, such that tri-grids will require multiple sections. The entries relevant for `ne30pg2_r05_IcoswISC30E3r5` are shown below.

```xml
    <gridmap atm_grid="ne30np4.pg2" ocn_grid="IcoswISC30E3r5">
      <map name="ATM2OCN_FMAPNAME">cpl/gridmaps/ne30pg2/map_ne30pg2_to_IcoswISC30E3r5_traave.20231121.nc</map>
      <map name="ATM2OCN_VMAPNAME">cpl/gridmaps/ne30pg2/map_ne30pg2_to_IcoswISC30E3r5_trbilin.20231121.nc</map>
      <map name="ATM2OCN_SMAPNAME">cpl/gridmaps/ne30pg2/map_ne30pg2_to_IcoswISC30E3r5-nomask_trbilin.20231121.nc</map>
      <map name="OCN2ATM_FMAPNAME">cpl/gridmaps/IcoswISC30E3r5/map_IcoswISC30E3r5_to_ne30pg2_traave.20231121.nc</map>
      <map name="OCN2ATM_SMAPNAME">cpl/gridmaps/IcoswISC30E3r5/map_IcoswISC30E3r5_to_ne30pg2_traave.20231121.nc</map>
      <map name="ATM2ICE_FMAPNAME_NONLINEAR">cpl/gridmaps/ne30pg2/map_ne30pg2_to_IcoswISC30E3r5_trfvnp2.20231121.nc</map>
      <map name="ATM2OCN_FMAPNAME_NONLINEAR">cpl/gridmaps/ne30pg2/map_ne30pg2_to_IcoswISC30E3r5_trfvnp2.20231121.nc</map>
    </gridmap>
```

```xml
    <gridmap atm_grid="ne30np4.pg2" lnd_grid="r05">
      <map name="ATM2LND_FMAPNAME">cpl/gridmaps/ne30pg2/map_ne30pg2_to_r05_traave.20231130.nc</map>
      <map name="ATM2LND_FMAPNAME_NONLINEAR">cpl/gridmaps/ne30pg2/map_ne30pg2_to_r05_trfvnp2.230516.nc</map>
      <map name="ATM2LND_SMAPNAME">cpl/gridmaps/ne30pg2/map_ne30pg2_to_r05_trbilin.20231130.nc</map>
      <map name="LND2ATM_FMAPNAME">cpl/gridmaps/ne30pg2/map_r05_to_ne30pg2_traave.20231130.nc</map>
      <map name="LND2ATM_SMAPNAME">cpl/gridmaps/ne30pg2/map_r05_to_ne30pg2_traave.20231130.nc</map>
    </gridmap>
```

```xml
    <gridmap atm_grid="ne30np4.pg2" rof_grid="r05">
      <map name="ATM2ROF_FMAPNAME">cpl/gridmaps/ne30pg2/map_ne30pg2_to_r05_traave.20231130.nc</map>
      <map name="ATM2ROF_FMAPNAME_NONLINEAR">cpl/gridmaps/ne30pg2/map_ne30pg2_to_r05_trfvnp2.230516.nc</map>
      <map name="ATM2ROF_SMAPNAME">cpl/gridmaps/ne30pg2/map_ne30pg2_to_r05_trbilin.20231130.nc</map>
    </gridmap>
```

Note that all of these paths are relative to the input data path set as `DIN_LOC_ROOT` which has a default for each machine.

### Defining a New Atmosphere Grid

When defining a new atmosphere grid, information needs to be provided on how the grid is constructed.

To define a new atmosphere grid a line must be added to `components/eam/bld/config_files/horiz_grid.xml` that indicates the numebr of elements and physics columns. In the lines below for `ne30np4` (without the physgrid) and `ne30pg2` (with the physgrid) you can see the value of `ne` is the same (number of elements along a cube edge), but the number of physics columns is different.

```xml
<horiz_grid dyn="se"    hgrid="ne30np4"      ncol="48602"   csne="30"  csnp="4" npg="0" />
<horiz_grid dyn="se"    hgrid="ne30np4.pg2"  ncol="21600"   csne="30"  csnp="4" npg="2" />
```

An explanation of how to calculate the number of physics columns can be found here: [Atmosphere Grid Overview](../../../EAM/tech-guide/atmosphere-grid-overview.md).

For a grid with regional refinement, follow the conventions of other grids in this file. There is no formula to calculate the number of columns for RRM grids, but the value can be obtained from the grid files used for mapping.

```xml
<horiz_grid dyn="se" hgrid="ne0np4_antarcticax4v1.pg2"        ncol="48836"  csne="0" csnp="4" npg="2" />
```

### Defining a New Land Grid

If you are creating a new grid that will be used by the land model the grid name needs to be added to the list `valid_values` associated with the `res` entry in the file `components/elm/bld/namelist_files/namelist_definition.xml` that holds the definition of namelist variables used by the land model.

```xml
<entry id="res" type="char*30" category="default_settings"
       group="default_settings"
       valid_values=
"512x1024,360x720cru,128x256,64x128,...">
Horizontal resolutions
Note: 0.1x0.1, 0.5x0.5, 5x5min, 10x10min, 3x3min and 0.33x0.33 are only used for CLM tools
</entry>
```

Simply add the name of your new grid to the list of `valid_values`.

## Namelist Variable Defaults

Each new grid will likely need various new default parameter values to be specified. These parameters can be set for individual simulations by editing the `user_nl_*` files in the case directory, but for these to become defaults any time the grid is used then new defaults need to be specified.

The lists below show namelist parameters that may need to be specified for a new grid. The creator of a new grid is responsible for understanding these parameters and deciding when new defaults are appropriate.

### Atmosphere Namelist Parameters

- `drydep_srf_file`     - Data file for surface aerosol deposition
- `bnd_topo`            - Surface topography (smoothed for target grid)
- `mesh_file`           - HOMME np4 mesh file (exodus format)
- `se_tstep`            - HOMME time step [seconds]
- `dt_remap_factor`     - HOMME vertical remap factor
- `dt_tracer_factor`    - HOMME tracer advection factor
- `hypervis_subcycle_q` - HOMME tracer hyperviscosity factor

### Land Namelist Parameters

- `fsurdat`             - Surface data file
- `finidat`             - Land model initial condition file
- `flanduse_timeseries` - Time-evolving land-use data file

Back to step-by-step guide for [Adding Support for New Grids](../adding-grid-support-step-by-step-guide.md)
