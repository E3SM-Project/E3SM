# Model output

EAMxx allows the user to configure the desired model output via [YAML](https://yaml.org/) files,
with each YAML file associated to a different output file. In order to add an output stream,
one needs to run `atmchange output_yaml_files+=/path/to/my/output/yaml` (more information on how
to use `atmchange` can be found [here](./model_input.md#changing-model-inputs-atmchange)).
During the `buildnml` phase of the case management system, a copy of these YAML files will be copied
into the RUNDIR/data folder. During this process, the files will be parsed, and any CIME-related
variable will be resolved accordingly. Therefore, it is not advised to put the original YAML files in RUNDIR/data,
since upon `buildnml` execution, all the CIME vars will no longer be available in the YAML file,
making it harder to tweak it, and even harder to share with other users/cases. Another consequence
of this is that the user should not modify the yaml files in RUNDIR/data, since any modification will
be lost on the next run of `buildnml`.

## Basic output

The following is a basic example of an output request.

```yaml
%YAML 1.1
---
filename_prefix: my_output
Averaging Type: Average
Max Snapshots Per File: 10
Fields:
  Physics:
    Field Names:
      - T_mid
      - qv
  Dynamics:
    Field Names:
      - dp3d_dyn
      - omega_dyn
output_control:
  Frequency: 6
  frequency_units: nhours
```

Notice that lists can be equivalently specified in YAML as `Field Names: [f1, f2, f3]`.
The user can specify fields to be outputted from any of the grids used in the simulation.
In the example above, we requested fields from both the Physics and Dynamics grid.
The meaning of the other parameters is as follows:

- `Averaging Type`: how the fields are integrated in time before being saved. Valid
  options are

    - Instant: no integration, each time frame saved corresponds to instantaneous values
      of the fields
    - Average/Max/Min: the fields undergo the corresponding operation over the time
      interval specified in the `output_control` section. In the case above, each snapshot
      saved to file corresponds to an average of the output fields over 6h windows.

- `filename_prefix`: the prefix of the output file, which will be created in the run
  directory. The full filename will be `$prefix.$avgtype.$frequnits_x$freq.$timestamp.nc`,
  where $timestamp corresponds to the first snapshot saved in the file for Instant output,
  or the beginning of the first averaging window for the other averaging types
- `Max Snapshots Per File`: specifies how many time snapshots can be put in a file. Once
  this number is reached, EAMxx will close the file and open a new one.
- `Frequency`: how many units of time are between two consecutive writes to file. For
  Instant output the fields are "sampled" at this frequency, while for other averaging
  types the fields are "integrated" in time over this window
- `frequency_units`: units of the output frequency. Valid options are `nsteps` (the
  number of atmosphere time steps), `nsecs`, `nmins`, `nhours`, `ndays`, `nmonths`,
  `nyears`.

## Diagnostic output

In addition to the fields computed by EAMxx as part of the timestep, the user can
request to output derived quantities, which will be computed on the fly by the
I/O interface of EAMxx. There are two types of diagnostic outputs:

- quantities computed as a function of EAMxx fields. These are simply physical quantities
  that EAMxx does not keep in persistent storage. As of May 2024, the available
  derived quantities are (case sensitive):

    - `PotentialTemperature`
    - `AtmosphereDensity`
    - `Exner`
    - `VirtualTemperature`
    - `z_int`
    - `z_mid`
    - `geopotential_int`
    - `geopotential_mid`
    - `dz`
    - `DryStaticEnergy`
    - `SeaLevelPressure`
    - `LiqWaterPath`
    - `IceWaterPath`
    - `VapWaterPath`
    - `RainWaterPath`
    - `RimeWaterPath`
    - `ShortwaveCloudForcing`
    - `LongwaveCloudForcing`
    - `RelativeHumidity`
    - `ZonalVapFlux`
    - `MeridionalVapFlux`
    - `precip_liq_surf_mass_flux`
    - `precip_ice_surf_mass_flux`
    - `precip_total_surf_mass_flux`
    - `surface_upward_latent_heat_flux`
    - `wind_speed`
    - `AerosolOpticalDepth550nm`
    - `NumberPath`
    - `AeroComCld`

  TODO: add some information about what each diagnostic is, perhaps a formula

- lower-dimensional slices of a field. These are hyperslices of an existing field or of
  another diagnostic output. As of August 2023, given a field X, the available options
  are:

    - `X_at_lev_N`: slice the field `X` at the N-th vertical level index. Recall that
      in EAMxx N=0 corresponds to the model top.
    - `X_at_model_bot`, `X_at_model_top`: special case for top and bottom of the model.
    - `X_at_Ymb`, `X_at_YPa`, `X_at_YhPa`: interpolates the field `X` at a vertical position
      specified by the give pressure `Y`. Available units are `mb` (millibar), `Pa`, and `hPa`.
    - `X_at_Ym`: interpolates the field `X` at a vertical height of `Y` meters.

## Remapped output

The following options can be used to to save fields on a different grid from the one
they are computed on.

- `horiz_remap_file`: a path to a map file (as produced by `ncremap`) between the grid
  where the fields are defined and a coarser grid. EAMxx will use this to remap fields
  on the fly, allowing to reduce the size of the output file. Note: with this feature,
  the user can only specify fields from a single grid.
- `vertical_remap_file`: similar to the previous option, this map file is used to
  refine/coarsen fields in the vertical direction.
- `IOGrid`: this parameter can be specified inside one of the grids sections, and will
  denote the grid (which must exist in the simulation) where the fields must be remapped
  before being saved to file. This feature is really only used to save fields on the
  dynamics grid without saving twice the DOFs at the interface of two spectral elements.
  In fact, native output from the Dynamics grid would produce `6*num_elems*ngp*ngp`,
  where `ngp` is the number of Gauss points along each axis in the 2d spectral element.
  Note: this feature cannot be used along with the horizontal/vertical remapper.

## Tendencies output

EXPLAIN how to request tendencies on a per-process basis

## Additional options

LIST ADDITIONAL OPTIONS, like flush_frequency, file_max_storage_type, etc.

