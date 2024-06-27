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
of this is that the user should not modify the YAML files in RUNDIR/data, since any modification will
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
    - `X_at_Ym_above_Z`: interpolates the field `X` at a vertical height of `Y` meters above
      `Z`, with `Z=surface` or `Z=sealevel`.

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
  E.g., for a scalar quantity defined only in the horizontal direction, native output
  from the Dynamics grid would produce arrays of length `nelems*ngp*ngp`, where `ngp` 
  is the number of Gauss points along each axis in the 2d spectral element, and `nelems`
  is the number of horizontal elements. However, due to continuity, the values on the
  Gauss points on the element boundary must match the values on the neighboring element,
  resulting in duplicated data. By remapping to a "unique" version of the dynamics grid
  (which in EAMxx is referred to as "Physics GLL"), we can save roughly 45% of storage.
  Note: this feature cannot be used along with the horizontal/vertical remap.

## Tendencies output

It is also possible to request tendencies of fields that are updated by atmosphere processes,
on a per-process basis (here, "updated" means that the field is both an input as well as an output
of the atmosphere process). Since the term "tendency" can be used with slightly different connotations,
we clarify what we mean by that when it comes to model output: if process P updates field A,
by the tendency of A from process P we mean `(A_after_P - A_before_P) / dt`, where `dt` is
the atmosphere timestep.

As of May 2024, the user needs two things in order to get tendencies from a process. E.g., to get
the tendencies of `T_mid` and `horiz_winds` from the process `shoc`, one needs:

- `atmchange shoc::compute_tendencies=T_mid,horiz_winds`;
- add `shoc_T_mid_tend` and `shoc_horiz_winds_tend` to the list of fields in the desired output YAML file.

## Additional options

The YAML file shown at the top of this section, together with the remap options in the following
section, covers most of the options used in a typical run. There are however particular use cases
that require some less common options, which we list here (in parentheses, the location in the YAML
file and the type of the parameter value).

- `flush_frequency` (toplevel list, integer): this parameter  can be used to specify how often the IO
  library should sync the in-memory data to file. If not specified, the IO library is free to decide
  when it should flush the data. This option can be helpful for debugging, in case a crash is occurring
  after a certain number of steps, but before the IO library would automatically flush to file.
- `Floating Point Precision` (toplevel list, string): this parameter specifies the precision to be used for floating
  point variables in the output file. By default, EAMxx uses single precision. Valid values are
  `single`, `float`, `double`, and `real`. The first two are synonyms, while the latter resolves
  to `single` or `double` depending on EAMxx cmake configuration parameter `SCREAM_DOUBLE_PRECISION`.
- `file_max_storage_type` (toplevel list, string): this parameter determines how the capacity of the file is specified.
  By default, it is set to `num_snapshots`, which makes EAMxx read `Max Snapshots Per File` (explained
  in the first section). However, the user can specify `one_year` or `one_month`, which will make
  EAMxx create one output file per year/month of simulation, fitting however many snapshots are needed
  in each file (depending on the output frequency). If `one_year` or `one_month` are used, the option
  `Max Snapshots Per File` is ignored.
- `MPI Ranks in Filename` (toplevel list, boolean): this option specifies whether the number of MPI ranks in the atm
  communicator should be put inside the output file name. By default, this is `false`, since it is
  usually not important. This option is mostly important for standalone unit testing, where several
  versions of the same test (corresponding to different numbers of MPI ranks) are running concurrently,
  so that different file names are needed to avoid resource contention.
- `save_grid_data` (`output_control` sublist, boolean): this option allows to specify whether grid data
  (such as `lat`/`lon`) should be added to the output stream. By default, it is `true`.
- `iotype` (toplevel list, string): this option allows the user to request a particular format for the output
  file. The possible values are `default`, `netcdf`, `pnetcdf, `adios`, `hdf5`, where `default` means
  "whatever is the PIO type from the case settings".
- `skip_t0_output` (`output_control` sublist, boolean): this option is relevant only for `Instant` output,
  where fields are also outputed at the case start time (i.e., after initialization but before the beginning
  of the first timestep). By default it is set to `false`.
- restart options: when performing a restart, EAMxx attempts to restart every output stream listed in
  the `output_yaml_files` atm option (which can be queried via `atmquery output_yaml_files`). The user
  can specify a few options, in order to tweak the restart behavior:
  - `Perform Restart` (`Restart` sublist, boolean): this parameter is `true` by default, but can be set
    to `false` to force the model to ignore any history restart files, and start the output stream from
    scratch, as if this was an initial run.
  - `filename_prefix` (`Restart` sublist, string): by default, this parameter is set to match the value
    of `filename_prefix` from the toplevel list. It can be set to something else in case we want to
    restart a previous simulation that was using a different filename prefix.
  - `force_new_file` (`Restart` sublist, boolean): this parameter allows to start a fresh new output file
    upon restarts. By default, is is set to `false`, so that EAMxx will attempt to resume filling the
    last produced output file (if any, and if it can accommodate more snapshots).



