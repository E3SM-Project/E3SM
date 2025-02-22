# Model Configuration

## Model Inputs

This section explains how input parameters are passed to EAMxx, and how the
user can change their value.
<!-- The full list of the currently configurable runtime parameters for
EAMxx can be found [here](../common/eamxx_params.md). -->

The infrastructure for reading inputs into EAMxx involves a few scripts/files:

1. `atmchange` and `atmquery`: these scripts are located in
`SRCDIR/components/eamxx/scripts`, and are soft-linked in the case folder.
As their names suggest, they can be used to query and change the runtime
configuration parameter of EAMxx.
Since these two scripts are the only scripts that the average user needs to
know and interact with, in the next sections we give a brief overview of how
they can be used and how their output can be interpreted.
Additionally, for both of the scripts a short help can also be obtained
using the `-h` flag.
2. `buildnml`: this script, located in `SRCDIR/components/eamxx/cime_config`,
is called by CIME's case management scripts (`case.setup`, `case.build`, and
`case.submit`), and is responsible for creating the input files that EAMxx will
read to load its runtime parameters.
Users should not have to modify this script,
nor should they have to manually call it, but it is useful to know what it does.
When `buildnml` runs, it creates a few files, containing EAMxx input parameters:
    - `eamxx_input.yaml`: this [YAML](https://yaml.org/spec/1.2.2)
    file is located in the `RUNDIR/data` folder, and will be read by EAMxx at
    runtime to load all of its configuration parameters.
    More precisely, this file contains parameters that need to be used inside
    the EAMxx interfaces.
    - `namelist.nl`: this
    [namelist](https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnad/index.html)
    file is located in the `RUNDIR/data` folder, and will be parsed at runtime
    to get all the parameters for the HOMME dycore (ADD REF).
    This file only contains dycore-specific parameters that are only recognized
    inside HOMME, and does not contain any parameter pertaining EAMxx infrastructure.
    - `namelist_eamxx.xml`: this XML file is located in the case directory,
    and contains all the runtime parameters that EAMxx will read in at runtime.
    `buildnml` uses this XML file as an intermediate file during the generation
    of `eamxx_input.yaml` and `namelist.nl`.
    More specifically, `buildnml` generates this file using case information to
    select the proper configurations from the file
    `namelist_defaults_eamxx.xml`, located in
    `SRCDIR/components/eamxx/cime_config`.
    Despite the fact that the only files that are needed at runtime are
    `eamxx_input.yaml` and `namelist.nl`, we generate and keep this XML file
    around to make the implementation of `atmquery` easier.

    Since these files are automatically generated when `buildnml` runs, users
    should not manually modify them.
    Any manual modification will be lost the next time `buildnml` runs
    (e.g., at `case.submit` time).

## Querying model inputs: `atmquery`

This script is the simplest way for the user to check the value and properties
of EAMxx input parameters.
A basic usage of the script is

``` {.shell .copy}
./atmquery my_param
```

which will retrieve the value of the parameter called `my_param`, by locating
the XML node "my_param" in the file `namelist_eamxx.xml` in the RUNDIR folder.
Obviously, an XML file can have multiple nodes with the same tag,
and the script is implemented to error out if multiple matches are found.
In such a scenario, the user needs to provide also the parents nodes names,
using enough parents to uniquely identify the node (in most cases,
one parent is enough).
To specify a parent, the user can prepend the parent name and `::` to the node name:

``` {.shell .copy}
./atmquery foo::my_param
```

The output will contain the fully scoped parameter name, along with the value. E.g.,

``` {.shell .copy}
$ ./atmquery foo::my_param
    namelist_defaults::node1::node2::foo::my_param:   10
```

If the input parameter is not a leaf but a node with sub-elements, the output
will recursively print all sub-elements, properly indented:

``` {.shell .copy}
$ ./atmquery homme
    homme
        Moisture: moist
        BfbHash: 18
        number_of_subcycles: 1
        enable_precondition_checks: true
        enable_postcondition_checks: true
        repair_log_level: trace
        internal_diagnostics_level: 0
        compute_tendencies: None
```

It is sometimes desirable to query _all_ the nodes that have a particular name,
or that contain a particular string.
We can do that by using the `--grep` flag:

``` {.shell .copy}
$ ./atmquery --grep sub
    iop_options::iop_dosubsidence: false
    ctl_nl::hypervis_subcycle: 1
    ctl_nl::hypervis_subcycle_tom: 1
    ctl_nl::hypervis_subcycle_q: 6
    atmosphere_processes::number_of_subcycles: 1
    sc_import::number_of_subcycles: 1
    homme::number_of_subcycles: 1
    physics::number_of_subcycles: 1
```

Similarly to the CIME utility `xmlquery`, the options `--value`, `--type`,
`--valid-values`, and `--full` can be used to respectively retrieve just the
parameter value (useful for shell scripting), the parameter's type,
a list of valid values for the parameter (when applicable), or all of the above:

``` {.shell .copy}
$ ./atmquery atm_log_level --value
    info
$ ./atmquery atm_log_level --type
    namelist_defaults::driver_options::atm_log_level: string
$ ./atmquery atm_log_level --valid-values
    namelist_defaults::driver_options::atm_log_level: ['trace', 'debug', 'info',
                                                       'warn', 'error']
$ ./atmquery atm_log_level --full
    namelist_defaults::driver_options::atm_log_level
        value: info
        type: string
        valid values: ['trace', 'debug', 'info', 'warn', 'error']
```

Finally, the option `--listall` can be used to list the whole content of the
XML file, which will be displayed with each node indented in its parent scope:

``` {.shell .copy}
$ ./atmquery --listall
    namelist_defaults
        grids_manager
            Type: Homme
            physics_grid_type: PG2
            physics_grid_rebalance: None
            dynamics_namelist_file_name: ./data/namelist.nl
            vertical_coordinate_filename: /some/path/to/coords/file.nc
        initial_conditions
            Filename: /some/path/to/ic/file.nc
            topography_filename: /some/path/to/topo/file.nc
    [...]
```

## Changing model inputs: `atmchange`

When `buildnml` runs, the model inputs are deduced from the case configuration
settings (e.g., the grid, the compset, etc.) and the
`namelist_eamxx_defaults.xml` file, located in the eamxx source tree.
The user can change any of these parameters using the `atmchange` script.
A basic usage of the script is

``` {.shell .copy}
./atmchange my_param=10
```

As for `atmquery`, if there are multiple matches for a given parameter name,
the user must specify a unique scoped name, which allows `atmchange` to
uniquely identify the XML node to modify:

``` {.shell .copy}
$ ./atmquery homme::number_of_subcycles
    namelist_defaults::atmosphere_processes::homme::number_of_subcycles: 1
$ ./atmchange number_of_subcycles=10
ERROR: internal_diagnostics_level is ambiguous. Use ANY in the node path to
allow multiple matches. Matches:
    namelist_defaults::atmosphere_processes::number_of_subcycles
    namelist_defaults::atmosphere_processes::sc_import::number_of_subcycles
    namelist_defaults::atmosphere_processes::homme::number_of_subcycles
    namelist_defaults::atmosphere_processes::physics::number_of_subcycles
    namelist_defaults::atmosphere_processes::physics::mac_aero_mic::number_of_subcycles
    namelist_defaults::atmosphere_processes::physics::mac_aero_mic::tms::number_of_subcycles
    namelist_defaults::atmosphere_processes::physics::mac_aero_mic::shoc::number_of_subcycles
    namelist_defaults::atmosphere_processes::physics::mac_aero_mic::cldFraction::number_of_subcycles
    namelist_defaults::atmosphere_processes::physics::mac_aero_mic::spa::internal_diagnostics_level
    namelist_defaults::atmosphere_processes::physics::mac_aero_mic::p3::number_of_subcycles
    namelist_defaults::atmosphere_processes::physics::rrtmgp::number_of_subcycles
    namelist_defaults::atmosphere_processes::sc_export::number_of_subcycles
$ ./atmchange homme::number_of_subcycles=10
Regenerating /path/to/namelist_eamxx.xml. Manual edits will be lost.
$ ./atmquery homme::number_of_subcycles
    namelist_defaults::atmosphere_processes::homme::number_of_subcycles: 10
```

In some cases, the user may be interested in changing _all_ nodes with a
given name.
In that case, you can use 'ANY' as a node name:

``` {.shell .copy}
$ ./atmquery --grep number_of_subcycles
    atmosphere_processes::number_of_subcycles: 1
    sc_import::number_of_subcycles: 1
    homme::number_of_subcycles: 1
    physics::number_of_subcycles: 1
    mac_aero_mic::number_of_subcycles: 24
    tms::number_of_subcycles: 1
    shoc::number_of_subcycles: 1
    cldFraction::number_of_subcycles: 1
    spa::number_of_subcycles: 1
    p3::number_of_subcycles: 1
    rrtmgp::number_of_subcycles: 1
    sc_export::number_of_subcycles: 1
$ ./atmchange ANY::number_of_subcycles=3
Regenerating /path/to/namelist_eamxx.xml. Manual edits will be lost.
$ ./atmquery --grep number_of_subcycles
    atmosphere_processes::number_of_subcycles: 3
    sc_import::number_of_subcycles: 3
    homme::number_of_subcycles: 3
    physics::number_of_subcycles: 3
    mac_aero_mic::number_of_subcycles: 3
    tms::number_of_subcycles: 3
    shoc::number_of_subcycles: 3
    cldFraction::number_of_subcycles: 3
    spa::number_of_subcycles: 3
    p3::number_of_subcycles: 3
    rrtmgp::number_of_subcycles: 3
    sc_export::number_of_subcycles: 3
```

In addition, "ANY" can be used in a "scoped" string, to limit the set of matches:

``` {.shell .copy}
$ ./atmchange mac_aero_mic::ANY::number_of_subcycles=1
Regenerating /path/to/namelist_eamxx.xml. Manual edits will be lost.
$ ./atmquery --grep number_of_subcycles
    atmosphere_processes::number_of_subcycles: 3
    sc_import::number_of_subcycles: 3
    homme::number_of_subcycles: 3
    physics::number_of_subcycles: 3
    mac_aero_mic::number_of_subcycles: 1
    tms::number_of_subcycles: 1
    shoc::number_of_subcycles: 1
    cldFraction::number_of_subcycles: 1
    spa::number_of_subcycles: 1
    p3::number_of_subcycles: 1
    rrtmgp::number_of_subcycles: 3
    sc_export::number_of_subcycles: 3
```

Since the XML file stores constraints on the parameter value
(like its type or valid values), attempting to use the wrong type will cause an error:

``` {.shell .copy}
$ ./atmquery --type se_ne
    namelist_defaults::ctl_nl::se_ne: integer
$ ./atmchange se_ne=hello
ERROR: Could not refine 'hello' as type 'integer':
```

There are three main types supported: integer, float, string, logical.
When passing a string to `atmchange`, the script will try to interpret it
according to the parameter type, and throw an error if that's not possible:

- For "string", anything works.
- For "integer", only digits are allowed, possibly with a negative sign in front;
- For "float", only digits are allowed, possibly with a negative sign in front
and a decimal point;
- For "logical", only the strings "true" and "false" are allowed (case insensitive).
There are two additional types supported: "file" and "array(T)", where "T" is
any of the other supported types (but not another array):
      - "file" is used to inform CIME of the input files that have to be download
      from E3SM data servers, like initial conditions files, or certain lookup tables.
      - "array(T)" allows to specify a list of items (of the same type),
      which will be parsed inside EAMxx as a `std::vector<T>`.

For type "string" and "array(T)", it is also possible to _append_ to the
currently stored value.

``` {.shell .copy}
$ ./atmquery homme::compute_tendencies
    namelist_defaults::atmosphere_processes::homme::compute_tendencies:
        value: a, b
        type: array(string)
        valid values: []
$ ./atmchange homme::compute_tendencies+=c
$ ./atmquery homme::compute_tendencies --full
    namelist_defaults::atmosphere_processes::homme::compute_tendencies
        value: a, b, c
        type: array(string)
        valid values: []
```

### Modifying the list of atmosphere processes

The `atmchange` script can be used to change any of the runtime parameters of EAMxx.
In particular, it can
be used to add, remove, or reorder atmosphere processes.
When adding an atmosphere process, we must first make sure that the defaults
for that process are present in `namelist_defaults_eamxx.xml`.
For instance, the default settings for the "physics" atmosphere process group
include the following:

``` {.shell .copy}
$ ./atmquery physics::atm_procs_list
    namelist_defaults::atmosphere_processes::physics::atm_procs_list: mac_aero_mic,rrtmgp
```

where "mac_aero_mic" is itself an atmosphere process group, consisting of
macrophysics, aerosols, and microphysics processes.
If one wanted to add the "cosp" atmosphere process to this list, and change the
number of its subcycles, it could do so via

``` {.shell .copy}
$ ./atmchange physics::atm_procs_list+=cosp
[...]
$ ./atmchange cosp::number_of_subcycles=3
```

Notice that if we swapped the two commands, we would get an error,
since the node "cosp" is not present in the XML generated from the defaults
until we decide to add it.

It is also possible to declare a new (empty) atmosphere process group,
which can then be filled with valid atmosphere processes via subsequent calls
to `atmchange`.
The syntax to trigger this behavior consists in specifying a process name that
begins and ends with an underscore:

``` {.shell .copy}
./atmchange physics::atm_procs_list+=_my_group_
```

This adds a new process to the list of processes in "physics",
called "\_my_group\_", and which is itself an atmosphere process group.
Hence, we can then do

``` {.shell .copy}
./atmchange _my_group_::atm_procs_list+=A,B
```

where A and B must be valid atmosphere process names (i.e., present in `namelist_defaults_eamxx.xml`)
or be themselves new atmosphere process groups (i.e., beginning/ending with an underscore)

`atmchange` can also be used to completely change a list of atmosphere processes:

``` {.shell .copy}
./atmchange physics::atm_procs_list=A,B,C
```

Notice that we used "=" instead of "+=", which means we will be overwriting the
value, rather than appending.
Any atmosphere process that was previously in the list but is no longer in it
will be removed from the generated `namelist_defaults.xml`
(and `eamxx_input.yaml`) files, along with all their nested parameters.

## Model Output

EAMxx allows the user to configure the desired model output via
[YAML](https://yaml.org/) files,
with each YAML file associated to a different output stream (i.e., a file).
In order to add an output stream,
one needs to run `atmchange output_yaml_files+=/path/to/my/output/yaml`
(more information on how to use `atmchange` can be found [here](#changing-model-inputs-atmchange)).
During the `buildnml` phase of the case management system, these YAML files
will be copied into the RUNDIR/data folder.
During this process, the files will be parsed, and any CIME-related
variable will be resolved accordingly.
Therefore, it is not advised to put the original YAML files
in RUNDIR/data, since upon `buildnml` execution, all the CIME vars will no
longer be present in the YAML file (replaced by their values),
making it harder to tweak it, and even harder to share with other users/cases.
Another consequence of this is that, much like `eamxx_input.yaml`,
the user should not modify the generated YAML files in RUNDIR/data,
since any modification will be lost on the next run of `buildnml`
(e.g., during `case.submit`).

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

Notice that lists can be equivalently specified in YAML as
`Field Names: [f1, f2, f3]`.
The user can specify fields to be outputted from any of the grids on which
they are available (although most fields are only available on ONE grid).
In the example above, we requested fields from both the Physics and Dynamics grid.
The meaning of the other parameters is as follows:

- `Averaging Type`: how the fields are integrated in time before being saved.
Valid options are:
      - `Instant`: no integration, each time frame saved corresponds to
      instantaneous values of the fields.
      - `Average`/`Max`/`Min`: the fields undergo the corresponding operation
      over the time interval since the last output step (as specified in the
      `output_control` section).
      In the case above, each snapshot saved to file corresponds to an average
      of the output fields over 6h windows.
- `filename_prefix`: the prefix of the output file, which will be created
in the run directory.
      - The full filename will be `$prefix.$avgtype.$frequnits_x$freq.$timestamp.nc`,
      where `$timestamp` corresponds to the first snapshot saved in the file for
      `Instant` output, or the beginning of the first averaging window for the other
      averaging types.
      - If not set, it defaults to `$casename.eamxx.h`.
- `Max Snapshots Per File`: specifies how many time snapshots can be put in a file.
      - Once this number is reached, EAMxx will close the file and open a new one.
      - If not set, it defaults to `-1`, signaling "unlimited storage".
- `Frequency`: how many units of time are between two consecutive writes to file.
      - For `Instant` output the fields are "sampled" at this frequency,
      while for other averaging types the fields are "integrated"
      in time over this window
- `frequency_units`: units of the output frequency.
      - Valid options are `nsteps` (the number of atmosphere time steps),
      `nsecs`, `nmins`, `nhours`, `ndays`, `nmonths`, `nyears`.

## Diagnostic output

In addition to the fields computed by EAMxx as part of the timestep, the user can
request to output derived quantities, which will be computed on the fly by the
I/O interface of EAMxx. There are two types of diagnostic outputs:

- **Quantities computed as a function of EAMxx fields**
      - These are simply physical quantities that EAMxx does not keep in
      persistent storage.
      - As of May 2024, the available derived quantities are (case sensitive):
        - `Exner`
        - `PotentialTemperature`
        - `LiqPotentialTemperature`
        - `dz`
        - `z_int`
        - `z_mid`
        - `geopotential_int`
        - `geopotential_mid`
        - `AtmosphereDensity`
        - `VirtualTemperature`
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

- **Lower-dimensional slices of a field**
      - These are hyperslices of an existing field or of another diagnostic output.
      - As of August 2023, given a field `X`, the available options are:
          - `X_at_lev_N`: slice the field `X` at the `N`-th vertical level index.
              - Recall that in EAMxx `N = 0` corresponds to the model top.
          - `X_at_model_bot`, `X_at_model_top`:
          special case for top and bottom of the model.
          - `X_at_Ymb`, `X_at_YPa`, `X_at_YhPa`:
          interpolates the field `X` at a vertical position specified by the
          given pressure `Y`.
              - Available units are `mb` (millibar), `Pa`, and `hPa`.
          - `X_at_Ym_above_Z`:
          interpolates the field `X` at a vertical height of `Y` meters above
          `Z`, with `Z = surface` or `Z = sealevel`.

## Remapped output

The following options can be used to to save fields on a different grid from
the one they are computed on.

- `horiz_remap_file`: a path to a map file (as produced by `ncremap`) between
the grid where the fields are defined and a coarser grid.
      - EAMxx will use this to remap fields on the fly, allowing to reduce
      the size of the output file.
          - **Note:** with this feature, the user can only specify fields
          from a single grid.
- `vertical_remap_file`: similar to the previous option, this map file is used to
refine/coarsen fields in the vertical direction.
- `IOGrid`: this parameter can be specified inside one of the grids sections,
and will denote the grid (which must exist in the simulation) where the fields
must be remapped before being saved to file.
      - This feature is really only used to save fields on the dynamics grid
      without saving twice the DOFs at the interface of two spectral elements.
      - E.g., for a scalar quantity defined only in the horizontal direction,
      native output from the Dynamics grid would produce arrays of length
      `nelems * ngp * ngp`, where `ngp` is the number of Gauss points along
      each axis in the 2D spectral element, and `nelems` is the number of
      horizontal elements.
      - However, due to continuity, the values on the Gauss points on the
      element boundary must match the values on the neighboring element,
      resulting in duplicated data.
      - By remapping to a "unique" version of the dynamics grid
      (which in EAMxx is referred to as "Physics GLL"), we can save roughly
      45% of storage.
      - **Note:** this feature cannot be used along with the
      horizontal/vertical remap.

## Tendencies output

It is also possible to request tendencies of fields that are updated by
atmosphere processes, on a per-process basis (here, "updated" means that the
field is both an input as well as an output of the atmosphere process).
Since the term "tendency" can be used with slightly different connotations,
we clarify what we mean by that when it comes to model output:
if process P updates field A, by the tendency of A from process P we mean
`(A_after_P - A_before_P) / dt`, where `dt` is the atmosphere timestep.

As of May 2024, the user needs two things in order to get tendencies from a
process. E.g., to get the tendencies of `T_mid` and `horiz_winds` from the
process `shoc`, one needs to:

- Set `atmchange shoc::compute_tendencies=T_mid,horiz_winds`.
- Add `shoc_T_mid_tend` and `shoc_horiz_winds_tend` to the list of fields in
the desired output YAML file.

## Additional options

The YAML file shown at the top of this section, together with the remap options
in the following section, covers most of the options used in a typical run.
There are, however, particular use cases that require some less common options,
which we list here (in parentheses, the location in the YAML file and the type
of the parameter value).

- `flush_frequency` (top-level list, `integer`)
      - This parameter can be used to specify how often the IO library
      should sync the in-memory data to file.
      - If not specified, the IO library is free to decide when it should flush
      the data.
      - This option can be helpful for debugging, in case a crash is occurring
      after a certain number of steps, but before the IO library would
      automatically flush to file.
- `Floating Point Precision` (top-level list, `string`):
      - This parameter specifies the precision to be used for floating point
      variables in the output file.
      - By default, EAMxx uses single precision.
      - Valid values are `single`, `float`, `double`, and `real`.
          - The first two are synonyms, while the latter resolves to `single`
          or `double` depending on EAMxx CMake configuration parameter `EAMXX_DOUBLE_PRECISION`.
- `file_max_storage_type` (top-level list, `string`):
      - This parameter determines how the capacity of the file is specified.
        - By default, it is set to `num_snapshots`, which makes EAMxx read
        `Max Snapshots Per File` (explained in the first section).
        - However, the user can specify `one_year` or `one_month`,
        which will make EAMxx create one output file per year/month of simulation,
        fitting however many snapshots are needed in each file
        (depending on the output frequency).
            - If `one_year` or `one_month` are used, the option
            `Max Snapshots Per File` is ignored.
- `MPI Ranks in Filename` (top-level list, `boolean`):
      - This option specifies whether the number of MPI ranks in the atmosphere
      communicator should be put inside the output file name.
      - By default, this is `false`, since it is usually not important.
      - This option is mostly important for standalone unit testing, where several
      versions of the same test (corresponding to different numbers of MPI ranks)
      are running concurrently, so that different file names are needed to
      avoid resource contention.
- `iotype` (top-level list, `string`):
      - This option allows the user to request a particular format for the
      output file.
      - The possible values are 'default', 'netcdf', 'pnetcdf' 'adios',
      'hdf5', where 'default' means "whatever is the PIO type from the case settings".
- `save_grid_data` (`output_control` sub-list, `boolean`):
      - This option allows to specify whether grid data (such as `lat`/`lon`)
      should be added to the output stream.
      - By default, it is `true`.
- `skip_t0_output` (`output_control` sub-list, `boolean`):
      - This option is relevant only for `Instant` output, for which fields are
      also output at the case start time
          - I.e., after initialization but before the beginning of the first timestep).
      - By default it is set to `false`.
- **Restart Options:**
      - When performing a restart, EAMxx attempts to restart every output
      stream listed in the `output_yaml_files` atmosphere option
      (which can be queried via `atmquery output_yaml_files`).
      - The user can specify a few options, in order to tweak the restart behavior:
          - `Perform Restart` (`Restart` sub-list, `boolean`):
              - This parameter is `true` by default, but can be set
              to `false` to force the model to ignore any history restart files,
              and start the output stream from
            scratch, as if this was an initial run.
          - `filename_prefix` (`Restart` sub-list, `string`):
              - By default, this parameter is set to match the value of
              `filename_prefix` from the top-level list.
              - This can be set to a different value in case we want to restart
              a previous simulation that was using a different filename prefix.
          - `force_new_file` (`Restart` sub-list, `boolean`):
              - This parameter allows to start a fresh new output file upon restarts.
              - By default, is is set to `false`, so that EAMxx will attempt to
              resume filling the last produced output file (if any, and if it
              can accommodate more snapshots).
