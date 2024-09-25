# Model inputs

This section explains how input parameters are passed to EAMxx, and how the user can
change their value. The full list of the currently configuraable runtime parameters for
EAMxx can be found [here](../common/eamxx_params.md).

The infrastructure for reading inputs into EAMxx involves a few scripts/files:

1. `atmchange` and `atmquery`: these scripts are located in `SRCDIR/components/eamxx/scripts`,
    and are soft-linked in the case folder. As their names suggest, they can be used to query
    and change the runtime configuration parameter of EAMxx. Since these two scripts are the
    only scripts that the average user needs to know and interact with, in the next sections
    we give a brief overview of how they can be used and how their output can be interpreted.
    Additionally, for both of the scripts a short help can also be obtained using the `-h` flag.
2. `buildnml`: this script, located in `SRCDIR/components/eamxx/cime_config`, is called by CIME's
    case management scripts (`case.setup`, `case.build`, and `case.submit`), and is responsible for creating
    the input files that EAMxx will read to load its runtime parameters. Users should not have to modify this script,
    nor should they have to manually call it, but it is useful to know what it does.
    When `buildnml` runs, it creates a few files, containing EAMxx input parameters:
    -   `scream_input.yaml`: this [YAML](https://yaml.org/spec/1.2.2) file is located in the `RUNDIR/data` folder,
        and will be read by EAMxx at runtime to load all of its configuration parameters. More precisely,
        this file contains parameters that need to be used inside the EAMxx interfaces.
    -   `namelist.nl`: this [namelist](https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnad/index.html) file is
        located in the `RUNDIR/data` folder, and will be parsed at runtime to get all the parameters for the
        HOMME dycore (ADD REF). This file only contains dycore-specific parameters that are only recognized
        inside HOMME, and does not contain any parameter pertaining EAMxx infrastructure.
    -   `namelist_scream.xml`: this XML file is located in the case directory, and contains all the runtime parameters that EAMxx
        will read in at runtime. `buildnml` uses this XML file as an intermediate file during the generation of
        `scream_input.yaml` and `namelist.nl`. More specifically, `buildnml` generates this file using case information to select the
        proper configurations from the file `namelist_defaults_scream.xml`, located in `SRCDIR/components/eamxx/cime_config`.
        Despite the fact that the only files that are needed at runtime are `scream_input.yaml` and `namelist.nl`,
        we generate and keep this XML file around to make the implementation of `atmquery` easier.

    Since these files are automatically generated when `buildnml` runs, users should not manually modify them.
    Any manual modification will be lost the next time `buildnml` runs (e.g., at `case.submit` time).

## Querying model inputs: atmquery

This script is the simplest way for the user to check the value and properties of EAMxx input parameters.
A basic usage of the script is

```shell
$ ./atmquery my_param
```

which will retrieve the value of the parameter called `my_param`, by locating the XML node "my_param" in the
file `namelist_scream.xml` in the RUNDIR folder. Obviously, an XML file can have multiple nodes with the same tag,
and the script is implemented to error out if multiple matches are found. In such a scenario, the user needs
to provide also the parents nodes names, using enough parents to uniquely identify the node (in most cases,
one parent is enough). To specify a parent, the user can prepend the parent name and `::` to the node name:

```shell
$ ./atmquery foo::my_param
```

The output will contain the fully scoped parameter name, along with the value. E.g.,

```shell
$ ./atmquery foo::my_param
    namelist_defaults::node1::node2::foo::my_param:   10
```

It is sometimes desirable to query _all_ the nodes that have a particular name, or that contain a particular
string. We can do that by using the `--grep` flag:

```shell
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

TODO: This difference between basic and `--grep` is not really intuitive: as pointed out in [this
issue](https://github.com/E3SM-Project/scream/issues/2413), we should change this. If we do, don't forget
to update this following part of the docs.
Using the `--grep` option has another effect: if the match is not a leaf of the XML tree, all its subelements
are printed:

```shell
$ ./atmquery --grep homme
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

Similarly to the CIME utility `xmlchange`, the options `--value`, `--type`, `--valid-values`, and `--full` can be
used to respectively retrieve just the parameter value (useful for shell scripting), the parameter's type,
a list of valid values for parameter (when applicable), or all of the above:

```shell
$ ./atmquery atm_log_level --value
    info
$ ./atmquery atm_log_level --type
    namelist_defaults::driver_options::atm_log_level: string
$ ./atmquery atm_log_level --valid-values
    namelist_defaults::driver_options::atm_log_level: ['trace', 'debug', 'info', 'warn', 'error']
$ ./atmquery atm_log_level --full
      namelist_defaults::driver_options::atm_log_level
        value: info
        type: string
        valid values: ['trace', 'debug', 'info', 'warn', 'error']
```

Finally, the option `--listall` can be used to list the whole content of the XML file, which will be displayed
with each node indented in its parent scope:

```shell
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
    ...
```

## Changing model inputs: atmchange

When `buildnml` runs, the model inputs are deduced from the case configuration settings (e.g., the grid,
the compset, etc.) and the `namelist_scream_defaults.xml` file, located in the eamxx source tree.
The user can change any of these parameters using the `atmchange` script.
A basic usage of the script is

```shell
$ ./atmchange my_param=10
```

As for `atmquery`, if there are multiple matches for a given parameter name, the user must specify a unique
scoped name, which allows `atmchange` to uniquely identify the XML node to modify:

```shell
$ ./atmquery homme::number_of_subcycles
    namelist_defaults::atmosphere_processes::homme::number_of_subcycles: 1
$ ./atmchange number_of_subcycles=10
ERROR: internal_diagnostics_level is ambiguous. Use ANY in the node path to allow multiple matches. Matches:
  namelist_defaults::atmosphere_processes::number_of_subcycles
  namelist_defaults::atmosphere_processes::sc_import::number_of_subcycles
  namelist_defaults::atmosphere_processes::homme::number_of_subcycles
  namelist_defaults::atmosphere_processes::physics::number_of_subcycles
  namelist_defaults::atmosphere_processes::physics::mac_aero_mic::number_of_subcycles
  namelist_defaults::atmosphere_processes::physics::mac_aero_mic::tms::number_of_subcycles
  namelist_defaults::atmosphere_processes::physics::mac_aero_mic::shoc::number_of_subcycles
  namelist_defaults::atmosphere_processes::physics::mac_aero_mic::cldFraction::number_of_subcycles
  namelist_defaults::atmosphere_processes::physics::mac_aero_mic::p3::number_of_subcycles
  namelist_defaults::atmosphere_processes::physics::rrtmgp::number_of_subcycles
  namelist_defaults::atmosphere_processes::sc_export::number_of_subcycles
$ ./atmchange homme::number_of_subcycles=10
Regenerating /path/to/namelist_scream.xml. Manual edits will be lost.
$ ./atmquery homme::number_of_subcycles
    namelist_defaults::atmosphere_processes::homme::number_of_subcycles: 10
```

In some cases, the user may be interested in changing _all_ nodes with a given name. In that case,
you can use 'ANY' as a node name:

```shell
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
Regenerating /path/to/namelist_scream.xml. Manual edits will be lost.
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

Since the XML file stores constraints on the parameter value (like its type or valid values), attempting to use the
wrong type will cause an error:

```shell
$ ./atmquery --type se_ne
    namelist_defaults::ctl_nl::se_ne: integer
$ ./atmchange se_ne=hello
ERROR: Could not refine 'hello' as type 'integer':
```

There are three main types supported: integer, float, string, logical. When passing a string to `atmchange`,
the script will try to interpret it acoording to the parameter type, and throw an error if that's not possible:
for "string", anything works; for "integer", only digits are allowed, possibly with a negative sign in front;
for "float", only digits are allowed, possibly with a negative sign in front and a decimal point; for "logical",
only the strings "true" and "false" are allowed (case insensitive). There are two additional types supported:
"file" and "array(T)", where "T" is any of the other supported types (but not another array):
 - "file" is used to inform CIME of the input files that have to be download from E3SM data servers, like initial conditions files,
or certain lookup tables.
 - "array(T)" allows to specify a list of items (of the same type), which will be parsed inside EAMxx as
   a `std::vector<T>`.

For type "string" and "array(T)", it is also possible to _append_ to the currently stored value

```shell
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

The `atmchange` script can be used to change any of the runtime parameters of EAMxx. In particular, it can
be used to add, remove, or reorder atmosphere processes. When adding an atmosphere process, we must first
make sure that the defaults for that process are present in `namelist_defaults_scream.xml`. For instance,
the default settings for the "physics" atmosphere process group include the following:

```shell
$ ./atmquery physics::atm_procs_list
    namelist_defaults::atmosphere_processes::physics::atm_procs_list: mac_aero_mic,rrtmgp
```

where "mac_aero_mic" is itself an atmosphere process group, consisting of macrophysics, aerosols, and microphysics
processes. If one wanted to add the "cosp" atmosphere process to this list, and change the number of its
subcycles, it could do so via

```shell
$ ./atmchange physics::atm_procs_list+=cosp
$ ./atmchange cosp::number_of_subcycles=3
```

Notice that if we swapped the two commands, we would get an error, since the node "cosp" is not present in
the XML generated from the defaults until we decide to add it.

It is also possible to declare a new (empty) atmosphere process group, which can then be filled with valid
atmosphere processes via subsequent calls to `atmchange`. The syntax to trigger this behavior consists
in specifying a process name that begins and ends with an underscore:

```shell
$ ./atmchange physics::atm_procs_list+=_my_group_
```

This adds a new process to the list of processes in "physics", called "\_my_group\_", and which is itself
an atmosphere process group. Hence, we can then do

```shell
$ ./atmchange _my_group_::atm_procs_list+=A,B
```

where A and B must be valid atmosphere process names (i.e., present in `namelist_defaults_scream.xml`)
or be themselves new atmosphere process groups (i.e., beginning/ending with an underscore)

`atmchange` can also be used to completely change a list of atmosphere processes:

```shell
$ ./atmchange physics::atm_procs_list=A,B,C
```

Notice that we used "=" instead of "+=", which means we will be overwriting the value, rather than appending.
Any atmosphere process that was previously in the list but is no longer in it will be removed from
the generated `namelist_defaults.xml` (and `scream_input.yaml`) files, along with all their nested parameters.
