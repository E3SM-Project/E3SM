
# MPAS-Ocean Quick Start

This MPAS-Ocean Quick Start Guide describes how to set up and run MPAS-Ocean within E3SM. More details can be found in the [MPAS-Ocean User's Guide](https://zenodo.org/records/11098080), as well as instructions on running the stand-alone ocean model.

## Steps to build and run MPAS-Ocean

Step-by-step instructions on how to run E3SM can be found at [https://docs.e3sm.org/running-e3sm-guide](https://docs.e3sm.org/running-e3sm-guide).

This MPAS-Ocean Quick Start guide provides the additional information to run configurations that are not fully-coupled (e.g. C-case: active ocean only; G-case: active ocean and sea ice) within E3SM. This is done by changing the compset. Certain parameters, including the mesh, namelist parameters, input data files, and output file specifcations can also be modified. These are described below as ways to customize runs.

Templates of 1-month example E3SM run-scripts for a [G-case](../gcase.template.v3LR.chrysalis.sh) and [C-case](../ccase.template.v3LR.anvil.sh) are provided. Key information for the user to modify includes:

- `MACHINE` and `PROJECT` if applicable
- Paths for code repository in `CODE_ROOT` and case directory in `CASE_ROOT`, which includes the run directory for output.
- Simulation compsets, resolution and name (see below)
- Wallclock duration in `WALLTIME` and simulation duration with `STOP_OPTION` and `STOP_N`.

Additional runscript examples can be found [here](https://github.com/E3SM-Project/SimulationScripts/tree/master/archive/CoupledGroup/v3.LR) for v3.LR.

## Scientifically supported compsets and meshes

### Compsets

The compsets below are typical ocean and sea ice-focused compsets supported by E3SM:

`GMPAS-JRA1p5` - Active ocean-sea ice configuration forced by data atmosphere based on JRA55 v1.5 (covers 63 years, 1958-2020)

`GMPAS-IAF` - Active ocean-sea ice configuration forced by data atmosphere based on CORE-II (covers 62 years, 1948-2009)

`GMPAS-JRA1p5-DIB-DISMF` - Active ocean-sea ice configuration forced by JRA55 v1.5 atmosphere (as above), with data iceberg and data ice-shelf melt

`GMPAS-JRA1p5-DIB-PISMF` - Active ocean-sea ice configuration forced by JRA55 v1.5 atmosphere (as above), with data iceberg and prognostic ice-shelf melt

`CMPASO-JRA1p4` - Active ocean configuration forced by data atmosphere based on JRA55 v1.4 (covers 61 years, 1958-2018)

`CMPASO-IAF` - Active ocean configuration forced by data atmosphere based on CORE-II (covers 62 years, 1948-2009)

Additional compsets can be found in the [mpas-ocean `config_compsets.xml`](https://github.com/E3SM-Project/E3SM/blob/master/components/mpas-ocean/cime_config/config_compsets.xml). Note that the fully coupled compsets and their aliases can be found in the [cime allactive `config_compsets.xml`](https://github.com/E3SM-Project/E3SM/blob/master/cime_config/allactive/config_compsets.xml).
For more information on the schemes used within MPAS-Ocean, refer to the [MPAS-Ocean User's Guide](https://zenodo.org/records/11098080).

A full list of Compsets in the current repository can be listed using

```text
cd cime/scripts
./query_config --compsets
```

### Meshes

Some supported meshes for G- and C-cases include:

`TL319_IcoswISC30E3r5` - Icosahedral 30 km mesh with ice shelves cavities (wISC), E3SMv3 (E3) revision r5, TL319 is the grid for JRA.

`T62_IcoswISC30E3r5` - Icosahedral 30 km mesh with ice shelves cavities (wISC), E3SMv3 (E3) revision r5, T62 is the grid for CORE-II.

`T62_oQU240` - Quasi uniform 2-degree ocean mesh, to be used with CORE-II only. Good for rapid testing, used in nightly testing **not** production runs.  This grid is not scientifically validated .

Note: the mesh should be consistent with the compset (e.g. JRA vs CORE). Additional mesh information can be found [here](https://github.com/E3SM-Project/E3SM/blob/master/cime_config/config_grids.xml).

A full list of Meshes in the current repository can be listed using

```text
cd cime/scripts
./query_config --grids
```

## Customizing runs

### Namelist changes

Without additional input, E3SM will generate the namelist file `mpaso_in` in the run directory using the default values for the compset and mesh requested.
Namelist parameters can be changed from default values by modifying the `user_nl_mpaso` file, found in the `case_scripts` directory. This is done by entering ``[namelist option] = [changed value]`` as separate lines in the ``user_nl_mpaso`` file. All other options will remain defaults. These changes can be added at run time and will take effect in the next submission.

Refer to the [MPAS-Ocean User's Guide](https://zenodo.org/records/11098080) (Chapter 10) for a comprehensive description of the namelist parameters and the options that they correspond to. Namelist options may also be found in the code repository in the file `components/mpas-ocean/src/Registry.xml` for general flags, `components/mpas-ocean/src/tracer_groups/Registry_*.xml` for specific tracer group flags, and `components/mpas-ocean/src/analysis_members/Registry_*.xml` for analysis member flags.

#### Example of a namelist change via `user_nl_mpaso`

```text
 config_GM_closure = 'constant'
 config_gm_constant_kappa = 900
 config_time_integrator = 'split_explicit'
 config_am_timeseriesstatsmonthly_compute_interval = '00-00-01_00:00:00'
```

In this example, the namelist changes include changing the eddy closure (first 2 options for the type of closure and a parameter value), switching the time integration scheme, and modifying an analysis member (in this case, the interval from which the monthly analysis member is computed).

Reminder: `user_nl_mpaso` can be empty. All options not specified are defaults (given the compset and mesh). Some options (like interior restoring) require extra fields to be present in the input file.  

### Configuring input and output for MPAS-Ocean

The reading and writing of model fields in MPAS is handled by user-configurable streams. A stream
represents a fixed set of model fields, together with dimensions and attributes, that are all written
or read together to or from the same file or set of files. They are used for reading initial conditions, for writing and reading restart fields, and for writing additional model history fields.
Streams are defined in XML configuration files that are created at build time for each model core. The name of this XML file for the ocean core is `streams.ocean` (the sea ice has a similar `streams.seaice`). Importantly, the stream file is generated in the run directory during the ``./case.setup`` step, but **changes made into the run directory will not take effect**. To make changes to the output fields, **copy the `streams.ocean` file from the run directory into the``case_scripts/SourceMods/src.mpaso/`` directory**. Changes to the stream file made into the ``SourceMods`` sub-directory will take effect on the next case submission (there is no need to re-compile after making modifications to the XML file). Alternatively, changes to the streams file can be made directly in the code in ``components/mpas-ocean/cime_config/buildnml``.

#### Checking initial conditions

Key information to check regarding the input data typically include:

```text
<immutable_stream name="mesh"
                  type="none"
                  io_type="pnetcdf"
                  filename_template="/lcrc/group/e3sm/data/inputdata/ocn/mpas-o/IcoswISC30E3r5/mpaso.IcoswISC30E3r5.20231120.nc"
/>
<immutable_stream name="input"
                  type="input"
                  io_type="pnetcdf"
                  input_interval="initial_only"
                  filename_template="/lcrc/group/e3sm/data/inputdata/ocn/mpas-o/IcoswISC30E3r5/mpaso.IcoswISC30E3r5.20231120.nc"
/>
```

The `mesh` filename points to the mesh file used. The `input` filename points to the file containing the ocean initial conditions (if the run type is `initial`).
The streams file can be large, it is often useful to rely on the search function to navigate it. For larger meshes (millions of horizontal cells) the flag `io_type="pnetcdf"` must be changed to `io_type="pnetcdf,cdf5"`.

#### Checking and modifying the output data

By default, MPAS-Ocean will output a set of monthly-averaged variables. The streams file can be modified to include additional variables in the existing output files, produce additional output files, or change the output frequency (e.g. high frequency files shifting between daily or 5-daily frequencies).

The XML file is organized into blocks describing each stream. Typical streams for output include:

`timeSeriesStatsMonthlyOutput` - monthly averaged output

`highFrequencyOutput` - high frequency snapshots (not averaged) output with frequency `output_interval`

Under each block header is the list of variables (individual variables, variable structure, or variable arrays) that will be output within the relevant stream.

##### Example workflow for modifying the output fields

- copy the `streams.ocean` file from the run directory to the `SourceMods` directory (see above)
- identify the variable name for the variable of interest. You can find the variable name by searching the ocean Registry.xml (in the `src` directory) or Registry_package.xml in the `tracer` and `analysis_member` sub-directories. Note whether the variable of interest is included within a `var_array` or a `var_struct`.
- identify the output stream of interest (e.g. monthly averages, high frequency, others). You can search for a stream name, known output filename, or output interval.
- check whether the variable of interest is included in the `streams.ocean` file. Search for the `var name`, or the `var_array` and `var_struct` if applicable. If it is, copy the variable line from other streams into the stream of interest. If it is not included, copy it from the Registry.xml.
- check whether the relevant stream is turned on. This includes checking that `output_interval` in the stream header is not `None`.
- make further modifications: e.g. you can modify the `output_interval` for the high-frequency stream. If you are turning on a new stream, remove unnecessary variables from the stream.

##### Excerpts from a `streams.ocean` file

```text
<stream name="timeSeriesStatsMonthlyOutput"
        type="output"
        precision="single"
        io_type="pnetcdf"
        useMissingValMask="true"
        filename_template="GMPAS-JRA.TL319_IcoswISC30E3r5.anvil.mpaso.hist.am.timeSeriesStatsMonthly.$Y-$M-$D.nc"
        filename_interval="00-01-00_00:00:00"
        reference_time="01-01-01_00:00:00"
        output_interval="00-01-00_00:00:00"
        clobber_mode="truncate"
        packages="timeSeriesStatsMonthlyAMPKG"
        runtime_format="single_file">

    <var name="daysSinceStartOfSim"/>
    <var name="binBoundaryMerHeatTrans"/>
    <var name="binBoundaryZonalMean"/>
    <var name="ssh"/>
    <var_struct name="tracers"/>
    <var name="velocityMeridional"/>
    <var name="velocityZonal"/>
    <var name="layerThickness"/>
...
</stream>
```

```text
<stream name="highFrequencyOutput"
        type="output"
        precision="single"
        io_type="pnetcdf"
        filename_template="GMPAS-JRA.TL319_IcoswISC30E3r5.anvil.mpaso.highFrequencyOutput.$Y-$M-$D_$h.$m.$s.nc"
        filename_interval="00-01-00_00:00:00"
        reference_time="01-01-01_00:00:00"
        output_interval="00-00-05_00:00:00"
        clobber_mode="truncate"
        packages="highFrequencyOutputAMPKG">

    <var name="xtime"/>
    <var name="daysSinceStartOfSim"/>
    <var_array name="activeTracersAtSurface"/>
    <var_array name="activeTracersAt250m"/>
    <var_array name="activeTracersAtBottom"/>
    <var name="kineticEnergyAtSurface"/>
    <var name="kineticEnergyAt250m"/>
    <var name="relativeVorticityAt250m"/>
    <var name="ssh"/>
    <var name="pressureAdjustedSSH"/>
    <var name="boundaryLayerDepth"/>
    <var name="dThreshMLD"/>
    <var name="tThreshMLD"/>
    <var name="columnIntegratedSpeed"/>
    <var name="barotropicSpeed"/>
    <var name="landIceFreshwaterFlux"/>
    <var name="pressureAdjustedSSH"/>
    <var name="atmosphericPressure"/>
</stream>
```

A more comprehensive description of the streams options can be found in Chapter 6 of the [MPAS-Ocean User's Guide](https://zenodo.org/records/11098080).
