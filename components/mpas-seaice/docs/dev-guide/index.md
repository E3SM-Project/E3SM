# Developer's Guide

Development of the MPAS-seaice component should follow the general procedures outlined by the E3SM project.

[Development Guide for E3SM Code](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1868455/Development+Getting+Started+Guide)
[Development Guide for E3SM Documentation](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/3924787306/Developing+Documentation)

## Configuration Controls

MPAS-seaice is controlled using namelist options.  Namelist files for E3SM runs are found in ``E3SM/components/mpas-seaice/bld/namelist_files/``.  However, the values in these files are drawn from the Registry, following the convention of all MPAS components. Registry files are used directly for stand-alone MPAS-seaice runs, and E3SM scripts pass information from them into E3SM's namelist files when a PR is merged.  E3SM's namelist files need to be changed for development purposes.  It's easiest to change all of them when needed, to keep them consistent, taking care not to unintentionally change standalone MPAS-seaice configurations.

## MPAS Framework

MPAS-seaice is built on the MPAS Framework.

The MPAS Framework provides the foundation for a generalized geophysical fluid dynamics model on unstructured spherical and planar meshes. On top of the framework, implementations specific to the modeling of a particular physical system (e.g., sea ice, ocean) are created as MPAS cores. The MPAS design philosophy is to leverage the efforts of developers from the various MPAS cores to provide common framework functionality with minimal effort, allowing MPAS core developers to focus on development of the physics and features relevant to their application.

The framework code includes shared modules for fundamental model operation. Significant capabilities include:

- **Description of model data types.** MPAS uses a handful of fundamental Fortran derived types for basic model functionality. Core-specific model variables are handled through custom groupings of model fields called pools, for which custom access routines exist. Core-specific variables are defined in XML syntax in a Registry, and the framework parses the Registry, defines variables, and allocates memory as needed.
- **Mesh specification.** MPAS requires 36 fields to fully describe the mesh used in a simulation. These include the position, area, orientation, and connectivity of all cells, edges, and vertices in the mesh. The mesh specification can flexibly describe both spherical and planar meshes. For more information about the meshes, see the [Users Guide](../user-guide/index.md).
- **Distributed memory parallelization and domain decomposition.** The MPAS Framework provides needed routines for exchanging information between processors in a parallel environment using Message Passing Interface (MPI). This includes halo updates, global reductions, and global broadcasts. MPAS also supports decomposing multiple domain blocks on each processor to optimize model performance by minimizing transfer of data from disk to memory. Shared memory parallelization through OpenMP is also supported, but the implementation is left up to each core.
- **Parallel input and output capabilities.** MPAS performs parallel input and output of data from and to disk through the commonly used libraries of NetCDF, Parallel NetCDF (pnetcdf), and Parallel Input/Output (PIO). The Registry definitions control which fields can be input and/or output, and a framework "streams" functionality provides run-time configuration of what fields are to be written to what file name and at what frequency through an XML streams file. The MPAS framework includes additional functionality specific to providing a flexible model restart capability.
- **Advanced timekeeping.** MPAS uses a customized version of the timekeeping functionality of the Earth System Modeling Framework (ESMF), which includes a robust set of time and calendar tools used by many Earth System Models (ESMs). This allows explicit definition of model epochs in terms of years, months, days, hours, minutes, seconds, and fractional seconds and can be set to three different calendar types: Gregorian, Gregorian no leap, and 360 day. This flexibility helps enable multi-scale physics and simplifies coupling to ESMs. To manage the complex date/time types that ensue, MPAS framework provides routines for arithmetic of time intervals and the definition of alarm objects for handling events (e.g., when to write output, when the simulation should end).
- **Run-time configurable control of model options.** Model options are configured through namelist files that use standard Fortran namelist file format, and input/output are configured through streams files that use XML format. Both are completely adjustable at run time.
- **Online, run-time analysis framework.** A system for defining analysis of model states during run time, reducing the need for post-processing and model output.

Additionally, a number of shared operators exist to perform common operations on model data. These include geometric operations (e.g., length, area, and angle operations on the sphere or the plane), interpolation (linear, barycentric, Wachspress, radial basis functions, spline), vector and tensor operations (e.g., cross products, divergence), and vector reconstruction (e.g., interpolating from cell edges to cell centers). Most operators work on both spherical and planar meshes.

## Coupling of MPAS-seaice within E3SM

Several files within MPAS-SI control the coupling of fields between MPAS-SI and other E3SM components via the coupler. A list of all variables being passed to and from MPAS-SI can be found in ``components/mpas-seaice/driver/mpassi_cpl_indices.F``.
Variables named with the prefix ``x2i`` are passed from the coupler to MPAS-SI, variables with the prefix ``i2x`` are passed from MPAS-SI to the coupler.
For example, all coupled variables are named with the following convention:
``index_i2x_Si_ithick`` refers to a state variable of ice (``Si``) that gets passed from MPAS-SI to the coupler (``i2x``) for the ice thickness (``ithick``). ``index_x2i_So_t`` refers to an ocean state variable (``So``) passed from the the coupler to MPAS-SI (``x2i``) for ocean temperature (``t``).

## Icepack

For changes to Icepack, please consult the [CICE Consortium's recommendations for code contributions](https://github.com/CICE-Consortium/About-Us/wiki/Contributing).

To access the column physics in Icepack, MPAS-seaice uses methods defined in ``icepack_intfc.F90``. The 'init' and 'query' methods are used to set and retrieve Icepack values. A 'write' method is also available for documenting these values. MPAS-seaice follows the 'icepack_warnings' methodology where ``icepack_warnings_aborted`` is checked and ``icepack_warnings_flush`` is called after every call to an Icepack method. It does not directly “use” Icepack data, accessing Icepack data only through these interfaces.

### Physical constants

MPAS-Seaice centralizes its constants in ``src/shared/mpas_seaice_constants.F``. Values that represent shared physical references are imported from E3SM's PhysicalConstantsDictionary (PCD) and mapped to SeaIce names there. SeaIce-specific parameters, such as snow density, and values not currently available in PCD remain defined locally; PCD is not intended to replace every model parameter.

The SeaIce constants module passes the applicable values into Icepack through ``icepack_init_parameters`` in ``mpas_seaice_icepack.F``. Icepack does not depend directly on PCD: it receives values through its public initialization interface, and callers that omit optional inputs continue to use Icepack's defaults. For example, MPAS-Seaice passes PCD's latent heat of fusion explicitly; otherwise Icepack derives it as the latent heat of sublimation minus the latent heat of vaporization. Some quantities are converted or derived while mapping them, so follow the units in the receiving code (for example, carbon and iron molar masses are stored in PCD in kg/mol and converted to g/mol for SeaIce).

The iron molar-mass entry is based on NIST's standard atomic weight for naturally occurring iron, 55.845(2), and is recorded in PCD as 0.055845 kg/mol with an absolute uncertainty of 0.000002 kg/mol. SeaIce converts this to g/mol for its biogeochemistry calculations.

PCD's ``seconds_per_day`` defines the fixed day unit as exactly 86400 SI seconds, following [NIST SP 330 (2019), Section 4, Table 8](https://doi.org/10.6028/NIST.SP.330-2019). SeaIce maps this to ``seaiceSecondsPerDay`` and passes it to Icepack through the existing optional ``secday_in`` argument. Icepack retains its 86400-second default for callers that omit this input. This conversion is not a sidereal day; historical UTC days containing leap seconds differ from this duration. [CGPM Resolution 5 (2026)](https://doi.org/10.59161/CGPM2026RES5E) specifies that continuous UTC takes effect on 20 May 2027, ending future leap-second adjustments without changing the exact day-unit conversion or erasing historical leap seconds. The reciprocal remains derived from it, and ``seaiceSecondsPerYear`` retains its existing fixed 365-day convention; calendar handling is unchanged.

Other sea-ice quantities remain at their existing model values rather than being added to PCD. Snow density and thermal properties are parameterized model choices, seawater freezing temperature depends on state (including salinity and pressure), and brine viscosity is state- and formulation-dependent. The von Kármán constant is empirical and can depend on the turbulence formulation. These quantities therefore do not have a single generally applicable reference value to standardize here; this work leaves their existing SeaIce and Icepack values unchanged.

The CIME build compiles PCD through E3SM's shared ``share/util`` sources. The separate MPAS-Seaice Make build compiles the same PCD source as a prerequisite of ``mpas_seaice_constants.F``; it should not add another PCD source to the CIME component source list.

### Spherical Earth radius

SeaIce's ``seaiceEarthRadius`` and the shared ``SHR_CONST_REARTH`` use PCD's
GRS80 ``mean_radius``, approximately 6371008.771366666 m. This is the arithmetic
mean of the reference ellipsoid's two equatorial axes and polar axis, not the
equal-area or equal-volume sphere radius. It matches
[Omega's Earth-radius convention](https://github.com/E3SM-Project/Omega/blob/develop/components/omega/src/ocn/GlobalConstants.h).
Both SeaIce configuration defaults must track this PCD value.

E3SM spherical SeaIce runs reject a ``config_earth_radius`` that differs from
this reference by more than a relative tolerance of 1e-8. The tolerance is
increased to machine epsilon in single-precision MPAS builds, because MPAS
stores configuration and mesh metadata in its working precision. Standalone
idealized SeaIce runs may still override the radius; planar meshes are not
rescaled and do not undergo the spherical-radius check.

On non-restart initialization, SeaIce retains its existing mesh rescaling:
coordinates and edge lengths scale with radius, and areas scale with radius
squared. The mesh-pool and domain radius metadata are updated consistently.
Restart runs retain the exact comparison between the configured radius and the
stored mesh radius. Existing restarts at the legacy 6371229 m radius are therefore
not compatible with the new E3SM spherical convention. Restart migration is not
implemented; do not merely edit the radius attribute.

E3SM MPAS-Ocean forward and analysis runs also reject spherical input meshes
whose radius differs from the shared PCD reference. Standalone ocean builds
retain their existing behavior. Ocean input geometry must be generated or
migrated consistently: coordinates, lengths, areas, and radius metadata must
agree. When enabled, the existing MPAS-Ocean init-mode sphere expansion uses the
PCD radius in E3SM builds; standalone init-mode expansion retains its legacy
radius. Forward and analysis modes do not automatically migrate input geometry.
Omega likewise checks its input mesh against PCD. Changing metadata alone
does not migrate geometry. Ocean and SeaIce datasets and any area-dependent
coupling inputs must be reviewed together.

SeaIce and MPAS-Ocean export cell areas normalized by their mesh radius squared.
Coupler diagnostics convert these areas back to square meters using
``SHR_CONST_REARTH``. Mapping the shared radius to PCD aligns this convention, but
also changes other shared-radius consumers across E3SM. This is a non-bit-for-bit
geometry change, not just a SeaIce diagnostic change.

Before production use, validate migrated ocean/ice input geometry and normalized
coupler areas, then run existing coupled restart, conservation, MPI-layout, and
thread-count tests on a supported machine. The earlier SeaIce restart result
does not validate this radius change.

Basic Icepack development can be done in standalone mode using Icepack's testing scripts, directly in the submodule branch in MPAS-seaice. **We recommend that Icepack developments be thoroughly tested within E3SM's coupled framework throughout the development process, including fully coupled simulations.**

## E3SM-Polar-Developer Script

To accelerate early development stages, a script is available for configuring and testing MPAS-seaice (including the Icepack submodule) in D compsets, which have the sea ice component active and data models for the other components.

### View helpful information, including default values for duration, configuration, etc

```text
git clone git@github.com:E3SM-Project/SimulationScripts.git
cd SimulationScripts/archive/PolarGroup
./E3SM-Polar-Developer.sh -h
```

For debugging E3SM, search the script for 'debug' and follow the instructions.

The following examples describe how to use the script for development in Icepack.  Similar procedures could be used for any MPAS-SI physics development.

### Set up and run baselines

Create a file containing modified namelist options. The file ``nset01.nlk`` in this example turns off the ``snicar_ad`` radiation scheme.

```text
$ less nset01.nlk
[mpassi]
config_use_snicar_ad = {.false.}
```

Notes:

- A .nlk file without any config settings will create a baseline using default settings.

Fetch E3SM (choose any name for the directory baselines01):

```text
./E3SM-Polar-Developer.sh -s baselines01 -f git@github.com:E3SM-Project/E3SM
```

Set up a new case and build it:

```text
./E3SM-Polar-Developer.sh -s baselines01 -k nset01.nlk -e -n -b
```

Submit:

```text
./E3SM-Polar-Developer.sh -s baselines01 -k nset01.nlk -e -q
```

Examine the diagnostic output:

```text
./E3SM-Polar-Developer.sh -s baselines01 -k nset01.nlk -e -a -v
```

### Set up a sandbox for model development, to be compared with the baselines

Fetch E3SM (choose any name for the directory newdev01):

```text
./E3SM-Polar-Developer.sh -s newdev01 -f git@github.com:E3SM-Project/E3SM
```

Create a new development branch:

```text
cd ~/E3SM-Polar/code/newdev01
git branch newbranch
git checkout newbranch
```

Set up a new case and build it:

```text
./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -n -b
```

Develop and test...
Build/compile:

```text
./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -b
```

Submit:

```text
./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -q
```

Examine the diagnostic output:

```text
./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -a -v
```

Compare with the baselines case directory (use your D3 baselines directory):

```text
./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -a D3.nset01.baselines01.master.E3SM-Project.anvil -v
```

### Make changes in Icepack and PR to the Consortium

We recommend PR’ing Icepack changes first to the Consortium then to E3SM’s icepack fork, in order to keep the repositories in sync and to ensure the changes are robust outside of E3SM.  Some changes to Icepack require extensive changes to the driver code (e.g. MPAS-seaice or CICE), making this process challenging.  Contact the [CICE Consortium](https://github.com/CICE-Consortium/About-Us/wiki/Contributing) to discuss and identify a collaborative path forward.

First, create a baseline (standalone) Icepack test suite using the E3SM icepack fork or, if the Consortium code is different, using Consortium icepack main (see [Consortium documentation](https://cice-consortium-icepack.readthedocs.io/en/main/user_guide/ug_testing.html).)

Similarly test your branch of Icepack within E3SM and compare with the baseline.
When satisfied with E3SM testing, PR to Consortium icepack main:

```text
git remote add consortium git@github.com:cice-consortium/icepack.git
git pull consortium main
```

Fix conflicts if needed, then

```text
git add ...
git commit -m "update from cice-consortium main"
```

Continue testing. When satisfied,

```text
git push origin branch
```

Create a PR from branch to cice-consortium/icepack -b main.

Once the PR has been tested and merged into the main Icepack codebase, a new PR is submitted to E3SM.

[An extensive example of this workflow](https://e3sm.atlassian.net/wiki/spaces/pd/pages/3450339435/Project+Workflow) is available.

## CICE-QC Quality Control Testing

Example to run a CICE-QC comparison between two E3SM simulations with changes to the sea ice component.

### Set up and run simulations to be compared

```text
cd ~/SimulationScripts/archive/PolarGroup/
```

Create a `.nlk` file with namelist changes to include the thickness analysis member. Include changes to namelist values needed in both the baseline and the test here, if desired (append the last 3 lines here to the end of your standard D-case test .nlk).

```text
$ less qcbase.nlk
[mpassi]
config_AM_thicknesses_enable = {.true.}
config_AM_thicknesses_compute_on_startup = {.true.}
config_AM_thicknesses_write_on_startup = {.true.}
```

Use test script to clone E3SM, and create a sandbox

```text
./E3SM-Polar-Developer.sh -s qcbaseline -f git@github.com:E3SM-Project/E3SM
```

Edit ``~/E3SM-Polar/code/qcbaseline/components/mpas-seaice/cime_config/buildnml`` to change:

```text
lines.append('        output_interval="none">')
```

to

```text
lines.append('        output_interval="00-00-01_00:00:00">')
```

for ``stream name=“output”`` and add

```text
lines.append('    <var name="iceThicknessCell"/>')
```

a few lines below that:

```text
            lines.append('<stream name="output"')
            lines.append('        type="output"')
            lines.append('        io_type="{}"'.format(ice_pio_typename))
            lines.append('        filename_template="{}.mpassi{}.hist.$Y-$M-$D_$
S.nc"'.format(casename, inst_string))
            lines.append('        filename_interval="00-01-00_00:00:00"')
            lines.append('        clobber_mode="truncate"')
            lines.append('        reference_time="01-01-01_00:00:00"')
-           lines.append('        output_interval="none">')
+           lines.append('        output_interval="00-00-01_00:00:00">')
            lines.append('')
            lines.append('    <stream name="mesh"/>')
            lines.append('    <var name="xtime"/>')
            lines.append('    <var name="daysSinceStartOfSim"/>')
            lines.append('    <var name="iceAreaCell"/>')
            lines.append('    <var name="iceVolumeCell"/>')
            lines.append('    <var name="snowVolumeCell"/>')
+           lines.append('    <var name="iceThicknessCell"/>')
            lines.append('    <var name="uVelocityGeo"/>')
            lines.append('    <var name="vVelocityGeo"/>')
            lines.append('')
            lines.append('</stream>')
```

Build and run baseline case for 5 years (60 months):

```text
./E3SM-Polar-Developer.sh -s qcbaseline -k qcbase.nlk -e -d60 -nb
./E3SM-Polar-Developer.sh -s qcbaseline -k qcbase.nlk -e -d60 -q
```

Copy the thickness analysis member changes into your development directory:

```text
cd ~/E3SM-Polar/code/newdev01/components/mpas-seaice/cime_config/
cp ~/E3SM-Polar/code/qcbaseline/components/mpas-seaice/cime_config/buildnml .
```

If your development case adds namelist parameters, add the thickness analysis member to your .nlk file as above. This example uses the default configuration.

Build and run the development case:

```text
cd ~/SimulationScripts/archive/PolarGroup/
./E3SM-Polar-Developer.sh -s newdev01 -k qcbase.nlk -e -d60 -nb
./E3SM-Polar-Developer.sh -s newdev01 -k qcbase.nlk -e -d60 -q
```

### Run QC comparison

```text
cd ~/E3SM-Polar/code/newdev01/components/mpas-seaice/testing/cice-qc
```

See README.md.  This example is for chrysalis.

Edit ``job_script.cice-qc.chrysalis`` to export (insert your username)

```text
BASE = /lcrc/group/e3sm/[username]/E3SM-Polar/D12.qcbase.emc.qcbaseline.master.E3SM-Project.chrysalis/run.k000/
TEST = /lcrc/group/e3sm/[username]/E3SM-Polar/D12.qcbase.emc.newdev01.branch.E3SM-Project.chrysalis/run.k000
```

Submit QC test. Test results will be in the file ``qc_log.txt``.

```text
sbatch job_script.qc-testing-mpassi.chrysalis
less qc_log.txt
```

Example of desired result:

```text
Running QC test on the following directories:
  /lcrc/group/e3sm/ac.eclare/E3SM-Polar/D12.qcbase.emc.qcbaseline.master.E3SM-Project.chrysalis/run.k000/
  /lcrc/group/e3sm/ac.eclare/E3SM-Polar/D12.qcbase.emc.newdev01.branch.E3SM-Project.chrysalis/run.k000
Number of files: 61
2 Stage Test Passed
Quadratic Skill Test Passed for Northern Hemisphere
Quadratic Skill Test Passed for Southern Hemisphere
```

### Generate statistics from the CICE-QC runs

This only works if the .nlk filename is the same for both cases.  If comparing only namelist changes within MPAS-seaice, use the ``./E3SM-Polar-Developer.sh`` script with a single .nlk file that includes each option.

```text
cd ~/SimulationScripts/archive/PolarGroup/
$ ./E3SM-Polar-Developer.sh -s qcbaseline -k qcbase.nlk -e -d60 -a D12.qcbase.emc.newdev01.branch.E3SM-Project.chrysalis -v
```

### Create comparison plots

To generate MPAS-Analysis plots from the CICE-QC runs and compare:

```text
cd ~/SimulationScripts/archive/PolarGroup/QC_MPAS-SI
```

Edit each script for your run names, directories, etc (search for 'echmod' to find settings used for a QC comparison).

Edit and submit (on chrysalis) the job script 3 times, once for the baseline run, once for the test run, and finally to generate comparison plots.

```text
sbatch job_script
```

Browse the html output by navigating to the location indicated by ``htmlSubdirectory`` in the comparison script.
