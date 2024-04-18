Developer's Guide
=================

Development of the MPAS-seaice component should follow the general procedures outlined by the E3SM project.

[Development Guide for E3SM Code](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1868455/Development+Getting+Started+Guide)
[Development Guide for E3SM Documentation](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/3924787306/Developing+Documentation)

**Configuration Controls**
--------------------------

MPAS-seaice is controlled using namelist options.  Namelist files for E3SM runs are found in ``E3SM/components/mpas-seaice/bld/namelist_files/``.  However, the values in these files are drawn from the Registry, following the convention of all MPAS components. Registry files are used directly for stand-alone MPAS-seaice runs, and E3SM scripts pass information from them into E3SM's namelist files when a PR is merged.  E3SM's namelist files need to be changed for development purposes.  It's easiest to change all of them when needed, to keep them consistent, taking care not to unintentionally change standalone MPAS-seaice configurations.

**MPAS Framework**
------------------

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

**Icepack**
-----------

For changes to Icepack, please consult the [CICE Consortium's recommendations for code contributions](https://github.com/CICE-Consortium/About-Us/wiki/Contributing).

To access the column physics in Icepack, MPAS-seaice uses methods defined in ``icepack_intfc.F90``. The 'init' and 'query' methods are used to set and retrieve Icepack values. A 'write' method is also available for documenting these values. MPAS-seaice follows the 'icepack_warnings' methodology where ``icepack_warnings_aborted`` is checked and ``icepack_warnings_flush`` is called after every call to an Icepack method. It does not directly “use” Icepack data, accessing Icepack data only through these interfaces.

Basic Icepack development can be done in standalone mode using Icepack's testing scripts, directly in the submodule branch in MPAS-seaice. **We recommend that Icepack developments be thoroughly tested within E3SM's coupled framework throughout the development process, including fully coupled simulations.**

**E3SM-Polar-Developer Script**
-------------------------------

To accelerate early development stages, a script is available for configuring and testing MPAS-seaice (including the Icepack submodule) in D compsets, which have the sea ice component active and data models for the other components.

### View helpful information, including default values for duration, configuration, etc.

```text
git clone git@github.com:E3SM-Project/SimulationScripts.git
cd SimulationScripts/archive/PolarGroup
./E3SM-Polar-Developer.sh -h
```

For debugging E3SM, search the script for 'debug' and follow the instructions.

The following examples describe how to use the script for development in Icepack.  Similar procedures could be used for any MPAS-SI physics development.

### Set up and run baselines

Create a file containing modified namelist options. The file ``nset01.nlk`` in this example creates baselines for two types of column physics and turns off the ``snicar_ad`` radiation scheme.

```text
$ less nset01.nlk
[mpassi]
config_column_physics_type = {'column_package','icepack'}
config_use_snicar_ad = {.false.}
```

Notes:

- A .nlk file without any config settings will create a baseline using default settings.
- The ``column_package`` option is still available but is no longer being supported in MPAS-seaice.

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

Examine the diagnostic output (compares the icepack run with the column_package run in this example):

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

More extensive documentation of this workflow tool used for the Icepack merge project is available [here](https://acme-climate.atlassian.net/wiki/spaces/ICE/pages/3450339435/Project+Workflow).

**CICE-QC Quality Control Testing**
-----------------------------------

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

See README.md.  This example is for anvil.

Edit ``job_script.cice-qc.anvil`` to export (insert your username)

```text
BASE = /lcrc/group/e3sm/[username]/E3SM-Polar/D12.qcbase.emc.qcbaseline.master.E3SM-Project.anvil/run.k000/
TEST = /lcrc/group/e3sm/[username]/E3SM-Polar/D12.qcbase.emc.newdev01.branch.E3SM-Project.anvil/run.k000
```

Submit QC test. Test results will be in the file ``qc_log.txt``.

```text
sbatch job_script.qc-testing-mpassi.anvil
less qc_log.txt
```

Example of desired result:

```text
Running QC test on the following directories:
  /lcrc/group/e3sm/ac.eclare/E3SM-Polar/D12.qcbase.emc.qcbaseline.master.E3SM-Project.anvil/run.k000/
  /lcrc/group/e3sm/ac.eclare/E3SM-Polar/D12.qcbase.emc.newdev01.branch.E3SM-Project.anvil/run.k000
Number of files: 61
2 Stage Test Passed
Quadratic Skill Test Passed for Northern Hemisphere
Quadratic Skill Test Passed for Southern Hemisphere
```

### Generate statistics from the CICE-QC runs

This only works if the .nlk filename is the same for both cases.  If comparing only namelist changes within MPAS-seaice, use the ``./E3SM-Polar-Developer.sh`` script with a single .nlk file that includes each option.

```text
cd ~/SimulationScripts/archive/PolarGroup/
$ ./E3SM-Polar-Developer.sh -s qcbaseline -k qcbase.nlk -e -d60 -a D12.qcbase.emc.newdev01.branch.E3SM-Project.anvil -v
```

### Create comparison plots

To generate MPAS-Analysis plots from the CICE-QC runs and compare:

Copy the scripts in the file above to anvil or chrysalis - PROVIDE FILE

Edit each script for your run names, directories, etc (search for 'echmod' to find settings used for the qcPR19 comparison above)

Edit and submit (on chrysalis) the job script 3 times, once for icepack, once for column, and finally for the comparison.

Browse the html output, e.g. navigate to
``https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.eclare/icepack-testing/D12.qcPR19.emc.qcPR19.snicar_active.eclare108213.anvil/mpas_analysis_output/``
