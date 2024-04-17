Development of the MPAS-seaice component should follow the general procedures outlined by the E3SM project.

[Development Guide for E3SM Code](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1868455/Development+Getting+Started+Guide)    
[Development Guide for E3SM Documentation](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/3924787306/Developing+Documentation)

**Icepack**
-----------

For changes to Icepack, please consult the [CICE Consortium's recommendations for code contributions](https://github.com/CICE-Consortium/About-Us/wiki/Contributing).

To access the column physics in Icepack, MPAS-seaice uses methods defined in ``icepack_intfc.F90``. The 'init' and 'query' methods are used to set and retrieve Icepack values. A 'write' method is also available for documenting these values. MPAS-seaice follows the 'icepack_warnings' methodology where ``icepack_warnings_aborted`` is checked and ``icepack_warnings_flush`` is called after every call to an Icepack method. It does not directly “use” Icepack data, accessing Icepack data only through these interfaces.

Basic Icepack development can be done in standalone mode using Icepack's testing scripts, directly in the submodule branch in MPAS-seaice. **We recommend that Icepack developments be thoroughly tested within E3SM's coupled framework throughout the development process, including fully coupled simulations.**

**E3SM-Polar-Developer Script**
-----------------------------------

To accelerate early development stages, a script is available for configuring and testing MPAS-seaice (including the Icepack submodule) in D compsets, which have the sea ice component active and data models for the other components.

**View helpful information, including default values for duration, configuration, etc.**

```
git clone git@github.com:E3SM-Project/SimulationScripts.git
cd SimulationScripts/archive/PolarGroup
./E3SM-Polar-Developer.sh -h
```

For debugging E3SM, search the script for 'debug' and follow the instructions.

The following examples describe how to use the script for development in Icepack.  Similar procedures could be used for any MPAS-SI physics development.

**Set up and run baselines**

Create a file containing modified namelist options. The file ``nset01.nlk`` in this example creates baselines for two types of column physics and turns off the ``snicar_ad`` radiation scheme.
    
```
$ less nset01.nlk
[mpassi]
config_column_physics_type = {'column_package','icepack'}
config_use_snicar_ad = {.false.}
```

Notes:

 - A .nlk file without any config settings will create a baseline using default settings.
 - The ``column_package`` option is still available but is no longer being supported in MPAS-seaice.

Fetch E3SM (choose any name for the directory baselines01):

```
./E3SM-Polar-Developer.sh -s baselines01 -f git@github.com:E3SM-Project/E3SM
```

Set up a new case and build it:

```
./E3SM-Polar-Developer.sh -s baselines01 -k nset01.nlk -e -n -b
```

Submit:

```
./E3SM-Polar-Developer.sh -s baselines01 -k nset01.nlk -e -q
```

Examine the diagnostic output (compares the icepack run with the column_package run in this example):

```
./E3SM-Polar-Developer.sh -s baselines01 -k nset01.nlk -e -a -v
```


**Set up a sandbox for model development, to be compared with the baselines**

Fetch E3SM (choose any name for the directory newdev01):

```
./E3SM-Polar-Developer.sh -s newdev01 -f git@github.com:E3SM-Project/E3SM
```

Create a new development branch:

```
cd ~/E3SM-Polar/code/newdev01
git branch newbranch
git checkout newbranch
```

Set up a new case and build it:


```
./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -n -b
```

Develop and test...     
Build/compile:

```
./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -b
```

Submit:

```
./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -q
```

Examine the diagnostic output:

```
./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -a -v
```

Compare with the baselines case directory (use your D3 baselines directory):

```
./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -a D3.nset01.baselines01.master.E3SM-Project.anvil -v
```

**Make changes in Icepack and PR to the Consortium**

We recommend PR’ing Icepack changes first to the Consortium then to E3SM’s icepack fork, in order to keep the repositories in sync and to ensure the changes are robust outside of E3SM.  Some changes to Icepack require extensive changes to the driver code (e.g. MPAS-seaice or CICE), making this process challenging.  Contact the [CICE Consortium](https://github.com/CICE-Consortium/About-Us/wiki/Contributing) to discuss and identify a collaborative path forward.

First, create a baseline (standalone) Icepack test suite using the E3SM icepack fork or, if the Consortium code is different, using Consortium icepack main
(see [Consortium documentation](https://cice-consortium-icepack.readthedocs.io/en/main/user_guide/ug_testing.html).)

Similarly test your branch of Icepack within E3SM and compare with the baseline. 
When satisfied with E3SM testing, PR to Consortium icepack main:

```
git remote add consortium git@github.com:cice-consortium/icepack.git
git pull consortium main
```

Fix conflicts if needed, then

```
git add ...
git commit -m "update from cice-consortium main"
```

Continue testing. When satisfied,

```
git push origin branch
```

Create a PR from branch to cice-consortium/icepack -b main.

Once the PR has been tested and merged into the main Icepack codebase, a new PR is submitted to E3SM.

More extensive documentation of this workflow tool used for the Icepack merge project is available [here](https://acme-climate.atlassian.net/wiki/spaces/ICE/pages/3450339435/Project+Workflow).

**CICE-QC Quality Control Testing**
-----------------------------------

Example to run a CICE-QC comparison between two E3SM simulations with changes to the sea ice component.

**Set up and run simulations to be compared**

```
cd ~/SimulationScripts/archive/PolarGroup/
```

Create a `.nlk` file with namelist changes to include the thickness analysis member. Include changes to namelist values needed in both the baseline and the test here, if desired (append the last 3 lines here to the end of your standard D-case test .nlk).

```
$ less qcbase.nlk
[mpassi]
config_AM_thicknesses_enable = {.true.}
config_AM_thicknesses_compute_on_startup = {.true.}
config_AM_thicknesses_write_on_startup = {.true.}
```

Use test script to clone E3SM, and create a sandbox

```
./E3SM-Polar-Developer.sh -s qcbaseline -f git@github.com:E3SM-Project/E3SM
```

Edit ``~/E3SM-Polar/code/qcbaseline/components/mpas-seaice/cime_config/buildnml`` to change:

```
lines.append('        output_interval="none">')
```

to 

```
lines.append('        output_interval="00-00-01_00:00:00">')
```

for ``stream name=“output”`` and add

```
lines.append('    <var name="iceThicknessCell"/>')
```

a few lines below that:
```
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
```
./E3SM-Polar-Developer.sh -s qcbaseline -k qcbase.nlk -e -d60 -nb
./E3SM-Polar-Developer.sh -s qcbaseline -k qcbase.nlk -e -d60 -q
```

Copy the thickness analysis member changes into your development directory:

```
cd ~/E3SM-Polar/code/newdev01/components/mpas-seaice/cime_config/
cp ~/E3SM-Polar/code/qcbaseline/components/mpas-seaice/cime_config/buildnml .
```

If your development case adds namelist parameters, add the thickness analysis member to your .nlk file as above. This example uses the default configuration.

Build and run the development case:

```
cd ~/SimulationScripts/archive/PolarGroup/
./E3SM-Polar-Developer.sh -s newdev01 -k qcbase.nlk -e -d60 -nb
./E3SM-Polar-Developer.sh -s newdev01 -k qcbase.nlk -e -d60 -q
```

**Run QC comparison**

```
cd ~/E3SM-Polar/code/newdev01/components/mpas-seaice/testing/cice-qc
```

See README.md.  This example is for anvil.


Edit ``job_script.cice-qc.anvil`` to export (insert your username)

```
BASE = /lcrc/group/e3sm/[username]/E3SM-Polar/D12.qcbase.emc.qcbaseline.master.E3SM-Project.anvil/run.k000/
TEST = /lcrc/group/e3sm/[username]/E3SM-Polar/D12.qcbase.emc.newdev01.branch.E3SM-Project.anvil/run.k000
```

Submit QC test. Test results will be in the file ``qc_log.txt``.

```
sbatch job_script.qc-testing-mpassi.anvil
less qc_log.txt
```

Example of desired result:

```
Running QC test on the following directories:
  /lcrc/group/e3sm/ac.eclare/E3SM-Polar/D12.qcbase.emc.qcbaseline.master.E3SM-Project.anvil/run.k000/
  /lcrc/group/e3sm/ac.eclare/E3SM-Polar/D12.qcbase.emc.newdev01.branch.E3SM-Project.anvil/run.k000
Number of files: 61
2 Stage Test Passed
Quadratic Skill Test Passed for Northern Hemisphere
Quadratic Skill Test Passed for Southern Hemisphere
```

**Generate statistics from the CICE-QC runs**

This only works if the .nlk filename is the same for both cases.  If comparing only namelist changes within MPAS-seaice, use the ``./E3SM-Polar-Developer.sh`` script with a single .nlk file that includes each option.

```
cd ~/SimulationScripts/archive/PolarGroup/
$ ./E3SM-Polar-Developer.sh -s qcbaseline -k qcbase.nlk -e -d60 -a D12.qcbase.emc.newdev01.branch.E3SM-Project.anvil -v
```

**Create comparison plots**

To generate MPAS-Analysis plots from the CICE-QC runs and compare:

Copy the scripts in the file above to anvil or chrysalis - PROVIDE FILE

Edit each script for your run names, directories, etc (search for 'echmod' to find settings used for the qcPR19 comparison above)

Edit and submit (on chrysalis) the job script 3 times, once for icepack, once for column, and finally for the comparison.

Browse the html output, e.g. navigate to

    https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.eclare/icepack-testing/D12.qcPR19.emc.qcPR19.snicar_active.eclare108213.anvil/mpas_analysis_output/


