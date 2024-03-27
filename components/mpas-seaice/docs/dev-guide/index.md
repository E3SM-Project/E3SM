Development of the MPAS-seaice component should follow the general procedures outlined by the E3SM project.  Please refer to LINK for instructions. For changes to Icepack, please consult the [CICE Consortium's recommendations for code contributions](https://github.com/CICE-Consortium/About-Us/wiki/Contributing).

To access the column physics in Icepack, MPAS-seaice uses methods defined in ``icepack_intfc.F90``. The 'init' and 'query' methods are used to set and retrieve Icepack values. A 'write' method is also available for documenting these values. MPAS-seaice follows the 'icepack_warnings' methodology where ``icepack_warnings_aborted`` is checked and ``icepack_warnings_flush`` is called after every call to an Icepack method. It does not directly “use” Icepack data, accessing Icepack data only through these interfaces.

Basic Icepack development can be done in standalone mode using Icepack's testing scripts, directly in the submodule branch in MPAS-seaice.  We recommend that Icepack developments be thoroughly tested within E3SM's coupled framework throughout the development process, including fully coupled simulations.

**E3SM-Polar-Developer Script**
-----------------------------------

To accelerate early development stages, a script is available for configuring and testing MPAS-seaice (including the Icepack submodule) in D compsets, which have the sea ice component active and data models for the other components.

**View helpful information, including default values for duration, configuration, etc.**

    git clone git@github.com:E3SM-Project/SimulationScripts.git
    cd SimulationScripts/archive/PolarGroup
    ./E3SM-Polar-Developer.sh -h

For debugging E3SM, search the script for 'debug' and follow the instructions.

The following examples describe how to use the script for development in Icepack.  Similar procedures could be used for any MPAS-SI physics development.

**Set up and run baselines**

Create a file containing modified namelist options. The file ``nset01.nlk`` in this example creates baselines for two types of column physics and turns off the ``snicar_ad`` radiation scheme. (A file without ``config`` settings will create a single, default baseline.)

    $ less nset01.nlk
    [mpassi]
    config_column_physics_type = {'column_package','icepack'}
    config_use_snicar_ad = {.false.}

Fetch E3SM (choose any name for the directory baselines01):

    ./E3SM-Polar-Developer.sh -s baselines01 -f git@github.com:E3SM-Project/E3SM

Set up a new case and build it:

    ./E3SM-Polar-Developer.sh -s baselines01 -k nset01.nlk -e -n -b

Submit:

    ./E3SM-Polar-Developer.sh -s baselines01 -k nset01.nlk -e -q

Examine the diagnostic output (compares the icepack run with the column_package run in this example):

    ./E3SM-Polar-Developer.sh -s baselines01 -k nset01.nlk -e -a -v


**Set up a sandbox for model development, to be compared with the baselines**

Fetch E3SM (choose any name for the directory newdev01):

    ./E3SM-Polar-Developer.sh -s newdev01 -f git@github.com:E3SM-Project/E3SM

Create a new development branch:

    cd ~/E3SM-Polar/code/newdev01
    git branch newbranch
    git checkout newbranch

Set up a new case and build it:

    ./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -n -b

Develop and test...     
Build/compile:

    ./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -b

Submit:

    ./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -q

Examine the diagnostic output:

    ./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -e -a -v

Compare with the baselines case directory (use your D3 baselines directory):

    ./E3SM-Polar-Developer.sh -s newdev01 -k nset01.nlk -a D3.nset01.baselines01.master.E3SM-Project.anvil -v


**Make changes in Icepack and PR to the Consortium**

We recommend PR’ing Icepack changes first to the Consortium then to E3SM’s icepack fork, in order to keep the repositories in sync and to ensure the changes are robust outside of E3SM.  Some changes to Icepack require extensive changes to the driver code (e.g. MPAS-seaice or CICE), making this process challenging.  Contact the [CICE Consortium](https://github.com/CICE-Consortium/About-Us/wiki/Contributing) to discuss and identify a collaborative path forward.

First, create a baseline (standalone) Icepack test suite using the E3SM icepack fork or, if the Consortium code is different, using Consortium icepack main
(see [Consortium documentation](https://cice-consortium-icepack.readthedocs.io/en/main/user_guide/ug_testing.html).)

Similarly test your branch of Icepack within E3SM and compare with the baseline. 
When satisfied with E3SM testing, PR to Consortium icepack main:

    git remote add cice git@github.com:cice-consortium/icepack.git
    git pull cice main

Fix conflicts if needed, then

    git add ...
    git commit -m "update from cice-consortium main"

Continue testing. When satisfied,

    git push origin branch

Create a PR from branch to cice-consortium/icepack -b main.

Once the PR has been tested and merged into the main Icepack codebase, a new PR is submitted to E3SM.

More extensive documentation of this workflow tool used for the Icepack merge project is available [here](https://acme-climate.atlassian.net/wiki/spaces/ICE/pages/3450339435/Project+Workflow).

**CICE-QC Quality Control Testing**
-----------------------------------

Example to run a CICE-QC comparison between icepack and column_package on a branch

**Set up and run simulations to be compared**

    cd ~/SimulationScripts/archive/PolarGroup/

Create a `.nlk` file with namelist changes to include the thickness analysis member (append the last 3 lines here to the end of your standard D-case test .nlk)

    $ less qcPR19.nlk
    [mpassi]
    config_column_physics_type = {'column_package','icepack'}
    config_AM_thicknesses_enable = {.true.}
    config_AM_thicknesses_compute_on_startup = {.true.}
    config_AM_thicknesses_write_on_startup = {.true.}

Use test script to clone E3SM, and create a sandbox

    ./E3SM-Polar-Developer.sh -s qcPR19 -f git@github.com:eclare108213/E3SM snicar_active

Edit ``~/E3SM-Polar/code/qcPR19/components/mpas-seaice/cime_config/buildnml`` to change:

    lines.append('        output_interval="none"

to 

    lines.append('        output_interval="00-00-01_00:00:00">')

for ``stream name=“output”`` (line 451) and add

    lines.append('    <var name="iceThicknessCell"/>')

at line 458.

Create and build both cases to run 5 years, then submit:

    ./E3SM-Polar-Developer.sh -s qcPR19 -k qcPR19.nlk -e -d60 -nb
    ./E3SM-Polar-Developer.sh -s qcPR19 -k qcPR19.nlk -e -d60 -q

**Run QC comparison**

----UPDATE THIS SINCE THE PR HAS NOW BEEN MERGED ------

See ``README.md`` at [https://github.com/E3SM-Seaice-Discussion/E3SM/pull/8/files](https://github.com/E3SM-Seaice-Discussion/E3SM/pull/8/files).

    mkdir CICE-QC         --- if it doesn’t already exist
    cd CICE-QC
    git clone git@github.com:darincomeau/E3SM.git -b add-mpas-seaice-qc-testing
    cd E3SM/components/mpas-seaice/testing_and_setup/qc_testing/

[the following might not be necessary with the latest E3SM-Polar-Developer script]

Edit ``job_script.qc-testing-mpassi.anvil``, e.g.

    export BASE=/lcrc/group/acme/ac.eclare/E3SM-Polar/D12.qcPR19.emc.qcPR19.snicar_active.eclare108213.anvil/run.k000
    export TEST=/lcrc/group/acme/ac.eclare/E3SM-Polar/D12.qcPR19.emc.qcPR19.snicar_active.eclare108213.anvil/run.k001

Edit ``mpas-seaice.t-test.py`` to comment out the lines that insist on the path ending with ``run/`` but still define

    path_a = base_dir
    path_b = test_dir

Submit QC test:

    sbatch job_script.qc-testing-mpassi.anvil

Test results will be in the file ``qc_log.txt``.

**Create comparison plots**

To generate MPAS-Analysis plots from the CICE-QC runs and compare:

Copy the scripts in the file above to anvil or chrysalis - PROVIDE FILE

Edit each script for your run names, directories, etc (search for 'echmod' to find settings used for the qcPR19 comparison above)

Edit and submit (on chrysalis) the job script 3 times, once for icepack, once for column, and finally for the comparison.

Browse the html output, e.g. navigate to

    https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.eclare/icepack-testing/D12.qcPR19.emc.qcPR19.snicar_active.eclare108213.anvil/mpas_analysis_output/


