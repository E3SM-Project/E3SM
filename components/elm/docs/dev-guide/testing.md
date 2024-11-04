When a developer is ready to issue a PR to integrate their new developments into E3SM's `master` branch, they need to test their code development to ensure the new development does not introduce unexpected bugs. A land developer needs to atleast run the `e3sm_land_developer` test suite that involves following two steps: 

- **Step-1**: Generate the baselines using the hash on `master` that was the starting point of the 
  development branch, and
- **Step-2**: Switch to the development branch and re-run the test suite to compare against 
  the baselines generated in Step-1.

For the purpose of these document, let's assume the developer branch is `bishtgautam/lnd/emi-example`.

## Step-1: Generate Baselines

Find the starting hash on `master` by switching to the developer branch and looking at the Git
graph.

```bash
cd <e3sm-dir>

# Switch to your branch
git checkout bishtgautam/lnd/emi-example

# Let's look at the graph
git log --oneline --decorate --graph
* 9a2cd8459 (HEAD -> bishtgautam/lnd/emi-example, origin/bishtgautam/lnd/emi-example) Updates the Stub EM code
* d11fb8c52 Adds a readme and makefile for EMI Demo code
* 15a0e4b8d Minor update to the EMI demo code
* 63e2f6dcd Minor fix to ELM configure script
* 551fc79cc Updates ELM stub and demo code to use CNCarbonStateType
* 8c9d41f87 Adds code to exchange CNCarbonStateType
* 0000ded51 Minor fix
* 2baa4ebec Adds fixes to cmake
* f77ae2cea Rename few cmake files
* b48b66dd6 Adds a demo for a Stub EM
* fa94aea0d Fixes length of variable to store name of EMI data
* cca0e0be1 Adds capability to print all EMI Data in a list
* beb1fe9ee Adds .gitignore files for EMI
* 301acfe7a Renames ExternalModelInterfaceDataMod.F90
* f0eb5b40b Renames ExternalModelInterfaceDataDimensionMod.F90
* 076100882 Restructures ELM's EMI directory
*   cfb7fc2b2 Merge branch 'qzhu-lbl/lnd/ch4_inundation_bugfix' (PR #2814)
|\
| * ec99b1cb4 (origin/qzhu-lbl/lnd/ch4_inundation_bugfix)  - fix ch4 inundatation parameter read in bug
* |   c59c98116 Merge branch 'darincomeau/mpaso/enable_eddystats' (PR #2821)
|\ \
| * | c21bb46a4 Turn on eddy stats for oEC60to30v3 and oEC60to30v3wLI ocn grids by default
* | |   a2dd2fce5 Merge pull request #2811 from E3SM-Project/jayeshkrishna/pio2_cime_changes
```

The above Git graph tells us that `bishtgautam/lnd/emi-example` started with the `cfb7fc2b2` hash on the `master`. So, we will generate the baselines using `cfb7fc2b2`. Checkout `cfb7fc2b2` and update the submodels via:

```bash
git checkout cfb7fc2b2
Note: checking out 'cfb7fc2b2'.

You are in 'detached HEAD' state. You can look around, make experimental
changes and commit them, and you can discard any commits you make in this
state without impacting any branches by performing another checkout.

If you want to create a new branch to retain commits you create, you may
do so (now or later) by using -b with the checkout command again. Example:

  git checkout -b <new-branch-name>

HEAD is now at cfb7fc2b2... Merge branch 'qzhu-lbl/lnd/ch4_inundation_bugfix' (PR #2814)

git submodule update --init
```

Now we will use `cime/scripts/create_test` to generate baseline for the `e3sm_land_developer`
test suite.

```bash
cd cime/scripts
```

A few things that you need to decide at this stage include the directory location where the baselines will be saved, the name and ID for the baselines, the project allocation that you will use for running the simulation, if you would like to receive email notifications about tests, etc.

```bash
# Define the directory to hold the baseline 
export MY_BASELINE_DIR=/global/cscratch1/sd/gbisht/e3sm_baselines

# Do you want to use a name for the baseline? 
# One choice could be the git hash that is being used to generate the baselines.
export BASELINE_NAME=cfb7fc2b2

# Let's set TEST_ID to be same as BASELINE_NAME
export TEST_ID=${BASELINE_NAME}

# If you are a member of E3SM, you could use 'e3sm' project allocation
export PROJECT=e3sm

# Set your email
export MAIL_USER=<your-mail@something>

```

Use `./create_test --help` to get a complete list of arguments. Below are some additional
useful arguments for `./create_test`:

```bash
# Other arguments
# -v         : Verbose option
# -g         : Generate the baseline
# -q         : If you want to use a particular job queue (e.g. 'debug' queue on NERSC)
# --walltime : Specify the wall time for jobs (e.g. 30 min is max allowable for 'debug' queue on NERSC)
# --mail-user: If you want to receive emails about your jobs
# --mail-type: When to receive emails. Options are: never, all, begin, end, fail.
# -j         : Number of parallel jobs
```

Now run the `e3sm_land_developer`

```bash
./create_test e3sm_land_developer  \
--baseline-root ${MY_BASELINE_DIR} \
-b ${BASELINE_NAME}                \
-t ${TEST_ID}                      \
-q regular                         \
-p ${PROJECT}                      \
--walltime 00:30:00                \
--mail-user $MAIL_USER             \
--mail-type all                    \
-g                                 \
-v                                 \
-j 4
```

The cases would be named `*.G.*` to denote one is generating the baselines.
It will take a while to compile all the cases and submit the code.
It can take a long time (>30mins) for the test suite to run. To avoid interrupting the test suite
in the middle, one can run the test suite within a `screen` command.
After the cases have been successfully compiled and submitted, you can check that 
status of test by running the `cs.status.${TEST_ID}` file that was created in the scratch directory.



## Step-2: Compare against previously generated baselines

Now, switch to the development branch and be sure to update submodules.


```bash
cd <e3sm-dir>
git checkout bishtgautam/lnd/emi-example
Previous HEAD position was cfb7fc2b2... Merge branch 'qzhu-lbl/lnd/ch4_inundation_bugfix' (PR #2814)
Switched to branch 'bishtgautam/lnd/emi-example'
Your branch is up-to-date with 'origin/bishtgautam/lnd/emi-example'.

# Checkout the appropriate submodules
git submodule update --init
```

Again initialize few settings.

```bash
cd cime/scripts

# Let's use the settings as the last time
export MY_BASELINE_DIR=/global/cscratch1/sd/gbisht/e3sm_baselines

# IMPORTANT: One needs to use the same BASELINE_NAME as in Step-1 because one wants
#            to compare against the baselines previously generated in Step-1.
export BASELINE_NAME=cfb7fc2b2

# You can use the hash at the tip of your branch as the ID
export TEST_ID=9a2cd8459
```

Run the test suite and compare (via `-c`) against previously generated baselines.

```bash
./create_test e3sm_land_developer  \
--baseline-root ${MY_BASELINE_DIR} \
-b ${BASELINE_NAME}                \
-t ${TEST_ID}                      \
-q debug                           \
--walltime 00:30:00                \
--mail-user $MAIL_USER             \
--mail-type all                    \
-c                                 \
-v                                 \
-j 4
```

The cases would be named `*.C.*` to denote one is comparing against previously generated baselines. Similar to the last time, a new `cs.status.${TEST_ID}` would be created and
you can check the status of test by running it.
