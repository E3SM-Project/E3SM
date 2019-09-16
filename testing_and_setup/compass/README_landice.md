# Setting up COMPASS for landice test cases

## COMPASS conda environment

To set up and run landice test cases from COMPASS, you will need a conda
environment.  First, install Miniconda3 (if miniconda is not already
installed), then create a new conda environment as follows:
``` bash
conda create -n compass_py2.7 -c conda-forge python=2.7 geometric_features mpas_tools jigsaw metis pyflann scikit-image basemap pyamg pyqt
```
Each time you want to work with COMPASS, you will need to run:
```
conda activate compass_py2.7
```

Note: As of June 2019 Python 2.7 is still the supported version of python, but migration to Python 3.7 is underway.

## Setting config options

The file `general.config.landice` is a template containing a set of config
options that the COMPASS user must set in order to set up landice test cases.
Make a copy of this file (e.g. `general.config.landice.myEdits`) and set the options as follows.
* Set the appropriate paths in the `[namelists]` and `[streams]` sections
* Set the appropriate paths in the `[executables]` section.  The `model` path should be to the compiled `landice_model` file.  The other executables can currently be found in the `MPAS-Tools` repository at https://github.com/MPAS-Dev/MPAS-Tools.  Eventually they will be moved into the conda package, at which point they will no longer be required here.
* The [paths] section is currently not required.

## Setting up and running a specific test

* To see a list of available tests: `./list_testcases.py  -o landice`
* To set up a single test case, e.g.: `./setup_testcase.py -n 102 --work_dir ~/Documents/mpas-git/TESTS/ -f general.config.landice.myEdits`
  In this example test number 102 will be set up within the `~/Documents/mpas-git/TESTS/` directory using the information you have set up in the `general.config.landice.myEdits` file.  Note that within the `--work_dir` directory, a directory tree corresponding to the test will be created (e.g. `dome/2000m/halfar_analytic_test/`), even if a only a single test was requested.  See `/setup_testcase.py -h` for detailed usage options.
* Change directory to where the test was created.  Run the test with python driver script located there.  Note that individual steps can be run manually by cd'ing into the step directory and running the auto-generated python script for that step.

## Running a regression suite
Regression suites are defined that run a set of tests.  Regression suite are defined in the directory `landice/regression_suites/`.
* See `./manage_regression_suite.py -h` for detailed usage options
* A common usage is like: `./manage_regression_suite.py -t landice/regression_suites/combined_integration_test_suite.xml -f general.config.landice.myEdits -c --work_dir ~/Documents/mpas-git/TESTS/TRUSTED`
  This will run the regression suite without comparison to an older version of the code.
* To run the regression suite and compare to a previous run of the regression suite, the command can be called with the `-b` baseline option like: `./manage_regression_suite.py -t landice/regression_suites/combined_integration_test_suite.xml -f general.config.landice.myEdits -c -b ~/Documents/mpas-git/TESTS/TRUSTED --work_dir ~/Documents/mpas-git/TESTS/TESTING`

## Adding a new test
See `testing_and_setup/compass/README` and `testing_and_setup/compass/docs/*` for detailed information about the COMPASS system.
In general, COMPASS is organized like: `testing_and_setup/compass/<core_dir>/<configuration_dir>/<resolution_dir>/<test_dir>`
* `<core_dir>` is `landice`.
* `<configuration_dir>` is the parent directory for a test case configuration, such as the Halfar dome or Antarctica.  This directory level holds general purpose files for this configuration, such as initial condition setup scripts, pre/post-processing analysis scripts, namelist/streams template files, `albany_input.yaml` files, etc.
* `<resolution_dir>` is the grid resolution for a given configuration.  Note that in some cases this directory level is used for other purposes than resolution.
* `<test_dir>` is the specific type of test for a given configuration and resolution, e.g., a smoke test, a restart test, an analytic-solution test.  Note that in some cases this directory level is used for other purposes than a specific type of test.

To add a new test:
1. Create a new test `<configuration_dir>` directory and corresponding `<resolution_dir>` and `<test_dir>` levels in the `compass/landice` directory.
2. Add a namelist/streams template file, `albany_input.yaml` file, and any required pre/post-processing scripts (e.g. initial conditions, comparison to analytic solution) to the `<configuration_dir>` directory.
3. Add `config_*_step.xml` files in the `<test_dir>` for each step of the test required.  Each step is run in a separate directory.  A common layout is to create the mesh in one direcory with a `config_setup_mesh_step.xml` file and for the test to be run in a second directory with a `config_run_model_step.xml` file.
4. Add a `config_driver.xml` file to the `<test_dir>`.  This xml file references the individual step files.
5. Try out and debug your new test using the `setup_testcase.py` command described above.
6. A good example to use for setting up a new test is `testing_and_setup/compass/landice/dome/2000m/halfar_analytic_test`.
