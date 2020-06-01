# Setting up COMPASS for ocean test case

## COMPASS conda environment

To set up and run ocean test cases from COMPASS, you will need a conda
environment.  First, install Miniconda3 (if miniconda is not already
installed), then create a new conda environment as follows:
``` bash
conda create -n compass_0.1.6 -c conda-forge -c e3sm python=3.7 compass=0.1.6
```
Each time you want to work with COMPASS, you will need to run:
```
conda activate compass_0.1.6
```

An appropriate conda environment is already available on Los Alamos National
Laboratory's Institutional Computing (LANL IC) machines as well as Anvil, Compy
and Cori.  In each case, you will run:
```
source <base_path>/load_latest_compass.sh
```
Values of `<base_path>` are:
* grizzly and badger - `/usr/projects/climate/SHARED_CLIMATE/anaconda_envs`
* anvil (blues) - `/lcrc/soft/climate/e3sm-unified/`
* compy - `/share/apps/E3SM/conda_envs`
* cori - `/global/cfs/cdirs/acme/software/anaconda_envs`

## Setting config options

The file `general.config.ocean` is a template containing a set of config
options that the COMPASS user must set in order to set up ocean test cases.
Make a copy of this file (e.g. `config.ocean`) and set the options as follows.
In six places, replace `FULL_PATH_TO_MPAS_MODEL_REPO` with the path where you
have checked out (and built) the branch of MPAS-Model you are planning to use.
Five other paths are required, as explained below.

### mesh\_database, initial\_condition\_database and bathymetry\_database

These are directories for storing pre-generated mesh files, data sets for
creating initial conditions, and bathymetry data. These can be empty directories, in which case
meshes and other data sets will be downloaded as required during test-case
setup.  (If a test case appears to hang during setup, it is most likely
downloading mesh, initial-condition or bathymetry data.)

On LANL IC, the shared data bases can be found at:
```
mesh_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/mesh_database
initial_condition_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/initial_condition_database
bathymetry_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/bathymetry_database
```
