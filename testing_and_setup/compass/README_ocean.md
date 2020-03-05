# Setting up COMPASS for ocean test case

## COMPASS conda environment

To set up and run ocean test cases from COMPASS, you will need a conda
environment.  First, install Miniconda3 (if miniconda is not already
installed), then create a new conda environment as follows:
``` bash
conda create -n compass_py3.7 -c conda-forge -c xylar python=3.7 \
    geometric_features=0.1.4 mpas_tools=0.0.6 jigsaw=0.9.11 jigsawpy=0.0.2 \
    metis pyflann scikit-image basemap pyamg ffmpeg
```
Each time you want to work with COMPASS, you will need to run:
```
conda activate compass_py3.7
```

An appropriate conda environment is already available on Los Alamos National
Laboratory's Institutional Computing (LANL IC) machines:
```
source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base/etc/profile.d/conda.sh
conda activate compass_py3.7
```

## Setting config options

The file `general.config.ocean` is a template containing a set of config
options that the COMPASS user must set in order to set up ocean test cases.
Make a copy of this file (e.g. `config.ocean`) and set the options as follows.
In six places, replace `FULL_PATH_TO_MPAS_MODEL_REPO` with the path where you
have checked out (and built) the branch of MPAS-Model you are planning to use.
Five other paths are required, as explained below.

### geometric_data

```
geometric_data = FULL_PATH_TO_PATH_OF_CACHED_GEOMETRIC_FEATURES_DATA
```
In general, this can be any directory where geojson data files will be
automatically downloaded and cached as test cases require them.  However, if
working on compute nodes that cannot get to the internet (specifically GitHub),
it may be necessary to download the full set of `geometric_data` ahead of time
from:
[https://github.com/MPAS-Dev/geometric_features/tree/0.1/geometric_data](https://github.com/MPAS-Dev/geometric_features/tree/0.1/geometric_data)
The easiest way to do this is with:
```
git clone git@github.com:MPAS-Dev/geometric_features.git
cd geometric_features
git checkout 0.1
```
Then, point `geometric_data` in your config file to
`geometric_features/geometirc_data/` in the location you have just placed it.

On LANL IC, a full checkout is already available at:
```
geometric_data = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/geometric_data_v0.1
```

### mesh_database, initial_condition_database and bathymetry_database

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
