# Generating Domain Files

Domain files are needed at runtime by the coupler, data models, and land model.
The land model uses the mask to determine where to run and the coupler use the
land fraction to merge fluxes from multiple surface types to the atmosphere
above them.

Domain files are created from a conservative, monotone mapping file from the 
ocean grid (where the mask is defined) to the atmosphere grid.

a mapping files created in the 
previous step, using a tool provided with CIME in 
${e3sm_root}/cime/tools/mapping/gen_domain_files. 


## Environement

The new domain generation tool requires a few special packages, 
such as xarray, numba, and itertools.
These are all included in the E3SM unified environment:
https://e3sm.org/resources/tools/other-tools/e3sm-unified-environment/

Alternatively, a simple conda environment can be created with the following command:
```
conda create --name example_env --channel conda-forge xarray numpy numba scikit-learn netcdf4
```