# soil.py
# this module processes soil data for the landgen workflow
# the output is the complete set of soil data needed for landgen and elm
# currently these output data are static, but some source raw inputs may have timeseries

# run() function is the main entry point for this module, and will be called by landgen.py

import multiprocessing as mp
import importlib
from pathlib import Path


########## define helper functions for land_type run() here





########## run()

## arguments
## these first ones are module-specific parameters that are set in the config file
# active: true = module is run, false = module is skipped
# out_fname: output filename for the module
## these are general parameters for all modules that are set in the config file
# start_year: start year for processing
# end_year: end year for processing
# source_data_path: base path to the source data
# landgen_grid_path: path from source_data_path and the filename of the landgen grid
# out_path: base path for the output data
# landfrac: here this is the output landfrac data

## output

def run(active, out_fname, start_year, end_year, source_data_path, landgen_grid_path, out_path, landfrac):
    if active is False:
        print(f"Skipping soil module")
        return

    print(f"Processing soil module with parameters:")
    # todo: print the parameters here

    # processing code for soil

    ##### todo: need to figure out how to handle the data structure for soil_data and how to pass the landfrac data to the this module
    soil_data = None

    # soil data processing