# human.py
# this module processes human data for the landgen workflow
# the output is the complete set of human data needed for landgen and elm
# population denstiy data is a timeseries
# gdp data are static

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

def run(active, out_fname, com_config_dict, out_grid_data, manager, grid_manager):
    if active is False:
        print(f"Skipping human module")
        return

    # extract common parameters from shared config dict
    start_year        = com_config_dict['start_year']
    end_year          = com_config_dict['end_year']
    source_data_path  = com_config_dict['source_data_path']
    landgen_grid_path = com_config_dict['landgen_grid_path']
    out_path          = com_config_dict['out_path']
    landfrac          = out_grid_data.get_landfrac()

    print(f"Processing human module with parameters:")
    # todo: print the parameters here

    # processing code for human

    ##### todo: need to figure out how to handle the data structure for human_data and how to pass the landfrac data to the this module
    human_data = None

# human data processing