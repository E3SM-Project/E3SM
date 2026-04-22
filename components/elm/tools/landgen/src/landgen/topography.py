# topography.py
# this module processes topography data for the landgen workflow
# this needs to be run first, as the landfrac data are needed for each of the other modules
# if these data are already available, then this module will simply read in the landfrac data
# the output is the complete set of topography data needed for landgen and elm

# run() function is the main entry point for this module, and will be called by landgen.py

import multiprocessing as mp
import importlib
import sys
from pathlib import Path
from . import shared_data
from .shared_data import TopoData, TopoManager


########## define helper functions for land_type run() here





########## run()

## arguments
## these first ones are module-specific parameters that are set in the config file
# active: true = module is run, false = module is skipped
# out_fname: output filename for the module
# com_config_dict: the shared dictionary for the common parameters for all modules
# out_grid_data: the shared data structure for the landgen grid data

## output

def run(active, out_fname, com_config_dict, out_grid_data, manager, grid_manager):
    if active is False:
        print(f"Topography processing is not active, but reading in landfrac data for other active modules.")
        output_file = Path(com_config_dict['out_path']) / out_fname
        if output_file.exists():
            print(f"Output file {output_file} exists; reading in existing data.")
            ## todo: need to define read_landfrac to return the data in the correct format for the shared data structure
            ## may need to add elevation to GridData class and here
            out_grid_data.landfrac = read_landfrac(output_file)
        else:
            print(f"Error: Output file {output_file} does not exist; set active to True to process topography data.")
            sys.exit(1)

    # set up the topography module shared data structure
    topo_manager = TopoManager()
    topo_manager.start()
    topo_out_data = topo_manager.TopoData()
    topo_out_data.allocate()

    print(f"Processing topography module with parameters:")
    # todo: print the parameters here


    # topography data processing

    # mp pool and parallel processing code will happen in this module or its submoduels


    # free the module-specific shared data structure
    topo_out_data = None
    topo_manager.shutdown()
    return
