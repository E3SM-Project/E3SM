# the main module for landgen

# year processing can go forward or backward in time, so set the start and end years appropriately in the config file
# config_path is the full path the .json config file, including the file name, e.g. /path/to/config.json

import multiprocessing as mp
import importlib
import json
import sys
from pathlib import Path
from . import shared_data
from .shared_data import GridData, GridManager

def load_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

def main(config_path):
	# todo: need to deal with landfrac data structure
	landfrac = None
	config = load_config(config_path)
	modules = config.get('modules', [])

	# get the common parameters for all modules and store in a shared dictionary
	temp_dict = {
				'start_year': config.get('start_year', 2015),
				'end_year': config.get('end_year', 2015),
				'source_data_path': config.get('source_data_path', ''),
				'landgen_grid_path': config.get('landgen_grid_path', ''),
				'out_path': config.get('out_path', ''),
			}
	manager = mp.Manager()
	com_config_dict = manager.dict(temp_dict)

	# create the shared landgen out grid shared data structure
	grid_manager = GridManager()
	grid_manager.start()
	out_grid_data = grid_manager.GridData()
	out_grid_data.allocate()

	## todo: read in the grid file and set the values in out_grid_data

	for mod in modules:
		name = mod['name']
		params = mod.get('params', {})
		try:
			module = importlib.import_module(f'landgen.{name}')
			if hasattr(module, 'run'):
				print(f"Running module: {name}")
				run_list = [*params.values(), com_config_dict, out_grid_data, manager, grid_manager]
				module.run(*run_list)
			else:
				print(f"Module {name} does not have a 'run' function.")
		except ImportError as e:
			print(f"Could not import module {name}: {e}")

		
	# free the shared memory
	com_config_dict = None
	out_grid_data = None
	manager.shutdown()
	grid_manager.shutdown()
	return
