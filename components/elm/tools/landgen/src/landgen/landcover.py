# landcover.py
# this module processes land cover data for the landgen workflow
# the output is a complete land cover distribution
#    and includes some data associated with particular land covers

# run() function is the main entry point for this module, and will be called by process_single_year in land_type.py

import multiprocessing as mp
#import importlib
from pathlib import Path
import shared_data
import landcover_remote_sensing # not created yet
import transitions # not created yet
import normalize_cell # not created yet
import os

########## define helper functions for landcover run() here

##### landcover_process()

## arguments
# lc_data: land cover data structure that is passed between modules
# year: the year for which to process the land cover data
# prev_year: the previous year for which land cover data were processed
# source_data_path: base path to the source data
# landgen_grid_path: path from source_data_path and the filename of the landgen grid
# lc_rs_path: path from source_data_path and the filename of the land cover remote sensing data (if _use_lc_rs is True)
# out_path: base path for the output data; this is needed to read in the previous year's land cover data (if _use_lc_rs is False)
# prev_out_fname: the output filename for the previous year's land cover data (if _use_lc_rs is False)

## output

def landcover_process(lt_year_data, year, prev_year, prev_fname, lc_rs_path, lc_rs_name,
                            com_config_dict, out_grid_data, ll_limits, cell_ids,
                            man_lock, grid_lock, lt_lock):

    print(f"Processing landcover module year {year} with parameters:")
    # todo: print the parameters here

    # todo: use the with man_lock:, with grid_lock:, with lt_lock: syntax for accessing each managed data structure

    # todo: probably need to add lai path to land_type params and pass it through to here

    lc_rs_data = None
    prev_lt_data = None
    climate_data = None
    lai_data = None
    elm_data = None

    # todo: these data below will be read in based on ll_limits

    if lc_rs.use_lc_rs(year):
        # read modis cover data
        # reading both the igbp cover data and the veg continuous fields data
        lc_rs_data = lc_rs.read(year, com_config_dict['source_data_path'], lc_rs_path, lc_rs_name)
    else:
        # read and process previous year's landgen land type data and transitions to calculate this year's land cover distribution
        if prev_year is not None:
            # read previous year's landgen land type data
            prev_out_file = Path(com_config_dict['out_path']) / prev_fname
            if prev_out_file.exists():
                print(f"Reading previous year landgen land type data from {prev_out_file}")
                # todo: define this in a helper function
                prev_lt_data = read_prev_lt(prev_year, prev_out_file)
            else:
                print(f"Error: Previous year output file {prev_out_file} does not exist; cannot read previous year landgen land type data.")
                sys.exit(1)

            # Calculate this year's land cover distribution using the previous year's data and the transitions
            # these calculations are based on landgen land type outputs
            temp_lt_data = transitions.run(prev_lt_data, year, prev_year)
            # convert this year's land cover distribution to the lc rs classes
            lc_rs_data = lc_rs.convert_landgen_to_lc_rs(temp_lt_data, lc_rs_name)
        else:
            print(f"Error: No previous year data available for year {year}, and _use_lc_rs is False; cannot process land cover data.")
            sys.exit(1)

    # get climate data (1900-2020, four historical periods; and cmip 6 future scenarios, 1km) need to pick the correct period
    # if year < 1900, then use the 1900 climate data
    # todo: probably define this here becase these are specific data 
    climate_data = read_climate_data(year, com_config_dict['source_data_path'], lai_path)

    # get lai data for splitting tree/grass/shrub; this is based on li et al 1km lai data
    #    these data do have short timeseries? then need to select appropriate year
    #todo: this can be in a utils module because other modules need to read these source data 
    lai_data = read_lai_data(year, com_config_dict['source_data_path'], lai_path)

    # todo: use uraster to convert lc_rs_data and climate data and lai data to the landgen grid

    # convert lc_rs_data to the elm land types; this is igbp to generic elm land type mapping
    #    also use the veg continuous fields data; can set modis to elm mapping file name here and read it based on lc_rs_name
    elm_data = lc_rs.convert_lc_rs_to_elm(lc_rs_data, lc_rs_name)

    # split tree/grass/shrub pfts based on cliamte data and li et al 1km lai data
    elm_data = split_tree_grass_shrub(elm_data, climate_data, lai_data)

    # normalize cell by adjusting the land cover distribution to fill the cell land area and reconciling with ocean data (landfrac)
    elm_data = normalize_cell.fill_land(elm_data, landfrac)       # fill_land
    elm_data = normalize_cell.reconcile_ocean(elm_data, landfrac)  # reconcile_ocean

    # now put elm data into lt_year_data

    return





########## run()

## called by land_type.process_single_year() for each year, and this is where the multiprocessing happens for the landcover module
## this sets up the pool and calls the landcover_process() function for each chunk of data

def run(lt_year_data, year, prev_year, prev_fname, lc_rs_path, lc_rs_name,
                            com_config_dict, out_grid_data, ll_limits, cell_ids,
                            manager, grid_manager, lt_manager):



    print(f"Processing landcover module with parameters:")
    # todo: print the parameters here

    # number of available cpu cores (set by SBATCH during job submission)
    omp_threads_str = os.environ.get('OMP_NUM_THREADS')

    if omp_threads_str is not None:
        try:
            # Convert the string value to an integer
            omp_threads_int = int(omp_threads_str)
            print(f"OMP_NUM_THREADS is set to: {omp_threads_int}")
        except ValueError:
            print(f"OMP_NUM_THREADS is set to an invalid integer value: {omp_threads_str}")
    else:
        print("OMP_NUM_THREADS environment variable is not set.")
        # If not set, set to total cores on the node
        omp_threads_int = mp.cpu_count()
        print(f"Using total cores: {omp_threads_int}, but this may fail if
              SBATCH --cpus-per-task is set to a lower number or SBATCH --exclusive is not set")

    # set up the pool and call the landcover_process() function for each chunk of data
    # chunks are defined by the lat-lon limits and corresponding landgen grid cell ids for the chunk;
    #    these are created in land_type.process_single_year() and passed to this run() function as lists?
    # there are more chunks than cpus; the pool will manage this for efficiency because chunks vary in size
    # the results will be stored directly in the lt_year_data shared structure

    # get the manager locks for the shared data structures
    # should not need the managers in the worker process
    man_lock = manager.Lock()
    grid_lock = grid_manager.Lock()
    lt_lock = lt_manager.Lock()

## todo: figure out the data to pass here
# each chunk is a tuple of the arguments for landcover_process, residing in a list
# each tuple includes the lat/lon limits and cell ids for the chunk, and the static arguments that are repeated for each chunk
# e.g.: data_chunks = [(lt_year_data, year, prev_year, prev_fname, lc_rs_path, lc_rs_name,
#          com_config_dict, out_grid_data, ll_limits1, cell_ids1, man_lock, grid_lock, lt_lock),
#          (lt_year_data, year, prev_year, prev_fname, lc_rs_path, lc_rs_name,
#           com_config_dict, out_grid_data, ll_limits2, cell_ids2, man_lock, grid_lock, lt_lock), etc]


    with mp.Pool(processes=omp_threads_int) as pool:
        pool.starmap(landcover_process, data_chunks)

    return