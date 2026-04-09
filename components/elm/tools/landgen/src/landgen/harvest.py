# harvest.py
# this module processes harvest and grazing data for the landgen workflow
# the output is a complete harvest distribution
#    and includes some data associated with particular harvests

# run() function is the main entry point for this module, and will be called by process_single_year in land_type.py

import multiprocessing as mp
from pathlib import Path
import shared_data
import landgen_io
import normalize_cell # not created yet
import pandas as pd
import os

########## define helper functions for harvest run() here

##### harvest_process()

## arguments
# lc_data: land cover data structure that is passed between modules
# year: the year for which to process the land cover data
# source_data_path: base path to the source data
# landgen_grid_path: path from source_data_path and the filename of the landgen grid
# out_path: base path for the output data; this is needed to read in the previous year's harvest data

## output

def harvest_process(lt_year_data, year, harvest_path, harvest_name, grazing_path, grazing_names,
                    com_config_dict, out_grid_data, ll_limits, cell_ids,
                    global_mesh_df, man_lock, grid_lock, lt_lock):

    print(f"Processing harvest and grazing module year {year} with parameters:")

    # each worker writes its temp files to a unique subdirectory to avoid collisions
    min_lat, max_lat, min_lon, max_lon = ll_limits
    tmp_dir = (
        Path(com_config_dict['out_path'])
        / 'tmp'
        / f"harvest_{year}_{min_lat:.0f}_{min_lon:.0f}"
    )

    # --- read source data (full global arrays; slicing happens inside regrid) ---
    harvest_data = landgen_io.read_luh2_harvest(year, harvest_path, harvest_name)
    grazing_data = landgen_io.read_hyde_grazing(year, grazing_path, grazing_names)

    # --- regrid and store harvest variables into lt_year_data.harvest_frac ---
    # LUH2_HARVEST_VARS order matches the n_harvest=5 dimension in LtData:
    #   index 0: primf_harv, 1: primn_harv, 2: secmf_harv, 3: secyf_harv, 4: secnf_harv
    for i, varname in enumerate(landgen_io.LUH2_HARVEST_VARS):
        regridded = landgen_io.regrid_to_landgen_grid(
            harvest_data[varname],
            harvest_data['lat'],
            harvest_data['lon'],
            cell_ids, ll_limits,
            global_mesh_df,
            tmp_dir / varname,
            varname,
        )
        with lt_lock:
            lt_year_data.harvest_frac[cell_ids, i] = regridded

    # --- regrid and store grazing variables into lt_year_data.grazing_frac ---
    # grazing_names dict order matches n_grazing=2 dimension in LtData:
    #   index 0: pasture (grazing_land), index 1: rangeland
    for i, category in enumerate(grazing_names.keys()):
        regridded = landgen_io.regrid_to_landgen_grid(
            grazing_data[category],
            grazing_data['lat'],
            grazing_data['lon'],
            cell_ids, ll_limits,
            global_mesh_df,
            tmp_dir / category,
            category,
        )
        with lt_lock:
            lt_year_data.grazing_frac[cell_ids, i] = regridded

    return


########## run()

## called by land_type.process_single_year() for each year, and this is where the multiprocessing happens for the landcover module
## this sets up the pool and calls the harvest_process() function for each chunk of data

def run(lt_year_data, year, prev_year, harvest_path, harvest_name, grazing_path, grazing_names,
        com_config_dict, out_grid_data, manager, grid_manager, lt_manager):

    print(f"Processing harvest module with parameters:")
    # todo: print the parameters here

    # load the global HEALPix mesh parquet once here so worker processes
    # don't each re-read the 37 MB file; pass the DataFrame into each chunk tuple
    global_parquet_path = (
        Path(com_config_dict['source_data_path'])
        / Path(com_config_dict['landgen_grid_path']).parent
        / 'merged_land_cells.parquet'
    )
    global_mesh_df = pd.read_parquet(global_parquet_path)
    print(f"  Loaded HEALPix mesh: {len(global_mesh_df)} cells from {global_parquet_path}")

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

    # set up the pool and call the harvest_process() function for each chunk of data
    # chunks are defined by the lat-lon limits and corresponding landgen grid cell ids for the chunk;
    #    these are created in land_type.process_single_year() and passed to this run() function as lists?
    # there are more chunks than cpus; the pool will manage this for efficiency because chunks vary in size
    # the results will be stored directly in the lt_year_data shared structure

    # get the manager locks for the shared data structures
    # should not need the managers in the worker process
    man_lock = manager.Lock()
    grid_lock = grid_manager.Lock()
    lt_lock = lt_manager.Lock()

    # Build data_chunks: one tuple per spatial chunk covering the globe in 10x10 degree boxes.
    # For each chunk, filter the global mesh to cells whose centroid lat/lon falls within the box.
    import land_type as _lt
    decomp_box_size_degrees = 10
    chunk_ll_limits = _lt.calc_ll_limits(decomp_box_size_degrees)

    data_chunks = []
    for ll in chunk_ll_limits:
        min_lat, max_lat, min_lon, max_lon = ll
        mask = (
            (global_mesh_df['lat'] >= min_lat) & (global_mesh_df['lat'] < max_lat) &
            (global_mesh_df['lon'] >= min_lon) & (global_mesh_df['lon'] < max_lon)
        )
        chunk_cell_ids = global_mesh_df.loc[mask, 'cellid'].values
        if len(chunk_cell_ids) == 0:
            continue  # skip ocean-only or empty chunks
        data_chunks.append((
            lt_year_data, year, harvest_path, harvest_name, grazing_path, grazing_names,
            com_config_dict, out_grid_data, ll, chunk_cell_ids,
            global_mesh_df, man_lock, grid_lock, lt_lock,
        ))

    print(f"  Submitting {len(data_chunks)} harvest chunks to pool of {omp_threads_int} workers")

    with mp.Pool(processes=omp_threads_int) as pool:
        pool.starmap(harvest_process, data_chunks)

    return