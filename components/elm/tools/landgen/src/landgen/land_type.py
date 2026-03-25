# land_type.py
# this module processes land type data for the landgen workflow
# the output is a complete land type distribution
#    and includes some data associated with particular land types

# run() function is the main entry point for this module, and will be called by landgen.py

import multiprocessing as mp
import importlib
from pathlib import Path
import shared_data

########## define helper functions for land_type here

def calc_ll_limits(size_degrees):
    """Calculate and return a list of tuples (min_lat, max_lat, min_lon, max_lon)
    for each size_degrees x size_degrees chunk covering the globe.

    Latitude  spans -90 to  90 degrees.
    Longitude spans -180 to 180 degrees.

    Returns:
        list of tuples: [(min_lat, max_lat, min_lon, max_lon), ...]
    """
    ll_limits = []
    if 90 % size_degrees != 0 or 180 % size_degrees != 0:
        raise ValueError(
            f"size_degrees ({size_degrees}) must evenly divide "
            f"both 90 (latitude half-range) and 180 (longitude half-range)."
        )
    lat = -90.0
    while lat < 90.0:
        max_lat = min(lat + size_degrees, 90.0)
        lon = -180.0
        while lon < 180.0:
            max_lon = min(lon + size_degrees, 180.0)
            ll_limits.append((lat, max_lat, lon, max_lon))
            lon = max_lon
        lat = max_lat
    return ll_limits

##### process_single_year()
def _process_single_year(lt_year_data, year, prev_year, out_fname, lc_rs_path, lc_rs_name, crop_path, urban_path,
                         lake_path, ice_path, wetland_path, harvest_path, assoc_path, com_config_dict, out_grid_data,
                         manager, grid_manager, lt_manager):
    """Process land type data for a single year."""

    # arguments
    # lt_year_data: the shared data structure for the land type data for this year
    # year: the year for which to process the land type data
    # prev_year: the previous year for which land type data were processed

    # other arguments are described below for the run() function

    # data chunks are based on 10x10 degree lat-lon boxes (648 chunks)
    #    15x15 degree box gives 288 chunks, 30x30 box gives 72 chunks
    # todo: add this decomp box size to the config file and to com_config_dict
    decomp_box_size_degrees = 10
    # get the lat-lon limits for the landgen grid
    # ll_limits = list(float) of [(min_lat, max_lat, min_lon, max_lon),(min_lat, max_lat, min_lon, max_lon),... for each chunk]
    chunk_ll_limits = calc_ll_limits(decomp_box_size_degrees)


    # todo: set the multiprocessing chunk info here for all run functions
    # need the lat/lon limits for the chunk, and the corresponding landgen grid cell ids for the chunk
    # do this by lat-lon because all source raw data are on lat-lon grids, but at different resolutions 
    #    and since the source data are several and varied,
    #    it will be faster to read in just the corresponding landgren grid cell info,
    #    even though it is non-sequential
    # the chunks won't be equal size on the landgen grid, but will be consistent,
    #    and will be varied in size based on input source raw data res
    # so put the lat/lon limits and cell ids in corresponding chunks and use the process queue or pool for efficiency
    # arguments for each run function below will include data chunk list with the lat/lon limits and cell ids for the chunk
    #    or the entire list is created here
    # the chunked data are a list of tuples with each argument; 
    #    an example is that the static arguments here will be repeated in each tuple,
    #  and the lat/lon and cell ids will be different for each chunk
    # this info will have to be passed to the run functions below



    # Process landcover
    ## todo: determine prev_fname from out_fname and prev_year
    # each module's run function calls the multiple processes because these modules need to be done sequentially
    landcover = importlib.import_module('landcover')
    lc_data = landcover.run(lt_year_data, year, prev_year, prev_fname, lc_rs_path, lc_rs_name,
                            com_config_dict, out_grid_data, ll_limits, cell_ids, manager, grid_manager, lt_manager)

    # Process crop data - adjust lc crop area
    crop = importlib.import_module('crop')
    lc_data = crop.run(lt_year_data, year, prev_year, crop_path, com_config_dict, out_grid_data, ll_limits, cell_ids,
                       manager, grid_manager, lt_manager)

    # Process urban data - adjust lc urban area
    urban = importlib.import_module('urban')
    lc_data = urban.run(lt_year_data, year, prev_year, urban_path, com_config_dict, out_grid_data, ll_limits, cell_ids,
                        manager, grid_manager, lt_manager)

    # Process lake data - adjust lc lake area
    lake = importlib.import_module('lake')
    lc_data = lake.run(lt_year_data, year, prev_year, lake_path, com_config_dict, out_grid_data, ll_limits, cell_ids,
                       manager, grid_manager, lt_manager)

    # Process ice data - adjust lc ice area
    ice = importlib.import_module('ice')
    lc_data = ice.run(lt_year_data, year, prev_year, ice_path, com_config_dict, out_grid_data, ll_limits, cell_ids,
                      manager, grid_manager, lt_manager)

    # Process wetland data - adjust lc wetland area
    # (may not be needed as the main source is currently the modis cover data;
    #  can allow for this in the future)
    #wetland = importlib.import_module('wetland')
    #lc_data = wetland.run(lt_year_data, year, prev_year, wetland_path, com_config_dict, out_grid_data, ll_limits, cell_ids, manager, grid_manager, lt_manager)

    # Process harvest/grazing data - adjust harvest/grazing area
    harvest = importlib.import_module('harvest')
    lc_data = harvest.run(lt_year_data, year, prev_year, harvest_path, com_config_dict, out_grid_data, ll_limits, cell_ids,
                          manager, grid_manager, lt_manager)

    # Normalize cell
    normalize_cell = importlib.import_module('normalize_cell')
    lc_data = normalize_cell.fill_land(lt_year_data, out_grid_data, ll_limits, cell_ids, manager, grid_manager, lt_manager)       # fill_land
    lc_data = normalize_cell.reconcile_ocean(lt_year_data, out_grid_data, ll_limits, cell_ids, manager, grid_manager, lt_manager)  # reconcile_ocean

    # Process veg-associated data
    veg_assoc = importlib.import_module('veg_assoc')
    lc_data = veg_assoc.run(lt_year_data, year, prev_year, assoc_path, com_config_dict, out_grid_data, ll_limits, cell_ids,
                            manager, grid_manager, lt_manager)

    # Ensure consistency
    consistency = importlib.import_module('consistency')
    lc_data = consistency.run(lt_year_data, year, out_grid_data, ll_limits, cell_ids, manager, grid_manager, lt_manager)

    return

########## run()

## arguments
## these first ones are module-specific parameters that are set in the config file
# active: true = module is run, false = module is skipped
# out_fname: output filename for the module
# the rest of the params set in the config file
# com_config_dict: the shared dictionary for the common parameters for all modules
# out_grid_data: the shared data structure for the landgen grid data

## output

def run(active, out_fname, lc_rs_path, lc_rs_name, crop_path, urban_path, lake_path, ice_path,
        wetland_path, harvest_path, assoc_path, com_config_dict, out_grid_data, manager, grid_manager):
    if active is False:
        print(f"Skipping land_type module")
        return

    # set up the land_type module shared data structure
    # this holds only one year of data, so write it each year
    LtManager.register('LtData', LtData)
    lt_manager = LtManager()
	lt_manager.start()
	lt_year_data = lt_manager.LtData()
	lt_year_data.allocate()

    print(f"Processing land_type module with parameters:")
    # todo: print the parameters here

    # processing code for land_type
    years = np.arange(start_year, end_year + 1)
    output_file = Path(out_path) / out_fname

    prev_year = None

    # 1. Loop over desired years
    for year in years:
        # 2. Process single year
        print(f"  Processing year: {year}")
        _process_single_year(lt_year_data, year, prev_year, out_fname, lc_rs_path, lc_rs_name, crop_path, urban_path,
                             lake_path, ice_path, wetland_path, harvest_path, assoc_path, com_config_dict, out_grid_data,
                             manager, grid_manager, lt_manager)
        


        # append this year's data to the output file
        # todo: can we prepend data is going backwards in time, so that the output file is in chronological order?
        #    otherwise need to write it after all years are processed 

        prev_year = year



    # free the module-specific shared data structure
    lt_year_data = None
    lt_manager.shutdown()
    return
        
