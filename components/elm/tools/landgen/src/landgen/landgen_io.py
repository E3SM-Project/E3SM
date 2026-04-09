# landgen_io.py
# Utility functions for reading harvest data from LUH2 and grazing data from HYDE3.5

import xarray as xr
import numpy as np
from pathlib import Path
import pandas as pd
import rasterio
from rasterio.transform import from_bounds
import json
import shapely.wkb
from uraster.classes.uraster import uraster as URaster

# Default harvest variable names from LUH2 transitions.nc
LUH2_HARVEST_VARS = [
    'primf_harv',   # wood harvest area from primary forest land
    'primn_harv',   # wood harvest area from primary non forest land
    'secmf_harv',   # wood harvest area from secondary mature forest land
    'secyf_harv',   # wood harvest area from secondary young forest land
    'secnf_harv',   # wood harvest area from secondary non forest land
]

#--------------------------------------------------------------------------
def _get_year_idx(time_values, year, ncfile, time_units=None):
    """
    Find the time index in a NetCDF time axis for the requested calendar year.
    Handles two common CF-convention unit patterns:

      'years since <epoch>'  -- time_values are year offsets from epoch
      'days since <epoch>'   -- time_values are day offsets from epoch;
                                converted to fractional years via / 365.25

    If time_units is None, time_values are assumed to already be calendar years.
    Raises ValueError if the closest available year is more than 1 year away.

    Args:
        time_values (np.ndarray): 1D array of numeric time values.
        year (int):               Requested calendar year.
        ncfile (Path):            File path, used only in error messages.
        time_units (str|None):    CF-convention units string, e.g.
                                  'years since 850-01-01 0:0:0' or
                                  'days since 1-5-1 00:00:00'.

    Returns:
        int: Index into time_values.
    """
    units_lower = time_units.strip().lower() if time_units else ''

    if units_lower.startswith('years since'):
        # time_values are year offsets from epoch_year
        # e.g. 'years since 850-01-01' -> epoch_year=850; value 1160 -> year 2010
        epoch_str  = time_units.strip().split('since', 1)[1].strip()
        epoch_year = int(epoch_str.split('-')[0])
        search_value = year - epoch_year
        # convert matched offset back to calendar year for error message
        def _to_cal(offset): return offset + epoch_year

    elif units_lower.startswith('days since'):
        # time_values are day offsets from epoch_year
        # e.g. 'days since 1-5-1 00:00:00' -> epoch_year=1
        # convert both the time axis and the requested year to fractional years from epoch
        epoch_str  = time_units.strip().split('since', 1)[1].strip()
        epoch_year = int(epoch_str.split('-')[0])
        # convert time_values (days) to fractional calendar years
        time_as_years = time_values / 365.25 + epoch_year
        search_value  = float(year)           # search in calendar-year space
        # re-use time_as_years for the distance search below
        year_idx     = int(np.argmin(np.abs(time_as_years - search_value)))
        actual_year  = time_as_years[year_idx]
        if abs(actual_year - year) > 1:
            raise ValueError(
                f"Requested year {year} not found in {ncfile} "
                f"(closest available: {actual_year:.0f})"
            )
        return year_idx

    else:
        # assume time_values are already calendar years
        search_value = float(year)
        def _to_cal(offset): return offset

    year_idx    = int(np.argmin(np.abs(time_values - search_value)))
    actual_year = _to_cal(time_values[year_idx])

    if abs(actual_year - year) > 1:
        raise ValueError(
            f"Requested year {year} not found in {ncfile} "
            f"(closest available: {actual_year:.0f})"
        )
    return year_idx


#--------------------------------------------------------------------------
def read_luh2_harvest(year, harvest_path, harvest_name, variable_names=None):
    """
    Read LUH2 harvest variables for a given year.

    Args:
        year (int): Year to extract. LUH2 covers 850-2015.
        harvest_path (str or Path): Directory containing the LUH2 NetCDF file
        harvest_name (str): Filename of the LUH2 NetCDF file
        variable_names (list or None): Variables to extract. Defaults to all 5 LUH2_HARVEST_VARS if None.

    Returns:
        dict: {varname: 2D np.ndarray shape (lat, lon)} for the requested year.
              lat/lon coordinate arrays are included as 'lat' and 'lon' keys.
    """
    if variable_names is None:
        variable_names = LUH2_HARVEST_VARS

    ncfile = Path(harvest_path) / harvest_name
    if not ncfile.exists():
        raise FileNotFoundError(f"LUH2 harvest file not found: {ncfile}")

    ds = xr.open_dataset(ncfile, decode_times=False)

    # LUH2 time axis is 'years since 850-01-01'; values are offsets from 850, not calendar years
    time_units = ds['time'].attrs.get('units', None)
    year_idx = _get_year_idx(ds['time'].values, year, ncfile, time_units=time_units)
    print(f"  read_luh2_harvest: reading year {year} (time index {year_idx}) from {ncfile}")

    out = {'lat': ds['lat'].values, 'lon': ds['lon'].values}
    for v in variable_names:
        if v not in ds:
            raise KeyError(f"Variable '{v}' not found in {ncfile}. "
                           f"Available variables: {list(ds.data_vars)}")
        out[v] = ds[v].isel(time=year_idx).values  # shape: (lat, lon)

    ds.close()
    return out

#--------------------------------------------------------------------------
def read_hyde_grazing(year, grazing_path, grazing_names):
    """
    Read HYDE3.5 grazing data for a given year.

    grazing_names is a dict mapping a grazing category label to a NetCDF filename,
    as specified in config.json.  Each file is expected to contain a single variable
    whose name matches the file stem (e.g. 'pasture.nc' -> variable 'pasture').

    Args:
        year (int): Year to extract. HYDE3.5 baseline covers 10000 BCE - 2023 CE.
        grazing_path (str or Path): Directory containing the HYDE3.5 NetCDF files
        grazing_names (dict): Mapping of {category_label: filename}

    Returns:
        dict: {category_label: 2D np.ndarray shape (lat, lon)} for the requested year,
              plus 'lat' and 'lon' coordinate arrays (taken from the first file read).
    """
    if not isinstance(grazing_names, dict):
        raise TypeError(
            f"grazing_names must be a dict (e.g. {{'pasture': 'pasture.nc', "
            f"'rangeland': 'rangeland.nc'}}), got {type(grazing_names).__name__}. "
            f"Please update config.json accordingly."
        )

    grazing_path = Path(grazing_path)
    out = {}

    first = True
    for category, filename in grazing_names.items():
        ncfile = grazing_path / filename
        if not ncfile.exists():
            raise FileNotFoundError(f"HYDE3.5 grazing file not found: {ncfile}")
        # variable name is the file stem, e.g. 'pasture.nc' -> 'pasture'
        varname = Path(filename).stem
        ds = xr.open_dataset(ncfile, decode_times=False)
        time_units = ds['time'].attrs.get('units', None)
        year_idx = _get_year_idx(ds['time'].values, year, ncfile, time_units=time_units)
        print(f"  read_hyde_grazing: reading '{varname}' year {year} "
              f"(time index {year_idx}) from {ncfile}")
        if first:
            out['lat'] = ds['lat'].values
            out['lon'] = ds['lon'].values
            first = False
        out[category] = ds[varname].isel(time=year_idx).values  # shape: (lat, lon)
        ds.close()

    return out

#--------------------------------------------------------------------------
def write_latlon_to_geotiff(data_2d, lat, lon, ll_limits, tmp_path):
    """
    Slice a 2D (lat, lon) source array to a lat-lon bounding box and write
    it to a GeoTIFF file for use as uraster input.

    Args:
        data_2d (np.ndarray): 2D array shape (n_lat, n_lon), global coverage.
        lat (np.ndarray):     1D latitude coordinate array, degrees, south-to-north.
        lon (np.ndarray):     1D longitude coordinate array, degrees, west-to-east.
        ll_limits (tuple):    (min_lat, max_lat, min_lon, max_lon) for this chunk.
        tmp_path (str|Path):  Full path of the GeoTIFF file to write.

    Returns:
        Path: Path to the written GeoTIFF.
    """
    min_lat, max_lat, min_lon, max_lon = ll_limits

    # find row/col indices that fall within the bounding box
    # add a 1-cell buffer on each side to avoid edge interpolation artefacts
    lat_step = float(lat[1] - lat[0])
    lon_step = float(lon[1] - lon[0])
    lat_mask = (lat >= min_lat - abs(lat_step)) & (lat <= max_lat + abs(lat_step))
    lon_mask = (lon >= min_lon - abs(lon_step)) & (lon <= max_lon + abs(lon_step))

    lat_idx = np.where(lat_mask)[0]
    lon_idx = np.where(lon_mask)[0]

    if lat_idx.size == 0 or lon_idx.size == 0:
        raise ValueError(
            f"No source grid cells found within ll_limits {ll_limits}. "
            f"lat range: [{lat.min():.2f}, {lat.max():.2f}], "
            f"lon range: [{lon.min():.2f}, {lon.max():.2f}]"
        )

    # slice the data
    chunk_lat = lat[lat_idx]
    chunk_lon = lon[lon_idx]
    chunk_data = data_2d[np.ix_(lat_idx, lon_idx)].astype(np.float32)

    # rasterio uses (west, south, east, north) bounds
    west  = float(chunk_lon[0])  - abs(lon_step) / 2.0
    east  = float(chunk_lon[-1]) + abs(lon_step) / 2.0
    south = float(chunk_lat[0])  - abs(lat_step) / 2.0
    north = float(chunk_lat[-1]) + abs(lat_step) / 2.0

    n_rows, n_cols = chunk_data.shape
    transform = from_bounds(west, south, east, north, n_cols, n_rows)

    tmp_path = Path(tmp_path)
    with rasterio.open(
        tmp_path,
        'w',
        driver='GTiff',
        height=n_rows,
        width=n_cols,
        count=1,
        dtype=np.float32,
        crs='EPSG:4326',
        transform=transform,
        nodata=np.nan,
    ) as dst:
        # rasterio band 1 is row-ordered north-to-south; flip if lat is south-to-north
        if chunk_lat[0] < chunk_lat[-1]:
            dst.write(np.flipud(chunk_data), 1)
        else:
            dst.write(chunk_data, 1)

    print(f"  write_latlon_to_geotiff: wrote {n_rows}x{n_cols} chunk to {tmp_path}")
    return tmp_path


#--------------------------------------------------------------------------
def write_chunk_mesh_to_geojson(global_mesh_df, cell_ids, tmp_path):
    """
    Filter the global HEALPix mesh DataFrame to only the cells in cell_ids
    and write a chunk-sized GeoJSON file for use as uraster source mesh.

    GeoJSON is used (rather than Parquet) because it is always supported by
    GDAL/OGR without additional plugins (unlike the Parquet/DuckDB driver).

    The caller should load the global parquet once (e.g. in run()) and pass
    the resulting DataFrame here to avoid re-reading the 37 MB file for every
    variable and chunk.

    Args:
        global_mesh_df (pd.DataFrame): Full merged_land_cells DataFrame, already
                                       loaded. Must have 'cellid' (int) and
                                       'geometry' (WKB bytes) columns.
        cell_ids (array-like):         1D array of integer cellid values for chunk.
        tmp_path (str|Path):           Full path of the output GeoJSON file to write.

    Returns:
        Path: Path to the written GeoJSON.
    """
    chunk_df = global_mesh_df[global_mesh_df['cellid'].isin(cell_ids)]

    if chunk_df.empty:
        raise ValueError(
            f"No mesh cells found for the provided cell_ids. "
            f"First few cell_ids: {list(cell_ids[:5])}"
        )

    # build GeoJSON manually from WKB geometry column
    features = []
    for _, row in chunk_df.iterrows():
        geom = shapely.wkb.loads(row['geometry'])
        features.append({
            'type': 'Feature',
            'geometry': geom.__geo_interface__,
            'properties': {'cellid': int(row['cellid'])},
        })

    geojson = {'type': 'FeatureCollection', 'features': features}
    tmp_path = Path(tmp_path)
    tmp_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tmp_path, 'w') as f:
        json.dump(geojson, f)

    print(f"  write_chunk_mesh_to_geojson: wrote {len(features)} cells to {tmp_path}")
    return tmp_path


#--------------------------------------------------------------------------
def regrid_to_landgen_grid(data_2d, src_lat, src_lon, cell_ids, ll_limits,
                            global_mesh_df, tmp_dir, varname,
                            remap_method=3):
    """
    Regrid a single 2D source variable onto the landgen HEALPix grid cells
    for one spatial chunk, using uraster.

    Workflow:
        1. Slice source data to ll_limits and write to a temp GeoTIFF.
        2. Filter global_mesh_df to cell_ids and write a chunk parquet.
        3. Run uraster with iFlag_remap_method=3 (weighted average).
        4. Read the uraster output GeoJSON, extract 'mean' per cellid.
        5. Return a 1D array aligned to cell_ids order; cells with no overlap get 0.
        6. Clean up all temp files.

    Args:
        data_2d (np.ndarray):        2D source array shape (n_lat, n_lon).
        src_lat (np.ndarray):        1D source latitude array.
        src_lon (np.ndarray):        1D source longitude array.
        cell_ids (array-like):       1D array of integer cellid values for this chunk.
        ll_limits (tuple):           (min_lat, max_lat, min_lon, max_lon).
        global_mesh_df (pd.DataFrame): Pre-loaded merged_land_cells DataFrame.
                                     Load once in run() and pass here to avoid
                                     re-reading the 37 MB parquet for every variable.
        tmp_dir (str|Path):          Directory for temporary files (unique per worker).
        varname (str):               Variable name, used for temp file naming only.
        remap_method (int):          uraster iFlag_remap_method
                                     (1=nearest, 2=nearest, 3=weighted average).
                                     Default 3 is correct for area fraction data.

    Returns:
        np.ndarray: 1D float64 array of length len(cell_ids), regridded values
                    in the same order as cell_ids.
    """
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # unique suffix to avoid collisions between parallel workers
    suffix = f"{varname}_{ll_limits[0]:.0f}_{ll_limits[2]:.0f}"
    tmp_raster  = tmp_dir / f"src_{suffix}.tif"
    tmp_mesh_in = tmp_dir / f"mesh_in_{suffix}.geojson"   # GeoJSON: universal GDAL support
    tmp_mesh_out= tmp_dir / f"mesh_out_{suffix}.geojson"

    try:
        # 1. write source chunk as GeoTIFF
        write_latlon_to_geotiff(data_2d, src_lat, src_lon, ll_limits, tmp_raster)

        # 2. write filtered chunk mesh as GeoJSON (parquet requires libgdal-arrow-parquet)
        write_chunk_mesh_to_geojson(global_mesh_df, cell_ids, tmp_mesh_in)

        # 3. run uraster
        config = {
            'sFilename_source_mesh':    str(tmp_mesh_in),
            'aFilename_source_raster':  [str(tmp_raster)],
            'sFilename_target_mesh':    str(tmp_mesh_out),
            'iFlag_remap_method':       remap_method,
            'sField_unique_id':         'cellid',
            'iFlag_global':             0,   # chunk is regional, not global
            'iFlag_polar':              0,
        }
        processor = URaster(config)
        processor.setup()
        processor.run_remap()

        # 4. read output GeoJSON and extract 'mean' per cellid
        with open(tmp_mesh_out, 'r') as f:
            geojson = json.load(f)

        # build a cellid -> mean value lookup from the uraster output
        # cells with no raster overlap are missing the 'mean' key entirely;
        # default those to 0.0
        result_map = {}
        for feature in geojson['features']:
            props = feature['properties']
            cid = int(props['cellid'])
            val = props.get('mean', None)
            result_map[cid] = float(val) if (val is not None and not np.isnan(val)) else 0.0

        # 5. align to cell_ids order; missing cells default to 0
        out = np.array([result_map.get(int(cid), 0.0) for cid in cell_ids],
                       dtype=np.float64)

        print(f"  regrid_to_landgen_grid: '{varname}' -> {len(out)} cells, "
              f"non-zero: {np.count_nonzero(out)}")
        return out

    finally:
        # 6. clean up temp files regardless of success or failure
        for f in [tmp_raster, tmp_mesh_in, tmp_mesh_out]:
            try:
                Path(f).unlink(missing_ok=True)
            except Exception:
                pass
        # uraster also creates a '_fixed' copy of the input mesh; clean that up too
        try:
            Path(str(tmp_mesh_in).replace('.geojson', '_fixed.geojson')).unlink(missing_ok=True)
        except Exception:
            pass

#--------------------------------------------------------------------------