import json
from netCDF4 import Dataset
import numpy as np
import math
from hurricane_model import hurricane
from structures import geogrid
import datetime

def read_grid_file(grid_filename: str, grid_flag: int) -> (float, float):
    if grid_flag == 1:
        xll, yll, cellsize, numcells_lat, numcells_lon = read_raster_inputs(grid_filename)
        lon, lat = setup_regular_grid(xll, yll, cellsize, numcells_lat, numcells_lon)
    else:
        lon, lat = read_netcdf(grid_filename)

    return lon, lat

def read_input_file(filename: str) -> (str, int, str, float, float):
    try:
        f = open(filename, "r")
    except FileNotFoundError as fnf_error:
        raise fnf_error

    traj_filename = f.readline().rstrip('\n')
    grid_flag = f.readline().rstrip('\n').split()
    grid_flag = int(grid_flag[0])
    grid_filename = f.readline().rstrip('\n')
    ambient_pressure = f.readline().rstrip('\n').split()
    ambient_pressure = float(ambient_pressure[0])
    holland_b_param = f.readline().rstrip('\n').split()
    holland_b_param = float(holland_b_param[0])

    f.close()

    return traj_filename, grid_flag, grid_filename, ambient_pressure, holland_b_param


def setup_regular_grid(xll: float, yll: float, cellsize: float, numcells_lat: int, numcells_lon: int) -> (float, float):
    npoints = numcells_lat * numcells_lon
    lon = np.zeros((npoints, ))
    lat = np.zeros((npoints, ))
    k = 0
    for i in range(0, numcells_lon):
        for j in range(0, numcells_lat):
            lon[k] = xll + (float(i) + 0.5) * cellsize
            lat[k] = yll + (float(j) + 0.5) * cellsize
            k += 1

    lat = lat * math.pi / 180.  # Convert to radians
    lon = lon * math.pi / 180.  # Convert to radians

    return lon, lat


def read_raster_inputs(filename: str) -> (float, float, float, int, int):
    try:
        f = open(filename, "r")
    except FileNotFoundError as fnf_error:
        raise fnf_error

    # longitude of the south west corner in deg
    temp = f.readline().rstrip('\n').split()
    xll = float(temp[0])
    # latitude of the south west corner in deg
    temp = f.readline().rstrip('\n').split()
    yll = float(temp[0])
    # cell size in deg
    temp = f.readline().rstrip('\n').split()
    cellsize = float(temp[0])
    # number of cells for latitude
    temp = f.readline().rstrip('\n').split()
    numcells_lat = int(temp[0])
    # number of cells for longitude
    temp = f.readline().rstrip('\n').split()
    numcells_lon = int(temp[0])

    f.close()

    return xll, yll, cellsize, numcells_lat, numcells_lon


def read_json(filename: str):
    try:
        with open(filename) as json_data:
            json_raw = json.load(json_data)
            return json_raw

    except FileNotFoundError as fnf_error:
        raise fnf_error


def read_netcdf(filename: str) -> (float, float):
    # http://unidata.github.io/netcdf4-python/#section1
    # lat and lon from the netCDF file are assumed in radians
    try:
        nc = Dataset(filename)
        temp_lat = nc.variables['latCell'][:]
        temp_lon = nc.variables['lonCell'][:]

        # Convert to numpy array for subsequent processing
        lat = np.array(temp_lat)
        lon = np.array(temp_lon) - 2. * math.pi
        for i in range(0, len(lon)):
            if lon[i] <= -math.pi:
                lon[i] += 2. * math.pi

        return lon, lat

    except FileNotFoundError as fnf_error:
        raise fnf_error


def initialize_hurricane(traj_filename: str, ambient_pressure: float, holland_b_param: float) -> list:
    # JSON Specs
    # "timeUnits": "hours",
    # "distanceUnits": "miles",
    # "windspeedUnits": "knots",
    # "pressureUnits": "mb",

    json_raw = read_json(traj_filename)

    ref_date = datetime.datetime.strptime(json_raw['initialTime'],'%Y-%m-%d_%H:%M:%S')

    curr_hurricane = []
    traj = json_raw['stormTrack']['features']

    for it in range(0, len(traj)):
        coord = traj[it]['geometry']['coordinates']
        center_coord = [x * math.pi / 180. for x in coord]  # degree to rad
        extent = traj[it]['properties']['rMax'] * 1.60934   # miles to km
        pmin = traj[it]['properties']['minP']               # in mbar
        deltap = ambient_pressure - pmin                    # in mbar
        time = traj[it]['properties']['time']               # in hrs
        vmax = traj[it]['properties']['wMax'] * 1.852       # from knots to km/h

        curr_hurricane.append(hurricane.Hurricane(tuple(center_coord), extent, pmin, deltap, vmax,
                                                  holland_b_param, time, ref_date))

    # Compute the components of the forward velocity
    for it in range(0, len(traj) - 1):
        x1 = curr_hurricane[it].center[0]
        y1 = curr_hurricane[it].center[1]

        x2 = curr_hurricane[it + 1].center[0]
        y2 = curr_hurricane[it + 1].center[1]

        theta = math.atan2(y2 - y1, x2 - x1)
        vf = traj[it]['properties']['vf'] * 1.852
        curr_hurricane[it].set_vf((vf * math.cos(theta), vf * math.sin(theta)))

    return curr_hurricane


def initialize_grid(grid_filename: str, grid_flag: int) -> geogrid.GeoGrid:
    lon, lat = read_grid_file(grid_filename, grid_flag)
    return geogrid.GeoGrid(lon, lat)
