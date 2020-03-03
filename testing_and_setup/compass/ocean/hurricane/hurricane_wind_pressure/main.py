from winds_io import import_data
from winds_io import output_data
from structures import geogrid
import sys
import numpy as np
from winds import parameters
from winds import wind_model

def sim_hurricane():
    # Read in the input file to check which grid we are using
    print('Import user inputs')
    traj_filename, grid_flag, grid_filename, ambient_pressure, holland_b_param = \
        import_data.read_input_file('hurricane_inputs.txt')

    # Read grid-specific parameters and create grid
    print('Read-in grid')
    grid = import_data.initialize_grid(grid_filename, grid_flag)

    # Read hurricane trajectory and set hurricane parameters
    print('Initialize hurricane trajectory data')
    curr_hurricane = import_data.initialize_hurricane(traj_filename, ambient_pressure, holland_b_param)

    # Define parameters
    print('Define parameters')
    params = define_params(curr_hurricane)

    # Compute winds on grid
    print('Compute winds')
    winds = compute_winds(curr_hurricane, params, grid)

    # Output results
    print('Output results')
    output_data.write_netcdf('out.nc', curr_hurricane, grid, winds)

def compute_winds(curr_hurricane, params, grid: geogrid):
    ntimes = len(curr_hurricane) - 1
    mywinds = []
    for it in range(0, ntimes):
        print('Time iteration %d / %d' % (it + 1, len(curr_hurricane) - 1))
        mywinds.append(wind_model.WindModel(params, curr_hurricane[it], grid))

    return mywinds

def define_params(curr_hurricane):
    lat = []
    for i in range(0, len(curr_hurricane)):
        lat.append(curr_hurricane[i].center[1])
    return parameters.Parameters(np.mean(lat))


if __name__ == "__main__":
    sim_hurricane()

    print('Program executed succesfully')
    sys.exit(0)
    # # Read in the input file to check which grid we are using
    # traj_filename, grid_flag, grid_filename = import_data.read_input_file('hurricane_inputs.txt')
    #
    # # Read hurricane trajectory
    # traj = import_data.read_json(traj_filename)
    #
    # # Create trajectory object
    # curr_hurricane = initialize_hurricane(traj)
    #
    # # Read grid-specific parameters
    # if grid_flag == 1:
    #     xll, yll, cellsize, numcells_lat, numcells_lon = import_data.read_raster_inputs(grid_filename)
    # else:
    #     coord = import_data.read_netcdf(grid_filename)

    # Create the grid


