import netCDF4 as nc4

if __name__ == "__main__":

    mesh_file = sys.argv[1]
    nc_mesh = nc4.Dataset(mesh_file, 'r+')

    if 'cullCell' not in nc_vars:
        nc_mesh.createVariable('cullCell', 'i', ('nCells'))
    
    mesh.variables['cullCell'][:] = nc_mesh.variables['bathymetry'][:] > 20.0
    nc_mesh.close()
