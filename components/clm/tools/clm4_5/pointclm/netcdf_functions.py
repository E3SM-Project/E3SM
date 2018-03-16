#Python utilities for reading and writing variables to a netcdf file
#  using Scientific Python OR scipy, whichever available

def getvar(fname, varname):
    from netCDF4 import Dataset
    nffile = Dataset(fname,"r")
    var = nffile.variables[varname]
    print var[:]
    varvals = var[:]
    nffile.close()
    return varvals

def putvar(fname, varname, varvals):
    from netCDF4 import Dataset
    nffile = Dataset(fname,"a")
    var = nffile.variables[varname]
    var[:] = varvals
    nffile.close()
    ierr = 0
    return ierr
