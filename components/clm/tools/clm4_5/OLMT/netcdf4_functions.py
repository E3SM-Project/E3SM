#Python utilities for reading and writing variables to a netcdf file
#  using Scientific Python OR scipy, whichever available

def getvar(fname, varname):
    from netCDF4 import Dataset
    nffile = Dataset(fname,"r")
    if varname in nffile.variables:
      varvals = nffile.variables[varname][:]
    else:
      print('Warning: '+varname+' not in '+fname)
      varvals=[-1]
    nffile.close()
    return varvals

def putvar(fname, varname, varvals):
    from netCDF4 import Dataset
    import numpy as np
    nffile = Dataset(fname,"a")
    if (varname in nffile.variables):
      nffile.variables[varname][...] = varvals
    else:
      print('Warning: '+varname+' not in '+fname)
    nffile.close()
    ierr = 0
    return ierr
