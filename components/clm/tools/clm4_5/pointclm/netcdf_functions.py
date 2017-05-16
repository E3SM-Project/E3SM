#Python utilities for reading and writing variables to a netcdf file
#  using Scientific Python OR scipy, whichever available

def getvar(fname, varname):
    usescipy = False
    try:
        import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"r",mmap=False)
        var = nffile.variables[varname]
        varvals = var[:].copy()    #works for vector only?
        nffile.close()
    else:
        nffile = netcdf.NetCDFFile(fname,"r")
        var = nffile.variables[varname]
        varvals = var.getValue()
        nffile.close()
    return varvals

def putvar(fname, varname, varvals):
    usescipy = False
    try:
        import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"a",mmap=False)
        var = nffile.variables[varname]
        var[:] = varvals
        nffile.close()
    else:
        nffile = netcdf.NetCDFFile(fname,"a")
        var = nffile.variables[varname]
        var.assignValue(varvals)
        nffile.close()
    ierr = 0
    return ierr
