from utils import expect

from netCDF4 import Dataset

import pathlib

###########################################################################
def nc_file_init(output_file,ne=0,np=0,ncol=0,nlev=0,phys_grid="gll",overwrite=False):
###########################################################################

    expect (phys_grid=="gll", "Error! So far, only GLL physics grid supported.")

    ofile = pathlib.Path(output_file).resolve().absolute()
    expect ( (not ofile.exists()) or overwrite,
            "Error! Output file '{}' already exists. Use -o to overwrite.".format(ofile))

    ds = Dataset(ofile,"w",format="NETCDF3_64BIT_OFFSET",persist=True)

    expect (ncol >= 1 or (ne>1 and np>1),
            "Error! In order to create a dataset, do one of the following:\n"
            "   - specify ncol via -ncol flag (ncol>=1)\n"
            "   - specify ne via -ne flag, and possibly np via -np flag (ne>=2, np>=2)\n")
    expect (nlev > 1, "Error! Invalid value for nlev (must be >=2).")
    expect (ncol==0 or ne==0,
            "Error! You can specify either -ncol or -ne, but not both.")

    # Allow user to specify arbitrary ncols, which can be handy for small unit tests
    if ncol>0:
        ncols = ncol
    else:
        npm1 = np - 1
        nelems = ne*ne*6
        ncols = nelems*npm1*npm1 + 2
        ds.ne   = ne
        ds.np   = np

    ds.createDimension("time",1)
    ds.createDimension("ncol",ncols)
    ds.createDimension("lev",nlev)
    ds.createDimension("ilev",nlev+1)

    ds.ncol = ncols
    ds.nlev = nlev

    ds.close()

    return True
