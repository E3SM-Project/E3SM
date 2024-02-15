import os
import sys

# import shapely here (before any esmpy imports via cdms2) to prevent a
# segfault related to multiprocessing.  Do not import esmpy here to prevent
# issue with dask when using ESMF with system compilers.
import shapely

__version__ = "v2.11.0rc1"
INSTALL_PATH = os.path.join(sys.prefix, "share/e3sm_diags/")

# Disable MPI in cdms2, which is not currently supported by E3SM-unified
os.environ["CDMS_NO_MPI"] = "True"
# Must be done before any CDAT library is called.
os.environ["UVCDAT_ANONYMOUS_LOG"] = "no"
os.environ["CDAT_ANONYMOUS_LOG"] = "no"
# Needed for when using hdf5 >= 1.10.0,
# without this, errors are thrown on Edison compute nodes.
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
# Used by numpy, causes too many threads to spawn otherwise.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
