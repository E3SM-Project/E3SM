import os
import sys

__version__ = "v2.5.0rc1"
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
