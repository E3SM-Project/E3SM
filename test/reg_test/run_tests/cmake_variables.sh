# These variables are set by CMake
HOMME_DIR=@Homme_Build_DIR@
HOMME_TEST_RESULTS=@Homme_Results_DIR@

# The location of the baseline results
HOMME_BASELINE_DIR=@Homme_Baseline_DIR@

# The location of the baseline results
HOMME_DEFAULT_BASELINE_DIR=${HOMME_DIR}/tests/baseline

# The location of the tests directory
HOMME_TESTING_DIR=${HOMME_DIR}/tests

# Are we using queuing
HOMME_QUEUING=@HOMME_QUEUING@

# The "type" of submission (lsf, pbs, none mpi etc.) for creating the executable scripts 
HOMME_Submission_Type=@HOMME_SUBMISSION_TYPE@

# Account ID for charging
HOMME_ACCOUNT=@HOMME_PROJID@

# The cprnc Netcdf comparison tool
CPRNC_BINARY=@CPRNC_BINARY@

# The cprnc Netcdf comparison tool
PYTHON_EXECUTABLE=@PYTHON_EXECUTABLE@

# Whether to run the OpenMP tests
RUN_OPENMP=@ENABLE_OPENMP@

# User defined MPI executable
MPI_EXEC=@USE_MPIEXEC@

# User defined MPI options
MPI_OPTIONS="@USE_MPI_OPTIONS@"

