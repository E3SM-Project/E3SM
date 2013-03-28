# These variables are set by CMake
HOMME_DIR=@Homme_Build_DIR@
HOMME_TEST_RESULTS=@Homme_Results_DIR@

# The location of the baseline results
HOMME_BASELINE_DIR=@Homme_Baseline_DIR@

# The location of the tests directory
HOMME_TESTING_DIR=${HOMME_DIR}/tests

# The "type" of submission (lsf, pbs, standard mpi etc.) for creating the executable scripts 
HOMME_Submission_Type=@Homme_Submission_Type@

# Whether to use cprnc to diff the Netcdf files
USE_CPRNC=@TEST_USING_CPRNC@

# The cprnc Netcdf comparison tool
CPRNC_BINARY=@CPRNC_BINARY@

# The cprnc Netcdf comparison tool
PYTHON_EXECUTABLE=@PYTHON_EXECUTABLE@


