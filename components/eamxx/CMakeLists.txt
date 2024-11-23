####################################################################
#       Basic setup: version, cime debug, cmake paths,...          #
####################################################################
if (NOT DEFINED PROJECT_NAME)
  cmake_minimum_required(VERSION 3.14)
  cmake_policy(SET CMP0057 NEW)
  set(SCREAM_CIME_BUILD FALSE)

  # Print the sha of the last commit (useful to double check which version was tested on CDash)
  execute_process (COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE LAST_GIT_COMMIT_SHA
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(EAMXX_GIT_VERSION ${LAST_GIT_COMMIT_SHA} CACHE STRING "The sha of the last git commit." FORCE)
  message(STATUS "The sha of the last commit is ${EAMXX_GIT_VERSION}")
else()
  set(SCREAM_CIME_BUILD TRUE)
endif()

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  # The new behavior for cmp0074 makes cmake use (rather than ignore)
  # any <PackageName>_ROOT env/cmake variable previously set.
  cmake_policy(SET CMP0074 NEW)
endif()

set (EAMXX_VERSION_MAJOR 1)
set (EAMXX_VERSION_MINOR 0)
set (EAMXX_VERSION_PATCH 0)

if ($ENV{SCREAM_FORCE_CONFIG_FAIL})
  message(FATAL_ERROR "Failed, as instructed by environment")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_ci)
if (CMAKE_BUILD_TYPE_ci STREQUAL "debug")
  set (SCREAM_DEBUG TRUE)
else ()
  set (SCREAM_DEBUG FALSE)
endif()

# Add the ./cmake folder to cmake path. Also add EKAT's cmake folder
set (EKAT_CMAKE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../externals/ekat/cmake)
list(APPEND CMAKE_MODULE_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/cmake
     ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/modules
     ${EKAT_CMAKE_PATH}
     ${EKAT_CMAKE_PATH}/tpls
)
if (SCREAM_CIME_BUILD)
  list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cime)
endif ()

# We want to use C++17 in EAMxx
set(CMAKE_CXX_STANDARD 17)

if (NOT SCREAM_CIME_BUILD)
  project(SCREAM CXX C Fortran)

  if (SCREAM_CORI_HACK)
    list(APPEND CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "ifcore")
    list(REMOVE_ITEM CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "ifport")
  endif()

else()
  # Ensure our languages are all enabled
  enable_language(C CXX Fortran)
endif()

####################################################################
#                  Kokkos/YAKL-related settings                    #
####################################################################

if (Kokkos_ENABLE_CUDA)
  if (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
    if (Kokkos_ENABLE_DEBUG_BOUNDS_CHECK)
      string (CONCAT msg
          "Kokkos_ENALBE_CUDA_RELOCATABLE_DEVICE_CODE=ON, and Kokkos_ENALBE_DEBUG_BOUNDS_CHECK=ON.\n"
          "   -> Disabling bounds checks, to prevent internal compiler errors.")
      message(WARNING "${msg}")
      set (Kokkos_ENABLE_DEBUG_BOUNDS_CHECK OFF CACHE BOOL "" FORCE)
    else()
      string (CONCAT msg
          "Kokkos_ENALBE_CUDA_RELOCATABLE_DEVICE_CODE=ON.\n"
          "   -> Disabling bounds checks, to prevent internal compiler errors.")
      message(STATUS "${msg}")
      set (Kokkos_ENABLE_DEBUG_BOUNDS_CHECK OFF CACHE BOOL "" FORCE)
    endif()
  endif()
endif()

# Homme's composef90 library (built for f90 support) requires Kokkos::Serial
# to be defined, which in turn requires the CMake var Kokkos_ENABLE_SERIAL
# to be on. For now, simply ensure Kokkos Serial is enabled
option (Kokkos_ENABLE_SERIAL "" ON)

set (EAMXX_ENABLE_GPU FALSE CACHE BOOL "")
set (CUDA_BUILD FALSE CACHE BOOL "") #needed for yakl if kokkos vars are not visible there?
set (HIP_BUILD FALSE CACHE BOOL "") #needed for yakl if kokkos vars are not visible there?

# Determine if this is a Cuda build.
if (Kokkos_ENABLE_CUDA)
  # Add CUDA as a language for CUDA builds
  enable_language(CUDA)
  set (EAMXX_ENABLE_GPU TRUE  CACHE BOOL "" FORCE)
  set (CUDA_BUILD TRUE CACHE BOOL "" FORCE) #needed for yakl if kokkos vars are not visible there?
endif ()

# Determine if this is a Cuda build.
if (Kokkos_ENABLE_HIP)
  # Add CUDA as a language for CUDA builds
  enable_language(HIP)
  set (EAMXX_ENABLE_GPU TRUE CACHE BOOL "" FORCE)
  set (HIP_BUILD TRUE CACHE BOOL "" FORCE) #needed for yakl if kokkos vars are not visible there?
endif ()

if( NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "[Cc]lang" )
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
endif()

####################################################################
#                EAMxx main configuration options                  #
####################################################################

# First, Compute reasonable defaults. This needs to happen before the CACHE variables are set
set(DEFAULT_MAX_RANKS 4)
set(DEFAULT_MAX_THREADS 16)
set(DEFAULT_MIMIC_GPU FALSE)
set(DEFAULT_FPE FALSE)
set(DEFAULT_PACK_SIZE 16)
if (EAMXX_ENABLE_GPU)
  # On the GPU, the pack size must be 1
  set(DEFAULT_PACK_SIZE 1)
  set(DEFAULT_MAX_THREADS 1)
  # Limit to 1 rank, cause parallel builds testing might limit the number of available gpus
  set(DEFAULT_MAX_RANKS 1)
else()
  if (SCREAM_DEBUG)
    set(DEFAULT_MIMIC_GPU TRUE)
  endif()
endif ()
set(DEFAULT_FPMODEL "precise")
set(DEFAULT_LIB_ONLY FALSE)
if (SCREAM_CIME_BUILD)
  set(DEFAULT_LIB_ONLY TRUE)
endif()
set(DEFAULT_NUM_VERTICAL_LEV 72)

find_path(NF_CONFIG_SEARCH nf-config)
if (NF_CONFIG_SEARCH)
  execute_process(COMMAND ${NF_CONFIG_SEARCH}/nf-config --prefix
    RESULT_VARIABLE NF_STATUS
    OUTPUT_VARIABLE NF_CONFIG_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if (NF_STATUS EQUAL 0)
    set(DEFAULT_NetCDF_Fortran_PATH ${NF_CONFIG_OUTPUT})
  endif()
endif()

find_path(NC_CONFIG_SEARCH nc-config)
if (NC_CONFIG_SEARCH)
  execute_process(COMMAND ${NC_CONFIG_SEARCH}/nc-config --prefix
    RESULT_VARIABLE NC_STATUS
    OUTPUT_VARIABLE NC_CONFIG_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if (NC_STATUS EQUAL 0)
    set(DEFAULT_NetCDF_C_PATH ${NC_CONFIG_OUTPUT})
  endif()
endif()

if (DEFINED ENV{SCREAM_MACHINE})
  set(DEFAULT_SCREAM_MACHINE $ENV{SCREAM_MACHINE})
elseif(MACH)
  # MACH comes from CIME build system
  set(DEFAULT_SCREAM_MACHINE ${MACH})
endif()

set(DEFAULT_SMALL_KERNELS FALSE)
if (Kokkos_ENABLE_HIP)
  set(DEFAULT_SMALL_KERNELS TRUE)
endif()

if (SCREAM_DEBUG)
  set(DEFAULT_FPMODEL "strict")
  if (SCREAM_PACK_SIZE EQUAL 1 AND NOT EAMXX_ENABLE_GPU)
    set(DEFAULT_FPE TRUE)
  endif ()
endif()

### Now that reasonable defaults have been computed, set CACHE vars
set(SCREAM_MIMIC_GPU ${DEFAULT_MIMIC_GPU} CACHE BOOL "Mimic GPU to correctness-test inter-column parallelism on non-GPU platform")
set(SCREAM_PACK_CHECK_BOUNDS FALSE CACHE BOOL "If defined, scream::pack objects check indices against bounds")
set(SCREAM_LIB_ONLY ${DEFAULT_LIB_ONLY} CACHE BOOL "Only build libraries, no exes")
set(NetCDF_Fortran_PATH ${DEFAULT_NetCDF_Fortran_PATH} CACHE FILEPATH "Path to netcdf fortran installation")
set(NetCDF_C_PATH ${DEFAULT_NetCDF_C_PATH} CACHE FILEPATH "Path to netcdf C installation")
set(SCREAM_MACHINE ${DEFAULT_SCREAM_MACHINE} CACHE STRING "The CIME/SCREAM name for the current machine")
option(SCREAM_MPI_ON_DEVICE "Whether to use device pointers for MPI calls" ON)
option(SCREAM_ENABLE_MAM "Whether to enable MAM aerosol support" ON)
set(SCREAM_SMALL_KERNELS ${DEFAULT_SMALL_KERNELS} CACHE STRING "Use small, non-monolothic kokkos kernels for ALL components that support them")
set(SCREAM_P3_SMALL_KERNELS ${SCREAM_SMALL_KERNELS} CACHE STRING "Use small, non-monolothic kokkos kernels for P3 only")
set(SCREAM_SHOC_SMALL_KERNELS ${SCREAM_SMALL_KERNELS} CACHE STRING "Use small, non-monolothic kokkos kernels for SHOC only")
if (NOT SCREAM_P3_SMALL_KERNELS AND NOT SCREAM_SHOC_SMALL_KERNELS)
  set(EKAT_DISABLE_WORKSPACE_SHARING TRUE CACHE STRING "")
endif()

# Add RRTMGP settings. Note, we might consider also adding RRTMGP_EXPENSIVE_CHECKS
# to turn on the RRTMGP internal checks here as well, via
#     option (RRTMGP_EXPENSIVE_CHECKS "Turn on internal RRTMGP error checking" ${SCREAM_DEBUG})
# and then adding to scream_config.h:
#     #cmakedefine RRTMGP_EXPENSIVE_CHECKS
option (SCREAM_RRTMGP_DEBUG "Turn on extra debug checks in RRTMGP" ${SCREAM_DEBUG})

option(SCREAM_RRTMGP_ENABLE_YAKL "Use YAKL under rrtmgp" TRUE)
option(SCREAM_RRTMGP_ENABLE_KOKKOS "Use Kokkos under rrtmgp" FALSE)
if (SCREAM_RRTMGP_ENABLE_YAKL)
  add_definitions("-DRRTMGP_ENABLE_YAKL")
endif()

if (SCREAM_RRTMGP_ENABLE_KOKKOS)
  add_definitions("-DRRTMGP_ENABLE_KOKKOS")
endif()


set(SCREAM_DOUBLE_PRECISION TRUE CACHE BOOL "Set to double precision (default True)")

# For now, only used in share/grid/remap/refining_remapper_rma.*pp
option (EAMXX_ENABLE_EXPERIMENTAL_CODE "Compile one-sided MPI for refining remappers" OFF)

option (SCREAM_ENABLE_ML_CORRECTION "Whether to enable ML correction parametrization" OFF)

# Set number of vertical levels
set(SCREAM_NUM_VERTICAL_LEV ${DEFAULT_NUM_VERTICAL_LEV} CACHE STRING
    "The number of levels used in the vertical grid."
)
option(SCREAM_HAS_LEAP_YEAR "Whether scream uses leap years or not" ON)

set(SCREAM_FPMODEL ${DEFAULT_FPMODEL} CACHE STRING "Compiler floating point model")
set(SCREAM_FPE ${DEFAULT_FPE} CACHE BOOL "Enable floating point error exception")

# Whether to use XYZ as a method to detect memory usage.
option (SCREAM_ENABLE_GETRUSAGE "Whether getrusage can be used to get memory usage." OFF)
option (SCREAM_ENABLE_STATM "Whether /proc/self/statm can be used to get memory usage." OFF)

# Whether to disable warnings from tpls.
set (SCREAM_DISABLE_TPL_WARNINGS ON CACHE BOOL "")

# Dycore settings
set(DEFAULT_SCREAM_DYNAMICS_DYCORE "NONE")
if (SCREAM_CIME_BUILD AND SCREAM_DYN_TARGET STREQUAL "theta-l_kokkos")
  set (DEFAULT_SCREAM_DYNAMICS_DYCORE "Homme")
endif()

set(SCREAM_DYNAMICS_DYCORE ${DEFAULT_SCREAM_DYNAMICS_DYCORE} CACHE STRING
  "The name of the dycore to be used for dynamics. If NONE, then any code/test requiring dynamics is disabled.")

string(TOUPPER "${SCREAM_DYNAMICS_DYCORE}" SCREAM_DYNAMICS_DYCORE)
if (NOT ${SCREAM_DOUBLE_PRECISION})
  # Homme cannot handle single precision, for now. This causes tests to fail.
  # Fixing this requires adding a config parameter to homme, to switch between
  # single and double. That must be done in the upstream repo (E3SM), before
  # we can support it here.
  # So, for now, if Homme is the requested dyn dycore AND single precision is
  # requested, we disable dynamics, printing a warning.
  if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
    message("WARNING! Homme dycore cannot be used in a Single Precision build. Turning Homme off.")
    set(SCREAM_DYNAMICS_DYCORE "NONE")
  endif()
endif()

# Set the scream base and src directory, to be used across subfolders
set(SCREAM_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SCREAM_SRC_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SCREAM_BIN_DIR  ${CMAKE_CURRENT_BINARY_DIR})

option (EAMXX_ENABLE_PYSCREAM "Whether to enable python interface to eamxx" OFF)
if (EAMXX_ENABLE_PYSCREAM)
  # Pybind11 requires shared libraries
  set (BUILD_SHARED_LIBS ON)
endif()

####################################################################
#                    Packs-related settings                        #
####################################################################

# Determine the main pack size.
set(SCREAM_PACK_SIZE ${DEFAULT_PACK_SIZE} CACHE STRING
  "The number of scalars in a scream::pack::Pack and Mask. Larger packs have good performance on conditional-free loops due to improved caching.")

# Besides the the main pack size, we have a couple of other pack sizes used across EAMxx
# For some routines, SKX may have better performance with pack_size=1
set(SCREAM_SMALL_PACK_SIZE ${SCREAM_PACK_SIZE} CACHE STRING
  "The number of scalars in a scream::pack::SmallPack and SmallMask. Smaller packs can have better performance in loops with conditionals since more of the packs will have masks with uniform value.")
set(SCREAM_POSSIBLY_NO_PACK "${Kokkos_ARCH_SKX}" CACHE BOOL
  "Set possibly-no-pack to this value. You can set it to something else to restore packs on SKX for testing.")

if (SCREAM_POSSIBLY_NO_PACK)
  set (SCREAM_POSSIBLY_NO_PACK_SIZE 1)
else()
  set (SCREAM_POSSIBLY_NO_PACK_SIZE ${SCREAM_PACK_SIZE})
endif ()

function (check_pack_size master_pack_size pack_size name)
  math (EXPR PACK_MODULO "${master_pack_size} % ${pack_size}")
  if ((pack_size GREATER master_pack_size) OR (NOT PACK_MODULO EQUAL 0))
    message (FATAL_ERROR "Invalid '${name}' size of ${pack_size}. Needs to be <= ${master_pack_size} and be a factor of it")
  endif()
endfunction ()

# Checks on pack sizes relative to the master one:
check_pack_size(${SCREAM_PACK_SIZE} ${SCREAM_SMALL_PACK_SIZE} "small pack")
# This one is an internal check, as the user cannot set SCREAM_POSSIBLY_NO_PACK_SIZE now.
check_pack_size(${SCREAM_PACK_SIZE} ${SCREAM_POSSIBLY_NO_PACK_SIZE} "possibly no pack")

####################################################################
#                      Input-data locations                        #
####################################################################

# Handle input root
if (SCREAM_MACHINE AND NOT SCREAM_INPUT_ROOT)
    execute_process(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/query-cime ${SCREAM_MACHINE} DIN_LOC_ROOT
    RESULT_VARIABLE QC_STATUS
    OUTPUT_VARIABLE QC_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if (QC_STATUS EQUAL 0)
    set(DEFAULT_SCREAM_INPUT_ROOT ${QC_OUTPUT})
  endif()
endif()

set(SCREAM_INPUT_ROOT ${DEFAULT_SCREAM_INPUT_ROOT} CACHE PATH "Root of downloaded input files. Should match DIN_LOC_ROOT on CIME machines")
if (NOT SCREAM_INPUT_ROOT)
  string (CONCAT msg
    "  -> Value not found for SCREAM_INPUT_ROOT.\n"
    "     You have three ways for setting this folder. In order of precedence:\n"
    "      - define an env var SCREAM_MACHINE (must be a CIME-supported machine)\n"
    "      - define a cmake cache entry SCREAM_MACHINE (must be a CIME-supported machine)\n"
    "      - define a cmake cache entry SCREAM_INPUT_ROOT.\n"
    "     Ultimately, SCREAM_INPUT_ROOT is the folder where input data will be downloaded.\n")
  message ("${msg}")
  message (FATAL_ERROR "ERROR! Aborting...")
elseif(NOT IS_DIRECTORY ${SCREAM_INPUT_ROOT})
  string (CONCAT msg
    "  -> SCREAM_INPUT_ROOT=${SCREAM_INPUT_ROOT} is not a directory."
    "     You have three ways for setting this folder. In order of precedence:\n"
    "      - define an env var SCREAM_MACHINE (must be a CIME-supported machine)\n"
    "      - define a cmake cache entry SCREAM_MACHINE (must be a CIME-supported machine)\n"
    "      - define a cmake cache entry SCREAM_INPUT_ROOT.\n"
    "     Ultimately, SCREAM_INPUT_ROOT is the folder where input data will be downloaded.\n")
  message ("${msg}")
  message (FATAL_ERROR "ERROR! Aborting...")
else()
  execute_process(COMMAND test -w ${SCREAM_INPUT_ROOT}
    RESULT_VARIABLE res)
  if (NOT res EQUAL 0)
    string (CONCAT msg
      "  -> No write permissions on ${SCREAM_INPUT_ROOT}\n"
      "     If some input files need to be downloaded, this may result in an error."
      "     NOTE: You have three ways for setting this folder. In order of precedence:\n"
      "      - define an env var SCREAM_MACHINE (must be a CIME-supported machine)\n"
      "      - define a cmake cache entry SCREAM_MACHINE (must be a CIME-supported machine)\n"
      "      - define a cmake cache entry SCREAM_INPUT_ROOT.\n"
      "     Ultimately, SCREAM_INPUT_ROOT is the folder where input data will be downloaded.\n")
    message ("${msg}")
  endif()
endif()

set (SCREAM_DATA_DIR ${SCREAM_INPUT_ROOT}/atm/scream CACHE PATH "" FORCE)
set (TOPO_DATA_DIR ${SCREAM_INPUT_ROOT}/atm/cam/topo CACHE PATH "" FORCE)
set (IOP_DATA_DIR ${SCREAM_INPUT_ROOT}/atm/cam/scam/iop CACHE PATH "" FORCE)

####################################################################
#                    Tests-related settings                        #
####################################################################

if (NOT SCREAM_LIB_ONLY)

  # Assuming SCREAM_LIB_ONLY is FALSE (else, no exec is built at all), we need to decide
  # wether to build baseline-related execs, and wether we are generating baselines or
  # comparing against them. These options can help reducing a bit the code that is built
  # when generating baselines or when running memory-check tests (no baselines needed there)
  option(SCREAM_ONLY_GENERATE_BASELINES "Whether building only baselines-related executables" OFF)
  option(SCREAM_ENABLE_BASELINE_TESTS "Whether to run baselines-related tests" ON)
  if (SCREAM_ONLY_GENERATE_BASELINES AND NOT SCREAM_ENABLE_BASELINE_TESTS)
    message (FATAL_ERROR
      "Makes no sense to set SCREAM_ONLY_GENERATE_BASELINES=ON,\n"
      "but set SCREAM_ENABLE_BASELINE_TESTS=OFF.")
  endif()

  set(SCREAM_BASELINES_DIR "UNSET" CACHE PATH "Folder containing baselines data")
  if (SCREAM_ENABLE_BASELINE_TESTS)
    if (NOT EXISTS ${SCREAM_BASELINES_DIR}/data OR NOT IS_DIRECTORY ${SCREAM_BASELINES_DIR}/data)
      string (CONCAT msg
        "Error! Baselines tests enabled, but baseline dir is invalid.\n"
        "  SCREAM_BASELINES_DIR: ${SCREAM_BASELINES_DIR}")
      message ("${msg}")
      message (FATAL_ERROR "Aborting...")
    endif()
  endif()

  # All baselines tests will add a line to the baseline_list file,
  # specifyiing the full name of the baseline they generated.
  # When test-all-scream has to generate new baselines, it will run
  # ctest -L baseline_gen, and then read this file to find out all the
  # baseline files to copy into the baseline directory
  set (SCREAM_TEST_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/data)
  file(MAKE_DIRECTORY ${SCREAM_TEST_OUTPUT_DIR})
  file(TOUCH ${SCREAM_TEST_OUTPUT_DIR}/baseline_list)

  set(SCREAM_TEST_LEVEL_AT "0")
  set(SCREAM_TEST_LEVEL_NIGHTLY "1")
  set(SCREAM_TEST_LEVEL_EXPERIMENTAL "2")

  set(SCREAM_TEST_LEVEL "AT" CACHE STRING "The test level to run. Default is AT. NIGHTLY will run additional tests but is not guaranteed to PASS. EXPERIMENTAL will run even more tests with failures being more likely")

  if (SCREAM_TEST_LEVEL STREQUAL "AT")
    set(SCREAM_TEST_LEVEL ${SCREAM_TEST_LEVEL_AT})
  elseif (SCREAM_TEST_LEVEL STREQUAL "NIGHTLY")
    set(SCREAM_TEST_LEVEL ${SCREAM_TEST_LEVEL_NIGHTLY})
  elseif (SCREAM_TEST_LEVEL STREQUAL "EXPERIMENTAL")
    set(SCREAM_TEST_LEVEL ${SCREAM_TEST_LEVEL_EXPERIMENTAL})
  else()
    message(FATAL_ERROR "Unknown SCREAM_TEST_LEVEL '${SCREAM_TEST_LEVEL}'")
  endif()

  set(SCREAM_TEST_MAX_THREADS ${DEFAULT_MAX_THREADS} CACHE STRING "Upper limit on threads per rank for threaded tests")
  set(SCREAM_TEST_THREAD_INC 1 CACHE STRING "Thread count increment for threaded tests")
  set(SCREAM_TEST_MAX_RANKS ${DEFAULT_MAX_RANKS} CACHE STRING "Upper limit on ranks for mpi tests")
  math(EXPR DEFAULT_MAX_TOTAL_THREADS "${SCREAM_TEST_MAX_RANKS}*${SCREAM_TEST_MAX_THREADS}")
  set(SCREAM_TEST_MAX_TOTAL_THREADS ${DEFAULT_MAX_TOTAL_THREADS} CACHE STRING "Upper limit on nranks*threads for threaded tests")

  # Make sure SCREAM_TEST_MAX_RANKS and SCREAM_TEST_MAX_THREADS do not individually exceed SCREAM_TEST_MAX_TOTAL_THREADS
  if (SCREAM_TEST_MAX_THREADS GREATER ${SCREAM_TEST_MAX_TOTAL_THREADS})
    string(CONCAT msg
           "The requested number of max threads/rank (${SCREAM_TEST_MAX_THREADS}) is larger "
           "than the max total threads (${SCREAM_TEST_MAX_TOTAL_THREADS}). Setting "
           "SCREAM_TEST_MAX_THREADS=${SCREAM_TEST_MAX_TOTAL_THREADS}")
    message(STATUS "${msg}")
    set (SCREAM_TEST_MAX_THREADS ${SCREAM_TEST_MAX_TOTAL_THREADS})
  endif()
  if (SCREAM_TEST_MAX_RANKS GREATER ${SCREAM_TEST_MAX_TOTAL_THREADS})
    string(CONCAT msg
           "The requested number of max ranks (${SCREAM_TEST_MAX_RANKS}) is larger "
           "than the max total threads (${SCREAM_TEST_MAX_TOTAL_THREADS}). Setting "
           "SCREAM_TEST_MAX_RANKS=${SCREAM_TEST_MAX_TOTAL_THREADS}")
    message(STATUS "${msg}")
    set (SCREAM_TEST_MAX_RANKS ${SCREAM_TEST_MAX_TOTAL_THREADS})
  endif()

  # This is a meta-variable, which individual tests can use to set *different* degrees
  # of testing, in terms of resolutions. E.g., for SHORT use 3 timesteps, for MEDIUM use 10,
  # for LONG use 100. It is *completely* up to the test to decide what short, medium, and long mean.
  if (EKAT_ENABLE_COVERAGE OR EKAT_ENABLE_CUDA_MEMCHECK OR EKAT_ENABLE_VALGRIND OR EKAT_ENABLE_COMPUTE_SANITIZER)
    set (SCREAM_TEST_SIZE_DEFAULT SHORT)
  else()
    set (SCREAM_TEST_SIZE_DEFAULT MEDIUM)
  endif()

  set(SCREAM_TEST_SIZE ${SCREAM_TEST_SIZE_DEFAULT} CACHE STRING "The kind of testing to perform: SHORT, MEDIUM, LONG. Only applies to certain tests and is generally used to reduce test length when valgrind/cuda-memcheck are on.")
  set(SCREAM_TEST_VALID_SIZES "SHORT;MEDIUM;LONG" CACHE INTERNAL "List of valid values for SCREAM_TEST_SIZE")

  if (SCREAM_TEST_SIZE STREQUAL "SHORT")
    add_definitions(-DSCREAM_SHORT_TESTS)
  endif()

  enable_testing()
  include(CTest)
endif()

####################################################################
#               Configure all tpls and subfolders                  #
####################################################################

if (DEFINED ENV{SCREAM_FAKE_ONLY})
  # We don't really need to build ekat, but we do need to configure the test-launcher

  # Declare some vars that Ekat would have declared, and may be used later
  option (EKAT_ENABLE_MPI "Whether EKAT requires MPI." ON)
  option (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES "Whether test-launcher should try to manage thread distribution. Requires a ctest resource file to be effective." OFF)
  option (EKAT_ENABLE_VALGRIND "Whether to run tests with valgrind" OFF)
  set(EKAT_VALGRIND_SUPPRESSION_FILE "" CACHE FILEPATH "Use this valgrind suppression file if valgrind is enabled.")
  set (EKAT_ENABLE_GPU False)
  if (Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL)
    set (EKAT_ENABLE_GPU True)
  endif ()

  if (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES)
    set (TEST_LAUNCHER_MANAGE_RESOURCES True)
  else()
    set (TEST_LAUNCHER_MANAGE_RESOURCES False)
  endif()
  if (EKAT_ENABLE_GPU)
    set (TEST_LAUNCHER_ON_GPU True)
  else()
    set (TEST_LAUNCHER_ON_GPU False)
  endif()

  set (EKAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../externals/ekat)

  if (EKAT_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS C)

    # NOTE: may be an overkill, but depending on which FindMPI module is called,
    #       if _FIND_REQUIRED is not checked, we may not get a fatal error
    #       if the required components are not found. So check the _FOUND var.
    if (NOT MPI_C_FOUND)
      message (FATAL_ERROR "EKAT *requires* the C component of MPI to be found")
    endif()

    # We should avoid cxx bindings in mpi; they are already deprecated,
    # and can cause headaches at link time, cause they require -lmpi_cxx
    # (openpmi) or -lmpicxx (mpich) flags.
    include(EkatMpiUtils)
    DisableMpiCxxBindings()
    SetMpiRuntimeEnvVars()
  endif()

  # Configure/install test-launcher to build/install folder
  configure_file(${EKAT_SOURCE_DIR}/bin/test-launcher
                 ${CMAKE_BINARY_DIR}/bin/test-launcher)

  if (EKAT_ENABLE_VALGRIND)
    add_subdirectory(${EKAT_SOURCE_DIR}/valgrind_support
                     ${CMAKE_BINARY_DIR}/externals/ekat/valgrind_support)
  endif()

else()
  include(EkatBuildEkat)
  BuildEkat(PREFIX "SCREAM"
    ENABLE_TESTS OFF
    ENABLE_FPE   ON
    ENABLE_FPE_DEFAULT_MASK OFF)
endif()

# Set compiler-specific flags
include(EkatSetCompilerFlags)
ResetFlags()
SetCommonFlags()
SetProfilingFlags(PROFILER ${EKAT_PROFILING_TOOL} COVERAGE ${EKAT_ENABLE_COVERAGE})

include(EkatMpiUtils)
# We should avoid cxx bindings in mpi; they are already deprecated,
# and can cause headaches at link time, cause they require -lmpi_cxx
# (for openpmi; -lmpicxx for mpich) flag.
DisableMpiCxxBindings()

if (SCREAM_DOUBLE_PRECISION)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(SCREAM_Fortran_FLAGS -real-size 64)
  elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
    set(SCREAM_Fortran_FLAGS -s default32 -eZ)
  else()
    set(SCREAM_Fortran_FLAGS -fdefault-real-8 -fdefault-double-8)
  endif()
endif()

# This must be done using add_definitions because it is used to determine
# whether to include scream_config.h.
add_definitions(-DSCREAM_CONFIG_IS_CMAKE)

# Hooks for testing test-all-scream within test-scripts
if ($ENV{SCREAM_FORCE_BUILD_FAIL})
  add_definitions(-DSCREAM_FORCE_BUILD_FAIL)
endif()
if ($ENV{SCREAM_FORCE_RUN_FAIL})
  add_definitions(-DSCREAM_FORCE_RUN_FAIL)
endif()

# SCREAM_FAKE_ONLY is used to run some tests to ensure the enclosing testing
# scripts work correctly. We are not really interested in building/testing SCREAM.
if (NOT DEFINED ENV{SCREAM_FAKE_ONLY})

  add_subdirectory(tpls)
  add_subdirectory(src)

  if (NOT SCREAM_LIB_ONLY)
    add_subdirectory(tests)
  endif()

  # Generate scream_config.h and scream_config.f
  include (EkatUtils)
  EkatConfigFile(${CMAKE_CURRENT_SOURCE_DIR}/src/scream_config.h.in
                 ${CMAKE_CURRENT_BINARY_DIR}/src/scream_config.h
                 F90_FILE ${CMAKE_CURRENT_BINARY_DIR}/src/scream_config.f)

  # Build any tools in the scripts/ dir.
  add_subdirectory(scripts)

  # Generate scream_config.h and scream_config.f
  include (EkatUtils)
  EkatConfigFile(${CMAKE_CURRENT_SOURCE_DIR}/src/scream_config.h.in
                 ${CMAKE_CURRENT_BINARY_DIR}/src/scream_config.h
                 F90_FILE ${CMAKE_CURRENT_BINARY_DIR}/src/scream_config.f)
else()
  # This is a "fake" build of scream, to stress test the testing scripts.
  # We only need to build something in our tests folder
  add_subdirectory(tests/meta-tests)
endif()

####################################################################
#          Print some scream config settings to screen             #
####################################################################

message ("**************************************************")
message ("*       Summary of EAMxx config settings         *")
message ("**************************************************")

# Shortcut function, to print a variable
function (print_var var)
  message ("${var}: ${${var}}")
endfunction ()

print_var(EAMXX_ENABLE_GPU)
print_var(CUDA_BUILD)
print_var(HIP_BUILD)
print_var(SCREAM_MACHINE)
print_var(SCREAM_DYNAMICS_DYCORE)
print_var(SCREAM_DOUBLE_PRECISION)
print_var(SCREAM_MIMIC_GPU)
print_var(SCREAM_FPE)
print_var(SCREAM_NUM_VERTICAL_LEV)
print_var(SCREAM_PACK_SIZE)
print_var(SCREAM_SMALL_PACK_SIZE)
print_var(SCREAM_POSSIBLY_NO_PACK_SIZE)
print_var(SCREAM_LINK_FLAGS)
print_var(SCREAM_FPMODEL)
print_var(SCREAM_LIB_ONLY)
print_var(SCREAM_TPL_LIBRARIES)
print_var(SCREAM_TEST_MAX_THREADS)
print_var(SCREAM_TEST_THREAD_INC)
print_var(SCREAM_TEST_MAX_RANKS)

message ("**************************************************")
