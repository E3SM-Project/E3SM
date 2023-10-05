#ifndef SCREAM_CONFIG_H
#define SCREAM_CONFIG_H

// If defined, Real is double; if not, Real is float.
#define SCREAM_DOUBLE_PRECISION

// If defined, enable floating point exceptions.
/* #undef SCREAM_FPE */

// The number of scalars in a scream::pack::Pack and Mask.
#define SCREAM_PACK_SIZE 1

// The number of scalars in a scream::pack::SmallPack and SmallMask.
#define SCREAM_SMALL_PACK_SIZE 1

// The number of scalars in a possibly-no-pack. Use this packsize when a routine does better with pksize=1 on some architectures (SKX).
#define SCREAM_POSSIBLY_NO_PACK_SIZE FALSE

// How many levels to use for the vertical grid
#define SCREAM_NUM_VERTICAL_LEV 

// Whether this is a CUDA/HIP build
#define EAMXX_ENABLE_GPU

// Whether scream uses leap years or not
/* #undef SCREAM_HAS_LEAP_YEAR */

// What level of testing we are doing. 0=autotesting, 1=nightly, 2=experimental
#define SCREAM_TEST_LEVEL 

// Whether getrusage can be used to get memory usage
/* #undef SCREAM_ENABLE_GETRUSAGE */
// Whether /proc/self/statm can be used to get memory usage
/* #undef SCREAM_ENABLE_STATM */

#if defined(SCREAM_ENABLE_STATM) || defined(SCREAM_ENABLE_GETRUSAGE)
#define SCREAM_HAS_MEMORY_USAGE
#endif

#define SCREAM_MPI_ON_DEVICE 0

// Data directory for the scream project
#define SCREAM_DATA_DIR "/lus/gila/projects/CSC249ADSE15_CNDA/inputdata/atm/scream"

// Whether or not to run RRTMGP debug checks
/* #undef SCREAM_RRTMGP_DEBUG */

// Whether monolithic kernels are on
#define SCREAM_SMALL_KERNELS

// The sha of the last commit
#define EAMXX_GIT_VERSION ""

// The version of EAMxx
#define EAMXX_VERSION ".."

#endif
