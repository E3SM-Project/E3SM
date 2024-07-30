#ifndef HOMMEXX_CONFIG_H
#define HOMMEXX_CONFIG_H

// Whether the CUDA exec space has been selected
/* #undef HOMMEXX_CUDA_SPACE */

// Whether the OpenMP exec space has been selected
/* #undef HOMMEXX_OPENMP_SPACE */

// Whether the Threads exec space has been selected
/* #undef HOMMEXX_THREADS_SPACE */

// Whether the Serial exec space has been selected
/* #undef HOMMEXX_SERIAL_SPACE */

// Whether the Default Kokkos exec space has been selected
#define HOMMEXX_DEFAULT_SPACE

// Whether the debug parts of cxx code should be compiled or not
/* #undef HOMMEXX_DEBUG */

// Whether the MPI operations have to be performed directly on the device
#define HOMMEXX_MPI_ON_DEVICE 1

/* #undef HOMMEXX_CUDA_SHARE_BUFFER */

// Minimum and maximum number of warps to provide to a team
#define HOMMEXX_CUDA_MIN_WARP_PER_TEAM 8
#define HOMMEXX_CUDA_MAX_WARP_PER_TEAM 16

// User-defined VECTOR_SIZE
#define HOMMEXX_VECTOR_SIZE 8

#define HOMMEXX_SHA1 "d2ec49bd059355131679652011ad8f9c9745cd3c"

#endif // HOMMEXX_CONFIG_H
