/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_CONFIG_HPP
#define HOMMEXX_CONFIG_HPP

#ifdef HOMMEXX_CONFIG_IS_CMAKE
# include "Hommexx_config.h"
# ifdef HAVE_CONFIG_H
#  include "config.h.c"
# endif
#else
// Establish a good candidate vector size for eam builds
# ifdef CUDA_BUILD
#  define HOMMEXX_VECTOR_SIZE 1
# else
#  define HOMMEXX_VECTOR_SIZE 8
# endif
#endif

#if ! defined HOMMEXX_CUDA_SPACE && ! defined HOMMEXX_OPENMP_SPACE && ! defined HOMMEXX_THREADS_SPACE && ! defined HOMMEXX_SERIAL_SPACE
# define HOMMEXX_DEFAULT_SPACE
#endif

#ifndef HOMMEXX_MPI_ON_DEVICE
# define HOMMEXX_MPI_ON_DEVICE 1
#endif

#include <Kokkos_Core.hpp>

#ifdef KOKKOS_ENABLE_CUDA
# ifndef HOMMEXX_CUDA_MIN_WARP_PER_TEAM
#  define HOMMEXX_CUDA_MIN_WARP_PER_TEAM 8
# endif
# ifndef HOMMEXX_CUDA_MAX_WARP_PER_TEAM
#  define HOMMEXX_CUDA_MAX_WARP_PER_TEAM 16
# endif
#elif !defined(HOMMEXX_CONFIG_IS_CMAKE)
# define HOMMEXX_CUDA_MIN_WARP_PER_TEAM 1
# define HOMMEXX_CUDA_MAX_WARP_PER_TEAM 1
#endif

#endif // HOMMEXX_CONFIG_HPP
