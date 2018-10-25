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
#endif

#ifndef USE_KOKKOS_KERNELS
# define USE_KOKKOS_KERNELS
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
// AVX is not available on CUDA, so make certain this is 0
# ifdef HOMMEXX_AVX_VERSION
#  undef HOMMEXX_AVX_VERSION
# endif
# define HOMMEXX_AVX_VERSION 0
#elif ! defined HOMMEXX_AVX_VERSION
# define HOMMEXX_CUDA_MIN_WARP_PER_TEAM 1
# define HOMMEXX_CUDA_MAX_WARP_PER_TEAM 1
# if defined __AVX512F__
#  define HOMMEXX_AVX_VERSION 512
# elif defined __AVX2__
#  define HOMMEXX_AVX_VERSION 2
# elif defined __AVX__
#  define HOMMEXX_AVX_VERSION 1
# else
#  define HOMMEXX_AVX_VERSION 0
#  ifndef HOMMEXX_VECTOR_SIZE
#   define HOMMEXX_VECTOR_SIZE 8
#  endif
# endif
#endif

#endif
