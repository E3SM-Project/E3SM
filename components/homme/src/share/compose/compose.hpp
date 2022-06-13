#ifndef INCLUDE_COMPOSE_HPP
#define INCLUDE_COMPOSE_HPP

#ifdef HAVE_CONFIG_H
# include "config.h.c"
#endif

#include <Kokkos_Core.hpp>

#ifdef _OPENMP
# include <omp.h>
#endif

// Options

#ifdef NDEBUG
//# undef NDEBUG
#endif

#ifndef NDEBUG
# define COMPOSE_BOUNDS_CHECK
#endif

//#define COMPOSE_TIMERS
#ifdef COMPOSE_TIMERS
# include "gptl.h"
#endif

// Look for MPI-related memory leaks.
//#define COMPOSE_DEBUG_MPI

#if ! defined COMPOSE_PORT
# if defined HORIZ_OPENMP
#  define COMPOSE_HORIZ_OPENMP
# endif
# if defined COLUMN_OPENMP
#  define COMPOSE_COLUMN_OPENMP
# endif
#endif

#if defined COMPOSE_PORT
# if defined COMPOSE_HORIZ_OPENMP || defined COMPOSE_COLUMN_OPENMP
"This should not happen."
# endif
# ifndef KOKKOS_ENABLE_CUDA
// Mimic GPU threading on host to debug race conditions on a regular CPU.
//#  define COMPOSE_MIMIC_GPU
# endif
# if defined COMPOSE_MIMIC_GPU || defined KOKKOS_ENABLE_CUDA
// If defined, then certain buffers need explicit mirroring and copying.
#  define COMPOSE_PORT_SEPARATE_VIEWS
// If defined, do pass1 routines on host. This is for performance checking.
//#  define COMPOSE_PACK_NOSCAN
# endif
#endif

#if defined COMPOSE_BOUNDS_CHECK && defined NDEBUG
# pragma message "NDEBUG but COMPOSE_BOUNDS_CHECK"
#endif

#endif
