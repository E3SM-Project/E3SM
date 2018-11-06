/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/


#ifndef _PROFILING_HPP_
#define _PROFILING_HPP_

#include <Types.hpp>

#include "gptl.h"

#if defined(HOMMEXX_CUDA_SPACE) || \
  (defined(HOMMEXX_DEFAULT_SPACE) && defined(KOKKOS_ENABLE_CUDA)) // Can't use GPTL timers on CUDA
#define start_timer(name) {}
#define stop_timer(name) {}
#else
#define start_timer(name) { GPTLstart(name); }
#define stop_timer(name) { GPTLstop(name); }
#endif

#ifdef VTUNE_PROFILE
#include <ittnotify.h>

#define profiling_resume __itt_resume
#define profiling_pause __itt_pause

#elif defined(CUDA_PROFILE) // VTUNE_PROFILE
#include <cuda_profiler_api.h>

#define profiling_resume cuProfilerStart
#define profiling_pause cuProfilerStop

#elif defined(GPROF_PROFILE) // CUDA_PROFILE

/* Use undocumented GProf methods */
#pragma message "GProf support for this is undocumented and subject to change, this may fail as a result"

extern "C" void moncontrol(int);

#define profiling_resume() moncontrol(1)
#define profiling_pause() moncontrol(0)

#else // GPROF_PROFILE

#define profiling_resume() {}
#define profiling_pause() {}

#endif

#endif // _PROFILING_HPP_
