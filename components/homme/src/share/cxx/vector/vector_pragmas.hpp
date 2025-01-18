/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_VECTOR_PRAGMAS_HPP
#define HOMMEXX_VECTOR_PRAGMAS_HPP

#include "Config.hpp"

#if defined(__INTEL_COMPILER) || defined(__INTEL_CLANG_COMPILER) || defined(__INTEL_LLVM_COMPILER)

#define VECTOR_IVDEP_LOOP _Pragma("ivdep")
#define ALWAYS_VECTORIZE_LOOP _Pragma("vector always")
#define VECTOR_SIMD_LOOP _Pragma("omp simd")
#if HOMMEXX_VECTOR_SIZE == 1
# define VECTOR_SIMD_LOOP
#else
# define VECTOR_SIMD_LOOP _Pragma("omp simd")
#endif

#elif defined(__GNUG__) && !defined(__NVCC__)
#if(__GNUG__ == 4 && __GNUC_MINOR__ >= 9) || __GNUG__ > 4

#define VECTOR_IVDEP_LOOP _Pragma("GCC ivdep")
#define ALWAYS_VECTORIZE_LOOP _Pragma("GCC vector always")
#define VECTOR_SIMD_LOOP _Pragma("GCC ivdep")

#else // __GNUG__ ...
#pragma message( \
    "G++ <4.9 Does not support vectorization pragmas")
#define VECTOR_IVDEP_LOOP
#define ALWAYS_VECTORIZE_LOOP
#define VECTOR_SIMD_LOOP

#define HOMMEXX_NO_VECTOR_PRAGMAS

#endif // __GNUG__ ...

#else // defined(__INTEL_COMPILER) / defined(__GNUG__) ...
// fail gracefully when the compiler vectorization pragmas are unknown

#define VECTOR_IVDEP_LOOP
#define ALWAYS_VECTORIZE_LOOP
#define VECTOR_SIMD_LOOP

#define HOMMEXX_NO_VECTOR_PRAGMAS

#endif // defined(__INTEL_COMPILER) / defined(__GNUG__) ...

#endif // HOMMEXX_VECTOR_PRAGMAS_HPP
