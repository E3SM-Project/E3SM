/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_DIMENSIONS_HPP
#define HOMMEXX_DIMENSIONS_HPP

#include <Kokkos_Core.hpp>

#include "Config.hpp"

namespace Homme {

// Until whenever CUDA supports constexpr properly
#ifdef CUDA_BUILD

  #ifdef CAM
    #define QSIZE_D PCNST
  #endif

  #define VECTOR_SIZE         1
  #define VECTOR_END          (VECTOR_SIZE-1)

  #define NUM_PHYSICAL_LEV    PLEV
  #define NUM_TIME_LEVELS     3
  #define Q_NUM_TIME_LEVELS   2

  #define NUM_LEV             NUM_PHYSICAL_LEV
  #define NUM_LEV_P           (NUM_LEV + 1)
  #define NUM_INTERFACE_LEV   NUM_LEV_P

  // Note: on CUDA these are somewhat an overkill, since VECTOR_SIZE=1. However, if GPU change,
  //       and they start to support 2+ double precision vector units, then this would make sense
  #define LAST_LEV                (NUM_LEV-1)
  #define LAST_LEV_P              (NUM_LEV_P-1)
  #define LAST_MIDPOINT_VEC_IDX   ((NUM_PHYSICAL_LEV + VECTOR_SIZE - 1) % VECTOR_SIZE)
  #define LAST_INTERFACE_VEC_IDX  ((LAST_MIDPOINT_VEC_IDX+1) % VECTOR_SIZE)
#else

  #ifdef CAM
    static constexpr const int QSIZE_D = PCNST;
  #endif

  #if   (HOMMEXX_AVX_VERSION == 0)
    // Vector<VectorTag<SIMD<T, SpT>, l> > can use this for good results
    // on, e.g., Power9, where AVX doesn't exist.
    static constexpr const int VECTOR_SIZE = HOMMEXX_VECTOR_SIZE;
  #elif (HOMMEXX_AVX_VERSION == 1 || HOMMEXX_AVX_VERSION == 2)
    static constexpr const int VECTOR_SIZE = 4;
  #elif (HOMMEXX_AVX_VERSION == 512)
    static constexpr const int VECTOR_SIZE = 8;
  #endif
  static constexpr int VECTOR_END = VECTOR_SIZE-1;

  static_assert(VECTOR_SIZE>0, "Error: VECTOR_SIZE=0!");

  static constexpr const int NUM_PHYSICAL_LEV = PLEV;
  static constexpr const int NUM_LEV =
      (NUM_PHYSICAL_LEV + VECTOR_SIZE - 1) / VECTOR_SIZE;

  static constexpr const int NUM_INTERFACE_LEV = NUM_PHYSICAL_LEV + 1;
  static constexpr const int NUM_LEV_P =
      (NUM_INTERFACE_LEV + VECTOR_SIZE - 1) / VECTOR_SIZE;

  static constexpr const int NUM_TIME_LEVELS = 3;
  static constexpr const int Q_NUM_TIME_LEVELS = 2;

  // Note: the last midpoint is in the pack NUM_LEV, at position (NUM_PHYSICAL_LEV-1)%VECTOR_SIZE;
  //       the last interface is in the pack NUM_LEV_P, one position after the last midpoint (mod VECTOR_SIZE).
  //       This is true regardless of whether NUM_LEV_P>NUM_LEV or not, so one formula works for all cases:
  //         - if NUM_LEV=NUM_LEV_P, then last_midpoint_vec_idx<VECTOR_SIZE-1, so the formula clearly works
  //         - otherwise, last_midpoint_vec_idx=VECTOR_SIZE-1, so last_interface_vec_idx=VECTOR_SIZE%VECTOR_SIZE=0, which is correct.
  static constexpr int LAST_LEV               = NUM_LEV-1;
  static constexpr int LAST_LEV_P             = NUM_LEV_P-1;
  static constexpr int LAST_MIDPOINT_VEC_IDX  = (NUM_PHYSICAL_LEV + VECTOR_SIZE - 1) % VECTOR_SIZE;
  static constexpr int LAST_INTERFACE_VEC_IDX = (NUM_INTERFACE_LEV + VECTOR_SIZE - 1) % VECTOR_SIZE;
#endif // CUDA_BUILD


template<int ARRAY_LENGTH>
struct PackInfo {
  static constexpr int NumPacks       = (ARRAY_LENGTH + VECTOR_SIZE - 1) / VECTOR_SIZE;
  static constexpr int LastPackVecEnd = (ARRAY_LENGTH + VECTOR_SIZE - 1) % VECTOR_SIZE;
};

} // namespace TinMan

#endif // HOMMEXX_DIMENSIONS_HPP
