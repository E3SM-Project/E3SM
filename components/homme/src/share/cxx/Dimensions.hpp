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
#endif // CUDA_BUILD

template<int PHYSICAL_LENGTH>
struct ColInfo {
private:
  template<int LENGTH>
  struct Helper {
    static constexpr int NumPacks   = (LENGTH + VECTOR_SIZE - 1) / VECTOR_SIZE;
    static constexpr int LastVecLen = LENGTH - (NumPacks-1)*VECTOR_SIZE;

    static constexpr int LastPack   = NumPacks - 1;
    static constexpr int LastVecEnd = LastVecLen - 1;
  };
public:
  static constexpr int NumPacks   = Helper<PHYSICAL_LENGTH>::NumPacks;
  static constexpr int LastPack   = Helper<PHYSICAL_LENGTH>::LastPack;
  static constexpr int LastVecLen = Helper<PHYSICAL_LENGTH>::LastVecLen;
  static constexpr int LastVecEnd = Helper<PHYSICAL_LENGTH>::LastVecEnd;

  static constexpr int PenultLevPackIdx = Helper<PHYSICAL_LENGTH-1>::NumPacks;
  static constexpr int PenultLevVecIdx  = Helper<PHYSICAL_LENGTH-1>::LastVecEnd;
};

} // namespace TinMan

#endif // HOMMEXX_DIMENSIONS_HPP
