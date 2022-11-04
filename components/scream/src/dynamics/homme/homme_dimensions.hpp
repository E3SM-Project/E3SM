#ifndef SCREAM_HOMMEXX_DIMENSIONS_HPP
#define SCREAM_HOMMEXX_DIMENSIONS_HPP

// Get hommexx dimensions header
#include "Dimensions.hpp"

namespace scream {

// Until whenever CUDA supports constexpr properly
#ifdef EAMXX_ENABLE_GPU

#define HOMMEXX_QSIZE_D             QSIZE_D
#define HOMMEXX_NP                  NP

#define HOMMEXX_NUM_TIME_LEVELS     NUM_TIME_LEVELS
#define HOMMEXX_Q_NUM_TIME_LEVELS   Q_NUM_TIME_LEVELS

#define HOMMEXX_NUM_PHYSICAL_LEV    NUM_PHYSICAL_LEV
#define HOMMEXX_NUM_LEV             NUM_LEV
#define HOMMEXX_NUM_LEV_P           NUM_LEV_P
#define HOMMEXX_NUM_INTERFACE_LEV   NUM_INTERFACE_LEV

#define HOMMEXX_PACK_SIZE           VECTOR_SIZE

#else

constexpr int HOMMEXX_QSIZE_D            = QSIZE_D;
constexpr int HOMMEXX_NP                 = NP;

constexpr int HOMMEXX_NUM_PHYSICAL_LEV   = Homme::NUM_PHYSICAL_LEV;
constexpr int HOMMEXX_NUM_TIME_LEVELS    = Homme::NUM_TIME_LEVELS;
constexpr int HOMMEXX_Q_NUM_TIME_LEVELS  = Homme::Q_NUM_TIME_LEVELS;

constexpr int HOMMEXX_NUM_LEV            = Homme::NUM_LEV;
constexpr int HOMMEXX_NUM_LEV_P          = Homme::NUM_LEV_P;
constexpr int HOMMEXX_NUM_INTERFACE_LEV  = Homme::NUM_INTERFACE_LEV;

constexpr int HOMMEXX_PACK_SIZE          = Homme::VECTOR_SIZE;

#endif

} // namespace scream

#endif // SCREAM_HOMME_DIMENSIONS_HPP
