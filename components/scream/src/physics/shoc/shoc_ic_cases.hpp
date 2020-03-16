#ifndef INCLUDE_SCREAM_SHOC_IC_CASES_HPP
#define INCLUDE_SCREAM_SHOC_IC_CASES_HPP

#include "shoc_f90.hpp"

namespace scream {
namespace shoc {
namespace ic {

struct Factory {
  static FortranData::Ptr create(Int shcol = 1, Int nlev = 128);
};

} // namespace ic
} // namespace p3
} // namespace scream

#endif
