#ifndef INCLUDE_SCREAM_SHOC_IC_CASES_HPP
#define INCLUDE_SCREAM_SHOC_IC_CASES_HPP

#include "shoc_data.hpp"

namespace scream {
namespace shoc {
namespace ic {

struct Factory {
  enum IC { standard };
  static FortranData::Ptr create(IC ic, Int shcol = 1, Int nlev = 72,
                                 Int num_qtracers = 1);
};

} // namespace ic
} // namespace p3
} // namespace scream

#endif
