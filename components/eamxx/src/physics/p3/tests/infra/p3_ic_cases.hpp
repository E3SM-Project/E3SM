#ifndef INCLUDE_SCREAM_P3_IC_CASES_HPP
#define INCLUDE_SCREAM_P3_IC_CASES_HPP

#include "p3_data.hpp"

namespace scream {
namespace p3 {
namespace ic {

P3Data::Ptr make_mixed(Int ncol, Int nlev = 128);

struct Factory {
  enum IC { mixed };

  static P3Data::Ptr create(IC ic, Int ncol = 1, Int nlev = 128);
};

} // namespace ic
} // namespace p3
} // namespace scream

#endif
