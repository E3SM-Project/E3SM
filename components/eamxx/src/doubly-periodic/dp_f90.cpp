#include "dp_f90.hpp"
#include "physics_constants.hpp"

#include "ekat/ekat_assert.hpp"

using scream::Real;
using scream::Int;

extern "C" {
}

namespace scream {
namespace dp {

void dp_init(Int nlev, bool use_fortran, bool force_reinit) {
  static bool is_init = false;
  if (!is_init || force_reinit) {
    is_init = true;
  }
}

} // namespace dp
} // namespace scream
