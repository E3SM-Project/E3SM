#include "dp_f90.hpp"
#include "physics_constants.hpp"

#include "ekat/ekat_assert.hpp"

using scream::Int;

extern "C" {

void init_time_level_c (const int& nm1, const int& n0, const int& np1,
                        const int& nstep, const int& nstep0);

}

namespace scream {
namespace dp {

void dp_init(const bool force_reinit) {
  static bool is_init = false;
  if (!is_init || force_reinit) {
    init_time_level_c(10, 3, 11, 5, 4);
    is_init = true;
  }
}

} // namespace dp
} // namespace scream
