#ifndef SCREAM_DP_F90_HPP
#define SCREAM_DP_F90_HPP

#include "share/scream_types.hpp"

#include <memory>
#include <vector>

namespace scream {
namespace dp {

// Initialize DP with the given number of levels.
void dp_init(Int nlev, bool use_fortran=false, bool force_reinit=false);

}  // namespace dp
}  // namespace scream

#endif
