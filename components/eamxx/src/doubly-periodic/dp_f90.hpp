#ifndef SCREAM_DP_F90_HPP
#define SCREAM_DP_F90_HPP

#include "share/scream_types.hpp"

#include <memory>
#include <vector>

namespace scream {
namespace dp {

// Initialize DP. This is only for standalone DP testing.
void dp_init(const bool force_reinit=false);

}  // namespace dp
}  // namespace scream

#endif
