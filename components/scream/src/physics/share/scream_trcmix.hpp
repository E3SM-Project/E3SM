#ifndef SCREAM_TRCMIX_HPP
#define SCREAM_TRCMIX_HPP

#include "share/scream_types.hpp"

//
// Bridge function to call trcmix from CXX
//

namespace scream {
namespace physics {

void trcmix(const char* name, Int ncol, Int pcols, Int pver, Real* clat, Real* pmid, Real* q);

} // namespace physics
} // namespace scream

#endif
