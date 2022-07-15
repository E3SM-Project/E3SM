#include "scream_trcmix.hpp"

using scream::Real;
using scream::Int;

extern "C" {

void trcmix_cf2(const char* name, Int ncol, Int pcols, Int pver, Real* clat, Real* pmid, Real* q);

}

namespace scream {
namespace physics {

void trcmix(const char* name, Int ncol, Int pcols, Int pver, Real* clat, Real* pmid, Real* q)
{
  //d.transpose<ekat::TransposeDirection::c2f>();
  trcmix_cf2(name, ncol, pcols, pver, clat, pmid, q);
  //d.transpose<ekat::TransposeDirection::f2c>();
}

} // namespace physics
} // namespace scream
