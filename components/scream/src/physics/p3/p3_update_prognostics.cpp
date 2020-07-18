#include "p3_update_prognostics_impl.hpp"
#include "ekat/scream_types.hpp"

namespace scream {
namespace p3 {

  /*
   * Explicit instantiation for doing p3 update prognostics functions on Reals using the
   * default device.
   */

  template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
