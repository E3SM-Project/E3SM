#ifndef DP_CONSTANTS_HPP
#define DP_CONSTANTS_HPP

namespace scream {
namespace dp {

/*
 * Mathematical constants used by dp.
 */

template <typename Scalar>
struct Constants
{
  static constexpr Scalar iop_nudge_tq_low = 1050;
  static constexpr Scalar iop_nudge_tq_high = 0;
  static constexpr Scalar iop_nudge_tscale = 10800;
};

} // namespace dp
} // namespace scream

#endif
