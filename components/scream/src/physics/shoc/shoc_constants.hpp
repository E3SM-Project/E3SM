#ifndef SHOC_CONSTANTS_HPP
#define SHOC_CONSTANTS_HPP

namespace scream {
  namespace shoc {

    /*
     * Mathematical constants used by shoc.
     */

template <typename Scalar>
struct Constants
    {
      static constexpr Scalar mintke     = 0.0004;  // Minimum TKE [m2/s2]
      static constexpr Scalar maxtke     = 50.0;    // Maximum TKE [m2/s2]
      static constexpr Scalar maxlen     = 20000.0; // Upper limit for mixing length [m]
      static constexpr Scalar length_fac = 0.5;     // Mixing length scaling parameter
    };

  } // namespace shoc
} // namespace scream

#endif
