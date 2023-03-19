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
    static constexpr Scalar mintke         = 0.0004;       // Minimum TKE [m2/s2]
    static constexpr Scalar maxtke         = 50.0;         // Maximum TKE [m2/s2]
    static constexpr Scalar minlen         = 20.0;         // Lower limit for mixing length [m]
    static constexpr Scalar maxlen         = 20000.0;      // Upper limit for mixing length [m]
    static constexpr Scalar maxiso         = 20000.0;      // Upper limit for isotropy time scale [s]
    static constexpr Scalar length_fac     = 0.5;          // Mixing length scaling parameter
    static constexpr Scalar c_diag_3rd_mom = 7.0;          // Coefficient for diag third moment parameters
    static constexpr Scalar w3clip         = 1.2;          // Third moment of vertical velocity
    static constexpr Scalar ustar_min      = 0.01;         // Minimum surface friction velocity
    static constexpr Scalar largeneg       = -99999999.99; // Large negative value used for linear_interp threshold
    static constexpr bool   dothetal_skew  = false;        // Allow temperature skewness to be independent of moisture variance
    static constexpr Scalar pblmaxp        = 4e4;          // PBL max depth in pressure units
  };

  } // namespace shoc
} // namespace scream

#endif
