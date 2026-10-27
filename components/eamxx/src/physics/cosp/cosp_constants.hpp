#ifndef SCREAM_COSP_CONSTANTS_HPP
#define SCREAM_COSP_CONSTANTS_HPP

#include "share/core/eamxx_types.hpp"

namespace scream {
namespace cosp {

/*
 * Constants for the C++ port of COSP.
 *
 * A note on literals: the Fortran COSP library is compiled WITHOUT -fdefault-real-8
 * (the EAMxx Fortran flags only reach targets linking scream_share), so any literal
 * without a _wp suffix in COSP is single precision, even when it is assigned to, or
 * compared against, a real(wp) variable. To stay BFB with the Fortran, such literals
 * are written here as float literals cast to Scalar, e.g. Scalar(0.7490f).
 */

template <typename Scalar>
struct Constants
{
  // ---------------------------------------------------------------------------------------
  // Hydrometeor flags used in frac_out (subcolumn cloud type), see scops.F90
  static constexpr int I_LSC = 1; // Large-scale (stratiform) cloud
  static constexpr int I_CVC = 2; // Convective cloud

  // ---------------------------------------------------------------------------------------
  // MODIS simulator parameters, set in cosp_modis_init (cosp_modis_interface.F90)
  static constexpr int    num_trial_res = 15;
  static constexpr Scalar re_water_min  = 4.0;
  static constexpr Scalar re_water_max  = 30.0;
  static constexpr Scalar re_ice_min    = 5.0;
  static constexpr Scalar re_ice_max    = 90.0;
  static constexpr int phaseIsNone         = 0;
  static constexpr int phaseIsLiquid       = 1;
  static constexpr int phaseIsIce          = 2;
  static constexpr int phaseIsUndetermined = 3;
  static constexpr Scalar min_OpticalThickness          = 0.3;     // Minimum detectable tau
  static constexpr Scalar CO2Slicing_PressureLimit      = 70000.0; // Pa
  static constexpr Scalar CO2Slicing_TauLimit           = 1.0;
  static constexpr Scalar phase_TauLimit                = 1.0;
  static constexpr Scalar phaseDiscrimination_Threshold = 0.7;
  static constexpr Scalar re_fill                       = -999.0;

  // Ice density used in MODIS phase partitioning (cosp_optics.F90, 0.93_wp)
  static constexpr Scalar ice_density = 0.93;

  // ---------------------------------------------------------------------------------------
  // Missing value (cosp_config.F90; a single precision literal)
  static constexpr Scalar R_UNDEF = Scalar(-1.0E30f);

  // ---------------------------------------------------------------------------------------
  // Histogram bin boundaries (cosp_config.F90): real(wp) parameters initialized from
  // single precision literals
  static constexpr int ntau = 7;  // optical depth bins
  static constexpr int npres = 7; // cloud-top pressure bins
  static constexpr int nhgt = 16; // cloud-top height bins

  // tau_binBounds, ntau+1 values
  KOKKOS_INLINE_FUNCTION
  static Scalar tau_binBounds (const int i) {
    constexpr float v[ntau+1] = {0.0f, 0.3f, 1.3f, 3.6f, 9.4f, 23.f, 60.f, 10000.f};
    return Scalar(v[i]);
  }

  // pres_binBounds (hPa), npres+1 values
  KOKKOS_INLINE_FUNCTION
  static Scalar pres_binBounds (const int i) {
    constexpr float v[npres+1] = {0.f, 180.f, 310.f, 440.f, 560.f, 680.f, 800.f, 10000.f};
    return Scalar(v[i]);
  }

  // hgt_binBounds (km), nhgt+1 values
  KOKKOS_INLINE_FUNCTION
  static Scalar hgt_binBounds (const int i) {
    constexpr float v[nhgt+1] = {-.99f, 0.f, 0.5f, 1.f, 1.5f, 2.f, 2.5f, 3.f, 4.f, 5.f, 7.f,
                                 9.f, 11.f, 13.f, 15.f, 17.f, 99.f};
    return Scalar(v[i]);
  }

  // ---------------------------------------------------------------------------------------
  // MODIS simulator histogram sizes (cosp_config.F90)
  static constexpr int numMODISTauBins  = ntau;
  static constexpr int numMODISPresBins = npres;

  // ---------------------------------------------------------------------------------------
  // MISR simulator (MISR_simulator.F90)
  static constexpr int numMISRTauBins = ntau;
  static constexpr int numMISRHgtBins = nhgt;

  // tauchk = -1.*log(0.9999999), evaluated in single precision at compile time by the
  // Fortran compiler. Hard-coded, since a runtime logf can differ from it by one ulp.
  static constexpr Scalar misr_tauchk = Scalar(0x1.000002p-23);

  // ---------------------------------------------------------------------------------------
  // ISCCP simulator (icarus.F90)
  static constexpr int numISCCPTauBins  = ntau;
  static constexpr int numISCCPPresBins = npres;

  // tauchk = -1._wp*log(0.9999999_wp), as evaluated by the Fortran compiler
  static constexpr Scalar isccp_tauchk = Scalar(0x1.ad7f2b1049b9fp-24);
  static constexpr Scalar isccp_taumin = 0.3;
  static constexpr Scalar pstd         = 1013250.0; // Mean sea-level pressure (Pa)
  static constexpr Scalar isccp_t0     = 296.0;     // Mean surface temperature (K)

  // Physical constants used by ICARUS (cosp_constants.F90)
  static constexpr Scalar amw  = 18.01534;   // Water [g/mol]
  static constexpr Scalar amd  = 28.9644;    // Dry air [g/mol]
  static constexpr Scalar avo  = 6.023E23;   // Avogadro constant [1/mol]
  static constexpr Scalar grav = 9.806650;   // Gravitational acceleration [m/s2]
};

} // namespace cosp
} // namespace scream

#endif // SCREAM_COSP_CONSTANTS_HPP
