#ifndef SCREAM_COSP_CXX_FUNCTIONS_HPP
#define SCREAM_COSP_CXX_FUNCTIONS_HPP

#include "physics/cosp/cosp_constants.hpp"

#include "share/core/eamxx_types.hpp"

#include <cstdint>

namespace scream {
namespace cosp {

/*
 * Functions is a stateless struct encapsulating the C++/Kokkos port of the COSP
 * pieces used by EAMxx (subcolumn generator, optics, and the ISCCP, MISR and MODIS
 * simulators). Unlike P3/SHOC this port does not use packs: COSP's control flow is
 * per subcolumn and data dependent, and the priority is staying BFB with the Fortran.
 *
 * The functions are header-only (no ETI) for now. Each one mirrors a Fortran routine,
 * named in the comment above it, and is tested BFB against it in tests/.
 */

template <typename ScalarT, typename DeviceT>
struct Functions
{
  using Scalar = ScalarT;
  using Device = DeviceT;

  using KT         = ekat::KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;
  using C          = cosp::Constants<Scalar>;

  template <typename S> using view_1d = typename KT::template view_1d<S>;
  template <typename S> using view_2d = typename KT::template view_2d<S>;
  template <typename S> using view_3d = typename KT::template view_3d<S>;

  // ---------------------------------------------------------------------------------------
  // Random number generator (mo_rng.F90)

  // State of one random number stream (mo_rng.F90, type rng_state)
  using RngState = std::int32_t;

  // get_rng_1: return the next number in [0,1) and advance the stream
  KOKKOS_INLINE_FUNCTION
  static Scalar get_rng(RngState& seed);

  // ---------------------------------------------------------------------------------------
  // Optics (cosp_optics.F90, and the MODIS fits in modis_simulator.F90)

  // cosp_simulator_optics (2D and 3D): pick the stratiform or convective value
  // according to the subcolumn cloud type
  KOKKOS_INLINE_FUNCTION
  static Scalar simulator_optics(const Scalar flag, const Scalar var_cv, const Scalar var_ls);

  // modis_optics_partition: partition optical thickness into liquid and ice parts
  KOKKOS_INLINE_FUNCTION
  static void modis_optics_partition(const Scalar cloud_water, const Scalar cloud_ice,
                                     const Scalar water_size, const Scalar ice_size,
                                     const Scalar tau, Scalar& tau_liq, Scalar& tau_ice);

  // modis_optics: liquid fraction, asymmetry parameter and single scattering albedo
  KOKKOS_INLINE_FUNCTION
  static void modis_optics(const Scalar tau_liq, const Scalar size_liq,
                           const Scalar tau_ice, const Scalar size_ice,
                           Scalar& frac_liq, Scalar& g, Scalar& w0);

  // get_g_nir / get_ssa_nir: polynomial fits of g and w0 in MODIS band 7 vs size
  KOKKOS_INLINE_FUNCTION
  static Scalar get_g_nir(const int phase, const Scalar re);
  KOKKOS_INLINE_FUNCTION
  static Scalar get_ssa_nir(const int phase, const Scalar re);

  // ---------------------------------------------------------------------------------------
  // Statistics (cosp_stats.F90)

  // hist2D: joint histogram of npts samples; bin1/bin2 hold nbin1+1/nbin2+1 edges, and
  // hist is indexed (ibin1, ibin2). Runs serially, intended to be called once per column.
  template <typename V1, typename V2, typename B1, typename B2, typename H>
  KOKKOS_INLINE_FUNCTION
  static void hist2d(const int npts, const V1& var1, const V2& var2,
                     const int nbin1, const B1& bin1, const int nbin2, const B2& bin2,
                     const H& hist);
};

} // namespace cosp
} // namespace scream

#include "physics/cosp/impl/cosp_rng_impl.hpp"
#include "physics/cosp/impl/cosp_optics_impl.hpp"
#include "physics/cosp/impl/cosp_stats_impl.hpp"

#endif // SCREAM_COSP_CXX_FUNCTIONS_HPP
