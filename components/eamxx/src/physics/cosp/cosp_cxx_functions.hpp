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
  // Subcolumn generator (scops.F90)

  // RNG seed of one point from its surface pressure, as subsample_and_optics sets it
  KOKKOS_INLINE_FUNCTION
  static RngState scops_seed(const int npoints, const Scalar psfc, const Scalar psfc_min,
                             const Scalar psfc_max);

  // scops for one point; cc and conv by level, frac_out is (ncol, nlev), threshold is
  // scratch of size ncol
  template <typename VC, typename VV, typename VF, typename VT>
  KOKKOS_INLINE_FUNCTION
  static void scops(const int nlev, const int ncol, RngState seed, RngState seed_last,
                    const VC& cc, const VV& conv, int overlap, const VF& frac_out,
                    const VT& threshold);

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
  // ISCCP simulator (icarus.F90), for isccp_top_height = 1, isccp_top_height_direction = 2

  // Per-point quantities shared by all the subcolumns of a point
  struct IcarusPoint {
    Scalar ptrop, attrop, atmax;  // tropopause pressure and temperature, max T below it
    int    itrop;                 // tropopause level index (0-based)
    Scalar fluxtop_clrsky;        // clear-sky TOA flux
    Scalar meantbclr;             // clear-sky brightness temperature
    Scalar btcmin;                // flux of the minimum cloud brightness temperature
  };

  KOKKOS_INLINE_FUNCTION
  static Scalar icarus_dem_wv(const Scalar pfull, const Scalar phalf_top, const Scalar phalf_bot,
                              const Scalar qv, const Scalar at);
  KOKKOS_INLINE_FUNCTION
  static Scalar icarus_bb(const Scalar t);
  KOKKOS_INLINE_FUNCTION
  static Scalar icarus_tb(const Scalar flux);

  // Per-point part of icarus_subcolumn; qv, at, pfull by level, phalf by interface
  template <typename VQ, typename VA, typename VP, typename VPH>
  KOKKOS_INLINE_FUNCTION
  static IcarusPoint icarus_point(const int nlev, const Scalar skt, const Scalar emsfc_lw,
                                  const VQ& qv, const VA& at, const VP& pfull, const VPH& phalf);

  // icarus_subcolumn for one subcolumn of one point; dtau and dem by level
  template <typename VT, typename VE, typename VQ, typename VA, typename VP, typename VPH>
  KOKKOS_INLINE_FUNCTION
  static void icarus_subcolumn(const int nlev, const int sunlit, const IcarusPoint& pt,
                               const VT& dtau, const VE& dem, const Scalar skt,
                               const Scalar emsfc_lw, const VQ& qv, const VA& at,
                               const VP& pfull, const VPH& phalf,
                               Scalar& boxtau, Scalar& boxptop, Scalar& boxttop);

  // icarus_column for one point; boxtau, boxptop (Pa) by subcolumn, fq is (ntau, npres)
  template <typename VT, typename VP, typename H>
  KOKKOS_INLINE_FUNCTION
  static void icarus_column(const int ncol, const int sunlit, const VT& boxtau,
                            const VP& boxptop, const H& fq, Scalar& totalcldarea);

  // ---------------------------------------------------------------------------------------
  // MISR simulator (MISR_simulator.F90)

  // misr_subcolumn for one point and one subcolumn; dtau, zfull and at are indexed by level
  template <typename VT, typename VZ, typename VA>
  KOKKOS_INLINE_FUNCTION
  static void misr_subcolumn(const int nlev, const VT& dtau, const VZ& zfull, const VA& at,
                             const int sunlit, Scalar& tau_out, Scalar& ztop);

  // misr_column for one point; ztop and tau are indexed by subcolumn, fq is (ntau, nhgt)
  template <typename VZ, typename VT, typename H>
  KOKKOS_INLINE_FUNCTION
  static void misr_column(const int ncol, const VZ& ztop, const int sunlit, const VT& tau,
                          const H& fq);

  // ---------------------------------------------------------------------------------------
  // MODIS simulator (modis_simulator.F90), tau/CTP joint histogram only

  // two_stream_reflectance and two_stream (beam = 2) for one layer
  KOKKOS_INLINE_FUNCTION
  static Scalar two_stream_reflectance(const Scalar tau, const Scalar g, const Scalar w0);
  KOKKOS_INLINE_FUNCTION
  static void two_stream(const Scalar tau, const Scalar g, const Scalar w0, Scalar& ref,
                         Scalar& tra);

  // compute_toa_reflectace (with adding_doubling); tau, g, w0 by level
  template <typename VT, typename VG, typename VW>
  KOKKOS_INLINE_FUNCTION
  static Scalar compute_toa_reflectance(const int nlev, const VT& tau, const VG& g, const VW& w0);

  // cloud_top_pressure, as called from modis_subcolumn; phalf by interface
  template <typename VT, typename VPH>
  KOKKOS_INLINE_FUNCTION
  static Scalar modis_cloud_top_pressure(const int nlev, const VT& tau, const VPH& phalf,
                                         const Scalar tau_limit);

  // weight_by_extinction
  template <typename VT, typename VF>
  KOKKOS_INLINE_FUNCTION
  static Scalar modis_weight_by_extinction(const int nlev, const VT& tau, const VF& f,
                                           const Scalar tau_limit);

  // Size retrieval of modis_subcolumn (interpolate_to_min over the trial radii), in m
  KOKKOS_INLINE_FUNCTION
  static Scalar modis_retrieve_size(const int phase, const Scalar tau, const Scalar obs_refl);

  // modis_subcolumn for one subcolumn of a sunlit point; profiles by level, phalf by
  // interface, isccp_ctp is the ISCCP boxptop (Pa) of this subcolumn
  template <typename VPH, typename VT, typename VF, typename VG, typename VW>
  KOKKOS_INLINE_FUNCTION
  static void modis_subcolumn(const int nlev, const VPH& phalf, const VT& tau,
                              const VF& liq_frac, const VG& g, const VW& w0,
                              const Scalar isccp_ctp, int& phase, Scalar& ctp, Scalar& tau_ret,
                              Scalar& size_ret);

  // modis_column for one sunlit point; inputs by subcolumn, fq is (ntau, npres)
  template <typename VPh, typename VC, typename VT, typename VS, typename H>
  KOKKOS_INLINE_FUNCTION
  static void modis_column(const int ncol, const VPh& phase, const VC& ctp, const VT& tau,
                           const VS& size, const H& fq);

  // ---------------------------------------------------------------------------------------
  // Statistics (cosp_stats.F90)

  // hist2D: joint histogram of npts samples; bin1/bin2 hold nbin1+1/nbin2+1 edges, and
  // hist is indexed (ibin1, ibin2). Runs serially, intended to be called once per column.
  template <typename V1, typename V2, typename B1, typename B2, typename H>
  KOKKOS_INLINE_FUNCTION
  static void hist2d(const int npts, const V1& var1, const V2& var2,
                     const int nbin1, const B1& bin1, const int nbin2, const B2& bin2,
                     const H& hist);

  // ---------------------------------------------------------------------------------------
  // Driver (cosp_c2f_run in cosp_c2f.F90, and cosp_simulator in cosp.F90)

  // Scratch space for cosp_main: optical inputs by (column, subcolumn, level), and
  // simulator results by (column, subcolumn)
  struct Workspace {
    view_3d<Scalar> frac_out, tau, dem, liqfrac, g, w0;
    view_2d<Scalar> threshold, boxtau, boxptop, boxttop, misr_tau, misr_ztop, modis_ctp, modis_tau,
                    modis_size;
    view_2d<int>    modis_phase;

    Workspace() = default;
    Workspace(const int ncol, const int nsubcol, const int nlev);
  };

  // The whole COSP call for EAMxx: ISCCP total cloud area and tau/CTP histogram, MODIS
  // tau/CTP histogram (pressure from the surface up, as the Fortran returns them), and
  // MISR tau/CTH histogram, by column. Inputs are in-cloud values, levels go from the
  // top down. Night columns get R_UNDEF; the caller replaces them with its fill value.
  static void cosp_main(const int ncol, const int nsubcol, const int nlev,
                        const Scalar emsfc_lw,
                        const view_1d<const int>& sunlit, const view_1d<const Scalar>& skt,
                        const view_2d<const Scalar>& T_mid, const view_2d<const Scalar>& p_mid,
                        const view_2d<const Scalar>& p_int, const view_2d<const Scalar>& z_mid,
                        const view_2d<const Scalar>& qv, const view_2d<const Scalar>& qc,
                        const view_2d<const Scalar>& qi, const view_2d<const Scalar>& cldfrac,
                        const view_2d<const Scalar>& reff_qc,
                        const view_2d<const Scalar>& reff_qi,
                        const view_2d<const Scalar>& dtau067,
                        const view_2d<const Scalar>& dtau105,
                        const view_1d<Scalar>& isccp_cldtot,
                        const view_3d<Scalar>& isccp_ctptau,
                        const view_3d<Scalar>& modis_ctptau,
                        const view_3d<Scalar>& misr_cthtau, const Workspace& ws);
};

} // namespace cosp
} // namespace scream

#include "physics/cosp/impl/cosp_rng_impl.hpp"
#include "physics/cosp/impl/cosp_scops_impl.hpp"
#include "physics/cosp/impl/cosp_optics_impl.hpp"
#include "physics/cosp/impl/cosp_stats_impl.hpp"
#include "physics/cosp/impl/cosp_misr_impl.hpp"
#include "physics/cosp/impl/cosp_isccp_impl.hpp"
#include "physics/cosp/impl/cosp_modis_impl.hpp"
#include "physics/cosp/impl/cosp_main_impl.hpp"

#endif // SCREAM_COSP_CXX_FUNCTIONS_HPP
