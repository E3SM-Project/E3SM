#ifndef SPA_FUNCTIONS_IMPL_HPP
#define SPA_FUNCTIONS_IMPL_HPP

#include "physics/share/physics_constants.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "ekat/util/ekat_lin_interp.hpp"

/*-----------------------------------------------------------------*/
/* The main SPA routines used to convert SPA data into a format that
 * is usable by the rest of the atmosphere processes.
 *
 * SPA or Simple Prescribed Aerosols provides a way to prescribe
 * aerosols for an atmospheric simulation using pre-computed data.
 * 
 * The data is typically provided at a frequency of monthly, and
 * does not necessarily have to be on the same horizontal or vertical
 * domain as the atmospheric simulation.
 *
 * In order to accomodate coarse temporal resolution and a potentially
 * different spatial resolution it is necessary to perform a series
 * of interpolations, which make up the main body of the SPA routines.
 *
 * The interpolations can be broken into three categories.
 * 1. Horizontal Interpolation: TODO - not done yet.
 * The SPA data set does not have to be provided on the same grid as
 * the atmospheric simulation.  Whenever SPA data is loaded, it is
 * interpolated horizontally onto the simulation grid to provide
 * forcing at every location.  This can be done by,
 *   a) Preloaded remapping wieghts which are applied at every
 *      horizontal column.
 *   b) Online calculation of remapping weights given a set of lat/lon
 *      for the source data and comparing it with the lat/lon of each
 *      column in the simulation.  TODO: THIS HAS NOT BEEN IMPLEMENTED YET
 *
 * 2. Temporal Interpolation:
 * As noted above, the SPA data is provided at some fixed frequency.  Typically
 * as monthly data.  As a result, the data must be interpolated to the current
 * time of the simulation at each time step.  Temporal interpolation follows
 * a basic linear interpolation and is performed for all SPA data at all columns
 * and levels.
 * Note: There is also a temporal interpolation of the surface pressure for the SPA
 * data, which is used in the vertical reconstruction of the pressure profile.
 *
 * 3. Vertical Interpolation:
 * Given that the SPA data has been generated elsewhere it is very likely that
 * the vertical pressure profiles of the data won't match the simulation pressure
 * profiles.  The vertical SPA data structure must be remapped onto the simulation
 * pressure profile.
 * This is done using the EKAT linear interpolation code, see /externals/ekat/util/ekat_lin_interp.hpp
 * The SPA pressure profiles are calculated using the surface pressure which was
 * temporally interpolated in the last step and the set of hybrid coordinates (hyam and hybm)
 * that are used in EAM to construct the physics pressure profiles.
 * The SPA data is then projected onto the simulation pressure profile (pmid)
 * using EKAT linear interpolation. 
/*-----------------------------------------------------------------*/

namespace scream {
namespace spa {

// Helper function
template<typename ScalarT,typename ScalarS>
ScalarT linear_interp(const ScalarT& x0, const ScalarT& x1, const ScalarS& t_norm);
/*-----------------------------------------------------------------*/
// The main SPA routine which handles projecting SPA data onto the
// horizontal columns and vertical pressure profiles of the atmospheric
// state.
// Inputs:
//   time_state: A structure defined in spa_functions.hpp which handles
//     the current temporal state of the simulation.
//   pressure_state: A structure defined in spa_functions.hpp which handles
//     the vertical pressure profile for the atmospheric simulation state, and
//     all of the data needed to reconstruct the vertical pressure profile for
//     the SPA data.  See hybrid coordinate (hyam,hybm) and surface pressure (PS)
//   data_beg: A structure defined in spa_functions.hpp which handles the full
//     set of SPA data for the beginning of the month.
//   data_end: Similar to data_beg, but for SPA data for the end of the month.
//   data_out: A structure defined in spa_functions.hpp which handles the full
//     set of SPA data projected onto the pressure profile of the current atmosphere
//     state.  This is the data that will be passed to other processes.
//   ncols_atm, nlevs_atm: The number of columns and levels in the simulation grid.
//     (not to be confused with the number of columns and levels used for the SPA data, 
//      which can be different.)
//   nswbands, nlwbands: The number of shortwave (sw) and longwave (lw) aerosol bands 
//     for the data that will be passed to radiation.
template <typename S, typename D>
void SPAFunctions<S,D>
::spa_main(
  const SPATimeState& time_state,
  const SPAPressureState& pressure_state,
  const SPAData&   data_beg,
  const SPAData&   data_end,
  const SPAOutput& data_out,
  Int ncols_atm,
  Int nlevs_atm,
  Int nswbands,
  Int nlwbands)
{
  // Gather time stamp info
  auto& t_now = time_state.t_now;
  auto& t_beg = time_state.t_beg_month;
  auto& t_len = time_state.days_this_month;

  // For now we require that the Data in and the Data out have the same number of columns.
  EKAT_REQUIRE(ncols_atm==pressure_state.ncols);

  // Set up temporary arrays that will be used for the spa interpolation.
  view_2d<Spack> p_src("p_mid_src",ncols_atm,pressure_state.nlevs), 
                 ccn3_src("ccn3_src",ncols_atm,pressure_state.nlevs);
  view_3d<Spack> aer_g_sw_src("aer_g_sw_src",ncols_atm,nswbands,pressure_state.nlevs),
                 aer_ssa_sw_src("aer_ssa_sw_src",ncols_atm,nswbands,pressure_state.nlevs),
                 aer_tau_sw_src("aer_tau_sw_src",ncols_atm,nswbands,pressure_state.nlevs),
                 aer_tau_lw_src("aer_tau_lw_src",ncols_atm,nlwbands,pressure_state.nlevs);

  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nlevs_atm);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncols_atm, nk_pack);
  // SPA Main loop
  // Parallel loop order:
  // 1. Loop over all horizontal columns (i index)
  // 2. Loop over all aerosol bands (n index) - where applicable
  // 3. Loop over all vertical packs (k index)
  Kokkos::parallel_for(
    "spa main loop",
    policy,
    KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();  // SCREAM column index

    // Get single-column subviews of all 2D inputs, i.e. those that don't have aerosol bands
    const auto& ps_beg_sub                 = pressure_state.ps_this_month(i);
    const auto& ps_end_sub                 = pressure_state.ps_next_month(i);
    const auto& pmid_sub                   = ekat::subview(pressure_state.pmid, i);
    const auto& p_src_sub                  = ekat::subview(p_src, i);

    const auto& ccn3_beg_sub               = ekat::subview(data_beg.CCN3, i);
    const auto& ccn3_end_sub               = ekat::subview(data_end.CCN3, i);
    const auto& ccn3_src_sub               = ekat::subview(ccn3_src, i);

    // First Step: Horizontal Interpolation if needed - Skip for Now
  
    // Second Step: Temporal Interpolation
    // Use basic linear interpolation function y = b + mx
    auto t_norm = (t_now-t_beg)/t_len;
    /* Determine PS for the source data at this time */
    auto ps_src = linear_interp(ps_beg_sub,ps_end_sub,t_norm);
    {
    /* Reconstruct the vertical pressure profile for the data and time interpolation
     * of the data.
     * Note: CCN3 has the same dimensions as pressure so we handle that time interpolation
     *       in this loop as well.  */
    using C = scream::physics::Constants<Real>;
    static constexpr auto P0 = C::P0;
    const Int nk_pack = ekat::npack<Spack>(pressure_state.nlevs);
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
        // Reconstruct vertical temperature profile using the hybrid coordinate system.
        p_src_sub(k)    = ps_src * pressure_state.hybm(k) + P0 * pressure_state.hyam(k);
        // Time interpolation for CCN3
        ccn3_src_sub(k) = linear_interp(ccn3_beg_sub(k),ccn3_end_sub(k),t_norm);
    });
    team.team_barrier();
    }
    {
    /* Loop over all SW variables with nswbands */
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nswbands), [&] (int n) {
      const auto& aer_g_sw_beg_sub           = ekat::subview(data_beg.AER_G_SW, i, n);
      const auto& aer_g_sw_end_sub           = ekat::subview(data_end.AER_G_SW, i, n);
      const auto& aer_g_sw_src_sub           = ekat::subview(aer_g_sw_src, i, n);

      const auto& aer_ssa_sw_beg_sub         = ekat::subview(data_beg.AER_SSA_SW, i, n);
      const auto& aer_ssa_sw_end_sub         = ekat::subview(data_end.AER_SSA_SW, i, n);
      const auto& aer_ssa_sw_src_sub         = ekat::subview(aer_ssa_sw_src, i, n);
  
      const auto& aer_tau_sw_beg_sub         = ekat::subview(data_beg.AER_TAU_SW, i, n);
      const auto& aer_tau_sw_end_sub         = ekat::subview(data_end.AER_TAU_SW, i, n);
      const auto& aer_tau_sw_src_sub         = ekat::subview(aer_tau_sw_src, i, n);

      /* Now loop over fastest index, the number of vertical packs */
      Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(team, nk_pack), [&] (Int k) {
        aer_g_sw_src_sub(k)   = linear_interp(aer_g_sw_beg_sub(k),   aer_g_sw_end_sub(k),   t_norm);
        aer_ssa_sw_src_sub(k) = linear_interp(aer_ssa_sw_beg_sub(k), aer_ssa_sw_end_sub(k), t_norm);
        aer_tau_sw_src_sub(k) = linear_interp(aer_tau_sw_beg_sub(k), aer_tau_sw_end_sub(k), t_norm);
      });
    });
    team.team_barrier();
    }
    {
    /* Loop over all LW variables with nlwbands */
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nlwbands), [&] (int n) {
      const auto& aer_tau_lw_beg_sub         = ekat::subview(data_beg.AER_TAU_LW, i, n);
      const auto& aer_tau_lw_end_sub         = ekat::subview(data_end.AER_TAU_LW, i, n);
      const auto& aer_tau_lw_src_sub         = ekat::subview(aer_tau_lw_src, i, n);

      /* Now loop over fastest index, the number of vertical packs */
      Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(team, nk_pack), [&] (Int k) {
        aer_tau_lw_src_sub(k) = linear_interp(aer_tau_lw_beg_sub(k), aer_tau_lw_end_sub(k), t_norm);
      });
    });
    team.team_barrier();
    }
  });
  Kokkos::fence();

  // Third Step: Vertical interpolation, project the SPA data onto the pressure profile for this simulation.
  // This is done using the EKAT linear interpolation routine, see /externals/ekat/util/ekat_lin_interp.hpp
  // for more details. 
  using LIV = ekat::LinInterp<Real,Spack::n>;
  Real minthreshold = 0.0;  // Hard-code a minimum value for aerosol concentration to zero.

  LIV VertInterp(ncols_atm,pressure_state.nlevs,nlevs_atm,minthreshold);
  /* Parallel loop strategy:
   * 1. Loop over all simulation columns (i index)
   * 2. Where applicable, loop over all aerosol bands (n index)
   */ 
  const Int most_bands = std::max(nlwbands, nswbands);
  typename LIV::TeamPolicy band_policy(ncols_atm, ekat::OnGpu<typename LIV::ExeSpace>::value ? most_bands : 1, VertInterp.km2_pack());
  Kokkos::parallel_for("vertical-interp-spa",
    band_policy,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    const int i = team.league_rank();
    /* Setup the linear interpolater for this column. */
    VertInterp.setup(team,
                     ekat::subview(p_src,i),
                     ekat::subview(pressure_state.pmid,i));
    team.team_barrier();
    /* Conduct vertical interpolation for the 2D variable CCN3 */
    VertInterp.lin_interp(team,
                          ekat::subview(p_src,i),
                          ekat::subview(pressure_state.pmid,i),
                          ekat::subview(ccn3_src,i),
                          ekat::subview(data_out.CCN3,i));
    /* Conduct vertical interpolation for the LW banded data - nlwbands (n index) */
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nlwbands), [&] (int n) {
      const auto& tvr = Kokkos::ThreadVectorRange(team, VertInterp.km2_pack());
      VertInterp.lin_interp(team,
                            tvr,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_tau_lw_src,i,n),
                            ekat::subview(data_out.AER_TAU_LW,i,n));
    });
    /* Conduct vertical interpolation for the SW banded data - nswbands (n index) */
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nswbands), [&] (int n) {
      const auto& tvr = Kokkos::ThreadVectorRange(team, VertInterp.km2_pack());
      VertInterp.lin_interp(team,
                            tvr,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_g_sw_src,i,n),
                            ekat::subview(data_out.AER_G_SW,i,n));
      VertInterp.lin_interp(team,
                            tvr,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_ssa_sw_src,i,n),
                            ekat::subview(data_out.AER_SSA_SW,i,n));
      VertInterp.lin_interp(team,
                            tvr,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_tau_sw_src,i,n),
                            ekat::subview(data_out.AER_TAU_SW,i,n));
    });
  });
  Kokkos::fence();

}
/*-----------------------------------------------------------------*/
// A helper function to manage basic linear interpolation in time.
// The inputs of x0 and x1 represent the data to interpolate from at
// times t0 and t1, respectively.  To keep the signature of the function 
// simple we use
//    t_norm = (t-t0)/(t1-t0).
template<typename ScalarT, typename ScalarS>
inline ScalarT linear_interp(
  const ScalarT& x0,
  const ScalarT& x1,
  const ScalarS& t_norm)
{
  return (1.0 - t_norm) * x0 + t_norm * x1;
}
/*-----------------------------------------------------------------*/

} // namespace spa
} // namespace scream

#endif // SPA_FUNCTIONS_IMPL_HPP
