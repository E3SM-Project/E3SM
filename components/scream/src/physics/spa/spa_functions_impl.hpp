#ifndef SPA_FUNCTIONS_IMPL_HPP
#define SPA_FUNCTIONS_IMPL_HPP

#include "physics/share/physics_constants.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "ekat/util/ekat_lin_interp.hpp"

namespace scream {
namespace spa {

/*-----------------------------------------------------------------*/
template <typename S, typename D>
void SPAFunctions<S,D>
::spa_main(
  const SPATimeState& time_state,
  const SPAPressureState& pressure_state,
  const SPAData&   data_beg,
  const SPAData&   data_end,
  const SPAOutput& data_out,
  Int ncols_scream,
  Int nlevs_scream,
  Int nswbands,
  Int nlwbands)
{
  // Gather time stamp info
  auto& t_now = time_state.t_now;
  auto& t_beg = time_state.t_beg_month;
  auto& t_len = time_state.days_this_month;

  // For now we require that the Data in and the Data out have the same number of columns.
  EKAT_REQUIRE(ncols_scream==pressure_state.ncols);

  view_2d<Spack> p_src("p_mid_src",ncols_scream,pressure_state.nlevs), 
                 ccn3_src("ccn3_src",ncols_scream,pressure_state.nlevs);
  view_3d<Spack> aer_g_sw_src("aer_g_sw_src",ncols_scream,nswbands,pressure_state.nlevs),
                 aer_ssa_sw_src("aer_ssa_sw_src",ncols_scream,nswbands,pressure_state.nlevs),
                 aer_tau_sw_src("aer_tau_sw_src",ncols_scream,nswbands,pressure_state.nlevs),
                 aer_tau_lw_src("aer_tau_lw_src",ncols_scream,nlwbands,pressure_state.nlevs);

  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nlevs_scream);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncols_scream, nk_pack);
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
    auto slope = (t_now-t_beg)/t_len;
    /* Determine PS for the source data at this time */
    auto ps_src  =  ps_beg_sub + slope * (ps_end_sub-ps_beg_sub);
    /* Reconstruct the vertical pressure profile for the data and time interpolation
     * of the data */
    {
    using C = scream::physics::Constants<Real>;
    static constexpr auto P0 = C::P0;
    const Int nk_pack = ekat::npack<Spack>(pressure_state.nlevs);
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
        p_src_sub(k)    = ps_src * pressure_state.hybm(k) + P0 * pressure_state.hyam(k);
        ccn3_src_sub(k) = ccn3_beg_sub(k) + slope * (ccn3_end_sub(k) - ccn3_beg_sub(k));
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
        aer_g_sw_src_sub(k)   = aer_g_sw_beg_sub(k)   + slope * (aer_g_sw_end_sub(k)   - aer_g_sw_beg_sub(k));
        aer_ssa_sw_src_sub(k) = aer_ssa_sw_beg_sub(k) + slope * (aer_ssa_sw_end_sub(k) - aer_ssa_sw_beg_sub(k));
        aer_tau_sw_src_sub(k) = aer_tau_sw_beg_sub(k) + slope * (aer_tau_sw_end_sub(k) - aer_tau_sw_beg_sub(k));
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
        aer_tau_lw_src_sub(k) = aer_tau_lw_beg_sub(k) + slope * (aer_tau_lw_end_sub(k) - aer_tau_lw_beg_sub(k));
      });
    });
    team.team_barrier();
    }
  });
  Kokkos::fence();

  using LIV = ekat::LinInterp<Real,Spack::n>;
  Real minthreshold = 0.0;  // Hard-code a minimum value for aerosol concentration to zero.

  LIV VertInterp(ncols_scream,pressure_state.nlevs,nlevs_scream,minthreshold);
  Kokkos::parallel_for("vertical-interp-spa",
    VertInterp.m_policy,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    const int i = team.league_rank();
    VertInterp.setup(team,
                     ekat::subview(p_src,i),
                     ekat::subview(pressure_state.pmid,i));
    team.team_barrier();
    VertInterp.lin_interp(team,
                          ekat::subview(p_src,i),
                          ekat::subview(pressure_state.pmid,i),
                          ekat::subview(ccn3_src,i),
                          ekat::subview(data_out.CCN3,i));
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nlwbands), [&] (int n) {
      VertInterp.lin_interp(team,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_tau_lw_src,i,n),
                            ekat::subview(data_out.AER_TAU_LW,i,n));
    });
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nswbands), [&] (int n) {
      VertInterp.lin_interp(team,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_g_sw_src,i,n),
                            ekat::subview(data_out.AER_G_SW,i,n));
      VertInterp.lin_interp(team,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_ssa_sw_src,i,n),
                            ekat::subview(data_out.AER_SSA_SW,i,n));
      VertInterp.lin_interp(team,
                            ekat::subview(p_src,i),
                            ekat::subview(pressure_state.pmid,i),
                            ekat::subview(aer_tau_sw_src,i,n),
                            ekat::subview(data_out.AER_TAU_SW,i,n));
    });
  });
  Kokkos::fence();
}
/*-----------------------------------------------------------------*/

} // namespace spa
} // namespace scream

#endif // SPA_FUNCTIONS_IMPL_HPP
