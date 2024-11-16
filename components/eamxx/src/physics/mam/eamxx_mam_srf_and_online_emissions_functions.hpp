#ifndef EAMXX_MAM_SRF_AND_ONLINE_EMISSIONS_FUNCTIONS_HPP
#define EAMXX_MAM_SRF_AND_ONLINE_EMISSIONS_FUNCTIONS_HPP

namespace scream {

namespace {

using KT            = ekat::KokkosTypes<DefaultDevice>;
using view_1d       = typename KT::template view_1d<Real>;
using view_2d       = typename KT::template view_2d<Real>;
using const_view_1d = typename KT::template view_1d<const Real>;
using const_view_2d = typename KT::template view_2d<const Real>;

//-------- Inititlize gas and aerosol fluxes ------
void init_fluxes(const int &ncol,
                 view_2d &constituent_fluxes) {  // input-output

  constexpr int pcnst     = mam4::aero_model::pcnst;
  const int gas_start_ind = mam4::utils::gasses_start_ind();

  const auto policy =
      ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(
          ncol, pcnst - gas_start_ind);

  // Parallel loop over all the columns
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const KT::MemberType &team) {
        const int icol   = team.league_rank();
        view_1d flux_col = ekat::subview(constituent_fluxes, icol);

        // Zero out constituent fluxes only for gasses and aerosols
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, gas_start_ind, pcnst),
            [&](int icnst) { flux_col(icnst) = 0; });
      });
}  // init_fluxes ends

//-------- compute online emissions for dust, sea salt and marine organics -----
void compute_online_dust_nacl_emiss(
    const int &ncol, const int &nlev, const const_view_1d &ocnfrac,
    const const_view_1d &sst, const const_view_2d &u_wind,
    const const_view_2d &v_wind, const const_view_2d &dstflx,
    const const_view_1d &mpoly, const const_view_1d &mprot,
    const const_view_1d &mlip, const const_view_1d &soil_erodibility,
    const const_view_2d &z_mid,
    // output
    view_2d &constituent_fluxes) {
  const int surf_lev = nlev - 1;  // surface level

  Kokkos::parallel_for(
      "online_emis_fluxes", ncol, KOKKOS_LAMBDA(int icol) {
        // Input
        const const_view_1d dstflx_icol = ekat::subview(dstflx, icol);

        // Output
        view_1d fluxes_col = ekat::subview(constituent_fluxes, icol);

        // Compute online emissions
        // NOTE: mam4::aero_model_emissions calculates mass and number emission
        // fluxes in units of [kg/m2/s or #/m2/s] (MKS), so no need to convert
        mam4::aero_model_emissions::aero_model_emissions(
            sst(icol), ocnfrac(icol), u_wind(icol, surf_lev),
            v_wind(icol, surf_lev), z_mid(icol, surf_lev), dstflx_icol,
            soil_erodibility(icol), mpoly(icol), mprot(icol), mlip(icol),
            // out
            fluxes_col);
      });
}  // compute_online_dust_nacl_emiss ends

}  // namespace
}  // namespace scream

#endif  // EAMXX_MAM_SRF_AND_ONLINE_EMISSIONS_FUNCTIONS_HPP
