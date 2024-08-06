#ifndef EAMXX_MAM_CONSTITUTE_FLUXES_FUNCTIONS_HPP
#define EAMXX_MAM_CONSTITUTE_FLUXES_FUNCTIONS_HPP

#include <mam4xx/mam4.hpp>
#include <physics/mam/eamxx_mam_constituent_fluxes_interface.hpp>

namespace scream {

namespace {

void update_gas_aerosols_using_constituents(
    const int ncol, const int nlev, const double dt,
    const mam_coupling::DryAtmosphere &dry_atm,
    const MAMConstituentFluxes::const_view_2d &constituent_fluxes,
    // output
    const mam_coupling::AerosolState &wet_aero) {
  using C                      = physics::Constants<Real>;
  static constexpr auto gravit = C::gravit;  // Gravity [m/s2]
  static constexpr int pcnst   = mam4::aero_model::pcnst;

  // Declare local variables
  const int ncol_loc          = ncol;
  auto wet_aero_loc           = wet_aero;
  auto dry_atm_loc            = dry_atm;
  const int surface_lev       = nlev - 1;
  auto constituent_fluxes_loc = constituent_fluxes;

  // get the start index for gas species in the state_q array
  int istart = mam4::utils::gasses_start_ind();

  // number of constituents to update
  const int nconstituents = pcnst - istart;

  for(int icnst = 0; icnst < 6; ++icnst) {
    auto host_view = Kokkos::create_mirror_view(wet_aero.gas_mmr[icnst]);
    printf("BEFORE:::%e, %i\n", host_view(0, surface_lev), icnst + 9);
  }

  // Create a policy to loop over columns annd number of constituents
  // to update
  // const auto policy =
  //     ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol,
  //     nconstituents);
  const auto policy = ekat::ExeSpaceUtils<MAMConstituentFluxes::KT::ExeSpace>::
      get_default_team_policy(1, nconstituents);

  // Loop through all columns to update tracer mixing rations
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
        const int icol = team.league_rank();

        //----------------------------------------------------------------------
        // To form EAM like state%q array, we need prognostics (gas and aerosol
        // mmrs) atmosphere (qv, qc, nc, ni, etc.)
        //----------------------------------------------------------------------

        // Get prognostics
        mam4::Prognostics progs_at_col =
            mam_coupling::aerosols_for_column(wet_aero,  // output
                                              icol);     // input

        // Get atmospheric quantities
        const haero::Atmosphere haero_atm =
            atmosphere_for_column(dry_atm,  // output
                                  icol);    // input

        // Form state%q like array
        Real state_q_at_surf_lev[pcnst] = {};
        mam4::utils::extract_stateq_from_prognostics(
            progs_at_col, haero_atm,  // input
            state_q_at_surf_lev,      // output
            surface_lev);             // input

        // Compute the units conversion factor (kg/m2/s to kg/kg)
        EKAT_KERNEL_ASSERT_MSG(dry_atm.p_del(icol, surface_lev) != 0,
                               "Error! dry_atm.pdel must be non-zero!\n");
        const Real rpdel       = 1.0 / dry_atm.p_del(icol, surface_lev);
        const Real unit_factor = dt * gravit * rpdel;

        // Loop for
        auto pcnst_loc = pcnst;
        for(int icnst = 9; icnst < 15; ++icnst) {
          printf("BEFORE-state:%e,%e,%i, %i\n", state_q_at_surf_lev[icnst],
                 progs_at_col.q_gas[icnst - 9](surface_lev), icnst, icol);
        }
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, istart, pcnst_loc), [&](int icnst) {
              state_q_at_surf_lev[icnst] +=
                  constituent_fluxes(icol, icnst) * unit_factor;
            });  // pcsnt loop
        mam4::utils::inject_stateq_to_prognostics(state_q_at_surf_lev,  // input
                                                  progs_at_col,  // output
                                                  surface_lev);  // input

        for(int icnst = 9; icnst < 15; ++icnst) {
          printf("After-state:%e,%e,%i, %i\n", state_q_at_surf_lev[icnst],
                 progs_at_col.q_gas[icnst - 9](surface_lev), icnst, icol);
        }
      });  // icol loop
}

}  // namespace
}  // namespace scream

#endif  // ifdef EAMXX_MAM_CONSTITUTE_FLUXES_FUNCTIONS_HPP
