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
  const int surface_lev = nlev - 1;

  // get the start index for gas species in the state_q array
  int istart = mam4::utils::gasses_start_ind();

  // number of constituents to update (currently updating only MAM4xx
  // constituents)
  const int nconstituents = pcnst - istart;

  // Create a policy to loop over columns annd number of constituents
  // to update
  // FIXME: TODO:We don't need a team for "nconstituents", so we can make the
  // kookos_for simple by using just ncols
  const auto policy = ekat::ExeSpaceUtils<MAMConstituentFluxes::KT::ExeSpace>::
      get_default_team_policy(ncol, nconstituents);

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

        // Form state%q like array at surface level
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

        // Update state vector with constituent fluxes
        for(int icnst = istart; icnst < pcnst; ++icnst) {
          state_q_at_surf_lev[icnst] +=
              constituent_fluxes(icol, icnst) * unit_factor;
        }
        mam4::utils::inject_stateq_to_prognostics(state_q_at_surf_lev,  // input
                                                  progs_at_col,  // output
                                                  surface_lev);  // input
      });                                                        // icol loop
}

}  // namespace
}  // namespace scream

#endif  // ifdef EAMXX_MAM_CONSTITUTE_FLUXES_FUNCTIONS_HPP
