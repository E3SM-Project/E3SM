namespace scream::impl {

using View2D = DeviceType::view_2d<Real>;

KOKKOS_INLINE_FUNCTION
void compute_o3_column_density(
    const ThreadTeam &team, const haero::Atmosphere &atm,
    const mam4::Prognostics &progs, const View2D &invariants,
    const Real adv_mass_kg_per_moles[mam4::gas_chemistry::gas_pcnst],
    ColumnView o3_col_dens) {
  constexpr int gas_pcnst =
      mam4::gas_chemistry::gas_pcnst;            // number of gas phase species
  constexpr int nfs = mam4::gas_chemistry::nfs;  // number of "fixed species"
  constexpr int offset_aerosol = mam4::utils::gasses_start_ind();

  Real o3_col_deltas[mam4::nlev + 1] =
      {};  // o3 column density above model [1/cm^2]
  // NOTE: if we need o2 column densities, set_ub_col and setcol must be changed
  Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, mam4::nlev), [&](const int k) {
        Real pdel = atm.hydrostatic_dp(k);
        // extract aerosol state variables into "working arrays" (mass
        // mixing ratios) (in EAM, this is done in the gas_phase_chemdr
        // subroutine defined within
        //  mozart/mo_gas_phase_chemdr.F90)
        Real q[gas_pcnst]   = {};
        Real state_q[pcnst] = {};
        mam4::utils::extract_stateq_from_prognostics(progs, atm, state_q, k);

        for(int i = offset_aerosol; i < pcnst; ++i) {
          q[i - offset_aerosol] = state_q[i];
        }

        // convert mass mixing ratios to volume mixing ratios (VMR),
        // equivalent to tracer mixing ratios (TMR))
        Real vmr[gas_pcnst];
        mam_coupling::mmr2vmr(q, adv_mass_kg_per_moles, vmr);
        // ... compute invariants for this level
        Real invariants_k[nfs];
        for(int i = 0; i < nfs; ++i) {
          invariants_k[i] = invariants(k, i);
        }
        // compute the change in o3 density for this column above its neighbor
        mam4::mo_photo::set_ub_col(o3_col_deltas[k + 1], vmr, invariants_k,
                                   pdel);
      });
  // sum the o3 column deltas to densities
  mam4::mo_photo::setcol(o3_col_deltas, o3_col_dens);
}

}  // namespace scream::impl
