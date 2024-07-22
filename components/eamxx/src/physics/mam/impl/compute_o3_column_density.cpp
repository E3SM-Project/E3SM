namespace scream::impl {

using View2D = DeviceType::view_2d<Real>;

KOKKOS_INLINE_FUNCTION
void compute_o3_column_density(const ThreadTeam &team,
                               const haero::Atmosphere &atm,
                               const mam4::Prognostics &progs,
                               const View2D &invariants,
                               ColumnView o3_col_dens) {
  constexpr int gas_pcnst =
      mam4::gas_chemistry::gas_pcnst;            // number of gas phase species
  constexpr int nfs = mam4::gas_chemistry::nfs;  // number of "fixed species"

  Real o3_col_deltas[mam4::nlev + 1] =
      {};  // o3 column density above model [1/cm^2]
  // NOTE: if we need o2 column densities, set_ub_col and setcol must be changed
  Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, atm.num_levels()), [&](const int k) {
        Real temp = atm.temperature(k);
        Real pmid = atm.pressure(k);
        Real pdel = atm.hydrostatic_dp(k);

        // ... map incoming mass mixing ratios to working array
        Real q[gas_pcnst], qqcw[gas_pcnst];
        mam_coupling::transfer_prognostics_to_work_arrays(progs, k, q, qqcw);

        // ... Xform from mmr to vmr
        Real vmr[gas_pcnst], vmrcw[gas_pcnst];
        mam_coupling::convert_work_arrays_to_vmr(q, qqcw, vmr, vmrcw);

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
