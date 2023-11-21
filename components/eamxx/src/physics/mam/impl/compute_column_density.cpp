namespace scream::impl {

using KT      = ekat::KokkosTypes<DefaultDevice>;
using view_2d = typename KT::template view_2d<Real>;

KOKKOS_INLINE_FUNCTION
void compute_column_density(const ThreadTeam& team, const haero::Atmosphere& atm, view_2d& col_dens) {
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, atm.num_levels()), [&](const int k) {
    constexpr int nabscol = mam4::gas_chemistry::nabscol; // number of absorbing densities
    constexpr int gas_pcnst = mam4::gas_chemistry::gas_pcnst; // number of gas phase species
    constexpr int nfs = mam4::gas_chemistry::nfs;             // number of "fixed species"

    constexpr Real mwdry = 1.0/haero::Constants::molec_weight_dry_air;

    Real temp = atm.temperature(k);
    Real pmid = atm.pressure(k);
    Real pdel = atm.hydrostatic_dp(k);
    Real qv   = atm.vapor_mixing_ratio(k);

    // ... map incoming mass mixing ratios to working array
    Real mmr[gas_pcnst] = {};
    // FIXME: come back and fix this

    // ... set atmosphere mean mass to the molecular weight of dry air
    Real mbar = mwdry;

    // ... Xform from mmr to vmr
    Real vmr[gas_pcnst];
    for (int i = 0; i < gas_pcnst; ++i) {
      vmr[i] = mam4::conversions::vmr_from_mmr(mmr[i], mbar);
    }
    Real h2ovmr = mam4::conversions::vmr_from_mmr(qv, mbar);

    // ... compute invariants for this level
    Real invariants[nfs];
    // setinv(invariants, temp, h2ovmr, vmr, pmid); FIXME

    // compute the change in density for this column above its neighbor
    Real col_delta[nabscol]; // o2, o3 column density above model [1/cm^2]
    //set_ub_col(col_delta, vmr, invariants, pdel); FIXME
  });
  //setcol(col_delta, col_dens); // FIXME
}

} // namespace scream::impl
