#include "physics/p3/eamxx_p3_process_interface.hpp"

#ifdef EAMXX_HAS_PYTHON
#include "share/atm_process/atmosphere_process_pyhelpers.hpp"
#endif

#include <ekat_team_policy_utils.hpp>

namespace scream {

void P3Microphysics::run_impl (const double dt)
{
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  // Set the dt for p3 postprocessing
  p3_postproc.m_dt = dt;

  // Create policy for pre and post process pfor
  const auto nlev_packs  = ekat::npack<Pack>(m_num_levs);
  const auto policy = TPF::get_default_team_policy(m_num_cols, nlev_packs);

  // Assign values to local arrays used by P3, these are now stored in p3_loc.
  Kokkos::parallel_for(
    "p3_pre_process",
    policy,
    p3_preproc
  );
  Kokkos::fence();

  if (runtime_options.use_sdm_warm_emulator) {
    run_sdm_warm_emulator();
  }

  // Update the variables in the p3 input structures with local values.

  infrastructure.dt = dt;
  infrastructure.it++;

  // Reset internal WSM variables.
  workspace_mgr.reset_internals();

  // Run p3 main
  get_field_out("micro_liq_ice_exchange").deep_copy(0.0);
  get_field_out("micro_vap_liq_exchange").deep_copy(0.0);
  get_field_out("micro_vap_ice_exchange").deep_copy(0.0);

  // Optional extra p3 diags
  if (runtime_options.extra_p3_diags) {
    get_field_out("qr2qv_evap").deep_copy(0.0);
    get_field_out("qi2qv_sublim").deep_copy(0.0);
    get_field_out("qc2qr_accret").deep_copy(0.0);
    get_field_out("qc2qr_autoconv").deep_copy(0.0);
    get_field_out("qv2qi_vapdep").deep_copy(0.0);
    get_field_out("qc2qi_berg").deep_copy(0.0);
    get_field_out("qc2qr_ice_shed").deep_copy(0.0);
    get_field_out("qc2qi_collect").deep_copy(0.0);
    get_field_out("qr2qi_collect").deep_copy(0.0);
    get_field_out("qc2qi_hetero_freeze").deep_copy(0.0);
    get_field_out("qr2qi_immers_freeze").deep_copy(0.0);
    get_field_out("qi2qr_melt").deep_copy(0.0);
    get_field_out("qr_sed").deep_copy(0.0);
    get_field_out("qc_sed").deep_copy(0.0);
    get_field_out("qi_sed").deep_copy(0.0);
  }

  P3F::p3_main(runtime_options, prog_state, diag_inputs, diag_outputs, infrastructure,
               history_only, lookup_tables,
#ifdef SCREAM_P3_SMALL_KERNELS
               temporaries,
#endif
               workspace_mgr, m_num_cols, m_num_levs);

  // Conduct the post-processing of the p3_main output.
  Kokkos::parallel_for(
    "p3_post_process",
    policy,
    p3_postproc
  );
  Kokkos::fence();
}

void P3Microphysics::run_sdm_warm_emulator ()
{
#ifdef EAMXX_HAS_PYTHON
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  const auto nlev_packs = ekat::npack<Pack>(m_num_levs);
  const auto policy = TPF::get_default_team_policy(m_num_cols, nlev_packs);

  auto& emu_qc_field  = get_internal_field("sdm_warm_emulator_qc");
  auto& emu_nc_field  = get_internal_field("sdm_warm_emulator_nc");
  auto& emu_qr_field  = get_internal_field("sdm_warm_emulator_qr");
  auto& emu_nr_field  = get_internal_field("sdm_warm_emulator_nr");
  auto& emu_rho_field = get_internal_field("sdm_warm_emulator_rho");

  const auto emu_qc  = emu_qc_field.get_view<Real**>();
  const auto emu_nc  = emu_nc_field.get_view<Real**>();
  const auto emu_qr  = emu_qr_field.get_view<Real**>();
  const auto emu_nr  = emu_nr_field.get_view<Real**>();
  const auto emu_rho = emu_rho_field.get_view<Real**>();

  const auto qc        = prog_state.qc;
  const auto nc        = prog_state.nc;
  const auto qr        = prog_state.qr;
  const auto nr        = prog_state.nr;
  const auto dpres     = diag_inputs.dpres;
  const auto dz        = diag_inputs.dz;
  const auto cld_frac_l = diag_inputs.cld_frac_l;
  const auto nccn      = diag_inputs.nccn;

  const auto nlev = m_num_levs;
  const auto prescribed_ccn = infrastructure.prescribedCCN;
  const auto spa_ccn_to_nc_factor = runtime_options.spa_ccn_to_nc_factor;
  const auto spa_ccn_to_nc_exponent = runtime_options.spa_ccn_to_nc_exponent;
  constexpr Real g = PC::gravit.value;

  Kokkos::parallel_for(
    "p3_sdm_warm_emulator_prepare",
    policy,
    KOKKOS_LAMBDA(const KT::MemberType& team) {
      const int icol = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&] (const Int& ipack) {
        const auto qc_pack = qc(icol,ipack);
        const auto nc_pack = nc(icol,ipack);
        const auto qr_pack = qr(icol,ipack);
        const auto nr_pack = nr(icol,ipack);
        const auto dpres_pack = dpres(icol,ipack);
        const auto dz_pack = dz(icol,ipack);
        const auto cld_frac_l_pack = cld_frac_l(icol,ipack);

        for (Int lane = 0; lane < Pack::n; ++lane) {
          const Int ilev = ipack*Pack::n + lane;
          if (ilev < nlev) {
            Real nc_in = nc_pack[lane];
            if (prescribed_ccn) {
              const auto nccn_pack = nccn(icol,ipack);
              const Real nccn_scaled = nccn_pack[lane] * cld_frac_l_pack[lane];
              nc_in = Kokkos::max(nc_in,
                                   spa_ccn_to_nc_factor *
                                   Kokkos::pow(nccn_scaled, spa_ccn_to_nc_exponent));
            }

            emu_qc(icol,ilev)  = qc_pack[lane];
            emu_nc(icol,ilev)  = nc_in;
            emu_qr(icol,ilev)  = qr_pack[lane];
            emu_nr(icol,ilev)  = nr_pack[lane];
            emu_rho(icol,ilev) = dpres_pack[lane] / (g * dz_pack[lane]);
          }
        }
      });
    });
  Kokkos::fence();

  emu_qc_field.sync_to_host();
  emu_nc_field.sync_to_host();
  emu_qr_field.sync_to_host();
  emu_nr_field.sync_to_host();
  emu_rho_field.sync_to_host();

  py_module_call(
    "forward",
    get_py_field_host("sdm_warm_emulator_qc"),
    get_py_field_host("sdm_warm_emulator_nc"),
    get_py_field_host("sdm_warm_emulator_qr"),
    get_py_field_host("sdm_warm_emulator_nr"),
    get_py_field_host("sdm_warm_emulator_rho"),
    get_py_field_host("sdm_warm_emulator_qc2qr_autoconv_tend"),
    get_py_field_host("sdm_warm_emulator_qc2qr_accret_tend"),
    get_py_field_host("sdm_warm_emulator_ncautr"),
    get_py_field_host("sdm_warm_emulator_nc2nr_autoconv_tend"),
    get_py_field_host("sdm_warm_emulator_nc_accret_tend"),
    get_py_field_host("sdm_warm_emulator_nc_selfcollect_tend"),
    get_py_field_host("sdm_warm_emulator_nr_selfcollect_tend"),
    get_py_field_host("sdm_warm_emulator_use_cloud"),
    get_py_field_host("sdm_warm_emulator_use_rain"));

  get_internal_field("sdm_warm_emulator_qc2qr_autoconv_tend").sync_to_dev();
  get_internal_field("sdm_warm_emulator_qc2qr_accret_tend").sync_to_dev();
  get_internal_field("sdm_warm_emulator_ncautr").sync_to_dev();
  get_internal_field("sdm_warm_emulator_nc2nr_autoconv_tend").sync_to_dev();
  get_internal_field("sdm_warm_emulator_nc_accret_tend").sync_to_dev();
  get_internal_field("sdm_warm_emulator_nc_selfcollect_tend").sync_to_dev();
  get_internal_field("sdm_warm_emulator_nr_selfcollect_tend").sync_to_dev();
  get_internal_field("sdm_warm_emulator_use_cloud").sync_to_dev();
  get_internal_field("sdm_warm_emulator_use_rain").sync_to_dev();
#else
  EKAT_ERROR_MSG(
      "[P3Microphysics] Error! use_sdm_warm_emulator=true requires "
      "EAMXX_ENABLE_PYTHON=ON.\n");
#endif
}

} // namespace scream
