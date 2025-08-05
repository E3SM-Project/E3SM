#include "physics/p3/eamxx_p3_process_interface.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream {

void P3Microphysics::run_impl (const double dt)
{
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  // Set the dt for p3 postprocessing
  p3_postproc.m_dt = dt;

  // Create policy for pre and post process pfor
  const auto nlev_packs  = ekat::npack<Spack>(m_num_levs);
  const auto policy = TPF::get_default_team_policy(m_num_cols, nlev_packs);

  // Assign values to local arrays used by P3, these are now stored in p3_loc.
  Kokkos::parallel_for(
    "p3_pre_process",
    policy,
    p3_preproc
  );
  Kokkos::fence();

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

} // namespace scream
