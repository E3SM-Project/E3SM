#include "physics/p3/eamxx_p3_process_interface.hpp"

namespace scream {

void P3Microphysics::run_impl (const double dt)
{
  // Set the dt for p3 postprocessing
  p3_postproc.m_dt = dt;

  // Assign values to local arrays used by P3, these are now stored in p3_loc.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_preproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();

  // allow namelist flags to override sub-grid cloud fraction (set it to 1 everywhere)
  if (m_params.get<bool>("set_cld_frac_l_to_one", false)) {
    auto& cld_frac_l = p3_preproc.cld_frac_l;
    Kokkos::deep_copy(cld_frac_l,1.0);
  }
  if (m_params.get<bool>("set_cld_frac_r_to_one", false)) {
    auto& cld_frac_r = p3_preproc.cld_frac_r;
    Kokkos::deep_copy(cld_frac_r,1.0);
  }
  if (m_params.get<bool>("set_cld_frac_i_to_one", false)) {
    auto& cld_frac_i = p3_preproc.cld_frac_i;
    Kokkos::deep_copy(cld_frac_i,1.0);
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

  P3F::p3_main(runtime_options, prog_state, diag_inputs, diag_outputs, infrastructure,
               history_only, lookup_tables,
#ifdef SCREAM_P3_SMALL_KERNELS
               temporaries,
#endif
               workspace_mgr, m_num_cols, m_num_levs);

  // Conduct the post-processing of the p3_main output.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_postproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();
}

} // namespace scream
