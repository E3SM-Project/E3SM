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
