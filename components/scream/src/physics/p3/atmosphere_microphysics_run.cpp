#include "physics/p3/atmosphere_microphysics.hpp"

namespace scream {

void P3Microphysics::run_impl (const int dt)
{
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

  // Run p3 main
  P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
               history_only, workspace_mgr, m_num_cols, m_num_levs);

  // Reset internal WSM variables so it can be used in the next timestep.
  workspace_mgr.reset_internals();

  // Conduct the post-processing of the p3_main output.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_postproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();
}

} // namespace scream
