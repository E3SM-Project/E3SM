#include "diagnostics/exner.hpp"

namespace scream
{

// =========================================================================================
ExnerDiagnostic::ExnerDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void ExnerDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  Units nondim(0,0,0,0,0,0,0);

  const auto& grid_name = m_params.get<std::string>("Grid");
  auto grid  = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa, grid_name, ps);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d_layout_mid, nondim, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(ps);
  m_diagnostic_output.allocate_view();

}
// =========================================================================================
void ExnerDiagnostic::initialize_impl(const RunType /* run_type */)
{
  auto ts = timestamp(); 
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
void ExnerDiagnostic::run_impl(const int /* dt */)
{

  const auto& p_mid          = get_field_in("p_mid").get_view<const Pack**>();
  const auto& output         = m_diagnostic_output.get_view<Pack**>();

  const auto nk_pack  = ekat::npack<Spack>(m_num_levs);
  Kokkos::parallel_for("ExnerDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*nk_pack),
                       KOKKOS_LAMBDA(int idx) {
      const int icol  = idx / nk_pack;
      const int jpack = idx % nk_pack;
      output(icol,jpack) = PF::exner_function(p_mid(icol,jpack));
  });
  Kokkos::fence();

}
// =========================================================================================
void ExnerDiagnostic::finalize_impl()
{
  // Nothing to do
}
// =========================================================================================
} //namespace scream
