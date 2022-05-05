#include "diagnostics/sea_level_pressure.hpp"

namespace scream
{

// =========================================================================================
SeaLevelPressureDiagnostic::SeaLevelPressureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void SeaLevelPressureDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto m2 = m*m;
  const auto s2 = s*s;

  const auto& grid_name = m_params.get<std::string>("Grid");
  auto grid  = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar2d_layout_col { {COL}, {m_num_cols} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid",          scalar3d_layout_mid, K,  grid_name, ps);
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Required>("phis",           scalar2d_layout_col, m2/s2, grid_name, ps);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar2d_layout_col, m, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(ps);
  m_diagnostic_output.allocate_view();

}
// =========================================================================================
void SeaLevelPressureDiagnostic::initialize_impl(const RunType /* run_type */)
{
  const auto& T_mid              = get_field_in("T_mid").get_view<const Real**>();
  const auto& p_mid              = get_field_in("p_mid").get_view<const Real**>();
  const auto& phis               = get_field_in("phis").get_view<const Real*>();

  const auto& output             = m_diagnostic_output.get_view<Real*>();

  // Set surface geopotential for this diagnostic
  const Real surf_geopotential = 0.0;

  auto ts = timestamp(); 
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);

  run_diagnostic.set_variables(m_num_cols,m_num_levs,T_mid,p_mid,phis,output);
}
// =========================================================================================
void SeaLevelPressureDiagnostic::run_impl(const int /* dt */)
{

  Kokkos::parallel_for("SeaLevelPressureDiagnostic",
                       m_num_cols,
                       run_diagnostic
  );
  Kokkos::fence();

}
// =========================================================================================
void SeaLevelPressureDiagnostic::finalize_impl()
{
  // Nothing to do
}
// =========================================================================================
} //namespace scream
