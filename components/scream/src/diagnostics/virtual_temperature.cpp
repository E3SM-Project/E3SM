#include "diagnostics/virtual_temperature.hpp"

namespace scream
{

// =========================================================================================
VirtualTemperatureDiagnostic::VirtualTemperatureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void VirtualTemperatureDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto Q = kg/kg;
  Q.set_string("kg/kg");

  const auto& grid_name = m_params.get<std::string>("Grid");
  auto grid  = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid",       scalar3d_layout_mid, K,  grid_name, ps);
  add_field<Required>("qv",          scalar3d_layout_mid, Q,  grid_name, "tracers", ps);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d_layout_mid, K, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(ps);
  m_diagnostic_output.allocate_view();

}
// =========================================================================================
void VirtualTemperatureDiagnostic::initialize_impl(const RunType /* run_type */)
{
  const auto& T_mid          = get_field_in("T_mid").get_view<const Pack**>();
  const auto& qv_mid         = get_field_in("qv").get_view<const Pack**>();

  const auto& output         = m_diagnostic_output.get_view<Pack**>();

  auto ts = timestamp(); 
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);

  const auto nk_pack  = ekat::npack<Spack>(m_num_levs);

  run_diagnostic.set_variables(m_num_cols,nk_pack,qv_mid,T_mid,output);
}
// =========================================================================================
void VirtualTemperatureDiagnostic::run_impl(const int /* dt */)
{

  const auto nk_pack  = ekat::npack<Spack>(m_num_levs);
  Kokkos::parallel_for("VirtualTemperatureDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*nk_pack),
                       run_diagnostic
  );
  Kokkos::fence();

}
// =========================================================================================
void VirtualTemperatureDiagnostic::finalize_impl()
{
  // Nothing to do
}
// =========================================================================================
} //namespace scream
