#include "diagnostics/potential_temperature.hpp"

namespace scream
{

// =========================================================================================
PotentialTemperatureDiagnostic::PotentialTemperatureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void PotentialTemperatureDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
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
  add_field<Required>("T_mid",          scalar3d_layout_mid, K,  grid_name, ps);
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa, grid_name, ps);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d_layout_mid, K, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();

}
// =========================================================================================
void PotentialTemperatureDiagnostic::initialize_impl(const RunType /* run_type */)
{
  const auto& T_mid          = get_field_in("T_mid").get_view<const Pack**>();
  const auto& p_mid          = get_field_in("p_mid").get_view<const Pack**>();

  const auto& output         = m_diagnostic_output.get_view<Pack**>();

  const auto nk_pack  = ekat::npack<Spack>(m_num_levs);

  run_diagnostic.set_variables(m_num_cols,nk_pack,p_mid,T_mid,output);
}
// =========================================================================================
void PotentialTemperatureDiagnostic::run_impl(const int /* dt */)
{

  const auto nk_pack  = ekat::npack<Spack>(m_num_levs);
  Kokkos::parallel_for("PotentialTemperatureDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*nk_pack),
                       run_diagnostic
  );
  Kokkos::fence();

}
// =========================================================================================
void PotentialTemperatureDiagnostic::finalize_impl()
{
  // Nothing to do
}
// =========================================================================================
} //namespace scream
