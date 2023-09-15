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

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
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
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(ps);
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void PotentialTemperatureDiagnostic::compute_diagnostic_impl()
{

  const auto npacks  = ekat::npack<Pack>(m_num_levs);
  auto theta = m_diagnostic_output.get_view<Pack**>();
  auto T_mid = get_field_in("T_mid").get_view<const Pack**>();
  auto p_mid = get_field_in("p_mid").get_view<const Pack**>();
  Kokkos::parallel_for("PotentialTemperatureDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*npacks),
                       KOKKOS_LAMBDA (const int& idx) {
      const int icol  = idx / npacks;
      const int jpack = idx % npacks;
      theta(icol,jpack) = PF::calculate_theta_from_T(T_mid(icol,jpack),p_mid(icol,jpack));
  });
  Kokkos::fence();
}
// =========================================================================================

} //namespace scream
