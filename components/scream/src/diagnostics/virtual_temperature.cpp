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

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
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
void VirtualTemperatureDiagnostic::compute_diagnostic_impl()
{

  const auto npacks  = ekat::npack<Pack>(m_num_levs);
  const auto& virtualT = m_diagnostic_output.get_view<Pack**>();
  const auto& T_mid    = get_field_in("T_mid").get_view<const Pack**>();
  const auto& qv_mid   = get_field_in("qv").get_view<const Pack**>();

  Kokkos::parallel_for("VirtualTemperatureDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*npacks),
                       KOKKOS_LAMBDA(const int& idx) {
      const int icol  = idx / npacks;
      const int jpack = idx % npacks;
      virtualT(icol,jpack) = PF::calculate_virtual_temperature(T_mid(icol,jpack),qv_mid(icol,jpack));
  });
  Kokkos::fence();

  const auto ts = get_field_in("qv").get_header().get_tracking().get_time_stamp();
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
} //namespace scream
