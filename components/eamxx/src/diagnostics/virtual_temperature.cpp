#include "diagnostics/virtual_temperature.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"

namespace scream
{

VirtualTemperatureDiagnostic::
VirtualTemperatureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

void VirtualTemperatureDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  auto scalar3d = grid->get_3d_scalar_layout(true);

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid", scalar3d, K,     grid_name);
  add_field<Required>("qv",    scalar3d, kg/kg, grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d, K, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void VirtualTemperatureDiagnostic::compute_diagnostic_impl()
{
  using PF = scream::PhysicsFunctions<DefaultDevice>;

  const auto& virtualT = m_diagnostic_output.get_view<Real**>();
  const auto& T_mid    = get_field_in("T_mid").get_view<const Real**>();
  const auto& qv_mid   = get_field_in("qv").get_view<const Real**>();

  int nlevs = m_num_levs;
  Kokkos::parallel_for("VirtualTemperatureDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*nlevs),
                       KOKKOS_LAMBDA(const int& idx) {
      const int icol = idx / nlevs;
      const int ilev = idx % nlevs;
      virtualT(icol,ilev) = PF::calculate_virtual_temperature(T_mid(icol,ilev),qv_mid(icol,ilev));
  });
  Kokkos::fence();
}

} //namespace scream
