#include "diagnostics/atm_density.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"

namespace scream
{

AtmDensityDiagnostic::
AtmDensityDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

void AtmDensityDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  // Boiler Plate
  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  // Set Field Layouts
  auto scalar3d = grid->get_3d_scalar_layout(true);

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid",          scalar3d, K,     grid_name);
  add_field<Required>("pseudo_density", scalar3d, Pa,    grid_name);
  add_field<Required>("p_mid",          scalar3d, Pa,    grid_name);
  add_field<Required>("qv",             scalar3d, kg/kg, grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d, kg/(m*m*m), grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void AtmDensityDiagnostic::compute_diagnostic_impl()
{
  using PF = scream::PhysicsFunctions<DefaultDevice>;

  const auto atm_dens = m_diagnostic_output.get_view<Real**>();
  const auto T_mid    = get_field_in("T_mid").get_view<const Real**>();
  const auto p_mid    = get_field_in("p_mid").get_view<const Real**>();
  const auto qv_mid   = get_field_in("qv").get_view<const Real**>();
  const auto pseudo_density_mid = get_field_in("pseudo_density").get_view<const Real**>();

  int nlevs = m_num_levs;
  Kokkos::parallel_for("AtmosphereDensityDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*nlevs),
                       KOKKOS_LAMBDA(const int& idx) {
      const int icol  = idx / nlevs;
      const int ilev = idx % nlevs;
      auto dz = PF::calculate_dz(pseudo_density_mid(icol,ilev),p_mid(icol,ilev),T_mid(icol,ilev),qv_mid(icol,ilev));
      atm_dens(icol,ilev) = PF::calculate_density(pseudo_density_mid(icol,ilev),dz);
  });
  Kokkos::fence();
}

} //namespace scream
