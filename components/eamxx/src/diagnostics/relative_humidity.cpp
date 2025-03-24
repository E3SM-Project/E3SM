#include "diagnostics/relative_humidity.hpp"
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs
#include "physics/share/physics_saturation_impl.hpp"

namespace scream
{

RelativeHumidityDiagnostic::
RelativeHumidityDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

void RelativeHumidityDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto nondim = Units::nondimensional();

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  auto scalar3d = grid->get_3d_scalar_layout(true);

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid",              scalar3d, K,     grid_name, SCREAM_PACK_SIZE);
  add_field<Required>("p_dry_mid",          scalar3d, Pa,    grid_name, SCREAM_PACK_SIZE);
  add_field<Required>("qv",                 scalar3d, kg/kg, grid_name, SCREAM_PACK_SIZE);
  add_field<Required>("pseudo_density",     scalar3d, Pa,    grid_name, SCREAM_PACK_SIZE);
  add_field<Required>("pseudo_density_dry", scalar3d, Pa,    grid_name, SCREAM_PACK_SIZE);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d, nondim, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(SCREAM_PACK_SIZE);
  m_diagnostic_output.allocate_view();
}

void RelativeHumidityDiagnostic::compute_diagnostic_impl()
{
  using Pack          = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  const auto npacks  = ekat::npack<Pack>(m_num_levs);
  auto theta     = m_diagnostic_output.get_view<Pack**>();
  auto T_mid     = get_field_in("T_mid").get_view<const Pack**>();
  auto p_dry_mid = get_field_in("p_dry_mid").get_view<const Pack**>();
  auto dp_wet    = get_field_in("pseudo_density").get_view<const Pack**>();
  auto dp_dry    = get_field_in("pseudo_density_dry").get_view<const Pack**>();
  auto qv_mid    = get_field_in("qv").get_view<const Pack**>();
  const auto& RH = m_diagnostic_output.get_view<Pack**>();

  using physics = scream::physics::Functions<Real, DefaultDevice>;

  Int num_levs = m_num_levs;
  Kokkos::parallel_for("RelativeHumidityDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*npacks),
                       KOKKOS_LAMBDA (const int& idx) {
      const int icol  = idx / npacks;
      const int jpack = idx % npacks;
      const auto range_pack = ekat::range<Pack>(jpack*Pack::n);
      const auto range_mask = range_pack < num_levs;
      auto qv_sat_l = physics::qv_sat_wet(T_mid(icol,jpack),  p_dry_mid(icol,jpack), true, range_mask, dp_wet(icol,jpack), dp_dry(icol,jpack),
                                           physics::MurphyKoop, "RelativeHumidityDiagnostic::compute_diagnostic_impl"); 
      RH(icol,jpack) = qv_mid(icol,jpack)/qv_sat_l;

  });
  Kokkos::fence();
}

} //namespace scream
