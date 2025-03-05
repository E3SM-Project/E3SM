#include "diagnostics/shortwave_cloud_forcing.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

namespace scream
{

ShortwaveCloudForcingDiagnostic::
ShortwaveCloudForcingDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

void ShortwaveCloudForcingDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  Units m2 (m*m,"m2");

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  auto scalar2d = grid->get_2d_scalar_layout();
  auto scalar3d = grid->get_3d_scalar_layout(true);

  // The fields required for this diagnostic to be computed
  add_field<Required>("SW_flux_dn",        scalar3d, W/m2, grid_name);
  add_field<Required>("SW_flux_up",        scalar3d, W/m2, grid_name);
  add_field<Required>("SW_clrsky_flux_dn", scalar3d, W/m2, grid_name);
  add_field<Required>("SW_clrsky_flux_up", scalar3d, W/m2, grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar2d, W/m2, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void ShortwaveCloudForcingDiagnostic::compute_diagnostic_impl()
{
  using KT         = KokkosTypes<DefaultDevice>;
  using ESU        = ekat::ExeSpaceUtils<KT::ExeSpace>;
  using MemberType = typename KT::MemberType;

  const auto default_policy = ESU::get_default_team_policy(m_num_cols,1);

  const auto& SWCF              = m_diagnostic_output.get_view<Real*>();
  const auto& SW_flux_dn        = get_field_in("SW_flux_dn").get_view<const Real**>();
  const auto& SW_flux_up        = get_field_in("SW_flux_up").get_view<const Real**>();
  const auto& SW_clrsky_flux_dn = get_field_in("SW_clrsky_flux_dn").get_view<const Real**>();
  const auto& SW_clrsky_flux_up = get_field_in("SW_clrsky_flux_up").get_view<const Real**>();

  Kokkos::parallel_for("ShortwaveCloudForcingDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
    SWCF(icol) = (SW_flux_dn(icol,0) - SW_flux_up(icol,0)) - (SW_clrsky_flux_dn(icol,0) - SW_clrsky_flux_up(icol,0));
  });
  Kokkos::fence();
}

} //namespace scream
