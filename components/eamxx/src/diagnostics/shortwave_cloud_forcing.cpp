#include "diagnostics/shortwave_cloud_forcing.hpp"

namespace scream
{

// =========================================================================================
ShortwaveCloudForcingDiagnostic::ShortwaveCloudForcingDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void ShortwaveCloudForcingDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto m2 = m*m;

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar2d_layout_col{ {COL}, {m_num_cols} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("SW_flux_dn",          scalar3d_layout_mid, W/m2,  grid_name, ps);
  add_field<Required>("SW_flux_up",          scalar3d_layout_mid, W/m2,  grid_name, ps);
  add_field<Required>("SW_clrsky_flux_dn",          scalar3d_layout_mid, W/m2,  grid_name, ps);
  add_field<Required>("SW_clrsky_flux_up",          scalar3d_layout_mid, W/m2,  grid_name, ps);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar2d_layout_col, W/m2, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation();
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void ShortwaveCloudForcingDiagnostic::compute_diagnostic_impl()
{
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols,1);

  const auto& SWCF              = m_diagnostic_output.get_view<Real*>();
  const auto& SW_flux_dn        = get_field_in("SW_flux_dn").get_view<const Pack**>();
  const auto& SW_flux_up        = get_field_in("SW_flux_up").get_view<const Pack**>();
  const auto& SW_clrsky_flux_dn = get_field_in("SW_clrsky_flux_dn").get_view<const Pack**>();
  const auto& SW_clrsky_flux_up = get_field_in("SW_clrsky_flux_up").get_view<const Pack**>();

  Kokkos::parallel_for("ShortwaveCloudForcingDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
    SWCF(icol) = (SW_flux_dn(icol,0)[0] - SW_flux_up(icol,0)[0]) - (SW_clrsky_flux_dn(icol,0)[0] - SW_clrsky_flux_up(icol,0)[0]);
  });
  Kokkos::fence();

  const auto ts = get_field_in("SW_flux_dn").get_header().get_tracking().get_time_stamp();
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
} //namespace scream
