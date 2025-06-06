#include "diagnostics/longwave_cloud_forcing.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

namespace scream
{

// =========================================================================================
LongwaveCloudForcingDiagnostic::LongwaveCloudForcingDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void LongwaveCloudForcingDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  Units m2 (m*m,"m2");

  auto grid  = grids_manager->get_grid("physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar2d_layout_col{ {COL}, {m_num_cols} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  // The fields required for this diagnostic to be computed
  add_field<Required>("LW_flux_up",        scalar3d_layout_mid, W/m2,  grid_name);
  add_field<Required>("LW_clrsky_flux_up", scalar3d_layout_mid, W/m2,  grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar2d_layout_col, W/m2, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation();
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void LongwaveCloudForcingDiagnostic::compute_diagnostic_impl()
{
  using KT         = KokkosTypes<DefaultDevice>;
  using ESU        = ekat::ExeSpaceUtils<KT::ExeSpace>;
  using MemberType = typename KT::MemberType;

  const auto default_policy = ESU::get_default_team_policy(m_num_cols,1);

  const auto& LWCF              = m_diagnostic_output.get_view<Real*>();
  const auto& LW_flux_up        = get_field_in("LW_flux_up").get_view<const Real**>();
  const auto& LW_clrsky_flux_up = get_field_in("LW_clrsky_flux_up").get_view<const Real**>();

  // NOTE: as part of fixing https://github.com/E3SM-Project/E3SM/issues/6803, this hack
  //       may have to be revised. Namely, the "radiation_ran" extra data should be removed,
  //       in favor of a more structured and uniform approach
  const auto radiation_ran = get_field_in("LW_flux_up").get_header().get_extra_data<bool>("radiation_ran");
  if (not radiation_ran) {
    m_diagnostic_output.deep_copy(constants::DefaultFillValue<Real>().value);
    return;
  }

  Kokkos::parallel_for("LongwaveCloudForcingDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
    LWCF(icol) =  LW_clrsky_flux_up(icol,0) -  LW_flux_up(icol,0) ;
  });
  Kokkos::fence();
}

} //namespace scream
