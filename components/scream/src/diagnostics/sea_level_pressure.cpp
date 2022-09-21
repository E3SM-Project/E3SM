#include "diagnostics/sea_level_pressure.hpp"

namespace scream
{

// =========================================================================================
SeaLevelPressureDiagnostic::SeaLevelPressureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void SeaLevelPressureDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto Q = kg/kg;
  Q.set_string("kg/kg");
  const auto m2 = m*m;
  const auto s2 = s*s;

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar2d_layout_col{ {COL}, {m_num_cols} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid",          scalar3d_layout_mid, K,  grid_name, ps);
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Required>("phis",           scalar2d_layout_col, m2/s2, grid_name, ps);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar2d_layout_col, Pa, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation();
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void SeaLevelPressureDiagnostic::compute_diagnostic_impl()
{

  const auto npacks     = ekat::npack<Pack>(m_num_levs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols,1);

  const auto& p_sealevel         = m_diagnostic_output.get_view<Real*>();
  const auto& T_mid              = get_field_in("T_mid").get_view<const Pack**>();
  const auto& p_mid              = get_field_in("p_mid").get_view<const Pack**>();
  const auto& phis               = get_field_in("phis").get_view<const Real*>();

  const int pack_surf = std::min(m_num_levs / Pack::n, npacks-1);
  const int idx_surf  = m_num_levs % Pack::n;
  Kokkos::parallel_for("SeaLevelPressureDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
    p_sealevel(icol) = PF::calculate_psl(T_mid(icol,pack_surf)[idx_surf],p_mid(icol,pack_surf)[idx_surf],phis(icol));
  });
  Kokkos::fence();

  const auto ts = get_field_in("T_mid").get_header().get_tracking().get_time_stamp();
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
} //namespace scream
