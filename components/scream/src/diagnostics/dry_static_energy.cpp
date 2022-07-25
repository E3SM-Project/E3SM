#include "diagnostics/dry_static_energy.hpp"

namespace scream
{

// =========================================================================================
DryStaticEnergyDiagnostic::DryStaticEnergyDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void DryStaticEnergyDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto Q = kg/kg;
  Q.set_string("kg/kg");
  const auto m2 = m*m;
  const auto s2 = s*s;

  const auto& grid_name = m_params.get<std::string>("Grid");
  auto grid  = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar2d_layout_col{ {COL}, {m_num_cols} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid",          scalar3d_layout_mid, K,  grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Required>("qv",             scalar3d_layout_mid, Q,  grid_name, "tracers", ps);
  add_field<Required>("phis",           scalar2d_layout_col, m2/s2, grid_name, ps);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d_layout_mid, m2/s2, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(ps);
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void DryStaticEnergyDiagnostic::initialize_impl(const RunType /* run_type */)
{

  auto ts = timestamp(); 
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);

}
// =========================================================================================
void DryStaticEnergyDiagnostic::compute_diagnostic_impl()
{

  const auto npacks     = ekat::npack<Pack>(m_num_levs);
  const auto npacks_p1  = ekat::npack<Pack>(m_num_levs+1);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, npacks);

  const auto& dse                = m_diagnostic_output.get_view<Pack**>();
  const auto& T_mid              = get_field_in("T_mid").get_view<const Pack**>();
  const auto& p_mid              = get_field_in("p_mid").get_view<const Pack**>();
  const auto& qv_mid             = get_field_in("qv").get_view<const Pack**>();
  const auto& pseudo_density_mid = get_field_in("pseudo_density").get_view<const Pack**>();
  const auto& phis               = get_field_in("phis").get_view<const Real*>();

  // Set surface geopotential for this diagnostic
  const Real surf_geopotential = 0.0;

  const int num_levs = m_num_levs;
  view_1d dz("",npacks);
  view_1d z_int("",npacks_p1);
  view_1d z_mid("",npacks);
  Kokkos::deep_copy(dz,0.0);
  Kokkos::deep_copy(z_int,0.0);
  Kokkos::deep_copy(z_mid,0.0);
  Kokkos::parallel_for("DryStaticEnergyDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, npacks), [&] (const Int& jpack) {
      dz(jpack) = PF::calculate_dz(pseudo_density_mid(icol,jpack), p_mid(icol,jpack), T_mid(icol,jpack), qv_mid(icol,jpack));
    });
    team.team_barrier();
    PF::calculate_z_int(team,num_levs,dz,surf_geopotential,z_int);
    team.team_barrier();
    PF::calculate_z_mid(team,num_levs,z_int,z_mid);
    team.team_barrier();
    const auto& dse_s = ekat::subview(dse,icol);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, npacks), [&] (const Int& jpack) {
      dse_s(jpack) = PF::calculate_dse(T_mid(icol,jpack),z_mid(jpack),phis(icol));
    });
    team.team_barrier();
  });
  Kokkos::fence();

}
// =========================================================================================
} //namespace scream
