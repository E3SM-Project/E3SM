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

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
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

  // Initialize a 2d view of dz to be used in compute_diagnostic
  const auto npacks     = ekat::npack<Pack>(m_num_levs);
  const auto npacks_p1  = ekat::npack<Pack>(m_num_levs+1);
  m_tmp_mid = view_2d("",m_num_cols,npacks);
  m_tmp_int = view_2d("",m_num_cols,npacks_p1);
}
// =========================================================================================
void DryStaticEnergyDiagnostic::compute_diagnostic_impl()
{

  const auto npacks     = ekat::npack<Pack>(m_num_levs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_num_cols, npacks);

  const auto& dse                = m_diagnostic_output.get_view<Pack**>();
  const auto& T_mid              = get_field_in("T_mid").get_view<const Pack**>();
  const auto& p_mid              = get_field_in("p_mid").get_view<const Pack**>();
  const auto& qv_mid             = get_field_in("qv").get_view<const Pack**>();
  const auto& pseudo_density_mid = get_field_in("pseudo_density").get_view<const Pack**>();
  const auto& phis               = get_field_in("phis").get_view<const Real*>();

  // Set surface geopotential for this diagnostic
  const Real surf_geopotential = 0.0;

  const int num_levs = m_num_levs;
  auto      tmp_mid  = m_tmp_mid;
  auto      tmp_int  = m_tmp_int;
  Kokkos::deep_copy(tmp_mid,0.0);
  Kokkos::deep_copy(tmp_int,0.0);
  Kokkos::parallel_for("DryStaticEnergyDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
    const auto& dz_s    = ekat::subview(tmp_mid,icol);
    const auto& z_int_s = ekat::subview(tmp_int,icol);
    const auto& z_mid_s = dz_s; // Reuse the memory for z_mid, but set a new variable for code readability.
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, npacks), [&] (const Int& jpack) {
      dz_s(jpack) = PF::calculate_dz(pseudo_density_mid(icol,jpack), p_mid(icol,jpack), T_mid(icol,jpack), qv_mid(icol,jpack));
    });
    team.team_barrier();
    PF::calculate_z_int(team,num_levs,dz_s,surf_geopotential,z_int_s);
    team.team_barrier();
    PF::calculate_z_mid(team,num_levs,z_int_s,z_mid_s);
    team.team_barrier();
    const auto& dse_s = ekat::subview(dse,icol);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, npacks), [&] (const Int& jpack) {
      dse_s(jpack) = PF::calculate_dse(T_mid(icol,jpack),z_mid_s(jpack),phis(icol));
    });
    team.team_barrier();
  });
  Kokkos::fence();

  const auto ts = get_field_in("T_mid").get_header().get_tracking().get_time_stamp();
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
} //namespace scream
