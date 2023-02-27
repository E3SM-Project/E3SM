#include "diagnostics/vertical_layer_midpoint.hpp"

namespace scream
{

// =========================================================================================
VerticalLayerMidpointDiagnostic::
VerticalLayerMidpointDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void VerticalLayerMidpointDiagnostic::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
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
  add_field<Required>("T_mid",          scalar3d_layout_mid, K,  grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Required>("qv",             scalar3d_layout_mid, Q,  grid_name, "tracers", ps);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d_layout_mid, m, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(ps);
  m_diagnostic_output.allocate_view();

  // Initialize 2d view of dz to be used in compute_diagnostic
  const auto npacks_p1  = ekat::npack<Pack>(m_num_levs+1);
  m_z_int = view_2d("",m_num_cols,npacks_p1);
}
// =========================================================================================
void VerticalLayerMidpointDiagnostic::compute_diagnostic_impl()
{

  const auto npacks     = ekat::npack<Pack>(m_num_levs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_num_cols, npacks);
  const auto& z_mid              = m_diagnostic_output.get_view<Pack**>();
  const auto& T_mid              = get_field_in("T_mid").get_view<const Pack**>();
  const auto& p_mid              = get_field_in("p_mid").get_view<const Pack**>();
  const auto& qv_mid             = get_field_in("qv").get_view<const Pack**>();
  const auto& pseudo_density_mid = get_field_in("pseudo_density").get_view<const Pack**>();

  // Set surface geopotential for this diagnostic
  const Real surf_geopotential = 0.0;

  const int num_levs = m_num_levs;
  auto z_int = m_z_int;
  Kokkos::parallel_for("VerticalLayerMidpointDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
    const auto& z_mid_s = ekat::subview(z_mid, icol);
    const auto& dz_s    = z_mid_s; // Use the memory in z_mid for dz, since we don't set z_mid until after dz is no longer needed.
    const auto& z_int_s = ekat::subview(z_int, icol);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, npacks), [&] (const Int& jpack) {
      dz_s(jpack) = PF::calculate_dz(pseudo_density_mid(icol,jpack), p_mid(icol,jpack), T_mid(icol,jpack), qv_mid(icol,jpack));
    });
    team.team_barrier();
    PF::calculate_z_int(team,num_levs,dz_s,surf_geopotential,z_int_s);
    PF::calculate_z_mid(team,num_levs,z_int_s,z_mid_s);
  });

  const auto ts = get_field_in("qv").get_header().get_tracking().get_time_stamp();
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
} //namespace scream
