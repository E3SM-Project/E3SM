#include "diagnostics/vertical_layer_midpoint.hpp"

namespace scream
{

// =========================================================================================
VerticalLayerMidpointDiagnostic::VerticalLayerMidpointDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void VerticalLayerMidpointDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto Q = kg/kg;
  Q.set_string("kg/kg");

  const auto& grid_name = m_params.get<std::string>("Grid");
  auto grid  = grids_manager->get_grid(grid_name);
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

}
// =========================================================================================
void VerticalLayerMidpointDiagnostic::initialize_impl(const RunType /* run_type */)
{

  auto ts = timestamp(); 
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);

}
// =========================================================================================
void VerticalLayerMidpointDiagnostic::run_impl(const int /* dt */)
{

  const auto& T_mid              = get_field_in("T_mid").get_view<const Pack**>();
  const auto& p_mid              = get_field_in("p_mid").get_view<const Pack**>();
  const auto& qv_mid             = get_field_in("qv").get_view<const Pack**>();
  const auto& pseudo_density_mid = get_field_in("pseudo_density").get_view<const Pack**>();

  const auto& output             = m_diagnostic_output.get_view<Pack**>();

  // Set surface geopotential for this diagnostic
  const Real surf_geopotential = 0.0;

  const auto nlev_packs     = ekat::npack<Spack>(m_num_levs);
  const auto nlev_packs_p1  = ekat::npack<Spack>(m_num_levs+1);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);
  view_1d dz("",nlev_packs);
  view_1d z_int("",nlev_packs_p1);
  Kokkos::parallel_for("VerticalLayerInterfaceDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
      dz(k) = PF::calculate_dz(pseudo_density_mid(i,k), p_mid(i,k), T_mid(i,k), qv_mid(i,k));
    });
    const auto& output_s = ekat::subview(output, i);
    PF::calculate_z_int(team,m_num_levs,dz,surf_geopotential,z_int);
    PF::calculate_z_mid(team,m_num_levs,z_int,output_s);
  });

}
// =========================================================================================
void VerticalLayerMidpointDiagnostic::finalize_impl()
{
  // Nothing to do
}
// =========================================================================================
} //namespace scream
