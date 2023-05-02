#include "diagnostics/zonal_vapor_flux.hpp"

namespace scream
{

// =========================================================================================
ZonalVapFluxDiagnostic::ZonalVapFluxDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void ZonalVapFluxDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto Q = kg/kg;
  Q.set_string("kg/kg");

  auto vel = m/s;
  vel.set_string("m/s");

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  FieldLayout scalar2d_layout_mid { {COL},     {m_num_cols}            };
  FieldLayout horiz_wind_layout { {COL,CMP,LEV}, {m_num_cols,2,m_num_levs} };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,  grid_name, ps);
  add_field<Required>("qv",             scalar3d_layout_mid, Q,   grid_name, "tracers", ps);
  // Note both u and v are packaged into the single field horiz_winds.
  add_field<Required>("horiz_winds",    horiz_wind_layout,   m/s, grid_name, ps);


  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar2d_layout_mid, kg/m/s, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation();
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void ZonalVapFluxDiagnostic::compute_diagnostic_impl()
{

  using PC         = scream::physics::Constants<Real>;
  constexpr Real gravit = PC::gravit;
  const auto npacks         = ekat::npack<Pack>(m_num_levs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, npacks);
  const auto& qv_vert_integrated_flux_u = m_diagnostic_output.get_view<Real*>();
  const auto& qv_mid             = get_field_in("qv").get_view<const Pack**>();
  const auto& pseudo_density_mid = get_field_in("pseudo_density").get_view<const Pack**>();
  const auto& horiz_winds        = get_field_in("horiz_winds").get_view<const Pack***>();

  const auto num_levs = m_num_levs;
  Kokkos::parallel_for("ZonalVapFluxDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, num_levs), [&] (const Int& idx, Real& lsum) {
      const int jpack = idx / Pack::n;
      const int klev  = idx % Pack::n;
      // Note, horiz_winds contains u (index 0) and v (index 1).  Here we want u
      lsum += horiz_winds(icol,0,jpack)[klev] * qv_mid(icol,jpack)[klev] * pseudo_density_mid(icol,jpack)[klev]/gravit;
    },qv_vert_integrated_flux_u(icol));
    team.team_barrier();
  });

  const auto ts = get_field_in("qv").get_header().get_tracking().get_time_stamp();
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
} //namespace scream
