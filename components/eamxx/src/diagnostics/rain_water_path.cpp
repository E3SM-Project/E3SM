#include "diagnostics/rain_water_path.hpp"

namespace scream
{

// =========================================================================================
RainWaterPathDiagnostic::RainWaterPathDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void RainWaterPathDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
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
  FieldLayout scalar2d_layout_mid { {COL},     {m_num_cols}            };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Required>("qr",             scalar3d_layout_mid, Q,  grid_name, "tracers", ps);

  // Construct and allocate the diagnostic field
  const auto m2 = m*m;
  FieldIdentifier fid (name(), scalar2d_layout_mid, kg/m2, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation();
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void RainWaterPathDiagnostic::compute_diagnostic_impl()
{

  using PC         = scream::physics::Constants<Real>;
  constexpr Real gravit = PC::gravit;
  const auto npacks         = ekat::npack<Pack>(m_num_levs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, npacks);
  const auto& lwp                = m_diagnostic_output.get_view<Real*>();
  const auto& qr_mid             = get_field_in("qr").get_view<const Pack**>();
  const auto& pseudo_density_mid = get_field_in("pseudo_density").get_view<const Pack**>();

  const auto num_levs = m_num_levs;
  Kokkos::parallel_for("RainWaterPathDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, num_levs), [&] (const Int& idx, Real& lsum) {
      const int jpack = idx / Pack::n;
      const int klev  = idx % Pack::n;
      lsum += qr_mid(icol,jpack)[klev] * pseudo_density_mid(icol,jpack)[klev]/gravit;
    },lwp(icol));
    team.team_barrier();
  });

  const auto ts = get_field_in("qr").get_header().get_tracking().get_time_stamp();
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
} //namespace scream
