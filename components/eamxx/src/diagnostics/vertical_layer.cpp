#include "diagnostics/vertical_layer.hpp"

namespace scream
{

// =========================================================================================
VerticalLayerDiagnostic::
VerticalLayerDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  m_diag_name = params.get<std::string>("diag_name");
  EKAT_REQUIRE_MSG(m_diag_name == "z_int"            or m_diag_name == "z_mid"            or
                   m_diag_name == "geopotential_int" or m_diag_name == "geopotential_mid" or
                   m_diag_name == "dz",
                   "Error! VerticalLayerDiagnostic has been given an unknown name: "+m_diag_name+".\n");

  m_only_compute_dz = (m_diag_name == "dz");
  m_is_interface_layout = m_diag_name.find("_int") != std::string::npos;

  // Whether or not diagnostic is computed from sea level depends on the name.
  // "z_" -> from sea level, "geopotential_" -> from topography data.
  // This boolean is irrelevant for vertical layer thickness (dz).
  m_from_sea_level = m_diag_name.find("z_") != std::string::npos;
}
// ========================================================================================
void VerticalLayerDiagnostic::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto Q = kg/kg;
  Q.set_string("kg/kg");
  auto m2 = m*m;
  auto s2 = s*s;

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar2d_layout     { {COL     }, {m_num_cols             } };
  FieldLayout scalar3d_layout_mid { {COL,LEV }, {m_num_cols,m_num_levs  } };
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {m_num_cols,m_num_levs+1} };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid",          scalar3d_layout_mid, K,     grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,    grid_name, ps);
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa,    grid_name, ps);
  add_field<Required>("qv",             scalar3d_layout_mid, Q,     grid_name, ps);

  // Only need phis if computing geopotential_*
  if (not m_only_compute_dz and not m_from_sea_level) {
    add_field<Required>("phis", scalar2d_layout, m2/s2, grid_name);
  }

  // Construct and allocate the diagnostic field based on the diagnostic name.
  const auto diag_layout = m_is_interface_layout ? scalar3d_layout_int : scalar3d_layout_mid;
  FieldIdentifier fid (name(), diag_layout, m, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(ps);
  m_diagnostic_output.allocate_view();

  // Initialize temporary views based on need.
  if (not m_only_compute_dz) {
    if (m_is_interface_layout) {
      const auto npacks = ekat::npack<Pack>(m_num_levs);
      m_tmp_midpoint_view = view_2d("tmp_mid",m_num_cols,npacks);
    } else {
      const auto npacks_p1 = ekat::npack<Pack>(m_num_levs+1);
      m_tmp_interface_view = view_2d("tmp_int",m_num_cols,npacks_p1);
    }
  }
}
// =========================================================================================
void VerticalLayerDiagnostic::compute_diagnostic_impl()
{
  const auto npacks         = ekat::npack<Pack>(m_num_levs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_num_cols, npacks);

  const auto& T_mid              = get_field_in("T_mid").get_view<const Pack**>();
  const auto& p_mid              = get_field_in("p_mid").get_view<const Pack**>();
  const auto& qv_mid             = get_field_in("qv").get_view<const Pack**>();
  const auto& pseudo_density_mid = get_field_in("pseudo_density").get_view<const Pack**>();

  view_1d_const phis;
  if (not m_only_compute_dz and not m_from_sea_level) {
    phis = get_field_in("phis").get_view<const Real*>();
  }

  const bool only_compute_dz     = m_only_compute_dz;
  const bool is_interface_layout = m_is_interface_layout;
  const bool from_sea_level      = m_from_sea_level;
  const int  num_levs            = m_num_levs;

  // Alias correct view for diagnostic output and for tmp class views
  view_2d interface_view;
  view_2d midpoint_view;
  if (only_compute_dz) {
    midpoint_view  = m_diagnostic_output.get_view<Pack**>();
  } else if (is_interface_layout) {
    interface_view = m_diagnostic_output.get_view<Pack**>();
    midpoint_view  = m_tmp_midpoint_view;
  } else {
    midpoint_view  = m_diagnostic_output.get_view<Pack**>();
    interface_view = m_tmp_interface_view;
  }

  Kokkos::parallel_for("VerticalLayerDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();

    // Calculate dz. Use the memory in tmp_mid_view for dz and z_mid,
    // since we don't set z_mid until after dz is no longer needed.
    const auto& dz_s = ekat::subview(midpoint_view, icol);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, npacks), [&] (const Int& jpack) {
      dz_s(jpack) = PF::calculate_dz(pseudo_density_mid(icol,jpack), p_mid(icol,jpack), T_mid(icol,jpack), qv_mid(icol,jpack));
    });
    team.team_barrier();

    if (not only_compute_dz) {
      // Calculate z_int if this diagnostic is not dz
      const auto& z_int_s = ekat::subview(interface_view, icol);
      const Real surf_geopotential = from_sea_level ? 0.0 : phis(icol);
      PF::calculate_z_int(team,num_levs,dz_s,surf_geopotential,z_int_s);

      if (not is_interface_layout) {
        // Calculate z_mid if this diagnostic is not dz or an interface value
        const auto& z_mid_s = ekat::subview(midpoint_view, icol);
        PF::calculate_z_mid(team,num_levs,z_int_s,z_mid_s);
      }
    }
  });
}
// =========================================================================================
} //namespace scream
