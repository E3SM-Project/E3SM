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
                   m_diag_name == "altitude_int"     or m_diag_name == "altitude_mid"     or
                   m_diag_name == "dz",
                   "Error! VerticalLayerDiagnostic has been given an unknown name: "+m_diag_name+".\n");

  m_is_interface_layout = m_diag_name.find("_int") != std::string::npos;

  m_geopotential = m_diag_name.substr(0,12)=="geopotential";
  m_from_sea_level = m_diag_name[0]=='z' or m_geopotential;
}
// ========================================================================================
void VerticalLayerDiagnostic::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto m2 = pow(m,2);
  auto s2 = pow(s,2);

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
  add_field<Required>("qv",             scalar3d_layout_mid, kg/kg, grid_name, ps);

  // Only need phis if computing geopotential_*
  if (not m_geopotential) {
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
  if (m_diag_name!="dz") {
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
  using column_ops  = ColumnOps<DefaultDevice,Real>;
  // To use in column_ops, since we integrate from surface
  constexpr bool FromTop = false;

  const auto npacks         = ekat::npack<Pack>(m_num_levs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_num_cols, npacks);

  const auto& T_mid              = get_field_in("T_mid").get_view<const Pack**>();
  const auto& p_mid              = get_field_in("p_mid").get_view<const Pack**>();
  const auto& qv_mid             = get_field_in("qv").get_view<const Pack**>();
  const auto& pseudo_density_mid = get_field_in("pseudo_density").get_view<const Pack**>();

  view_1d_const phis;
  if (m_from_sea_level) {
    phis = get_field_in("phis").get_view<const Real*>();
  }

  const bool only_compute_dz     = m_diag_name=="dz";
  const bool is_interface_layout = m_is_interface_layout;
  const bool from_sea_level      = m_from_sea_level;
  const bool do_geopotential     = m_geopotential;
  const int  num_levs            = m_num_levs;
  constexpr auto g = scream::physics::Constants<Real>::gravit;

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

    // Whatever the output needs, the first thing to compute is dz.
    const auto& dz = ekat::subview(midpoint_view, icol);
    PF::calculate_dz(team,ekat::subview(pseudo_density_mid,icol),
                          ekat::subview(p_mid,icol),
                          ekat::subview(T_mid,icol),
                          ekat::subview(qv_mid,icol),
                          dz);
    team.team_barrier();

    // If dz is all we need, we're done
    if (only_compute_dz) { return; }

    // Now integrate to compute quantity at interfaces
    const auto& v_int = ekat::subview(interface_view, icol);

    // phi and z are related by phi=z*g, so dphi=dz*g, and z_surf = phis/g
    if (do_geopotential) {
      auto dphi = [&](const int ilev) {
        return dz(ilev) * g;
      };
      column_ops::template column_scan<FromTop>(team,num_levs,dphi,v_int,phis(icol));
    } else {
      const Real surf_val = from_sea_level ? phis(icol)/g : 0;
      column_ops::template column_scan<FromTop>(team,num_levs,dz,v_int,surf_val);
    }

    // If we need quantity at midpoints, simply do int->mid averaging
    if (not is_interface_layout) {
      team.team_barrier();
      const auto& v_mid = ekat::subview(midpoint_view, icol);
      column_ops::compute_midpoint_values(team,num_levs,v_int,v_mid);
    }
  });
}

} //namespace scream
