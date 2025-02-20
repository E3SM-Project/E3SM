#include "diagnostics/vertical_layer.hpp"

#include "physics/share/physics_constants.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"
#include "share/util/eamxx_column_ops.hpp"

namespace scream
{

// =========================================================================================
VerticalLayerDiagnostic::
VerticalLayerDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  m_diag_name = params.get<std::string>("diag_name");
  std::vector<std::string> supported = {"z","geopotential","height","dz"};

  EKAT_REQUIRE_MSG(ekat::contains(supported,m_diag_name),
      "[VerticalLayerDiagnostic] Error! Invalid diag_name.\n"
      "  - diag_name  : " + m_diag_name + "\n"
      "  - valid names: " + ekat::join(supported,", ") + "\n");

  auto vert_pos = params.get<std::string>("vert_location");
  EKAT_REQUIRE_MSG (vert_pos=="mid" || vert_pos=="int" ||
                    vert_pos=="midpoints" || vert_pos=="interfaces",
      "[VerticalLayerDiagnostic] Error! Invalid 'vert_location'.\n"
      "  - input value: " + vert_pos + "\n"
      "  - valid names: mid, midpoints, int, interfaces\n");
  m_is_interface_layout = vert_pos=="int" || vert_pos=="interfaces";

  m_geopotential = m_diag_name=="geopotential";
  m_from_sea_level = m_diag_name=="z" or m_geopotential;

  if (m_diag_name!="dz") {
    m_diag_name += m_is_interface_layout ? "_int" : "_mid";
  }
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

  const auto scalar2d     = grid->get_2d_scalar_layout();
  const auto scalar3d_mid = grid->get_3d_scalar_layout(true);
  const auto scalar3d_int = grid->get_3d_scalar_layout(false);

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid",          scalar3d_mid, K,     grid_name);
  add_field<Required>("pseudo_density", scalar3d_mid, Pa,    grid_name);
  add_field<Required>("p_mid",          scalar3d_mid, Pa,    grid_name);
  add_field<Required>("qv",             scalar3d_mid, kg/kg, grid_name);

  // Only need phis if computing geopotential_*
  if (m_from_sea_level) {
    add_field<Required>("phis", scalar2d, m2/s2, grid_name);
  }
}

void VerticalLayerDiagnostic::
initialize_impl (const RunType /*run_type*/)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  auto m2 = pow(m,2);
  auto s2 = pow(s,2);

  const auto& T    = get_field_in("T_mid");
  const auto& rho  = get_field_in("pseudo_density");
  const auto& p    = get_field_in("p_mid");
  const auto& qv   = get_field_in("qv");
  const auto& phis = m_from_sea_level ? get_field_in("phis") : T; // unused if m_from_sea_level=false

  // Construct and allocate the diagnostic field.
  // Notes:
  //  - consider diag name to set long name
  //  - check input fields alloc props to set alloc props for output

  const auto& grid_name = T.get_header().get_identifier().get_grid_name();
  const auto VLEV  = m_is_interface_layout ? ILEV : LEV;
  const auto nlevs = m_is_interface_layout ? m_num_levs+1 : m_num_levs;
  FieldLayout diag_layout ({COL,VLEV},{m_num_cols,nlevs});
  FieldIdentifier fid (m_diag_name, diag_layout, m_geopotential ? m2/s2 : m, grid_name);

  m_diagnostic_output = Field(fid);
  auto& diag_fap = m_diagnostic_output.get_header().get_alloc_properties();
  int ps = SCREAM_PACK_SIZE;
  for (const auto& f : {T,rho,p,qv,phis}) {
    const auto& fap = f.get_header().get_alloc_properties();
    const auto& f_ps = fap.get_largest_pack_size();

    // We must use a pack size that works with all inputs, so pick the smallest
    ps = std::min(ps,f_ps);
  }
  diag_fap.request_allocation(ps);
  m_diagnostic_output.allocate_view();

  using stratts_t = std::map<std::string,std::string>;
  auto& io_atts = m_diagnostic_output.get_header().get_extra_data<stratts_t>("io: string attributes");
  auto& long_name = io_atts["long_name"];
  if (m_diag_name=="dz") {
    long_name = "level thickness";
  } else if (m_diag_name=="z_mid") {
    long_name = "elevation above sealevel at level midpoints";
  } else if (m_diag_name=="z_int") {
    long_name = "elevation above sealevel at level interfaces";
  } else if (m_diag_name=="height_mid") {
    long_name= "elevation above surface at level midpoints";
  } else if (m_diag_name=="height_int") {
    long_name = "elevation above surface at level interfaces";
  } else if (m_diag_name=="geopotential_mid") {
    long_name = "geopotential height relative to sealevel at level midpoints";
  } else {
    long_name = "geopotential height relative to sealevel at level interfaces";
  }

  // Initialize temporary views based on need. Can alias the diag if a temp is not needed
  auto create_temp = [&](const std::string& name, int levs) {
    auto u = Units::nondimensional();
    FieldLayout fl({COL,LEV},{m_num_cols,levs});
    FieldIdentifier fid (name,fl,u,grid_name);
    Field f = Field(fid);
    f.get_header().get_alloc_properties().request_allocation(ps);
    f.allocate_view();
    return f;
  };
  if (m_diag_name == "dz") {
    m_tmp_midpoint = m_diagnostic_output;
    m_tmp_interface = m_diagnostic_output; // Not really used
  } else if (m_is_interface_layout) {
    m_tmp_midpoint = create_temp("tmp_mid",m_num_levs);
    m_tmp_interface = m_diagnostic_output;
  } else {
    m_tmp_interface = create_temp("tmp_int",m_num_levs+1);
    m_tmp_midpoint = m_diagnostic_output;
  }
}
// =========================================================================================
void VerticalLayerDiagnostic::compute_diagnostic_impl()
{
  const auto& fap = m_diagnostic_output.get_header().get_alloc_properties();
  if (fap.get_largest_pack_size()==SCREAM_PACK_SIZE) {
    do_compute_diagnostic_impl<SCREAM_PACK_SIZE>();
  } else {
    do_compute_diagnostic_impl<1>();
  }
}

template<int PackSize>
void VerticalLayerDiagnostic::do_compute_diagnostic_impl()
{
  using column_ops  = ColumnOps<DefaultDevice,Real>;
  using PackT = ekat::Pack<Real,PackSize>;
  using KT    = KokkosTypes<DefaultDevice>;
  using MemberType = typename KT::MemberType;
  using PF    = PhysicsFunctions<DefaultDevice>;

  // To use in column_ops, since we integrate from surface
  constexpr bool FromTop = false;

  const auto npacks = ekat::npack<PackT>(m_num_levs);
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_num_cols, npacks);

  const auto& T    = get_field_in("T_mid").get_view<const PackT**>();
  const auto& p    = get_field_in("p_mid").get_view<const PackT**>();
  const auto& qv   = get_field_in("qv").get_view<const PackT**>();
  const auto& rho  = get_field_in("pseudo_density").get_view<const PackT**>();
  const auto  phis = m_from_sea_level ? get_field_in("phis").get_view<const Real*>() : typename KT::view_1d<const Real>();

  const bool only_compute_dz     = m_diag_name=="dz";
  const bool is_interface_layout = m_is_interface_layout;
  const bool from_sea_level      = m_from_sea_level;
  const bool geopotential        = m_geopotential;
  const int  num_levs            = m_num_levs;
  constexpr auto g = scream::physics::Constants<Real>::gravit;

  // Alias correct view for diagnostic output and for tmp class views
  auto tmp_mid = m_tmp_midpoint.get_view<PackT**>();
  auto tmp_int = m_tmp_interface.get_view<PackT**>();

  // Define the lambda, then dispatch the ||for
  auto lambda = KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();

    // Whatever the output needs, the first thing to compute is dz.
    const auto& dz = ekat::subview(tmp_mid, icol);
    PF::calculate_dz(team,ekat::subview(rho,icol),
                          ekat::subview(p,icol),
                          ekat::subview(T,icol),
                          ekat::subview(qv,icol),
                          dz);
    team.team_barrier();

    // If dz is all we need, we're done
    if (only_compute_dz) { return; }

    // Now integrate to compute quantity at interfaces
    const auto& v_int = ekat::subview(tmp_int, icol);

    // phi and z are related by phi=z*g, so dphi=dz*g, and z_surf = phis/g
    if (geopotential) {
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
      const auto& v_mid = ekat::subview(tmp_mid, icol);
      column_ops::compute_midpoint_values(team,num_levs,v_int,v_mid);
    }
  };
  Kokkos::parallel_for(m_diag_name, policy, lambda);
}

} //namespace scream
