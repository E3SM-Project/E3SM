#include "vertical_layer.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/util/eamxx_column_ops.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream
{

VerticalLayer::
VerticalLayer (const ekat::Comm& comm, const ekat::ParameterList& params,
               const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  m_diag_name = params.get<std::string>("diag_name");
  std::vector<std::string> supported = {"z","geopotential","height","dz"};

  EKAT_REQUIRE_MSG(ekat::contains(supported,m_diag_name),
      "[VerticalLayer] Error! Invalid diag_name.\n"
      "  - diag_name  : " + m_diag_name + "\n"
      "  - valid names: " + ekat::join(supported,", ") + "\n");

  auto vert_pos = params.get<std::string>("vert_location");
  EKAT_REQUIRE_MSG (vert_pos=="mid" || vert_pos=="int" ||
                    vert_pos=="midpoints" || vert_pos=="interfaces",
      "[VerticalLayer] Error! Invalid 'vert_location'.\n"
      "  - input value: " + vert_pos + "\n"
      "  - valid names: mid, midpoints, int, interfaces\n");
  m_is_interface_layout = vert_pos=="int" || vert_pos=="interfaces";

  m_geopotential = m_diag_name=="geopotential";
  m_from_sea_level = m_diag_name=="z" or m_geopotential;

  if (m_diag_name!="dz") {
    m_diag_name += m_is_interface_layout ? "_int" : "_mid";
  }

  m_field_in_names.push_back("T_mid");
  m_field_in_names.push_back("pseudo_density");
  m_field_in_names.push_back("p_mid");
  m_field_in_names.push_back("qv");

  if (m_from_sea_level) {
    m_field_in_names.push_back("phis");
  }
}

void VerticalLayer::initialize_impl ()
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // Construct and allocate the diagnostic field.
  auto m2 = pow(m,2);
  auto s2 = pow(s,2);
  const auto vtag  = m_is_interface_layout ? ILEV : LEV;
  auto diag_layout = m_grid->get_3d_scalar_layout(vtag);
  FieldIdentifier fid (m_diag_name, diag_layout, m_geopotential ? m2/s2 : m, m_grid->name());
  m_diagnostic_output = Field(fid);
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
  if (m_diag_name == "dz") {
    m_tmp_midpoint = m_diagnostic_output;
    m_tmp_interface = m_diagnostic_output; // Not really used
  } else if (m_is_interface_layout) {
    auto fl = m_grid->get_3d_scalar_layout(LEV);
    FieldIdentifier fid ("tmp_mid",fl,none,m_grid->name());
    m_tmp_midpoint = Field(fid,true);
    m_tmp_interface = m_diagnostic_output;
  } else {
    auto fl = m_grid->get_3d_scalar_layout(ILEV);
    FieldIdentifier fid ("tmp_int",fl,none,m_grid->name());
    m_tmp_interface = Field(fid,true);
    m_tmp_midpoint = m_diagnostic_output;
  }
}

void VerticalLayer::compute_diagnostic_impl()
{
  using column_ops = ColumnOps<DefaultDevice,Real>;
  using KT         = KokkosTypes<DefaultDevice>;
  using MemberType = typename KT::MemberType;
  using PF         = PhysicsFunctions<DefaultDevice>;
  using TPF        = ekat::TeamPolicyFactory<KT::ExeSpace>;

  // To use in column_ops, since we integrate from surface
  constexpr bool FromTop = false;

  const int ncols = m_grid->get_num_local_dofs();
  const int nlevs = m_grid->get_num_vertical_levels();
  const bool only_compute_dz     = m_diag_name=="dz";
  const bool is_interface_layout = m_is_interface_layout;
  const bool from_sea_level      = m_from_sea_level;
  const bool geopotential        = m_geopotential;

  const auto& T    = m_fields_in.at("T_mid").get_view<const Real**>();
  const auto& p    = m_fields_in.at("p_mid").get_view<const Real**>();
  const auto& qv   = m_fields_in.at("qv").get_view<const Real**>();
  const auto& rho  = m_fields_in.at("pseudo_density").get_view<const Real**>();
  const auto  phis = m_from_sea_level ? m_fields_in.at("phis").get_view<const Real*>() : typename KT::view_1d<const Real>();

  constexpr Real g = scream::physics::Constants<Real>::gravit.value;

  // Alias correct view for diagnostic output and for tmp class views
  auto tmp_mid = m_tmp_midpoint.get_view<Real**>();
  auto tmp_int = m_tmp_interface.get_view<Real**>();

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
    if (only_compute_dz)
      return;

    // Now integrate to compute quantity at interfaces
    const auto& v_int = ekat::subview(tmp_int, icol);

    // phi and z are related by phi=z*g, so dphi=dz*g, and z_surf = phis/g
    if (geopotential) {
      auto dphi = [&](const int ilev) {
        return dz(ilev) * g;
      };
      column_ops::template column_scan<FromTop>(team,nlevs,dphi,v_int,phis(icol));
    } else {
      const Real surf_val = from_sea_level ? phis(icol)/g : 0;
      column_ops::template column_scan<FromTop>(team,nlevs,dz,v_int,surf_val);
    }

    // If we need quantity at midpoints, simply do int->mid averaging
    if (not is_interface_layout) {
      team.team_barrier();
      const auto& v_mid = ekat::subview(tmp_mid, icol);
      column_ops::compute_midpoint_values(team,nlevs,v_int,v_mid);
    }
  };
  const auto policy = TPF::get_thread_range_parallel_scan_team_policy(ncols, nlevs);
  Kokkos::parallel_for(m_diag_name, policy, lambda);
}

} //namespace scream
