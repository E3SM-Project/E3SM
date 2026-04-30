#include "vapor_flux.hpp"
#include "share/physics/physics_constants.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream
{

VaporFlux::
VaporFlux (const ekat::Comm& comm, const ekat::ParameterList& params,
                      const std::shared_ptr<const AbstractGrid>& grid)
  : AbstractDiagnostic(comm,params,grid)
{
  EKAT_REQUIRE_MSG (params.isParameter("wind_component"),
      "Error! VaporFlux requires 'wind_component' in its input parameters.\n");

  const auto& comp = m_params.get<std::string>("wind_component");
  if (comp=="Zonal") {
    m_component = 0;
  } else if (comp=="Meridional") {
    m_component = 1;
  } else {
    EKAT_ERROR_MSG (
        "Error! Invalid choice for 'wind_component' in VaporFlux.\n"
        "  - input value: " + comp + "\n"
        "  - valid values: Zonal, Meridional\n");
  }
  m_name = comp + "VapFlux";

  m_field_in_names.push_back("pseudo_density");
  m_field_in_names.push_back("qv");
  m_field_in_names.push_back("horiz_winds");

  using namespace ekat::units;

  auto diag_layout = m_grid->get_2d_scalar_layout();
  FieldIdentifier fid (m_name, m_grid->get_2d_scalar_layout(), kg/m/s, m_grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void VaporFlux::compute_impl()
{
  using PC  = scream::physics::Constants<Real>;
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using TPF = ekat::TeamPolicyFactory<typename KT::ExeSpace>;

  constexpr Real g = PC::gravit.value;

  const auto diag  = m_diagnostic_output.get_view<Real*>();
  const auto qv    = m_fields_in.at("qv").get_view<const Real**>();
  const auto rho   = m_fields_in.at("pseudo_density").get_view<const Real**>();
  const auto wind  = m_fields_in.at("horiz_winds").get_component(m_component).get_view<const Real**>();

  const auto nlevs = m_grid->get_num_vertical_levels();
  const auto ncols = m_grid->get_num_local_dofs();
  const auto policy = TPF::get_default_team_policy(ncols, nlevs);
  Kokkos::parallel_for("Compute " + m_name, policy,
                       KOKKOS_LAMBDA(const MT& team) {
    const int icol = team.league_rank();

    auto qv_icol   = ekat::subview(qv,icol);
    auto rho_icol  = ekat::subview(rho,icol);
    auto wind_icol = ekat::subview(wind,icol);

    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, nlevs),
                            [&] (const int& ilev, Real& lsum) {
      lsum += wind_icol(ilev) * qv_icol(ilev) * rho_icol(ilev) / g;
    },diag(icol));
    team.team_barrier();
  });
}

} //namespace scream
