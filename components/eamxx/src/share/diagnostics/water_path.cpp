#include "water_path.hpp"
#include "share/physics/physics_constants.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream
{

WaterPath::
WaterPath (const ekat::Comm& comm, const ekat::ParameterList& params,
           const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  using namespace ekat::units;

  EKAT_REQUIRE_MSG (params.isParameter("water_kind"),
      "Error! WaterPath requires 'water_kind' in its input parameters.\n");

  m_kind = m_params.get<std::string>("water_kind");
  if (m_kind=="Liq") {
    m_qname = "qc";
  } else if (m_kind=="Ice") {
    m_qname = "qi";
  } else if (m_kind=="Rain") {
    m_qname = "qr";
  } else if (m_kind=="Rime") {
    m_qname = "qm";
  } else if (m_kind=="Vap") {
    m_qname = "qv";
  } else {
    EKAT_ERROR_MSG (
        "Error! Invalid choice for 'WaterKind' in WaterPath.\n"
        "  - input value: " + m_kind + "\n"
        "  - valid values: Liq, Ice, Rain, Rime, Vap\n");
  }

  m_field_in_names.push_back("pseudo_density");
  m_field_in_names.push_back(m_qname);

  auto m2 = pow (m,2);
  auto diag_layout = m_grid->get_2d_scalar_layout();
  FieldIdentifier fid (m_kind + name(), diag_layout, kg/m2, m_grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void WaterPath::compute_diagnostic_impl()
{
  using PC  = scream::physics::Constants<Real>;
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using TPF = ekat::TeamPolicyFactory<typename KT::ExeSpace>;

  constexpr Real g = PC::gravit.value;

  const auto wp     = m_diagnostic_output.get_view<Real*>();
  const auto q      = m_fields_in.at(m_qname).get_view<const Real**>();
  const auto rho    = m_fields_in.at("pseudo_density").get_view<const Real**>();

  const auto nlevs = m_grid->get_num_vertical_levels();
  const auto ncols = m_grid->get_num_local_dofs();
  const auto policy = TPF::get_default_team_policy(ncols, nlevs);

  auto lambda = KOKKOS_LAMBDA(const MT& team) {
    const int icol = team.league_rank();
    auto q_icol    = ekat::subview(q,icol);
    auto rho_icol  = ekat::subview(rho,icol);
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, nlevs),
                            [&] (const int& ilev, Real& lsum) {
      lsum += q_icol(ilev) * rho_icol(ilev) / g;
    },wp(icol));
    team.team_barrier();
  };
  Kokkos::parallel_for("Compute " + m_kind + name(), policy, lambda);
}

} //namespace scream
