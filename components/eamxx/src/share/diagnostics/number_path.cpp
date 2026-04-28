#include "number_path.hpp"

#include <ekat_team_policy_utils.hpp>

#include "share/physics/physics_constants.hpp"

namespace scream {

NumberPath::
NumberPath(const ekat::Comm &comm,
           const ekat::ParameterList &params,
           const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  using namespace ekat::units;

  EKAT_REQUIRE_MSG(params.isParameter("number_kind"),
      "Error! NumberPath requires 'number_kind' in its input parameters.\n");

  m_kind = m_params.get<std::string>("number_kind");
  if(m_kind == "Liq") {
    m_qname = "qc";
    m_nname = "nc";
  } else if(m_kind == "Ice") {
    m_qname = "qi";
    m_nname = "ni";
  } else if(m_kind == "Rain") {
    m_qname = "qr";
    m_nname = "nr";
  } else {
    EKAT_ERROR_MSG(
        "Error! Invalid choice for 'NumberKind' in NumberPath.\n"
        "  - input value: " +
        m_kind +
        "\n"
        "  - valid values: Liq, Ice, Rain\n");
  }

  // The fields required for this diagnostic to be computed
  m_field_in_names.push_back("pseudo_density");
  m_field_in_names.push_back(m_qname);
  m_field_in_names.push_back(m_nname);

  // Construct and allocate the diagnostic field
  auto m2 = pow(m,2);
  auto diag_layout = m_grid->get_2d_scalar_layout();
  FieldIdentifier fid(m_kind + name(), diag_layout, kg/(kg*m2), m_grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void NumberPath::compute_diagnostic_impl()
{
  using PC  = scream::physics::Constants<Real>;
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using TPF = ekat::TeamPolicyFactory<typename KT::ExeSpace>;

  constexpr Real g = PC::gravit.value;

  const auto np  = m_diagnostic_output.get_view<Real *>();
  const auto q   = m_fields_in.at(m_qname).get_view<const Real **>();
  const auto n   = m_fields_in.at(m_nname).get_view<const Real **>();
  const auto rho = m_fields_in.at("pseudo_density").get_view<const Real **>();

  const auto nlevs  = m_grid->get_num_vertical_levels();
  const auto ncols  = m_grid->get_num_local_dofs();
  const auto policy = TPF::get_default_team_policy(ncols, nlevs);

  auto lambda = KOKKOS_LAMBDA(const MT &team) {
    const int icol = team.league_rank();
    auto q_icol    = ekat::subview(q, icol);
    auto n_icol    = ekat::subview(n, icol);
    auto rho_icol  = ekat::subview(rho, icol);
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, nlevs),
        [&](const int &ilev, Real &lsum) {
          lsum += q_icol(ilev) * n_icol(ilev) * rho_icol(ilev) / g;
        },
        np(icol));
  };
  Kokkos::parallel_for("Compute " + m_kind + name(), policy, lambda);
}

}  // namespace scream
