#include "diagnostics/number_path.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

#include "physics/share/physics_constants.hpp"

namespace scream {

NumberPathDiagnostic::NumberPathDiagnostic(const ekat::Comm &comm,
                                           const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  EKAT_REQUIRE_MSG(params.isParameter("Number Kind"),
                   "Error! NumberPathDiagnostic requires 'Number Kind' in its "
                   "input parameters.\n");

  m_kind = m_params.get<std::string>("Number Kind");
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
        "Error! Invalid choice for 'NumberKind' in NumberPathDiagnostic.\n"
        "  - input value: " +
        m_kind +
        "\n"
        "  - valid values: Liq, Ice, Rain\n");
  }
}

void NumberPathDiagnostic::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  auto m2 = pow(m,2);

  auto grid             = grids_manager->get_grid("Physics");
  const auto &grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs();  // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  auto scalar2d = grid->get_2d_scalar_layout();
  auto scalar3d = grid->get_3d_scalar_layout(true);

  // The fields required for this diagnostic to be computed
  add_field<Required>("pseudo_density", scalar3d, Pa, grid_name);
  add_field<Required>(m_qname, scalar3d, kg / kg, grid_name);
  add_field<Required>(m_nname, scalar3d, 1 / kg, grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid(m_kind + name(), scalar2d, kg/(kg*m2), grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void NumberPathDiagnostic::compute_diagnostic_impl() {
  using PC  = scream::physics::Constants<Real>;
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  constexpr Real g = PC::gravit;

  const auto np  = m_diagnostic_output.get_view<Real *>();
  const auto q   = get_field_in(m_qname).get_view<const Real **>();
  const auto n   = get_field_in(m_nname).get_view<const Real **>();
  const auto rho = get_field_in("pseudo_density").get_view<const Real **>();

  const auto num_levs = m_num_levs;
  const auto policy   = ESU::get_default_team_policy(m_num_cols, m_num_levs);
  Kokkos::parallel_for(
      "Compute " + m_kind + name(), policy, KOKKOS_LAMBDA(const MT &team) {
        const int icol = team.league_rank();
        auto q_icol    = ekat::subview(q, icol);
        auto n_icol    = ekat::subview(n, icol);
        auto rho_icol  = ekat::subview(rho, icol);
        Kokkos::parallel_reduce(
            Kokkos::TeamVectorRange(team, num_levs),
            [&](const int &ilev, Real &lsum) {
              lsum += q_icol(ilev) * n_icol(ilev) * rho_icol(ilev) / g;
            },
            np(icol));
        team.team_barrier();
      });
}

}  // namespace scream
