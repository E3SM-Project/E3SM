#include "diagnostics/number_path.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

#include "physics/share/physics_constants.hpp"

namespace scream {

NumberPathDiagnostic::NumberPathDiagnostic(const ekat::Comm &comm,
                                           const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  EKAT_REQUIRE_MSG(params.isParameter("Number Kinds"),
                   "Error! NumberPathDiagnostic requires 'Number Kinds' in its "
                   "input parameters.\n");

  m_kinds = m_params.get<std::vector<std::string>>("Number Kinds");
  // check if "Liq" or "Ice" or "Rain" is in m_kinds
  EKAT_REQUIRE_MSG(
      std::find(m_kinds.begin(), m_kinds.end(), "Liq") != m_kinds.end() ||
          std::find(m_kinds.begin(), m_kinds.end(), "Ice") != m_kinds.end() ||
          std::find(m_kinds.begin(), m_kinds.end(), "Rain") != m_kinds.end(),
      "Error! NumberPathDiagnostic requires 'Number Kinds' to contain at "
      "least one of 'Liq', 'Ice', or 'Rain'.\n");
  // populate m_qnames and m_nnames
  for(const auto &kind : m_kinds) {
    if(kind == "Liq") {
      m_qnames[kind] = "qc";
      m_nnames[kind] = "nc";
    } else if(kind == "Ice") {
      m_qnames[kind] = "qi";
      m_nnames[kind] = "ni";
    } else if(kind == "Rain") {
      m_qnames[kind] = "qr";
      m_nnames[kind] = "nr";
    } else {
      EKAT_ERROR_MSG(
          "Error! Invalid choice for 'NumberKind' in NumberPathDiagnostic.\n"
          "  - input value: " +
          kind +
          "\n"
          "  - valid values: Liq, Ice, Rain\n");
    }
  }
}

void NumberPathDiagnostic::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  auto m2 = pow(m, 2);

  auto grid             = grids_manager->get_grid("Physics");
  const auto &grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs();  // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  auto scalar2d = grid->get_2d_scalar_layout();
  auto scalar3d = grid->get_3d_scalar_layout(true);

  // The fields required for this diagnostic to be computed
  add_field<Required>("pseudo_density", scalar3d, Pa, grid_name);
  for(const auto &kind : m_kinds) {
    add_field<Required>(m_qnames.at(kind), scalar3d, kg / kg, grid_name);
    add_field<Required>(m_nnames.at(kind), scalar3d, 1 / kg, grid_name);
  }

  // Construct and allocate the diagnostic field
  for(const auto &kind : m_kinds) {
    FieldIdentifier fid(kind + name(), scalar2d, kg / (kg * m2), grid_name);
    m_diagnostic_output = Field(fid);
    m_diagnostic_output.allocate_view();
    m_diagnostic_fields[kind + name()] = m_diagnostic_output;
  }
}

void NumberPathDiagnostic::compute_diagnostic_impl() {
  using PC  = scream::physics::Constants<Real>;
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  constexpr Real g = PC::gravit;
  const auto rho   = get_field_in("pseudo_density").get_view<const Real **>();

  auto compute_number_path = [this](const std::string &kind,
                                    const decltype(rho) &rho) {
    const auto np = m_diagnostic_output.get_view<Real *>();
    const auto q  = get_field_in(m_qnames[kind]).get_view<const Real **>();
    const auto n  = get_field_in(m_nnames[kind]).get_view<const Real **>();

    const auto num_levs = m_num_levs;
    const auto policy   = ESU::get_default_team_policy(m_num_cols, m_num_levs);
    Kokkos::parallel_for(
        "Compute " + kind + name(), policy, KOKKOS_LAMBDA(const MT &team) {
          const int icol = team.league_rank();
          auto q_icol    = ekat::subview(q, icol);
          auto n_icol    = ekat::subview(n, icol);
          auto rho_icol  = ekat::subview(rho, icol);
          Real lsum      = 0.0;
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(team, num_levs),
              [&](const int &ilev, Real &lsum) {
                lsum += q_icol(ilev) * n_icol(ilev) * rho_icol(ilev) / g;
              },
              lsum);
          np(icol) = lsum;
          // TODO: should we team_barrier here?
        });
  };

  // TODO: should we parallelize over the kinds?
  for(const auto &kind : m_kinds) {
    compute_number_path(kind, rho);
  }
}

}  // namespace scream
