#include "diagnostics/unary_ops.hpp"

#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_reduction_utils.hpp>

namespace scream {

// parse string to get the operator code
int get_unary_operator_code(const std::string& op) {
  if (op == "log") return 0;
  if (op == "exp") return 1;
  if (op == "sqrt") return 2;
  if (op == "abs") return 3;
  if (op == "square") return 4;
  if (op == "inverse") return 5;
  return -1; // invalid operator
}

UnaryOpsDiag::UnaryOpsDiag(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  m_fn = m_params.get<std::string>("field_name");
  m_op = m_params.get<std::string>("unary_op");
}

void UnaryOpsDiag::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  const auto &gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_fn, gname);
}

void UnaryOpsDiag::initialize_impl(const RunType /*run_type*/) {
  const auto &f   = get_field_in(m_fn);
  const auto &fid = f.get_header().get_identifier();
  const auto &gn  = fid.get_grid_name();
  const auto &layout = fid.get_layout();

  using namespace ekat::units;
  ekat::units::Units diag_units = fid.get_units();

  // log|exp|sqrt|abs|square|inverse
  if (m_op == "log") {
    diag_units = ekat::units::Units::nondimensional();
  } else if (m_op == "exp") {
    diag_units = ekat::units::Units::nondimensional();
  } else if (m_op == "sqrt") {
    diag_units = pow(fid.get_units(), 1/2);
  } else if (m_op == "abs") {
    diag_units = fid.get_units();
  } else if (m_op == "square") {
    diag_units = pow(fid.get_units(), 2);
  } else if (m_op == "inverse") {
    diag_units = pow(fid.get_units(), -1);
  } else {
    EKAT_ERROR_MSG("Error! Unsupported unary operation: " + m_op + "\n");
  }

  FieldIdentifier fid_out(m_op + "_" + m_fn, layout.clone(), diag_units, gn);
  m_diagnostic_output = Field(fid_out);
  m_diagnostic_output.allocate_view();

  // for log and sqrt, we need a mask
  if (m_op == "sqrt" || m_op == "log") {
    auto mask = m_diagnostic_output.clone("mask_" + m_op + "_" + m_fn);
    m_diagnostic_output.get_header().set_extra_data("mask_field", mask);
  }
}

void UnaryOpsDiag::compute_diagnostic_impl() {
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using TPF = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  using RU  = ekat::ReductionUtils<typename KT::ExeSpace>;

  constexpr auto fill_val = constants::fill_value<Real>;

  const auto f = get_field_in(m_fn).get_view<const Real **>();
  const auto d = m_diagnostic_output.get_view<Real **>();
  const auto m = (m_op == "sqrt" || m_op == "log") ? m_diagnostic_output.get_header().get_extra_data<Field>("mask_field").get_view<Real **>() : d;
  const auto& layout = get_field_in(m_fn).get_header().get_identifier().get_layout();
  const auto num_cols = layout.dims()[0];
  const auto num_levs = layout.dims()[1];
  const auto policy   = TPF::get_default_team_policy(num_cols, num_levs);

  auto op_code = get_unary_operator_code(m_op);
  Kokkos::parallel_for("UnaryOpsDiag", policy, KOKKOS_LAMBDA(const MT& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levs), [&] (const int ilev) {
      switch (op_code) {
        case 0: // log
          d(icol, ilev) = (f(icol, ilev) > 0) ? Kokkos::log(f(icol, ilev)) : fill_val;
          m(icol, ilev) = (f(icol, ilev) > 0) ? 1 : 0;
          break;
        case 1: // exp
          d(icol, ilev) = Kokkos::exp(f(icol, ilev));
          break;
        case 2: // sqrt
          d(icol, ilev) = (f(icol, ilev) > 0) ? Kokkos::sqrt(f(icol, ilev)) : fill_val;
          m(icol, ilev) = (f(icol, ilev) > 0) ? 1 : 0;
          break;
        case 3: // abs
          d(icol, ilev) = Kokkos::abs(f(icol, ilev));
          break;
        case 4: // square
          d(icol, ilev) = Kokkos::pow(f(icol, ilev), 2.0);
          break;
        case 5: // inverse
          d(icol, ilev) = Kokkos::pow(f(icol, ilev), -1.0);
          break;
        default:
          // Should not reach here if input validation is correct
          d(icol, ilev) = f(icol, ilev);
          break;
      }
    });
  });
}

}  // namespace scream
