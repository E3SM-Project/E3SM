#include "diagnostics/below_or_above_interface.hpp"

#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_reduction_utils.hpp>

namespace scream {

BelowOrAboveInterface::BelowOrAboveInterface(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  m_fn = m_params.get<std::string>("field_name");
  m_op = m_params.get<std::string>("unary_op");
}

void BelowOrAboveInterface::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  const auto &gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_name, gname);
}

void AODVis::initialize_impl(const RunType /*run_type*/) {
  const auto &f   = get_field_in(m_name);
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
    diag_units = std::pow(fid.get_units(), 0.5);
  } else if (m_op == "abs") {
    diag_units = fid.get_units();
  } else if (m_op == "square") {
    diag_units = std::pow(fid.get_units(), 2.0);
  } else if (m_op == "inverse") {
    diag_units = std::pow(fid.get_units(), -1.0);
  } else {
    EKAT_ERROR_MSG("Error! Unsupported unary operation: " + m_op + "\n");
  }

  FieldIdentifier fid(m_op + "_" + m_name, layout.clone(), diag_units, gn);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void AODVis::compute_diagnostic_impl() {
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using TPF = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  using RU  = ekat::ReductionUtils<typename KT::ExeSpace>;

  constexpr auto fill_val = constants::fill_value<Real>;

  const auto f = get_field_in(m_name).get_view<const Real **>();
  const auto d = m_diagnostic_output.get_view<Real **>();
  const auto num_cols = get_field_in(m_name).get_header().get_num_local_dofs();
  const auto num_levs = get_field_in(m_name).get_header().get_num_vertical_levels();
  const auto policy   = TPF::get_default_team_policy(num_cols, num_levs);

  auto d_op = m_op;
  Kokkos::parallel_for("BelowOrAboveInterface", policy, KOKKOS_LAMBDA(const MT& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levels), [&] (const int ilev) {
      if (d_op == "log") {
        d(icol, ilev) = (f(icol, ilev) > 0) ? Kokkos::log(f(icol, ilev)) : fill_val;
      } else if (d_op == "exp") {
        d(icol, ilev) = Kokkos::exp(f(icol, ilev));
      } else if (d_op == "sqrt") {
        d(icol, ilev) = (f(icol, ilev) > 0) ? Kokkos::sqrt(f(icol, ilev)) : fill_val;
      } else if (d_op == "abs") {
        d(icol, ilev) = Kokkos::abs(f(icol, ilev));
      } else if (d_op == "square") {
        d(icol, ilev) = Kokkos::pow(f(icol, ilev), 2.0);
      } else if (d_op == "inverse") {
        d(icol, ilev) = Kokkos::pow(f(icol, ilev), -1.0);
      }
    });
  });
}

}  // namespace scream
