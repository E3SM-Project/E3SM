#include "diagnostics/vert_derivative.hpp"

#include "share/util/eamxx_common_physics_functions.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream {

VertDerivativeDiag::VertDerivativeDiag(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {}

void VertDerivativeDiag::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto &fn = m_params.get<std::string>("field_name");
  const auto &gn = m_params.get<std::string>("grid_name");
  const auto g   = grids_manager->get_grid(gn);

  add_field<Required>(fn, gn);

  // we support p or z
  m_derivative_method = m_params.get<std::string>("derivative_method");
  EKAT_REQUIRE_MSG(m_derivative_method == "p" || m_derivative_method == "z",
                   "Error! VertDerivativeDiag only supports 'p' or 'z' as derivative_method.\n"
                   " - derivative_method: " +
                       m_derivative_method + "\n");
  m_diag_name = fn + "_" + m_derivative_method + "vert_derivative";

  auto scalar3d = g->get_3d_scalar_layout(true);
  if (m_derivative_method == "p") {
    add_field<Required>("pseudo_density", scalar3d, Pa, gn);
  } else if (m_derivative_method == "z") {
    add_field<Required>("pseudo_density", scalar3d, Pa, gn);
    add_field<Required>("qv", scalar3d, kg / kg, gn);
    add_field<Required>("p_mid", scalar3d, Pa, gn);
    add_field<Required>("T_mid", scalar3d, K, gn);
  }
}

void VertDerivativeDiag::initialize_impl(const RunType /*run_type*/) {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto &f      = get_fields_in().front();
  const auto &fid    = f.get_header().get_identifier();
  const auto &layout = fid.get_layout();

  // TODO: support higher-dimensioned input fields
  EKAT_REQUIRE_MSG(layout.rank() >= 2 && layout.rank() <= 2,
                   "Error! Field rank not supported by VertDerivativeDiag.\n"
                   " - field name: " +
                       fid.name() +
                       "\n"
                       " - field layout: " +
                       layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags().back() == LEV,
                   "Error! VertDerivativeDiag diagnostic expects a layout ending "
                   "with the 'LEV' tag.\n"
                   " - field name  : " +
                       fid.name() +
                       "\n"
                       " - field layout: " +
                       layout.to_string() + "\n");

  ekat::units::Units diag_units = fid.get_units();

  m_denominator = get_field_in("pseudo_density").clone("denominator");

  if (m_derivative_method == "p") {
    diag_units = fid.get_units() / Pa;
  } else if (m_derivative_method == "z") {
    diag_units = fid.get_units() / m;
  }

  FieldIdentifier d_fid(m_diag_name, layout, diag_units, fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();
}

void VertDerivativeDiag::compute_diagnostic_impl() {
  // f: input field; o: output field; d: denominator field
  const auto &f = get_fields_in().front();
  const auto &o = m_diagnostic_output;
  // keep dp around
  const auto &dp = get_field_in("pseudo_density");

  // TODO: support higher-dimensioned input fields
  auto f2d = f.get_view<const Real **>();
  auto o2d = o.get_view<Real **>();

  auto dp2d = dp.get_view<const Real **>();

  using KT          = KokkosTypes<DefaultDevice>;
  using MT          = typename KT::MemberType;
  using TPF         = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  const int ncols   = m_denominator.get_header().get_identifier().get_layout().dim(0);
  const int nlevs   = m_denominator.get_header().get_identifier().get_layout().dim(1);
  const auto policy = TPF::get_default_team_policy(ncols, nlevs);

  // get the denominator first
  if (m_derivative_method == "p") {
    m_denominator.update(dp, sp(1.0), sp(0.0));
  } else if (m_derivative_method == "dz") {
    // TODO: for some reason the z_mid field keeps getting set to 0
    // TODO: as a workaround, just calculate z_mid here (sigh...)
    // m_denominator.update(get_field_in("z_mid"), 1.0, 0.0);
    using PF  = scream::PhysicsFunctions<DefaultDevice>;
    auto zm_v = m_denominator.get_view<Real **>();
    auto pm_v = get_field_in("p_mid").get_view<const Real **>();
    auto tm_v = get_field_in("T_mid").get_view<const Real **>();
    auto qv_v = get_field_in("qv").get_view<const Real **>();

    Kokkos::parallel_for(
        "Compute dz for " + m_diagnostic_output.name(), policy, KOKKOS_LAMBDA(const MT &team) {
          const int icol = team.league_rank();
          auto zm_icol   = ekat::subview(zm_v, icol);
          auto dp_icol   = ekat::subview(dp2d, icol);
          auto pm_icol   = ekat::subview(pm_v, icol);
          auto tm_icol   = ekat::subview(tm_v, icol);
          auto qv_icol   = ekat::subview(qv_v, icol);
          PF::calculate_dz(team, dp_icol, pm_icol, tm_icol, qv_icol, zm_icol);
        });
  }

  auto d_v = m_denominator.get_view<Real **>();
  Kokkos::parallel_for(
      "Compute df / denominator for " + m_diagnostic_output.name(), policy,
      KOKKOS_LAMBDA(const MT &team) {
        const int icol = team.league_rank();
        auto f_icol    = ekat::subview(f2d, icol); // field at midpoint
        auto o_icol    = ekat::subview(o2d, icol); // output at midnpoint
        auto d_icol =
            ekat::subview(d_v, icol); // recall denominator is already a difference of interfaces
        auto dpicol =
            ekat::subview(dp2d, icol); // in case of z deriv, d_icol and dpicol are not the same

        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&](const int ilev) {
          // boundary points
          if (ilev == 0) {
            o_icol(0) = (f_icol(1) - f_icol(0)) / d_icol(0);
          } else if (ilev == nlevs - 1) {
            o_icol(nlevs - 1) = (f_icol(nlevs - 1) - f_icol(nlevs - 2)) / d_icol(nlevs - 1);
          }
          // interior points
          else {
            // general formula for interface (assuming all fields are weighted by dp):
            // f_int(k+1) = (f_m(k+1)dp(k) + f_m(k)dp(k+1)) / (dp(k) + dp(k+1))
            // if dp(k) << dp(k+1), the interface is closer to fm(k) than fm(k+1)
            auto f_int_kp1 = (f_icol(ilev + 1) * dpicol(ilev) + f_icol(ilev) * dpicol(ilev + 1)) /
                             (dpicol(ilev) + dpicol(ilev + 1));
            auto f_int_kp0 = (f_icol(ilev) * dpicol(ilev - 1) + f_icol(ilev - 1) * dpicol(ilev)) /
                             (dpicol(ilev - 1) + dpicol(ilev));
            o_icol(ilev) = (f_int_kp1 - f_int_kp0) / d_icol(ilev);
          }
        });
      });
}

} // namespace scream
