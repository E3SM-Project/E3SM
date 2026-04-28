#include "vert_derivative.hpp"

#include "share/physics/eamxx_common_physics_functions.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream {

VertDerivative::
VertDerivative(const ekat::Comm &comm, const ekat::ParameterList &params,
               const std::shared_ptr<const AbstractGrid> &grid)
 : AbstractDiagnostic(comm, params, grid)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  m_field_name = m_params.get<std::string>("field_name");

  m_field_in_names.push_back(m_field_name);

  // we support p or z
  m_derivative_method = m_params.get<std::string>("derivative_method");
  EKAT_REQUIRE_MSG(m_derivative_method == "p" || m_derivative_method == "z",
                   "Error! VertDerivative only supports 'p' or 'z' as derivative_method.\n"
                   " - derivative_method: " +
                       m_derivative_method + "\n");

  m_field_in_names.push_back("pseudo_density");
  if (m_derivative_method == "z") {
    m_field_in_names.push_back("dz");
  }
}

void VertDerivative::initialize_impl()
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto& f      = m_fields_in.at(m_field_name);
  const auto& fid    = f.get_header().get_identifier();
  const auto& layout = fid.get_layout();

  // TODO: support higher-dimensioned input fields
  EKAT_REQUIRE_MSG(layout.rank() >= 2 && layout.rank() <= 2,
      "Error! Field rank not supported by VertDerivative.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags().back() == LEV,
      "Error! VertDerivative diagnostic expects a layout ending with the 'LEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  auto diag_units = fid.get_units();

  if (m_derivative_method == "p") {
    diag_units = fid.get_units() / Pa;
  } else if (m_derivative_method == "z") {
    diag_units = fid.get_units() / m;
  }

  auto diag_name = m_field_name + "_" + m_derivative_method + "vert_derivative";
  auto d_fid = fid.clone(diag_name).reset_units(diag_units);
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  if (f.has_valid_mask()) {
    m_diagnostic_output.create_valid_mask();
    m_diagnostic_output.get_header().set_may_be_filled(true);
  }
}

void VertDerivative::compute_diagnostic_impl()
{
  // f: input field; o: output field; d: denominator field
  const auto& f  = m_fields_in.at(m_field_name);
  const auto& o  = m_diagnostic_output;
  const auto& dp = m_fields_in.at("pseudo_density");

  // TODO: support higher-dimensioned input fields
  auto f2d = f.get_view<const Real **>();
  auto o2d = o.get_view<Real **>();

  auto dp2d = dp.get_view<const Real **>();

  using KT          = KokkosTypes<DefaultDevice>;
  using MT          = typename KT::MemberType;
  using TPF         = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  using mview_t     = Field::view_dev_t<int**>;
  using cmview_t    = Field::view_dev_t<const int**>;
  const int ncols   = f.get_header().get_identifier().get_layout().dim(0);
  const int nlevs   = f.get_header().get_identifier().get_layout().dim(1);
  const auto policy = TPF::get_default_team_policy(ncols, nlevs);

  auto d_v = (m_derivative_method == "z") ? m_fields_in.at("dz").get_view<Real **>() : dp2d;

  bool masked = m_diagnostic_output.has_valid_mask();
  auto d_mask = masked ? m_diagnostic_output.get_valid_mask().get_view<int**>() : mview_t{};
  auto f_mask = masked ? f.get_valid_mask().get_view<const int**>() : cmview_t{};
  int last_lev = nlevs-1;
  constexpr auto fv = constants::fill_value<Real>;
  auto lambda = KOKKOS_LAMBDA(const MT &team) {
    const int icol = team.league_rank();
    auto f_icol    = ekat::subview(f2d, icol); // field at midpoint
    auto o_icol    = ekat::subview(o2d, icol); // output at midnpoint
    auto d_icol    = ekat::subview(d_v, icol); // recall denominator is already a difference of interfaces
    auto dpicol    = ekat::subview(dp2d, icol); // in case of z deriv, d_icol and dpicol are not the same

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&](const int ilev) {
      // general formula for interface (assuming all fields are weighted by dp):
      // f_int(k+1) = (f_m(k+1)dp(k) + f_m(k)dp(k+1)) / (dp(k) + dp(k+1))
      // if dp(k) << dp(k+1), the interface is closer to fm(k) than fm(k+1)
      // boundary points: assume constant extrapolation (i.e., f_int(0)=f_mid(0))

      if (masked) {
        d_mask(icol,ilev) = f_mask(icol,ilev) *
                            (ilev<last_lev ? f_mask(icol,ilev+1) : 1) *
                            (ilev>0 ? f_mask(icol,ilev-1) : 1);
      }

      if (not masked or d_mask(icol,ilev)!=0) {
        auto f_int_kp1 =
            (ilev < nlevs - 1)
                ? (f_icol(ilev + 1) * dpicol(ilev) + f_icol(ilev) * dpicol(ilev + 1)) /
                      (dpicol(ilev) + dpicol(ilev + 1))
                : f_icol(nlevs - 1);

        auto f_int_kp0 =
            (ilev > 0) 
                ? (f_icol(ilev) * dpicol(ilev - 1) + f_icol(ilev - 1) * dpicol(ilev)) /
                       (dpicol(ilev - 1) + dpicol(ilev))
                : f_icol(0);

        o_icol(ilev) = (f_int_kp1 - f_int_kp0) / d_icol(ilev);
      } else {
        // Remove this when IO solely relies on mask fields
        o_icol(ilev) = fv;
      }
    });
  };
  Kokkos::parallel_for("Compute df / denominator for " + m_diagnostic_output.name(), policy, lambda);
}

} // namespace scream
