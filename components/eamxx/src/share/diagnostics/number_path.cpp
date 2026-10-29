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

  // The fields required for this diagnostic to be computed
  m_field_in_names = {"pseudo_density","qc","nc","qi","ni","qr","nr"};

  // Construct and allocate the diagnostic fields
  auto m2 = pow(m,2);
  auto diag_layout = m_grid->get_2d_scalar_layout();
  auto create = [&](const std::string& kind) {
    FieldIdentifier fid(kind + name(), diag_layout, kg/(kg*m2), m_grid->name());
    Field f(fid);
    f.allocate_view();
    add_output(f);
    return f;
  };
  m_liq  = create("Liq");
  m_ice  = create("Ice");
  m_rain = create("Rain");
}

void NumberPath::compute_impl()
{
  using PC  = scream::physics::Constants<Real>;
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using TPF = ekat::TeamPolicyFactory<typename KT::ExeSpace>;

  constexpr Real g = PC::gravit.value;

  const auto lnp = m_liq.get_view<Real *>();
  const auto inp = m_ice.get_view<Real *>();
  const auto rnp = m_rain.get_view<Real *>();
  const auto qc  = m_fields_in.at("qc").get_view<const Real **>();
  const auto nc  = m_fields_in.at("nc").get_view<const Real **>();
  const auto qi  = m_fields_in.at("qi").get_view<const Real **>();
  const auto ni  = m_fields_in.at("ni").get_view<const Real **>();
  const auto qr  = m_fields_in.at("qr").get_view<const Real **>();
  const auto nr  = m_fields_in.at("nr").get_view<const Real **>();
  const auto rho = m_fields_in.at("pseudo_density").get_view<const Real **>();

  const auto nlevs  = m_grid->get_num_vertical_levels();
  const auto ncols  = m_grid->get_num_local_dofs();
  const auto policy = TPF::get_default_team_policy(ncols, nlevs);

  // One reduction per kind, as in separate kernels, so each sum is
  // accumulated in the same order as when each kind was computed alone
  auto lambda = KOKKOS_LAMBDA(const MT &team) {
    const int icol = team.league_rank();
    auto rho_icol  = ekat::subview(rho, icol);
    auto path = [&](const auto& q, const auto& n, Real& np) {
      auto q_icol = ekat::subview(q, icol);
      auto n_icol = ekat::subview(n, icol);
      Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team, nlevs),
          [&](const int &ilev, Real &lsum) {
            lsum += q_icol(ilev) * n_icol(ilev) * rho_icol(ilev) / g;
          },
          np);
    };
    path(qc, nc, lnp(icol));
    path(qi, ni, inp(icol));
    path(qr, nr, rnp(icol));
  };
  Kokkos::parallel_for("Compute " + name(), policy, lambda);
}

}  // namespace scream
