#include "diagnostics/vert_divergence.hpp"

#include "share/util/eamxx_common_physics_functions.hpp"

namespace scream {

VertDivergenceDiag::VertDivergenceDiag(const ekat::Comm &comm,
                                   const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {}

void VertDivergenceDiag::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto &fn = m_params.get<std::string>("field_name");
  const auto &gn = m_params.get<std::string>("grid_name");
  const auto g   = grids_manager->get_grid(gn);

  add_field<Required>(fn, gn);

  // we support dp or dz
  m_divergence_method = m_params.get<std::string>("divergence_method");
  EKAT_REQUIRE_MSG(
      m_divergence_method == "dp" || m_divergence_method == "dz",
      "Error! VertDivergenceDiag only supports 'dp' or 'dz' as divergence_method.\n"
      " - divergence_method: " + m_divergence_method + "\n");
  m_diag_name = "d" + fn + m_divergence_method + "_vert_divergence";

  auto scalar3d = g->get_3d_scalar_layout(true);
  if (m_divergence_method == "dp") {
    add_field<Required>("pseudo_density", scalar3d, Pa, gn); 
  }
  else if (m_divergence_method == "dz")
  {
    add_field<Required>("pseudo_density", scalar3d, Pa, gn);
    add_field<Required>("qv", scalar3d, kg / kg, gn);
    add_field<Required>("p_mid", scalar3d, Pa, gn);
    add_field<Required>("T_mid", scalar3d, K, gn);
  }
}

void VertDivergenceDiag::initialize_impl(const RunType /*run_type*/) {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto &f      = get_fields_in().front();
  const auto &fid    = f.get_header().get_identifier();
  const auto &layout = fid.get_layout();

  // TODO: support higher-dimensioned input fields
  EKAT_REQUIRE_MSG(layout.rank() >= 2 && layout.rank() <= 2,
                   "Error! Field rank not supported by VertDivergenceDiag.\n"
                   " - field name: " + fid.name() + "\n"
                   " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags().back() == LEV,
                   "Error! VertDivergenceDiag diagnostic expects a layout ending "
                   "with the 'LEV' tag.\n"
                   " - field name  : " + fid.name() + "\n"
                   " - field layout: " + layout.to_string() + "\n");

  ekat::units::Units diag_units = fid.get_units();

  m_denominator = get_field_in("pseudo_density").clone("denominator");

  if (m_divergence_method == "dp") {
    diag_units = fid.get_units() / Pa;
  }
  else if (m_divergence_method == "dz") {
    diag_units = fid.get_units() / m;
  }

  FieldIdentifier d_fid(m_diag_name, layout, diag_units, fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();
}


void VertDivergenceDiag::compute_diagnostic_impl() {
  const auto &f = get_fields_in().front();
  const auto &d = m_diagnostic_output;

  // TODO: support higher-dimensioned input fields
  auto f2d = f.get_view<const Real**>(); 
  auto d2d = d.get_view<Real**>();

  using CO = scream::ColumnOps<DefaultDevice,Real>;
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  const int ncols = m_denominator.get_header().get_identifier().get_layout().dim(0);
  const int nlevs = m_denominator.get_header().get_identifier().get_layout().dim(1);
  const auto policy = ESU::get_default_team_policy(ncols, nlevs);

  // get the denominator first
  if (m_divergence_method == "dp") {
    m_denominator.update(get_field_in("pseudo_density"), sp(1.0), sp(0.0));
  }
  else if (m_divergence_method == "dz") {
    // TODO: for some reason the dz field keeps getting set to 0
    // TODO: as a workaround, just calculate dz here (sigh...)
    // m_denominator.update(get_field_in("dz"), 1.0, 0.0);

    auto dz_v = m_denominator.get_view<Real**>();
    auto dp_v = get_field_in("pseudo_density").get_view<const Real**>();
    auto pm_v = get_field_in("p_mid").get_view<const Real**>();
    auto tm_v = get_field_in("T_mid").get_view<const Real**>();
    auto qv_v = get_field_in("qv").get_view<const Real**>();
    Kokkos::parallel_for(
      "Compute dz for " + m_diagnostic_output.name(), policy, KOKKOS_LAMBDA(const MT &team) {
        const int icol = team.league_rank();
        auto dz_icol = ekat::subview(dz_v, icol);
        auto dp_icol = ekat::subview(dp_v, icol);
        auto pm_icol = ekat::subview(pm_v, icol);
        auto tm_icol = ekat::subview(tm_v, icol);
        auto qv_icol = ekat::subview(qv_v, icol);
        PF::calculate_dz(team, dp_icol, pm_icol, tm_icol, qv_icol, dz_icol);
      });
  }

  auto d_v = m_denominator.get_view<Real**>();
  // calculate df / denominator (setting the first level to be zero)
  // TODO: perhaps we could calculate the values of the field at interfaces using col ops?
  // TODO: though not sure it will be worth it...
  Kokkos::parallel_for(
    "Compute df / denominator for " + m_diagnostic_output.name(), policy, KOKKOS_LAMBDA(const MT &team) {
      const int icol = team.league_rank();
      auto df_icol = ekat::subview(f2d, icol);
      auto dd_icol = ekat::subview(d2d, icol);
      auto dv_icol = ekat::subview(d_v, icol);
      
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&](const int ilev) {
          if (ilev == 0) {
            dd_icol(ilev) = 0;
          }
          else {
            dd_icol(ilev) = (df_icol(ilev) - df_icol(ilev-1)) / dv_icol(ilev);
          }
      });
    });
}

}  // namespace scream
