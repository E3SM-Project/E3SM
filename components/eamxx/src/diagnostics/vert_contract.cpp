#include "diagnostics/vert_contract.hpp"

#include "physics/share/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream {

VertContractDiag::VertContractDiag(const ekat::Comm &comm,
                                   const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {}

void VertContractDiag::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto &fn = m_params.get<std::string>("field_name");
  const auto &gn = m_params.get<std::string>("grid_name");
  const auto g   = grids_manager->get_grid(gn);

  add_field<Required>(fn, gn);

  // we support either sum or avg
  m_contract_method = m_params.get<std::string>("contract_method");
  EKAT_REQUIRE_MSG(
      m_contract_method == "avg" || m_contract_method == "sum",
      "Error! VertContractDiag only supports 'avg' or 'sum' as contract_method.\n"
      " - contract_method: " + m_contract_method + "\n");
  // we support either dp or dz weighting, or no weighting at all (none)
  m_weighting_method = m_params.get<std::string>("weighting_method", "none");
  EKAT_REQUIRE_MSG(
      m_weighting_method == "dp" || m_weighting_method == "dz" || m_weighting_method == "none",
      "Error! VertContractDiag only supports 'dp' or 'dz' or 'none' as weighting_method.\n"
      " - weighting_method: " + m_weighting_method + "\n");
  
  m_diag_name = fn + "_vert_" + m_contract_method;
  // append weighting_method to name if needed
  m_diag_name = (m_weighting_method == "none") ? m_diag_name : m_diag_name + "_" + m_weighting_method + "_weighted";

  auto scalar3d = g->get_3d_scalar_layout(true);
  if (m_weighting_method == "dp") {
    add_field<Required>("pseudo_density", scalar3d, Pa, gn);
  } else if (m_weighting_method == "dz") {
    // TODO: for some reason the dz field keeps getting set to 0
    // TODO: as a workaround, just calculate dz here (sigh...)
    // add_field<Required>("dz", scalar3d, m, gn);
    add_field<Required>("pseudo_density", scalar3d, Pa, gn);
    add_field<Required>("qv", scalar3d, kg / kg, gn);
    add_field<Required>("p_mid", scalar3d, Pa, gn);
    add_field<Required>("T_mid", scalar3d, K, gn);

  }
}

void VertContractDiag::initialize_impl(const RunType /*run_type*/) {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto &f      = get_fields_in().front();
  const auto &fid    = f.get_header().get_identifier();
  const auto &layout = fid.get_layout();

  // if dp or dz, gotta give us ncol,nlev as inputs (so rank>=2)
  // else, can give us just nlev
  auto min_rank = (m_weighting_method == "none") ? 1 : 2;
  EKAT_REQUIRE_MSG(layout.rank() >= min_rank && layout.rank() <= 3,
                   "Error! Field rank not supported by VertContractDiag.\n"
                   "For dp- and dz-weighted contractions, need at rank-2 fields.\n"
                   " - field name: " + fid.name() + "\n"
                   " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags().back() == LEV,
                   "Error! VertContractDiag diagnostic expects a layout ending "
                   "with the 'LEV' tag.\n"
                   " - field name  : " + fid.name() + "\n"
                   " - field layout: " + layout.to_string() + "\n");

  ekat::units::Units diag_units = fid.get_units();

  // set up the weighting fields
  if (m_weighting_method == "dp") {
    m_weighting = get_field_in("pseudo_density").clone("vert_contract_wts");
  } else if (m_weighting_method == "dz") {
    // TODO: for some reason the dz field keeps getting set to 0
    // TODO: as a workaround, just calculate dz here (sigh...)
    // m_weighting = get_field_in("dz").clone("vert_contract_wts");
    m_weighting = get_field_in("pseudo_density").clone("vert_contract_wts");
  } else {
    // no weighting needed, so we set it to 1 with layout of (lev)
    FieldLayout layout_wts = {{LEV}, {layout.dim(LEV)}};
    FieldIdentifier f_id("vert_contract_wts", layout_wts, ekat::units::Units::nondimensional(), fid.get_grid_name());
    m_weighting = Field(f_id);
    m_weighting.allocate_view();
    m_weighting.deep_copy(sp(1));
  }

  if (m_weighting_method == "dp" && m_contract_method == "sum") {
    // we scale by the weighting, so we use fid units * Pa (but we scale by 1/g for dp)
    diag_units = fid.get_units() * Pa / (m/(s*s));
  } else if (m_weighting_method == "dz" && m_contract_method == "sum") {
    // we scale by the weighting, so we use fid units * m
    diag_units = fid.get_units() * m;
  }

  if (m_contract_method == "avg") {
    auto wts_layout = m_weighting.get_header().get_identifier().get_layout();
    FieldIdentifier wts_sum_fid("vert_contract_wts_sum", wts_layout.clone().strip_dim(LEV), diag_units, fid.get_grid_name());
    m_weighting_sum = Field(wts_sum_fid);
    m_weighting_sum.allocate_view();
    m_weighting_one = m_weighting.clone("vert_contract_wts_one");
    m_weighting_one.deep_copy(sp(1));
    vert_contraction<Real>(m_weighting_sum, m_weighting, m_weighting_one);
    VertContractDiag::scale_wts(m_weighting, m_weighting_sum);
  }

  FieldIdentifier d_fid(m_diag_name, layout.clone().strip_dim(LEV), diag_units, fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  if (f.get_header().has_extra_data("mask_data")) {
    m_diagnostic_output.get_header().set_extra_data("mask_data", m_diagnostic_output.clone(m_diag_name+"_mask"));
    m_diagnostic_output.get_header().set_extra_data("mask_value", f.get_header().get_extra_data<Real>("mask_value"));
  }
}

// TODO: move this to field_utils.hpp
//       by allowing update fxns there to op
//       on fields of ranks \in ranks, e.g.,
//       f1.scale_inv(f2) should work for:
//       - f2 scalar (rank-0)
//       - f2 with same layout as f1
//       - f2 with layout that is a subset of f1
//         (ncol,lev) is subset of (ncol, dim, lev)
//         (ncol) is subset of (ncol, lev), etc.
void VertContractDiag::scale_wts(Field &wts, const Field &wts_sum) {
  using KT = KokkosTypes<DefaultDevice>;
  using RP = typename KT::RangePolicy;

  auto wts_l          = wts.get_header().get_identifier().get_layout();
  const auto wts_rank = wts_l.rank();

  if (wts_rank == 1) {
    // no ncols, just nlevs
    const int nlevs      = wts_l.dim(0);
    const auto wts_v     = wts.get_view<Real *>();
    const auto wts_sum_v = wts_sum.get_view<const Real>();

    Kokkos::parallel_for(
        "VertContractDiag::scale_wts" + m_diag_name, RP(0, nlevs), KOKKOS_LAMBDA(const int &ilev) {
          if (wts_sum_v() != 0) {
            wts_v(ilev) /= wts_sum_v();
          } else {
            wts_v(ilev) = 0; // Handle division by zero by setting to 0
          }
        });

  } else if (wts_rank == 2) {
    // we have both ncols and nlevs
    const int ncols      = wts_l.dim(0);
    const int nlevs      = wts_l.dim(1);
    const auto wts_v     = wts.get_view<Real **>();
    const auto wts_sum_v = wts_sum.get_view<const Real *>();
    Kokkos::parallel_for(
        "VertContractDiag::scale_wts" + m_diag_name, RP(0, nlevs * ncols),
        KOKKOS_LAMBDA(const int &idx) {
          const int icol = idx / nlevs;
          const int ilev = idx % nlevs;
          if (wts_sum_v(icol) != 0) {
            wts_v(icol, ilev) /= wts_sum_v(icol);
          } else {
            wts_v(icol, ilev) = 0; // Handle division by zero by setting to 0
          }
        });
  } else {
    // you shouldn't have arrived here, error out
    EKAT_ERROR_MSG("Error!! VertContractDiag::scale_wts, unexpected field layout.");
  }
}

void VertContractDiag::compute_diagnostic_impl() {
  const auto &f = get_fields_in().front();
  const auto &d = m_diagnostic_output;

  // update the weights; if weighting by dp, we need to scale by 1/g
  if (m_weighting_method == "dp") {
    auto g = scream::physics::Constants<Real>::gravit;
    m_weighting.update(get_field_in("pseudo_density"), 1 / g, sp(0.0));
  } else if (m_weighting_method == "dz") {
    // TODO: for some reason the dz field keeps getting set to 0
    // TODO: as a workaround, just calculate dz here (sigh...)
    // m_weighting.update(get_field_in("dz"), 1.0, 0.0);
    using KT          = KokkosTypes<DefaultDevice>;
    using MT          = typename KT::MemberType;
    using TPF         = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
    using PF          = scream::PhysicsFunctions<DefaultDevice>;
    const int ncols   = m_weighting.get_header().get_identifier().get_layout().dim(0);
    const int nlevs   = m_weighting.get_header().get_identifier().get_layout().dim(1);
    const auto policy = TPF::get_default_team_policy(ncols, nlevs);

    auto dz_v = m_weighting.get_view<Real **>();
    auto dp_v = get_field_in("pseudo_density").get_view<const Real **>();
    auto pm_v = get_field_in("p_mid").get_view<const Real **>();
    auto tm_v = get_field_in("T_mid").get_view<const Real **>();
    auto qv_v = get_field_in("qv").get_view<const Real **>();
    Kokkos::parallel_for(
        "Compute dz for " + m_diagnostic_output.name(), policy, KOKKOS_LAMBDA(const MT &team) {
          const int icol = team.league_rank();
          auto dz_icol   = ekat::subview(dz_v, icol);
          auto dp_icol   = ekat::subview(dp_v, icol);
          auto pm_icol   = ekat::subview(pm_v, icol);
          auto tm_icol   = ekat::subview(tm_v, icol);
          auto qv_icol   = ekat::subview(qv_v, icol);
          PF::calculate_dz(team, dp_icol, pm_icol, tm_icol, qv_icol, dz_icol);
        });
  }

  // if dp|dz_weighted and avg, we need to scale the weighting by its 1/sum
  if ((m_weighting_method == "dp" || m_weighting_method == "dz") && m_contract_method == "avg") {
    vert_contraction<Real>(m_weighting_sum, m_weighting, m_weighting_one);
    VertContractDiag::scale_wts(m_weighting, m_weighting_sum);
  }

  // call the vert_contraction impl that will take care of everything
  // if f has a mask and we are averaging, need to call the avg specialization
  if (m_contract_method == "avg" && f.get_header().has_extra_data("mask_data")) {
    vert_contraction<Real,1>(d, f, m_weighting);
  } else {
    vert_contraction<Real,0>(d, f, m_weighting);
  }
}

} // namespace scream
