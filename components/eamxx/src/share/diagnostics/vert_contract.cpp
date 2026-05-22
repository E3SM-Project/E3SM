#include "vert_contract.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream {

VertContract::
VertContract(const ekat::Comm &comm,
             const ekat::ParameterList &params,
             const std::shared_ptr<const AbstractGrid> &grid)
 : AbstractDiagnostic(comm, params, grid)
{
  m_field_name = m_params.get<std::string>("field_name");

  m_field_in_names.push_back(m_field_name);

  // we support either sum or avg
  m_contract_method = m_params.get<std::string>("contract_method");
  EKAT_REQUIRE_MSG(
      m_contract_method == "avg" || m_contract_method == "sum",
      "Error! VertContract only supports 'avg' or 'sum' as contract_method.\n"
      " - contract_method: " + m_contract_method + "\n");
  // we support either dp or dz weighting, or no weighting at all (none)
  m_weighting_method = m_params.get<std::string>("weighting_method", "none");
  EKAT_REQUIRE_MSG(
      m_weighting_method == "dp" || m_weighting_method == "dz" || m_weighting_method == "none",
      "Error! VertContract only supports 'dp' or 'dz' or 'none' as weighting_method.\n"
      " - weighting_method: " + m_weighting_method + "\n");

  if (m_weighting_method == "dp") {
    m_field_in_names.push_back("pseudo_density");
  } else if (m_weighting_method == "dz") {
    m_field_in_names.push_back("pseudo_density");
    m_field_in_names.push_back("qv");
    m_field_in_names.push_back("p_mid");
    m_field_in_names.push_back("T_mid");
  }
}

void VertContract::initialize_impl()
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto &f      = m_fields_in.at(m_field_name);
  const auto &fid    = f.get_header().get_identifier();
  const auto &layout = fid.get_layout();
  const auto vtag    = layout.tags().back();

  // if dp or dz, gotta give us ncol,nlev as inputs (so rank>=2)
  // else, can give us just nlev
  auto min_rank = (m_weighting_method == "none") ? 1 : 2;
  EKAT_REQUIRE_MSG(layout.rank() >= min_rank && layout.rank() <= 3,
                   "Error! Field rank not supported by VertContract.\n"
                   "For dp- and dz-weighted contractions, need at rank-2 fields.\n"
                   " - field name: " + fid.name() + "\n"
                   " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags().back() == LEV,
                   "Error! VertContract diagnostic expects a layout ending "
                   "with the 'LEV' tag.\n"
                   " - field name  : " + fid.name() + "\n"
                   " - field layout: " + layout.to_string() + "\n");

  auto diag_units = fid.get_units();

  auto w_units = none;
  // set up the weighting field
  if (m_weighting_method == "dp") {
    m_weight = m_fields_in.at("pseudo_density").clone("vert_contract_wts");
    constexpr auto g_units = physics::Constants<Real>::gravit.units;
    w_units = m_weight.get_header().get_identifier().get_units() / g_units;
  } else if (m_weighting_method == "dz") {
    // TODO: for some reason the dz field keeps getting set to 0
    // TODO: as a workaround, just calculate dz here (sigh...)
    // m_weight = get_field_in("dz").clone("vert_contract_wts");
    m_weight = m_fields_in.at("pseudo_density").clone("vert_contract_wts");
    w_units = m;
  } else {
    // no weighting needed, so we set it to 1 with layout of (lev)
    FieldLayout layout_wts = {{vtag}, {layout.dim(vtag)}};
    FieldIdentifier w_fid("vert_contract_wts", layout_wts, w_units, fid.get_grid_name());
    m_weight = Field(w_fid,true);
    m_weight.deep_copy(1);
  }

  diag_units = diag_units * w_units;
  if (m_contract_method == "avg") {
    diag_units = diag_units / w_units;
  }

  auto diag_name = m_field_name + "_vert_" + m_contract_method;
  if (m_weighting_method != "none")
    diag_name += "_" + m_weighting_method + "_weighted";

  auto d_fid = fid.clone(diag_name).reset_layout(layout.clone().strip_dim(LEV)).reset_units(diag_units);
  m_diagnostic_output = Field(d_fid,true);

  if (m_contract_method == "avg") {
    // m_weight_sum must have the same layout as d (f's layout without vtag)
    // so that scale_inv(m_weight_sum) is valid regardless of input rank.
    auto wsum_fid = fid.clone("weight_sum")
                       .reset_layout(layout.clone().strip_dim(vtag))
                       .reset_units(w_units);
    m_weight_sum = Field(wsum_fid,true);

    // m_ones is a ones-field with f's layout (not m_weight's layout).
    // Using f's layout ensures vert_contraction(m_weight_sum, m_ones, m_weight)
    // produces the correct (COL[xCMP]) denominator for any input rank.
    m_ones = f.clone("ones");
    m_ones.deep_copy(1);
    if (f.has_valid_mask()) {
      // Share the valid_mask so the denominator honours the same masking
      m_ones.set_valid_mask(f.get_valid_mask());
      // Output gets a valid_mask: 1 where at least one level was valid
      m_diagnostic_output.create_valid_mask(Field::MaskInit::None);
      m_diagnostic_output.get_header().set_may_be_filled(true);
    }
  }
}

void VertContract::compute_impl()
{
  const auto &f = m_fields_in.at(m_field_name);
  auto& d = m_diagnostic_output;

  // Update the weights; if weighting by dp, we need to scale by 1/g
  if (m_weighting_method == "dp") {
    auto g = scream::physics::Constants<Real>::gravit.value;
    m_weight.update(m_fields_in.at("pseudo_density"), 1 / g, sp(0.0));
  } else if (m_weighting_method == "dz") {
    // TODO: for some reason the dz field keeps getting set to 0
    // TODO: as a workaround, just calculate dz here (sigh...)
    // m_weight.update(get_field_in("dz"), 1.0, 0.0);
    using KT          = KokkosTypes<DefaultDevice>;
    using MT          = typename KT::MemberType;
    using TPF         = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
    using PF          = scream::PhysicsFunctions<DefaultDevice>;
    const int ncols   = m_weight.get_header().get_identifier().get_layout().dim(0);
    const int nlevs   = m_weight.get_header().get_identifier().get_layout().dim(1);
    const auto policy = TPF::get_default_team_policy(ncols, nlevs);

    auto dz_v = m_weight.get_view<Real **>();
    auto dp_v = m_fields_in.at("pseudo_density").get_view<const Real **>();
    auto pm_v = m_fields_in.at("p_mid").get_view<const Real **>();
    auto tm_v = m_fields_in.at("T_mid").get_view<const Real **>();
    auto qv_v = m_fields_in.at("qv").get_view<const Real **>();
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

  // Now do the contraction itself
  vert_contraction(d, f, m_weight);

  if (m_contract_method == "avg") {
    // Denominator: sum_lev(weight * 1 * mask); layout matches d
    vert_contraction(m_weight_sum, m_ones, m_weight);
    if (d.has_valid_mask()) {
      // Mark where denominator is nonzero, then divide, then fill elsewhere
      auto& valid_out = d.get_valid_mask();
      compute_mask(m_weight_sum, 0, Comparison::NE, valid_out);
      d.scale_inv(m_weight_sum, valid_out);
      // IO relies on fill_value for masked-out entries
      d.deep_copy(constants::fill_value<Real>, valid_out, true);
    } else {
      d.scale_inv(m_weight_sum);
    }
  }
}

} // namespace scream

