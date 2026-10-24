#include "aodvis.hpp"

#include "share/util/eamxx_universal_constants.hpp"
#include "share/util/eamxx_utils.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_reduction_utils.hpp>

namespace scream {

AODVis::
AODVis(const ekat::Comm &comm, const ekat::ParameterList &params,
       const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  m_swbands = eamxx_swbands();
  m_vis_bnd = eamxx_vis_swband_idx();

  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The fields required for this diagnostic to be computed
  m_field_in_names.push_back("aero_tau_sw");
  m_field_in_names.push_back("sunlit_mask");

  // Construct and allocate the aodvis field
  FieldIdentifier fid(name(), m_grid->get_2d_scalar_layout(), none, m_grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void AODVis::initialize_impl() {
  // The sunlit mask, under a name of our own: IO shares an avg count between diags
  // whose masks have the same name, and that is for diags masked the same way by design
  const auto& sunlit = m_fields_in.at("sunlit_mask");
  m_diagnostic_output.set_valid_mask(sunlit.alias(m_diagnostic_output.name()+"_valid_mask"));
  // Night columns hold the fill value, so averaging must skip them
  m_diagnostic_output.get_header().set_may_be_filled(true);
}

void AODVis::compute_impl()
{
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using TPF = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  using RU  = ekat::ReductionUtils<typename KT::ExeSpace>;

  constexpr auto fill_value = constants::fill_value<Real>;

  const auto aod     = m_diagnostic_output.get_view<Real *>();
  const auto tau_vis = m_fields_in.at("aero_tau_sw")
                           .subfield(1, m_vis_bnd)
                           .get_view<const Real **>();
  const auto sunlit = m_fields_in.at("sunlit_mask").get_view<const int *>();

  const auto nlevs = m_grid->get_num_vertical_levels();
  const auto ncols = m_grid->get_num_local_dofs();
  const auto policy   = TPF::get_default_team_policy(ncols, nlevs);
  Kokkos::parallel_for(
      "Compute " + m_diagnostic_output.name(), policy, KOKKOS_LAMBDA(const MT &team) {
        const int icol = team.league_rank();
        if(sunlit(icol)) {
          auto tau_icol = ekat::subview(tau_vis, icol);
          aod(icol)     = RU::view_reduction(team, 0, nlevs, tau_icol);
        } else {
          aod(icol) = fill_value;
        }
      });
}

}  // namespace scream
