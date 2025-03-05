#include "diagnostics/aerocom_cld.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <string>

#include "diagnostics/aerocom_cld_util.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"

namespace scream {

AeroComCld::AeroComCld(const ekat::Comm &comm,
                       const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  EKAT_REQUIRE_MSG(params.isParameter("AeroComCld Kind"),
                   "Error! AeroComCld requires 'AeroComCld Kind' in its "
                   "input parameters.\n");

  m_topbot = m_params.get<std::string>("AeroComCld Kind");
  // check if m_topbot is "Bot" or "Top", else error out
  EKAT_REQUIRE_MSG(m_topbot == "Bot" || m_topbot == "Top",
                   "Error! AeroComCld requires 'AeroComCld Kind' "
                   "to be 'Bot' or 'Top' in its input parameters.\n");
}

void AeroComCld::
set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  auto grid             = grids_manager->get_grid("Physics");
  const auto &grid_name = grid->name();

  const auto nondim = Units::nondimensional();
  const auto micron = m / 1000000;

  // Set the index map and units map
  AeroComCldDiagUtil aercom_util;
  m_index_map = aercom_util.index_map;
  m_units_map = aercom_util.units_map;
  m_ndiag     = aercom_util.size;

  // Ensure m_index_map and m_units_map match
  EKAT_REQUIRE_MSG(m_index_map.size() == m_units_map.size(),
                   "Error! Some inconsistency in AeroComCld: index and units "
                   "maps do not match!\n");
  // Ensure m_index_map and m_ndiag match
  EKAT_REQUIRE_MSG(static_cast<int>(m_index_map.size()) == m_ndiag,
                   "Error! Some inconsistency in AeroComCld: index and units "
                   "maps do not match!\n");

  m_ncols = grid->get_num_local_dofs();
  m_nlevs = grid->get_num_vertical_levels();

  // Define layouts we need (both inputs and outputs)
  auto scalar3d = grid->get_3d_scalar_layout(true);
  auto vector2d = grid->get_2d_vector_layout(m_ndiag);

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid",          scalar3d, K,       grid_name);
  add_field<Required>("pseudo_density", scalar3d, Pa,      grid_name);
  add_field<Required>("p_mid",          scalar3d, Pa,      grid_name);
  add_field<Required>("qv",             scalar3d, kg / kg, grid_name);
  add_field<Required>("qc",             scalar3d, kg / kg, grid_name);
  add_field<Required>("qi",             scalar3d, kg / kg, grid_name);
  add_field<Required>("eff_radius_qc",  scalar3d, micron,  grid_name);
  add_field<Required>("eff_radius_qi",  scalar3d, micron,  grid_name);
  add_field<Required>("cldfrac_tot",    scalar3d, nondim,  grid_name);
  add_field<Required>("nc",             scalar3d, 1 / kg,  grid_name);
  add_field<Required>("ni",             scalar3d, 1 / kg,  grid_name);

  // A field to store dz
  FieldIdentifier m_dz_fid("dz", scalar3d, m, grid_name);
  m_dz = Field(m_dz_fid);
  m_dz.allocate_view();

  // Construct and allocate the output field
  FieldIdentifier fid(name() + m_topbot, vector2d, nondim, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();

  // Self-document the outputs to parse in post-processing
  using stratt_t = std::map<std::string, std::string>;
  auto d         = get_diagnostic();
  auto &metadata =
      d.get_header().get_extra_data<stratt_t>("io: string attributes");
  for(const auto &it : m_index_map) {
    metadata[it.first] =
        std::to_string(it.second) + " (" + m_units_map[it.first] + ")";
  }
}

void AeroComCld::compute_diagnostic_impl() {
  /* The goal of this routine/impl is to calculate properties at cloud top
   * based on the AeroCom recommendation. See reference for routine
   * get_subcolumn_mask in rrtmpg, where equation 14 is used for the
   * maximum-random overlap assumption for subcolumn generation. We use
   * equation 13, the column counterpart. The order is reversed to calculate
   * cloud-bottom properties as well.
   */
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  using PF = scream::PhysicsFunctions<DefaultDevice>;

  // Get the input fields
  const auto tmid = get_field_in("T_mid").get_view<const Real **>();
  const auto pden = get_field_in("pseudo_density").get_view<const Real **>();
  const auto pmid = get_field_in("p_mid").get_view<const Real **>();
  const auto qv   = get_field_in("qv").get_view<const Real **>();
  const auto qc   = get_field_in("qc").get_view<const Real **>();
  const auto qi   = get_field_in("qi").get_view<const Real **>();
  const auto rel  = get_field_in("eff_radius_qc").get_view<const Real **>();
  const auto rei  = get_field_in("eff_radius_qi").get_view<const Real **>();
  const auto cld  = get_field_in("cldfrac_tot").get_view<const Real **>();
  const auto nc   = get_field_in("nc").get_view<const Real **>();
  const auto ni   = get_field_in("ni").get_view<const Real **>();

  auto dz = m_dz.get_view<Real **>();

  // Get gravity acceleration constant from constants
  using physconst = scream::physics::Constants<Real>;

  // Set two potentially tunable parameters
  auto q_threshold           = q_thresh_set();
  auto cldfrac_tot_threshold = cldfrac_tot_thresh_set();

  const auto policy = ESU::get_default_team_policy(m_ncols, m_nlevs);

  const auto out = m_diagnostic_output.get_view<Real **>();

  // zero out the outputs
  Kokkos::deep_copy(out, 0.0);

  const int nlevs = m_nlevs;
  bool is_top     = (m_topbot == "Top");

  // subview the out field for each variable
  auto o_tmid          = ekat::subview_1(out, m_index_map["T_mid"]);
  auto o_pmid          = ekat::subview_1(out, m_index_map["p_mid"]);
  auto o_cldfrac_ice   = ekat::subview_1(out, m_index_map["cldfrac_ice"]);
  auto o_cldfrac_liq   = ekat::subview_1(out, m_index_map["cldfrac_liq"]);
  auto o_cdnc          = ekat::subview_1(out, m_index_map["cdnc"]);
  auto o_eff_radius_qc = ekat::subview_1(out, m_index_map["eff_radius_qc"]);
  auto o_eff_radius_qi = ekat::subview_1(out, m_index_map["eff_radius_qi"]);
  auto o_nc            = ekat::subview_1(out, m_index_map["nc"]);
  auto o_ni            = ekat::subview_1(out, m_index_map["ni"]);
  auto o_cldfrac_tot   = ekat::subview_1(out, m_index_map["cldfrac_tot"]);

  Kokkos::parallel_for(
      "Compute " + m_diagnostic_output.name(), policy, KOKKOS_LAMBDA(const MT &team) {
        const int icol = team.league_rank();
        // Subview the inputs at icol
        auto tmid_icol = ekat::subview(tmid, icol);
        auto pden_icol = ekat::subview(pden, icol);
        auto pmid_icol = ekat::subview(pmid, icol);
        auto qv_icol   = ekat::subview(qv, icol);
        auto qc_icol   = ekat::subview(qc, icol);
        auto qi_icol   = ekat::subview(qi, icol);
        auto rel_icol  = ekat::subview(rel, icol);
        auto rei_icol  = ekat::subview(rei, icol);
        auto cld_icol  = ekat::subview(cld, icol);
        auto nc_icol   = ekat::subview(nc, icol);
        auto ni_icol   = ekat::subview(ni, icol);

        // We need dz too
        auto dz_icol = ekat::subview(dz, icol);
        PF::calculate_dz(team, pden_icol, pmid_icol, tmid_icol, qv_icol,
                         dz_icol);

        // Initialize the 1D "clear fraction" as 1 (totally clear)
        auto clr_icol = 1.0;

        // Loop over all layers in serial (due to accumulative
        // product), starting at 2 (second highest) layer because the
        // highest is assumed to hav no clouds for cldtop, but starting
        // at m_nlevs-1 for cldbot

        auto topbot_calcs = [&](int ilay) {
          // Only do the calculation if certain conditions are met
          if((qc_icol(ilay) + qi_icol(ilay)) > q_threshold &&
             (cld_icol(ilay) > cldfrac_tot_threshold)) {
            /* PART I: Probabilistically determining cloud top/bot */
            // Populate clr_tmp as the clear-sky fraction
            // probability of this level, where clr_icol is that of
            // the previous level
            auto clr_tmp =
                clr_icol *
                (1.0 - ekat::impl::max(cld_icol(ilay - 1), cld_icol(ilay))) /
                (1.0 - ekat::impl::min(cld_icol(ilay - 1),
                                       Real(1.0 - cldfrac_tot_threshold)));
            // Temporary variable for probability "weights"
            auto wts = clr_icol - clr_tmp;
            // Temporary variable for liquid "phase"
            auto phi = qc_icol(ilay) / (qc_icol(ilay) + qi_icol(ilay));
            /* PART II: The inferred properties */
            /* In general, converting a 3D property X to a 2D cloud-top
             * counterpart x follows: x(i) += X(i,k) * weights * Phase
             * but X and Phase are not always needed */
            o_tmid(icol) += tmid_icol(ilay) * wts;
            o_pmid(icol) += pmid_icol(ilay) * wts;
            o_cldfrac_ice(icol) += (1.0 - phi) * wts;
            o_cldfrac_liq(icol) += phi * wts;
            // cdnc
            /* We need to convert nc from 1/mass to 1/volume first, and
             * from grid-mean to in-cloud, but after that, the
             * calculation follows the general logic */
            auto cdnc = nc_icol(ilay) * pden_icol(ilay) / dz_icol(ilay) /
                        physconst::gravit / cld_icol(ilay);
            o_cdnc(icol) += cdnc * phi * wts;
            o_nc(icol) += nc_icol(ilay) * phi * wts;
            o_ni(icol) += ni_icol(ilay) * (1.0 - phi) * wts;
            o_eff_radius_qc(icol) += rel_icol(ilay) * phi * wts;
            o_eff_radius_qi(icol) += rei_icol(ilay) * (1.0 - phi) * wts;
            // Reset clr_icol to clr_tmp to accumulate
            clr_icol = clr_tmp;
          }
        };

        if(is_top) {
          for(int ilay = 1; ilay < nlevs; ++ilay) {
            topbot_calcs(ilay);
          }
        } else {
          for(int ilay = nlevs - 1; ilay > 0; --ilay) {
            topbot_calcs(ilay);
          }
        }

        // After the serial loop over levels, the cloudy fraction is
        // defined as (1 - clr_icol). This is true because
        // clr_icol is the result of accumulative probabilities
        // (their products)
        o_cldfrac_tot(icol) = 1.0 - clr_icol;
      });
}

}  // namespace scream
