#include "diagnostics/aerocom_cld.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

#include "share/util/scream_common_physics_functions.hpp"

namespace scream {

AeroComCld::AeroComCld(const ekat::Comm &comm,
                       const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  EKAT_REQUIRE_MSG(
      params.isParameter("AeroComCld Kind"),
      "Error! AeroComCld requires 'AeroComCld Kind' in its "
      "input parameters.\n");

  m_topbot = m_params.get<std::string>("AeroComCld Kind");
  // check if m_topbot is "Bot" or "Top", else error out
  EKAT_REQUIRE_MSG(m_topbot == "Bot" || m_topbot == "Top",
                    "Error! AeroComCld requires 'AeroComCld Kind' "
                    "to be 'Bot' or 'Top' in its input parameters.\n");
}

std::string AeroComCld::name() const { return "AeroComCld" + m_topbot; }

void AeroComCld::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto grid             = grids_manager->get_grid("Physics");
  const auto &grid_name = grid->name();

  const auto nondim = Units::nondimensional();
  const auto micron = m / 1000000;

  m_ncols = grid->get_num_local_dofs();
  m_nlevs = grid->get_num_vertical_levels();
  m_ndiag = 8;

  // Define layouts we need (both inputs and outputs)
  FieldLayout scalar2d_layout{{COL, LEV}, {m_ncols, m_nlevs}};
  FieldLayout vector1d_layout{{COL, CMP}, {m_ncols, m_ndiag}};

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid", scalar2d_layout, K, grid_name);
  add_field<Required>("pseudo_density", scalar2d_layout, Pa, grid_name);
  add_field<Required>("p_mid", scalar2d_layout, Pa, grid_name);
  add_field<Required>("qv", scalar2d_layout, kg / kg, grid_name);
  add_field<Required>("qc", scalar2d_layout, kg / kg, grid_name);
  add_field<Required>("qi", scalar2d_layout, kg / kg, grid_name);
  add_field<Required>("eff_radius_qc", scalar2d_layout, micron, grid_name);
  add_field<Required>("eff_radius_qi", scalar2d_layout, micron, grid_name);
  add_field<Required>("cldfrac_tot", scalar2d_layout, nondim, grid_name);
  add_field<Required>("nc", scalar2d_layout, kg / kg, grid_name);

  // Construct and allocate the aodvis field
  FieldIdentifier fid(name(), vector1d_layout, nondim, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();

  for(int i = 1; i < m_nlevs; ++i) {
    m_level_vector[i] = m_topbot == "Top" ? i : m_nlevs - 1 - i;
  }

  // Self-document the outputs to parse in post-processing
  using stratt_t = std::map<std::string, std::string>;
  auto d         = get_diagnostic();
  auto &metadata =
      d.get_header().get_extra_data<stratt_t>("io: string attributes");
  metadata["0"] = "T_mid";
  metadata["1"] = "p_mid";
  metadata["2"] = "cldfrac_ice";
  metadata["3"] = "cldfrac_liq";
  metadata["4"] = "cdnc";
  metadata["5"] = "eff_radius_qc";
  metadata["6"] = "eff_radius_qi";
  metadata["7"] = "cldfrac_tot";
}

void AeroComCld::compute_diagnostic_impl() {
  /* The goal of this routine/impl is to calculate properties at cloud top
   * based on the AeroCom recommendation. See reference for routine
   * get_subcolumn_mask in rrtmpg, where equation 14 is used for the
   * maximum-random overlap assumption for subcolumn generation. We use
   * equation 13, the column counterpart. The logic is reversed to calculate
   * cloud-bottom properties as well.
   */
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  using PF = scream::PhysicsFunctions<DefaultDevice>;

  const auto out = m_diagnostic_output.get_view<Real **>();

  // zero out the outputs
  Kokkos::deep_copy(out, 0.0);

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

  // Hack: phony field to store dz
  auto dz_f = get_field_in("pseudo_density").clone("dz");
  auto dz   = dz_f.get_view<Real **>();

  // Get gravity acceleration constant from constants
  using physconst = scream::physics::Constants<Real>;
  // TODO: move tunable constant to namelist
  constexpr auto q_threshold = 0.0;
  // TODO: move tunable constant to namelist
  constexpr auto cldfrac_tot_threshold = 0.001;

  const auto num_levs = m_nlevs;
  const auto policy   = ESU::get_default_team_policy(m_ncols, m_nlevs);
  Kokkos::parallel_for(
      "Compute " + name(), policy, KOKKOS_LAMBDA(const MT &team) {
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

        // We need dz too
        // TODO: deduce dz_icol type, etc. on the spot here w/o field cloning?
        auto dz_icol = ekat::subview(dz, icol);
        PF::calculate_dz(team, pden_icol, pmid_icol, tmid_icol, qv_icol,
                         dz_icol);

        // Initialize the 1D "clear fraction" as 1 (totally clear)
        auto clr_icol = 1.0;

        // Loop over all layers in serial (due to accumulative
        // product), starting at 2 (second highest) layer because the
        // highest is assumed to hav no clouds for cldtop, but starting
        // at m_nlevs-1 for cldbot

        for(int ilay : m_level_vector) {
          // Only do the calculation if certain conditions are met
          if((qc_icol(ilay) + qi_icol(ilay)) > q_threshold &&
             (cld_icol(ilay) > cldfrac_tot_threshold)) {
            /* PART I: Probabilistically determining cloud top */
            // Populate aerocom_tmp as the clear-sky fraction
            // probability of this level, where aerocom_clr is that of
            // the previous level
            auto aerocom_tmp =
                clr_icol *
                (1.0 - ekat::impl::max(cld_icol(ilay - 1), cld_icol(ilay))) /
                (1.0 - ekat::impl::min(cld_icol(ilay - 1),
                                       1.0 - cldfrac_tot_threshold));
            // Temporary variable for probability "weights"
            auto aerocom_wts = clr_icol - aerocom_tmp;
            // Temporary variable for liquid "phase"
            auto aerocom_phi = qc_icol(ilay) / (qc_icol(ilay) + qi_icol(ilay));
            /* PART II: The inferred properties */
            /* In general, converting a 3D property X to a 2D cloud-top
             * counterpart x follows: x(i) += X(i,k) * weights * Phase
             * but X and Phase are not always needed */
            // T_mid_at_cldtop
            out(icol, /* T_mid */ 0) += tmid_icol(ilay) * aerocom_wts;
            // p_mid_at_cldtop
            out(icol, /* p_mid */ 1) += pmid_icol(ilay) * aerocom_wts;
            // cldfrac_ice_at_cldtop
            out(icol, /* cldfrac_ice */ 2) += (1.0 - aerocom_phi) * aerocom_wts;
            // cldfrac_liq_at_cldtop
            out(icol, /* cldfrac_liq */ 3) += aerocom_phi * aerocom_wts;
            // cdnc_at_cldtop
            /* We need to convert nc from 1/mass to 1/volume first, and
             * from grid-mean to in-cloud, but after that, the
             * calculation follows the general logic */
            auto cdnc = nc_icol(ilay) * pden_icol(ilay) / dz_icol(ilay) /
                        physconst::gravit / cld_icol(ilay);
            out(icol, /* cdnc */ 4) += cdnc * aerocom_phi * aerocom_wts;
            // eff_radius_qc_at_cldtopbot
            out(icol, /* eff_radius_qc */ 5) +=
                rel_icol(ilay) * aerocom_phi * aerocom_wts;
            // eff_radius_qi_at_cldtopbot
            out(icol, /* eff_radius_qi */ 6) +=
                rei_icol(ilay) * (1.0 - aerocom_phi) * aerocom_wts;
            // Reset aerocom_clr to aerocom_tmp to accumulate
            clr_icol = aerocom_tmp;
          }
        }
        // After the serial loop over levels, the cloudy fraction is
        // defined as (1 - aerocom_clr). This is true because
        // aerocom_clr is the result of accumulative probabilities
        // (their products)
        out(icol, /* cldfrac_tot */ 7) = 1.0 - clr_icol;
      });
}

}  // namespace scream
