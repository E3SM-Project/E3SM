#include "diagnostics/pbl_entrainment_budget.hpp"

#include "diagnostics/pbl_entrainment_budget_util.hpp"
#include "ekat/ekat_workspace.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "share/util/scream_universal_constants.hpp"

namespace scream {

PBLEntrainmentBudget::PBLEntrainmentBudget(const ekat::Comm &comm,
                                           const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  // Nothing to do here
}

void PBLEntrainmentBudget::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto grid             = grids_manager->get_grid("Physics");
  const auto &grid_name = grid->name();

  const auto nondim = Units::nondimensional();

  // Set the index map and units map
  PBLEntrainmentBudgetDiagUtil eadu;
  m_index_map = eadu.index_map;
  m_units_map = eadu.units_map;
  m_ndiag     = eadu.size;

  if(eadu.pblinvalg == "temperature-inversion") {
    m_pblinvalg = 1;
  } else if(eadu.pblinvalg == "thetal-only") {
    m_pblinvalg = 2;
  } else if(eadu.pblinvalg == "qt_only") {
    m_pblinvalg = 3;
  } else {
    EKAT_ERROR_MSG(
        "Error! Invalid pblinvalg. Only temperature-inversion, thetal-only, "
        "and "
        "qt_only are currently supported.\n");
  }

  // Ensure m_index_map and m_units_map match
  EKAT_REQUIRE_MSG(
      m_index_map.size() == m_units_map.size(),
      "Error! Some inconsistency in PBLEntrainmentBudget: index and units "
      "maps do not match!\n");
  // Ensure m_index_map and m_ndiag match
  EKAT_REQUIRE_MSG(
      static_cast<int>(m_index_map.size()) == m_ndiag,
      "Error! Some inconsistency in PBLEntrainmentBudget: m_ndiag and index "
      "map do not match!\n");

  m_ncols = grid->get_num_local_dofs();
  m_nlevs = grid->get_num_vertical_levels();

  // Define layouts we need (both inputs and outputs)
  FieldLayout scalar2d_layout{{COL, LEV}, {m_ncols, m_nlevs}};
  FieldLayout vector1d_layout{{COL, CMP}, {m_ncols, m_ndiag}};

  // The fields required for this diagnostic to be computed
  // Get qc and qv
  add_field<Required>("qc", scalar2d_layout, kg / kg, grid_name);
  add_field<Required>("qv", scalar2d_layout, kg / kg, grid_name);
  // Get T_mid, p_mid
  add_field<Required>("T_mid", scalar2d_layout, K, grid_name);
  add_field<Required>("p_mid", scalar2d_layout, Pa, grid_name);
  // Get pseudo_density
  add_field<Required>("pseudo_density", scalar2d_layout, Pa, grid_name);
  // Get radiation up and down terms
  add_field<Required>("SW_flux_dn", scalar2d_layout, W / m * m, grid_name);
  add_field<Required>("SW_flux_up", scalar2d_layout, W / m * m, grid_name);
  add_field<Required>("LW_flux_dn", scalar2d_layout, W / m * m, grid_name);
  add_field<Required>("LW_flux_up", scalar2d_layout, W / m * m, grid_name);

  // Construct and allocate the output field
  FieldIdentifier fid("PBLEntrainmentBudget", vector1d_layout, nondim,
                      grid_name);
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

void PBLEntrainmentBudget::initialize_impl(const RunType /*run_type*/) {
  // Field qt will have units and layout similar to qc, qv
  const auto &qv   = get_field_in("qv");
  const auto &qvid = qv.get_header().get_identifier();
  const auto &qvgn = qvid.get_grid_name();
  const auto &qvlo = qvid.get_layout();
  FieldIdentifier qf_prev("qtot_prev", qvlo.clone(), qvid.get_units(), qvgn);
  m_prev_qt = Field(qf_prev);
  m_prev_qt.allocate_view();
  // Field tl will have units of and layout similar to T_mid
  const auto &tm   = get_field_in("T_mid");
  const auto &tmid = tm.get_header().get_identifier();
  const auto &tmgn = tmid.get_grid_name();
  const auto &tmlo = tmid.get_layout();
  FieldIdentifier tf_prev("tliq_prev", tmlo.clone(), tmid.get_units(), tmgn);
  m_prev_tl = Field(tf_prev);
  m_prev_tl.allocate_view();
}

void PBLEntrainmentBudget::calc_tl_qt(const view_2d &tm_v, const view_2d &pm_v,
                                      const view_2d &qv_v, const view_2d &qc_v,
                                      const view_2d &tl_v,
                                      const view_2d &qt_v) {
  int ncols = m_ncols;
  int nlevs = m_nlevs;
  Kokkos::parallel_for(
      Kokkos::RangePolicy<>(0, ncols * nlevs), KOKKOS_LAMBDA(const int &idx) {
        const int icol   = idx / nlevs;
        const int jlev   = idx % nlevs;
        qt_v(icol, jlev) = qc_v(icol, jlev) + qv_v(icol, jlev);
        tl_v(icol, jlev) = PF::calculate_thetal_from_theta(
            PF::calculate_theta_from_T(tm_v(icol, jlev), pm_v(icol, jlev)),
            tm_v(icol, jlev), qc_v(icol, jlev));
      });
}

void PBLEntrainmentBudget::init_timestep(const util::TimeStamp &start_of_step) {
  m_start_t = start_of_step;

  const auto &tm_v = get_field_in("T_mid").get_view<Real **>();
  const auto &pm_v = get_field_in("p_mid").get_view<Real **>();
  const auto &qv_v = get_field_in("qv").get_view<Real **>();
  const auto &qc_v = get_field_in("qc").get_view<Real **>();

  const auto &m_prev_qt_v = m_prev_qt.get_view<Real **>();
  const auto &m_prev_tl_v = m_prev_tl.get_view<Real **>();

  calc_tl_qt(tm_v, pm_v, qv_v, qc_v, m_prev_tl_v, m_prev_qt_v);
}

void PBLEntrainmentBudget::compute_diagnostic_impl() {
  using PC  = scream::physics::Constants<Real>;
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  constexpr Real g = PC::gravit;
  Real fill_value  = constants::DefaultFillValue<Real>().value;

  // Before doing anything, subview the out field for each variable
  auto out = m_diagnostic_output.get_view<Real **>();

  auto o_pm_hplus = ekat::subview_1(out, m_index_map["p+"]);
  auto o_tl_hplus = ekat::subview_1(out, m_index_map["tl+"]);
  auto o_tl_caret = ekat::subview_1(out, m_index_map["tl^"]);
  auto o_tl_ttend = ekat::subview_1(out, m_index_map["tl_ttend"]);
  auto o_qt_hplus = ekat::subview_1(out, m_index_map["qt+"]);
  auto o_qt_caret = ekat::subview_1(out, m_index_map["qt^"]);
  auto o_qt_ttend = ekat::subview_1(out, m_index_map["qt_ttend"]);
  auto o_df_inpbl = ekat::subview_1(out, m_index_map["dF"]);

  // Get the input views
  const auto &qc_v = get_field_in("qc").get_view<Real **>();
  const auto &qv_v = get_field_in("qv").get_view<Real **>();
  const auto &tm_v = get_field_in("T_mid").get_view<Real **>();
  const auto &pm_v = get_field_in("p_mid").get_view<Real **>();
  const auto &pd_v = get_field_in("pseudo_density").get_view<Real **>();
  const auto &sd_v = get_field_in("SW_flux_dn").get_view<Real **>();
  const auto &su_v = get_field_in("SW_flux_dn").get_view<Real **>();
  const auto &ld_v = get_field_in("LW_flux_dn").get_view<Real **>();
  const auto &lu_v = get_field_in("LW_flux_up").get_view<Real **>();

  // tracked stuff
  const auto &prev_qtot_v = m_prev_qt.get_view<Real **>();
  const auto &prev_tliq_v = m_prev_tl.get_view<Real **>();

  view_2d qt_v("qt_v", m_ncols, m_nlevs);
  view_2d tl_v("tl_v", m_ncols, m_nlevs);
  calc_tl_qt(tm_v, pm_v, qv_v, qc_v, tl_v, qt_v);

  const auto &curr_ts =
      get_field_in("qc").get_header().get_tracking().get_time_stamp();

  auto dt = curr_ts - m_start_t;

  const int num_levs  = m_nlevs;
  const int pblinvalg = m_pblinvalg;
  const auto policy   = ESU::get_default_team_policy(m_ncols, m_nlevs);

  constexpr int wsms = 1;
  using WSMgr        = ekat::WorkspaceManager<Real, DefaultDevice>;
  WSMgr wsm(num_levs, wsms, policy);

  Kokkos::parallel_for(
      "Compute " + name(), policy, KOKKOS_LAMBDA(const MT &team) {
        const int icol = team.league_rank();
        // inputs
        const auto qc_icol = ekat::subview(qc_v, icol);
        const auto qv_icol = ekat::subview(qv_v, icol);
        const auto tm_icol = ekat::subview(tm_v, icol);
        const auto pm_icol = ekat::subview(pm_v, icol);
        const auto pd_icol = ekat::subview(pd_v, icol);
        const auto sd_icol = ekat::subview(sd_v, icol);
        const auto su_icol = ekat::subview(su_v, icol);
        const auto ld_icol = ekat::subview(ld_v, icol);
        const auto lu_icol = ekat::subview(lu_v, icol);

        // tracked
        const auto qt_icol = ekat::subview(qt_v, icol);
        const auto tl_icol = ekat::subview(tl_v, icol);

        auto prev_qtot_icol = ekat::subview(prev_qtot_v, icol);
        auto prev_tliq_icol = ekat::subview(prev_tliq_v, icol);

        auto ws = wsm.get_workspace(team);
        ekat::Unmanaged<WSMgr::view_1d<Real>> tm_grad;
        ws.take_many_contiguous_unsafe<wsms>({"tm_grad"}, {&tm_grad});

        // We first want to find the PBL inversion. There are three methods to
        // do so. All our methods here rely on the *gradient* of state fields
        // (for now). First, we can simply find the first level from the surface
        // that has a a temperature "inversion" (temperature goes positive
        // instead of negative). Second, we can find the level which has the
        // biggest positive jump in theta_l. Third, we can find the level which
        // has the biggest negative jump in qt.

        int opt_tm_grad_lev = 1;
        // Find tm_grad (tm_grad is a catch-all for the 3 methods)
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, 1, num_levs), [&](int k) {
              auto pm_diff = pm_icol(k - 1) - pm_icol(k);
              if(pblinvalg == 1) {
                // pblinvalg = 1 ---> based solely on
                // d(T_mid)/d(p_mid); finding the min (because
                // d(p_mid) < 0), so keeping signs
                tm_grad(k) = (tm_icol(k - 1) - tm_icol(k)) / pm_diff;
              } else if(pblinvalg == 2) {
                // pblinvalg = 2 ---> based solely on
                // d(theta_l)/d(p_mid); finding the min
                // (because d(p_mid) < 0), so keeping signs
                tm_grad(k) = (tl_icol(k - 1) - tl_icol(k)) / pm_diff;
              } else if(pblinvalg == 3) {
                // pblinvalg = 3 ---> based solely on
                // d(q_t)/d(p_mid); finding the max (because
                // d(p_mid) < 0), so reversing signs
                tm_grad(k) = -(qt_icol(k - 1) - qt_icol(k)) / pm_diff;
              }
            });
        team.team_barrier();

        // Find minimum gradient, because d(p_mid) < 0 in definition above
        // Starting from the surface, and ensuring p_mid > 70000.0 Pa,
        // to avoid resolving to some odd place higher up in the atmosphere.
        using minloc_t       = Kokkos::MinLoc<Real, int>;
        using minloc_value_t = typename minloc_t::value_type;
        minloc_value_t minloc;
        Kokkos::parallel_reduce(
            Kokkos::TeamVectorRange(team, 1, num_levs),
            [&](const int &k, minloc_value_t &result) {
              if(tm_grad(k) < result.val && pm_icol(k) > 70000.0) {
                result.val = tm_grad(k);
                result.loc = k;
              }
            },
            minloc_t(minloc));
        team.team_barrier();
        opt_tm_grad_lev = minloc.loc;

        if(opt_tm_grad_lev < 2 || opt_tm_grad_lev > num_levs - 1) {
          // Weird stuff can happen near the top and bottom of atm, so fill_val
          o_pm_hplus(icol) = fill_value;
          o_tl_hplus(icol) = fill_value;
          o_qt_hplus(icol) = fill_value;
          o_df_inpbl(icol) = fill_value;
          o_tl_caret(icol) = fill_value;
          o_qt_caret(icol) = fill_value;
          o_tl_ttend(icol) = fill_value;
          o_qt_ttend(icol) = fill_value;
        } else {
          // Save some outputs just above the "mixed" PBL
          o_pm_hplus(icol) = pm_icol(opt_tm_grad_lev - 1);
          o_tl_hplus(icol) = tl_icol(opt_tm_grad_lev - 1);
          o_qt_hplus(icol) = qt_icol(opt_tm_grad_lev - 1);

          // Save the dF term (F(h) - F(0))
          o_df_inpbl(icol) = (sd_icol(opt_tm_grad_lev - 1) - sd_icol(0)) -
                             (su_icol(opt_tm_grad_lev - 1) - su_icol(0)) +
                             (ld_icol(opt_tm_grad_lev - 1) - ld_icol(0)) -
                             (lu_icol(opt_tm_grad_lev - 1) - lu_icol(0));

          // Now only need to compute below from opt_tm_grad_lev to num_levs
          // Integrate through the PBL, mass-weighted
          // TODO:
          // combine/refactor this once inner parallel_reduce
          // with multiple results/sums is supported...
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(team, opt_tm_grad_lev, num_levs),
              [&](const int &k, Real &result) {
                result += tl_icol(k) * pd_icol(k) / g;
              },
              o_tl_caret(icol));
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(team, opt_tm_grad_lev, num_levs),
              [&](const int &k, Real &result) {
                result += qt_icol(k) * pd_icol(k) / g;
              },
              o_qt_caret(icol));
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(team, opt_tm_grad_lev, num_levs),
              [&](const int &k, Real &result) {
                auto tl_tend = (tl_icol(k) - prev_tliq_icol(k)) / dt;
                result += tl_tend * pd_icol(k) / g;
              },
              o_tl_ttend(icol));
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(team, opt_tm_grad_lev, num_levs),
              [&](const int &k, Real &result) {
                auto qt_tend = (qt_icol(k) - prev_qtot_icol(k)) / dt;
                result += qt_tend * pd_icol(k) / g;
              },
              o_qt_ttend(icol));
        }
        // release stuff from wsm
        ws.release_many_contiguous<wsms>({&tm_grad});
      });
}

}  // namespace scream
