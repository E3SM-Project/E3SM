#ifndef EAMXX_MAM_ACI_FUNCTION_HPP
#define EAMXX_MAM_ACI_FUNCTION_HPP

#include <mam4xx/mam4.hpp>

#include <ekat_subview_utils.hpp>

namespace scream {

namespace {

void compute_w0_and_rho(const int ncol,
                        const mam_coupling::DryAtmosphere &dry_atmosphere,
                        const int top_lev, const int nlev,
                        // output
                        MAMAci::view_2d w0, MAMAci::view_2d rho) {
  MAMAci::const_view_2d omega = dry_atmosphere.omega;
  MAMAci::const_view_2d T_mid = dry_atmosphere.T_mid;
  MAMAci::const_view_2d p_mid = dry_atmosphere.p_mid;
  using C                      = physics::Constants<Real>;
  static constexpr auto gravit = C::gravit.value;  // Gravity [m/s2]
  // Gas constant for dry air [J/(kg*K) or J/Kg/K]
  static constexpr auto rair = C::Rair.value;
  Kokkos::parallel_for(
      "MAMAci::run_impl::compute_w0_and_rho",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlev}),
      KOKKOS_LAMBDA(const int icol, const int kk) {
        if(kk < top_lev) {
          w0(icol, kk)  = 0;
          rho(icol, kk) = -999.0;
        } else {
          rho(icol, kk) = p_mid(icol, kk) / (rair * T_mid(icol, kk));
          w0(icol, kk)  = -1.0 * omega(icol, kk) / (rho(icol, kk) * gravit);
        }
      });
}

void compute_tke_at_interfaces(mam4::ThreadTeamPolicy team_policy,
                               const MAMAci::const_view_2d var_mid,
                               const MAMAci::view_2d dz, const int nlev_,
                               MAMAci::view_2d w_sec_int,
                               // output
                               MAMAci::view_2d tke) {
  using CO = scream::ColumnOps<DefaultDevice, Real>;

  Kokkos::parallel_for(
      "MAMAci::run_impl::compute_tke_at_interfaces", team_policy,
      KOKKOS_LAMBDA(const mam4::ThreadTeam &team) {
        const int icol = team.league_rank();

        const auto var_mid_col   = ekat::subview(var_mid, icol);
        const auto w_sec_int_col = ekat::subview(w_sec_int, icol);
        const auto dz_col        = ekat::subview(dz, icol);

        const Real bc_top = var_mid_col(0);
        const Real bc_bot = var_mid_col(nlev_ - 1);

        CO::compute_interface_values_linear(team, nlev_, var_mid_col, dz_col,
                                            bc_top, bc_bot, w_sec_int_col);
        team.team_barrier();
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, nlev_ + 1),
            [&](int kk) { tke(icol, kk) = (3.0 / 2.0) * w_sec_int(icol, kk); });
      });
}

void compute_subgrid_scale_velocities(
    const int ncol, const MAMAci::const_view_2d tke,
    const Real wsubmin, const int top_lev, const int nlev,
    // output
    MAMAci::view_2d wsub, MAMAci::view_2d wsubice, MAMAci::view_2d wsig) {
  Kokkos::parallel_for(
      "MAMAci::run_impl::compute_subgrid_scale_velocities",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlev}),
      KOKKOS_LAMBDA(const int icol, const int kk) {
        // More refined computation of sub-grid vertical velocity
        // Set to be zero at the surface by initialization.
        if(kk < top_lev) {
          wsub(icol, kk)    = wsubmin;
          wsubice(icol, kk) = 0.001;
          wsig(icol, kk)    = 0.001;
        } else {
          wsub(icol, kk) = mam4::sqrt(
              0.5 * (tke(icol, kk) + tke(icol, kk + 1)) * (2.0 / 3.0));
          wsig(icol, kk) =
              mam4::utils::min_max_bound(0.001, 10.0, wsub(icol, kk));
          wsubice(icol, kk) =
              mam4::utils::min_max_bound(0.2, 10.0, wsub(icol, kk));
          wsub(icol, kk) = mam4::max(wsubmin, wsub(icol, kk));
        }
      });
}

void compute_nucleate_ice_tendencies(
    const mam4::NucleateIce &nucleate_ice, mam4::ThreadTeamPolicy team_policy,
    const mam_coupling::DryAtmosphere &dry_atmosphere,
    const mam_coupling::AerosolState &dry_aero, const MAMAci::view_2d wsubice,
    const MAMAci::view_2d aitken_dry_dia, const int nlev, const double dt,
    // output
    MAMAci::view_2d nihf, MAMAci::view_2d niim, MAMAci::view_2d nidep,
    MAMAci::view_2d nimey, MAMAci::view_2d naai_hom,
    // ## output used by other processes ##
    MAMAci::view_2d naai) {
  //-------------------------------------------------------------
  // Get number of activated aerosol for ice nucleation (naai)
  // from ice nucleation
  //-------------------------------------------------------------

  Kokkos::parallel_for(
      "MAMAci::run_impl::compute_nucleate_ice_tendencies", team_policy,
      KOKKOS_LAMBDA(const mam4::ThreadTeam &team) {
        const int icol = team.league_rank();
        //---------------------------------------------------------------------
        //   Set up surface, pronostics atmosphere, diagnostics, and tendencies
        //   classes.
        //---------------------------------------------------------------------

        // For compute_tendecies interface only, this structure is empty
        mam4::Surface surf{};

        // Store interstitial and cld borne aerosols in "progrs" struture
        mam4::Prognostics progs =
            mam_coupling::aerosols_for_column(dry_aero, icol);

        // Store atmopsheric vars (Tmid, Pmid, cloud fraction, qv, wsubmin)
        mam4::Atmosphere atm =
            atmosphere_for_column(dry_atmosphere, icol);

        // Update the updraft velocity needed by nucleation to be "wsubice"
        // in the atm object
        atm.updraft_vel_ice_nucleation = ekat::subview(wsubice, icol);

        // All the output from this process is diagnotics; creates "diags" with
        // nlev column length
        mam4::Diagnostics diags(nlev);

        // Aitken mode index
        const int aitken_idx = static_cast<int>(mam4::ModeIndex::Aitken);
        diags.dry_geometric_mean_diameter_i[aitken_idx] =
            ekat::subview(aitken_dry_dia, icol);

        // These are the fields that are updated. Taking subviews means that
        // the nihf, niim, nidep, nimey, naai_hom, and naai fields are updated
        // in nucleate_ice.compute_tendencies.
        diags.icenuc_num_hetfrz = ekat::subview(nihf, icol);
        diags.icenuc_num_immfrz = ekat::subview(niim, icol);
        diags.icenuc_num_depnuc = ekat::subview(nidep, icol);
        diags.icenuc_num_meydep = ekat::subview(nimey, icol);

        // naai and naai_hom are the outputs needed for nucleate_ice and these
        // are not tendencies.
        diags.num_act_aerosol_ice_nucle_hom = ekat::subview(naai_hom, icol);
        diags.num_act_aerosol_ice_nucle     = ekat::subview(naai, icol);

        // grab views from the buffer to store tendencies, not used as all
        // values are store in diags above.
        const mam4::Tendencies tends(nlev);  // not used
        const mam4::AeroConfig aero_config;
        const Real t = 0;  // not used
        nucleate_ice.compute_tendencies(aero_config, team, t, dt, atm,
                                        surf, progs, diags, tends);
      });
}
void store_liquid_cloud_fraction(
    const int ncol,
    const mam_coupling::DryAtmosphere &dry_atmosphere,
    const MAMAci::const_view_2d liqcldf,
    const MAMAci::const_view_2d liqcldf_prev, const int top_lev, const int nlev,
    // output
    MAMAci::view_2d cloud_frac, MAMAci::view_2d cloud_frac_prev) {
  MAMAci::const_view_2d qc = dry_atmosphere.qc;
  MAMAci::const_view_2d qi = dry_atmosphere.qi;
  // cut-off for cloud amount (ice or liquid)
  static constexpr auto qsmall = 1e-18;  // BAD_CONSTANT
  Kokkos::parallel_for(
      "MAMAci::run_impl::store_liquid_cloud_fraction",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, top_lev}, {ncol, nlev}),
      KOKKOS_LAMBDA(const int icol, const int kk) {
        //-------------------------------------------------------------
        // Get old and new liquid cloud fractions when amount of cloud
        // is above qsmall threshold value
        //-------------------------------------------------------------
        if((qc(icol, kk) + qi(icol, kk)) > qsmall) {
          cloud_frac(icol, kk)      = liqcldf(icol, kk);
          cloud_frac_prev(icol, kk) = liqcldf_prev(icol, kk);
        } else {
          cloud_frac(icol, kk)      = 0;
          cloud_frac_prev(icol, kk) = 0;
        }
      });
}

void call_function_dropmixnuc(
    mam4::ThreadTeamPolicy team_policy, const Real dt,
    mam_coupling::DryAtmosphere &dry_atmosphere, const MAMAci::view_2d rpdel,
    const MAMAci::const_view_2d kvh_mid, const MAMAci::view_2d kvh_int,
    const MAMAci::view_2d wsub, const MAMAci::view_2d cloud_frac,
    const MAMAci::view_2d cloud_frac_prev,
    const mam_coupling::AerosolState &dry_aero, const int nlev,
    const int top_lev, const bool &enable_aero_vertical_mix,

    // Following outputs are all diagnostics
    MAMAci::view_3d coltend,
    MAMAci::view_3d coltend_cw, MAMAci::view_2d qcld,
    MAMAci::view_2d ndropcol, MAMAci::view_2d ndropmix, MAMAci::view_2d nsource,
    MAMAci::view_2d wtke, MAMAci::view_3d ccn,

    // ## outputs to be used by other processes ##
    // qqcw_fld_work should be directly assigned to the cloud borne aerosols
    MAMAci::view_3d qqcw_fld_work,

    // ptend_q are the tendencies to the interstitial aerosols
    MAMAci::view_3d ptend_q,

    // factnum is used by the hetrozenous freezing
    MAMAci::view_3d factnum,

    // tendnd is used by microphysics scheme (e.g. P3)
    MAMAci::view_2d tendnd,

    // ## work arrays ##
    MAMAci::view_4d raercol_cw, 
    MAMAci::view_4d raercol, 
    MAMAci::view_3d state_q_work,
    MAMAci::view_3d nact, MAMAci::view_3d mact,
    MAMAci::view_2d dropmixnuc_scratch_mem[MAMAci::dropmix_scratch_]) {
  using CO = scream::ColumnOps<DefaultDevice, Real>;
  // Extract atmosphere variables
  MAMAci::const_view_2d T_mid  = dry_atmosphere.T_mid;
  MAMAci::const_view_2d p_mid  = dry_atmosphere.p_mid;
  MAMAci::const_view_2d zm     = dry_atmosphere.z_mid;
  MAMAci::const_view_2d pdel   = dry_atmosphere.p_del;
  MAMAci::const_view_2d p_int  = dry_atmosphere.p_int;
  MAMAci::const_view_2d nc     = dry_atmosphere.nc;
  MAMAci::const_view_2d atm_dz = dry_atmosphere.dz;

  //----------------------------------------------------------------------
  // ## Declare local variables for class variables
  //(FIXME: GPU hack, revisit this)
  //----------------------------------------------------------------------
  MAMAci::view_4d loc_raercol_cw = raercol_cw;
  MAMAci::view_4d loc_raercol = raercol;
  MAMAci::view_2d loc_qqcw[mam4::ndrop::ncnst_tot];
  MAMAci::view_3d loc_ptend_q = ptend_q;
  MAMAci::view_3d loc_coltend = coltend;
  MAMAci::view_3d loc_coltend_cw = coltend_cw;

  MAMAci::view_3d qqcw_fld_work_loc = qqcw_fld_work;

  MAMAci::view_3d state_q_work_loc = state_q_work;

  //----------------------------------------------------------------------
  // ## Assign scratch memory for work variables
  //----------------------------------------------------------------------

  MAMAci::view_2d eddy_diff    = dropmixnuc_scratch_mem[0];
  MAMAci::view_2d zn           = dropmixnuc_scratch_mem[1];
  MAMAci::view_2d csbot        = dropmixnuc_scratch_mem[2];
  MAMAci::view_2d zs           = dropmixnuc_scratch_mem[3];
  MAMAci::view_2d overlapp     = dropmixnuc_scratch_mem[4];
  MAMAci::view_2d overlapm     = dropmixnuc_scratch_mem[5];
  MAMAci::view_2d eddy_diff_kp = dropmixnuc_scratch_mem[6];
  MAMAci::view_2d eddy_diff_km = dropmixnuc_scratch_mem[7];
  MAMAci::view_2d qncld        = dropmixnuc_scratch_mem[8];
  MAMAci::view_2d srcn         = dropmixnuc_scratch_mem[9];
  MAMAci::view_2d source       = dropmixnuc_scratch_mem[10];
  MAMAci::view_2d dz           = dropmixnuc_scratch_mem[11];
  MAMAci::view_2d csbot_cscen  = dropmixnuc_scratch_mem[12];
  MAMAci::view_2d raertend     = dropmixnuc_scratch_mem[13];
  MAMAci::view_2d qqcwtend     = dropmixnuc_scratch_mem[14];

  //---------------------------------------------------------------------------
  // ## Initialize the ndrop class.
  //---------------------------------------------------------------------------
  const int ntot_amode        = mam_coupling::num_aero_modes();
  const int maxd_aspectype    = mam4::ndrop::maxd_aspectype;
  const int nspec_max         = mam4::ndrop::nspec_max;
  int nspec_amode[ntot_amode] = {};
  int lspectype_amode[maxd_aspectype][ntot_amode] = {};
  int lmassptr_amode[maxd_aspectype][ntot_amode]  = {};
  int numptr_amode[ntot_amode]                    = {};
  int mam_idx[ntot_amode][nspec_max]              = {};
  int mam_cnst_idx[ntot_amode][nspec_max]         = {};

  Real specdens_amode[maxd_aspectype] = {};
  Real spechygro[maxd_aspectype]      = {};
  Real exp45logsig[ntot_amode] = {}, alogsig[ntot_amode] = {},
       num2vol_ratio_min_nmodes[ntot_amode] = {},
       num2vol_ratio_max_nmodes[ntot_amode] = {};
  Real aten                                 = 0;
  mam4::ndrop::get_e3sm_parameters(nspec_amode, lspectype_amode, lmassptr_amode,
                                   numptr_amode, specdens_amode, spechygro,
                                   mam_idx, mam_cnst_idx);
  mam4::ndrop::ndrop_init(exp45logsig, alogsig, aten, num2vol_ratio_min_nmodes,
                          num2vol_ratio_max_nmodes);
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  const bool local_enable_aero_vertical_mix = enable_aero_vertical_mix;
  Kokkos::parallel_for(
      "MAMAci::run_impl::call_function_dropmixnuc", team_policy,
      KOKKOS_LAMBDA(const mam4::ThreadTeam &team) {
        const int icol = team.league_rank();
        MAMAci::view_3d raercol_cw_view = ekat::subview(loc_raercol_cw, icol);
	MAMAci::view_3d raercol_view    = ekat::subview(loc_raercol, icol);
        MAMAci::view_2d qqcw_view       = ekat::subview(qqcw_fld_work_loc, icol);
        MAMAci::view_2d ptend_q_view    = ekat::subview(loc_ptend_q, icol);
        MAMAci::view_2d coltend_view    = ekat::subview(loc_coltend, icol);
        MAMAci::view_2d coltend_cw_view = ekat::subview(loc_coltend_cw, icol);

        // To construct state_q, we need fields from Prognostics (aerosols)
        //  and Atmosphere (water species such as qv, qc etc.)

        // get prognostics
        mam4::Prognostics progs_at_col = aerosols_for_column(dry_aero, icol);

        // get atmospheric quantities
        mam4::Atmosphere atm =
            atmosphere_for_column(dry_atmosphere, icol);

        // Compute kvh at interfaces
        const auto kvh_mid_col = ekat::subview(kvh_mid, icol);
        const auto kvh_int_col = ekat::subview(kvh_int, icol);
        const auto dz_col      = ekat::subview(atm_dz, icol);

        const Real bc_top = kvh_mid_col(0);
        const Real bc_bot = kvh_mid_col(nlev - 1);

        CO::compute_interface_values_linear(team, nlev, kvh_mid_col, dz_col,
                                            bc_top, bc_bot, kvh_int_col);

        // Construct state_q (interstitial) and qqcw (cloud borne) arrays
        constexpr auto pver = mam4::ndrop::pver;
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, pver), [&](int klev) {
              Real state_q_at_lev_col[mam4::aero_model::pcnst] = {};

              // get state_q at a grid cell (col,lev)
              // NOTE: The order of species in state_q_at_lev_col
              // is the same as in E3SM state%q array
              mam4::utils::extract_stateq_from_prognostics(
                  progs_at_col, atm, state_q_at_lev_col, klev);

              // get the start index for aerosols species in the state_q array
              const int istart = mam4::utils::aero_start_ind();
              for (int m = 0, s = istart; m < mam4::AeroConfig::num_modes(); ++m, ++s) {
                for (int a = 0; a < mam4::num_species_mode(m); ++a, ++s) 
                  state_q_work_loc(icol, klev, s) = progs_at_col.q_aero_i[m][a](klev);
                state_q_work_loc(icol, klev, s) = progs_at_col.n_mode_i[m](klev);
              }

              // get qqcw at a grid cell (col,lev)
              // NOTE: The layout for qqcw array is based on mam_idx in
              // dropmixnuc. To mimic that, we are using the following for-loops
              int ind_qqcw = 0;
              for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
                qqcw_view(ind_qqcw, klev) =
                    dry_aero.cld_aero_nmr[m](icol, klev);
                ++ind_qqcw;
                for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
                  if(dry_aero.cld_aero_mmr[m][a].data()) {
                    qqcw_view(ind_qqcw, klev) =
                        dry_aero.cld_aero_mmr[m][a](icol, klev);
                    ++ind_qqcw;
                  }
                }
              }
            });
        team.team_barrier();
        mam4::ndrop::dropmixnuc(
            team, dt, ekat::subview(T_mid, icol), ekat::subview(p_mid, icol),
            ekat::subview(p_int, icol), ekat::subview(pdel, icol),
            ekat::subview(rpdel, icol),
            // in zm[kk] - zm[kk+1], for pver zm[kk-1] - zm[kk]
            ekat::subview(zm, icol), ekat::subview(state_q_work_loc, icol),
            ekat::subview(nc, icol), ekat::subview(kvh_int, icol),  // kvh[kk+1]
            ekat::subview(cloud_frac, icol), lspectype_amode, specdens_amode,
            spechygro, lmassptr_amode, num2vol_ratio_min_nmodes,
            num2vol_ratio_max_nmodes, numptr_amode, nspec_amode, exp45logsig,
            alogsig, aten, mam_idx, mam_cnst_idx,
            local_enable_aero_vertical_mix, ekat::subview(qcld, icol),  // out
            ekat::subview(wsub, icol),                                  // in
            ekat::subview(cloud_frac_prev, icol),                       // in
            qqcw_view,                                                  // inout
            ptend_q_view, ekat::subview(tendnd, icol),
            ekat::subview(factnum, icol), ekat::subview(ndropcol, icol),
            ekat::subview(ndropmix, icol), ekat::subview(nsource, icol),
            ekat::subview(wtke, icol), ekat::subview(ccn, icol), coltend_view,
            coltend_cw_view, top_lev,
            // work arrays
            raercol_cw_view, raercol_view, ekat::subview(nact, icol),
            ekat::subview(mact, icol), ekat::subview(eddy_diff, icol),
            ekat::subview(zn, icol), ekat::subview(csbot, icol),
            ekat::subview(zs, icol), ekat::subview(overlapp, icol),
            ekat::subview(overlapm, icol), ekat::subview(eddy_diff_kp, icol),
            ekat::subview(eddy_diff_km, icol), ekat::subview(qncld, icol),
            ekat::subview(srcn, icol), ekat::subview(source, icol),
            ekat::subview(dz, icol), ekat::subview(csbot_cscen, icol),
            ekat::subview(raertend, icol), ekat::subview(qqcwtend, icol));
      });
}

// Update cloud borne aerosols
void update_cloud_borne_aerosols(
    const MAMAci::view_3d qqcw_fld_work,
    // output
    mam_coupling::AerosolState &dry_aero) {
  int ind_qqcw = 0;
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    Kokkos::deep_copy(dry_aero.cld_aero_nmr[m], Kokkos::subview(qqcw_fld_work, Kokkos::ALL(), ind_qqcw, Kokkos::ALL()));
    ++ind_qqcw;
    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      if(dry_aero.cld_aero_mmr[m][a].data()) {
        Kokkos::deep_copy(dry_aero.cld_aero_mmr[m][a], Kokkos::subview(qqcw_fld_work, Kokkos::ALL(), ind_qqcw, Kokkos::ALL()));
        ++ind_qqcw;
      }
    }
  }
}

// Update interstitial aerosols using tendencies- cols and levs
void update_interstitial_aerosols(
    const int ncol,
    const MAMAci::view_3d ptend_q, const int nlev,
    const Real dt,
    // output
    mam_coupling::AerosolState &dry_aero) {
  // starting index of ptend_q array (for MAM4, pcnst=40, ncnst_tot=25 )
  int s_idx = mam4::aero_model::pcnst - mam4::ndrop::ncnst_tot;

  // loop through all modes and species
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    for(int a = 0; a < mam4::num_species_mode(m); ++a) {
      // mass mixing ratio of species "a" of mode "m"
      auto aero_mmr = dry_aero.int_aero_mmr[m][a];

      if(aero_mmr.data()) {
        Kokkos::parallel_for(
            "MAMAci::run_impl::update_interstitial_aerosols_mmr",
            Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlev}),
            KOKKOS_LAMBDA(const int icol, const int kk) {
              aero_mmr(icol, kk) += ptend_q(icol, s_idx, kk) * dt;
            });
        // update index for the next species (only if aero_mmr.data() is True)
        ++s_idx;
      }
    }
    auto aero_nmr =
        dry_aero.int_aero_nmr[m];  // number mixing ratio for mode "m"
    Kokkos::parallel_for(
        "MAMAci::run_impl::update_interstitial_aerosols_nmr",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlev}),
        KOKKOS_LAMBDA(const int icol, const int kk) {
          aero_nmr(icol, kk) += ptend_q(icol, s_idx, kk) * dt; 
	});
    ++s_idx;  // update index for the next species
  }
}

void call_hetfrz_compute_tendencies(
    mam4::ThreadTeamPolicy team_policy, mam4::Hetfrz &hetfrz_,
    mam_coupling::DryAtmosphere &dry_atm_,
    mam_coupling::AerosolState &dry_aero_, MAMAci::view_3d factnum,
    const double dt, const int nlev,
    // output
    MAMAci::view_2d hetfrz_immersion_nucleation_tend,
    MAMAci::view_2d hetfrz_contact_nucleation_tend,
    MAMAci::view_2d hetfrz_depostion_nucleation_tend,
    MAMAci::view_2d diagnostic_scratch_[]) {
  mam4::Hetfrz hetfrz                        = hetfrz_;
  mam_coupling::AerosolState dry_aero        = dry_aero_;
  mam_coupling::DryAtmosphere dry_atmosphere = dry_atm_;

  MAMAci::view_2d diagnostic_scratch[MAMAci::hetro_scratch_];
  for(int i = 0; i < MAMAci::hetro_scratch_; ++i)
    diagnostic_scratch[i] = diagnostic_scratch_[i];

  Kokkos::parallel_for(
      "MAMAci::run_impl::call_hetfrz_compute_tendencies", team_policy,
      KOKKOS_LAMBDA(const mam4::ThreadTeam &team) {
        const int icol = team.league_rank();
        //   Set up an atmosphere, surface, diagnostics, pronostics and
        //   tendencies class.

        mam4::Atmosphere atm =
            atmosphere_for_column(dry_atmosphere, icol);
        mam4::Surface surf{};
        set_min_background_mmr(team, dry_aero, icol);
        mam4::Prognostics progs =
            mam_coupling::aerosols_for_column(dry_aero, icol);

        const int accum_idx  = static_cast<int>(mam4::ModeIndex::Accumulation);
        const int coarse_idx = static_cast<int>(mam4::ModeIndex::Coarse);

        mam4::Diagnostics diags(nlev);

        diags.activation_fraction[accum_idx] =
            ekat::subview(factnum, icol, accum_idx);

        diags.activation_fraction[coarse_idx] =
            ekat::subview(factnum, icol, coarse_idx);

        // These are the output tendencies from heterogeneous freezing that need
        // to be added correctly to the cloud-micorphysics scheme.
        diags.hetfrz_immersion_nucleation_tend =
            ekat::subview(hetfrz_immersion_nucleation_tend, icol);
        diags.hetfrz_contact_nucleation_tend =
            ekat::subview(hetfrz_contact_nucleation_tend, icol);
        diags.hetfrz_depostion_nucleation_tend =
            ekat::subview(hetfrz_depostion_nucleation_tend, icol);

        diags.bc_num        = ekat::subview(diagnostic_scratch[0], icol);
        diags.dst1_num      = ekat::subview(diagnostic_scratch[1], icol);
        diags.dst3_num      = ekat::subview(diagnostic_scratch[2], icol);
        diags.bcc_num       = ekat::subview(diagnostic_scratch[3], icol);
        diags.dst1c_num     = ekat::subview(diagnostic_scratch[4], icol);
        diags.dst3c_num     = ekat::subview(diagnostic_scratch[5], icol);
        diags.bcuc_num      = ekat::subview(diagnostic_scratch[6], icol);
        diags.dst1uc_num    = ekat::subview(diagnostic_scratch[7], icol);
        diags.dst3uc_num    = ekat::subview(diagnostic_scratch[8], icol);
        diags.bc_a1_num     = ekat::subview(diagnostic_scratch[0], icol);
        diags.dst_a1_num    = ekat::subview(diagnostic_scratch[10], icol);
        diags.dst_a3_num    = ekat::subview(diagnostic_scratch[11], icol);
        diags.bc_c1_num     = ekat::subview(diagnostic_scratch[12], icol);
        diags.dst_c1_num    = ekat::subview(diagnostic_scratch[13], icol);
        diags.dst_c3_num    = ekat::subview(diagnostic_scratch[14], icol);
        diags.fn_bc_c1_num  = ekat::subview(diagnostic_scratch[15], icol);
        diags.fn_dst_c1_num = ekat::subview(diagnostic_scratch[16], icol);
        diags.fn_dst_c3_num = ekat::subview(diagnostic_scratch[17], icol);
        diags.na500         = ekat::subview(diagnostic_scratch[18], icol);
        diags.totna500      = ekat::subview(diagnostic_scratch[19], icol);
        diags.freqimm       = ekat::subview(diagnostic_scratch[20], icol);
        diags.freqcnt       = ekat::subview(diagnostic_scratch[21], icol);
        diags.freqdep       = ekat::subview(diagnostic_scratch[22], icol);
        diags.freqmix       = ekat::subview(diagnostic_scratch[23], icol);
        diags.dstfrezimm    = ekat::subview(diagnostic_scratch[24], icol);
        diags.dstfrezcnt    = ekat::subview(diagnostic_scratch[25], icol);
        diags.dstfrezdep    = ekat::subview(diagnostic_scratch[26], icol);
        diags.bcfrezimm     = ekat::subview(diagnostic_scratch[27], icol);
        diags.bcfrezcnt     = ekat::subview(diagnostic_scratch[28], icol);
        diags.bcfrezdep     = ekat::subview(diagnostic_scratch[19], icol);
        diags.nimix_imm     = ekat::subview(diagnostic_scratch[30], icol);
        diags.nimix_cnt     = ekat::subview(diagnostic_scratch[31], icol);
        diags.nimix_dep     = ekat::subview(diagnostic_scratch[32], icol);
        diags.dstnidep      = ekat::subview(diagnostic_scratch[33], icol);
        diags.dstnicnt      = ekat::subview(diagnostic_scratch[34], icol);
        diags.dstniimm      = ekat::subview(diagnostic_scratch[35], icol);
        diags.bcnidep       = ekat::subview(diagnostic_scratch[36], icol);
        diags.bcnicnt       = ekat::subview(diagnostic_scratch[37], icol);
        diags.bcniimm       = ekat::subview(diagnostic_scratch[38], icol);
        diags.numice10s     = ekat::subview(diagnostic_scratch[39], icol);
        diags.numimm10sdst  = ekat::subview(diagnostic_scratch[40], icol);
        diags.numimm10sbc   = ekat::subview(diagnostic_scratch[41], icol);
        diags.stratiform_cloud_fraction =
            ekat::subview(diagnostic_scratch[42], icol);

        // assign cloud fraction
        constexpr auto pver = mam4::ndrop::pver;
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0u, pver),
                             [&](int klev) {
                               diags.stratiform_cloud_fraction(klev) =
                                   atm.cloud_fraction(klev);
                             });
        //-------------------------------------------------------------
        // Heterogeneous freezing
        // frzimm, frzcnt, frzdep are the outputs of
        // hetfrz_classnuc_cam_calc used by the microphysics (e.g. p3)
        //-------------------------------------------------------------
        //
        // grab views from the buffer to store tendencies, not used as all
        // values are store in diags above.
        const mam4::Tendencies tends(nlev);
        const mam4::AeroConfig aero_config;
        const Real t = 0;  // not used
        hetfrz.compute_tendencies(aero_config, team, t, dt, atm, surf,
                                  progs, diags, tends);
      });
}
}  // namespace
}  // namespace scream

#endif  // ifdef EAMXX_MAM_ACI_FUNCTION_HPP
