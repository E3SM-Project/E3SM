#ifndef ZM_ZM_CONV_MAIN_IMPL_HPP
#define ZM_ZM_CONV_MAIN_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU
#include <ekat_subview_utils.hpp>

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_conv_main. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 *
 * zm_conv_main is a host function that launches Kokkos kernels over columns.
 * Temporary per-column arrays are allocated as 2D device views; per-column
 * scalars are stored as 1D device views.  Sub-functions (compute_dilute_cape,
 * zm_cloud_properties, zm_closure, zm_calc_output_tend) are KOKKOS_FUNCTIONs
 * called from within team-policy kernels, sharing a single WorkspaceManager.
 */

template<typename S, typename D>
typename Functions<S,D>::view_1d<bool> Functions<S,D>::zm_conv_main(
  // Inputs
  const ZmRuntimeOpt& runtime_opt,
  const Int& ncol,
  const Int& pver,
  const Int& pverp,
  const bool& is_first_step,
  const Real& time_step,
  const uview_2d<const Real>& t_mid,
  const uview_2d<const Real>& q_mid_in,
  const uview_2d<const Real>& omega,
  const uview_2d<const Real>& p_mid_in,
  const uview_2d<const Real>& p_int_in,
  const uview_2d<const Real>& p_del_in,
  const uview_1d<const Real>& geos,
  const uview_2d<const Real>& z_mid_in,
  const uview_2d<const Real>& z_int_in,
  const uview_1d<const Real>& pbl_hgt,
  const uview_1d<const Real>& tpert,
  const uview_1d<const Real>& landfrac,
  const uview_2d<const Real>& t_star,
  const uview_2d<const Real>& q_star,
  // Outputs
  const uview_1d<Int>& msemax_klev,
  const uview_1d<Int>& jctop,
  const uview_1d<Int>& jcbot,
  const uview_1d<Int>& jt,
  const uview_1d<Real>& prec,
  const uview_2d<Real>& heat,
  const uview_2d<Real>& qtnd,
  const uview_1d<Real>& cape,
  const uview_1d<Real>& dcape,
  const uview_2d<Real>& mcon,
  const uview_2d<Real>& pflx,
  const uview_2d<Real>& zdu,
  const uview_2d<Real>& mflx_up,
  const uview_2d<Real>& entr_up,
  const uview_2d<Real>& detr_up,
  const uview_2d<Real>& mflx_dn,
  const uview_2d<Real>& entr_dn,
  const uview_2d<Real>& p_del,
  const uview_1d<Real>& dsubcld,
  const uview_2d<Real>& ql,
  const uview_1d<Real>& rliq,
  const uview_2d<Real>& rprd,
  const uview_2d<Real>& dlf)
{
  //----------------------------------------------------------------------------
  // Purpose: Main driver for Zhang-McFarlane convection scheme
  //----------------------------------------------------------------------------
  using ExeSpace   = typename KT::ExeSpace;
  using RangePolicy = Kokkos::RangePolicy<ExeSpace>;

  constexpr Real gravit = PC::gravit.value;
  constexpr Real cpair  = PC::Cpair.value;

  const Int msg = runtime_opt.limcnv - 1;

  //----------------------------------------------------------------------------
  // Allocate temporary 2D device views  [ncol, pver] or [ncol, pverp]
  //----------------------------------------------------------------------------
  view_2d<Real>
    s_mid    ("s_mid",     ncol, pver),
    q_mid    ("q_mid",     ncol, pver),
    p_mid    ("p_mid",     ncol, pver),
    z_mid    ("z_mid",     ncol, pver),
    s_int    ("s_int",     ncol, pver),
    q_int    ("q_int",     ncol, pver),
    t_pcl    ("t_pcl",     ncol, pver),
    q_pcl_sat("q_pcl_sat", ncol, pver),
    t_pcl_m1 ("t_pcl_m1",  ncol, pver),
    q_pcl_sat_m1("q_pcl_sat_m1", ncol, pver),
    dqdt     ("dqdt",      ncol, pver),
    dsdt     ("dsdt",      ncol, pver),
    mflx_net ("mflx_net",  ncol, pver),
    q_upd    ("q_upd",     ncol, pver),
    s_upd    ("s_upd",     ncol, pver),
    q_dnd    ("q_dnd",     ncol, pver),
    s_dnd    ("s_dnd",     ncol, pver),
    qst      ("qst",       ncol, pver),
    cu       ("cu",        ncol, pver),
    evp      ("evp",       ncol, pver),
    dl       ("dl",        ncol, pver);
  view_2d<Real>
    p_int    ("p_int",     ncol, pverp),
    z_int    ("z_int",     ncol, pverp);

  // Per-column scalar temporaries [ncol]
  view_1d<Real>
    z_srf          ("z_srf",           ncol),
    t_pcl_lcl      ("t_pcl_lcl",       ncol),
    t_pcl_lcl_m1   ("t_pcl_lcl_m1",    ncol),
    cape_m1        ("cape_m1",          ncol),
    cld_base_mass_flux("cld_base_mass_flux", ncol),
    q_mx           ("q_mx",            ncol),
    t_mx           ("t_mx",            ncol),
    q_mx_m1        ("q_mx_m1",         ncol),
    t_mx_m1        ("t_mx_m1",         ncol);
  view_1d<Int>
    pbl_top        ("pbl_top",         ncol),
    lcl            ("lcl",             ncol),
    lel            ("lel",             ncol),
    lcl_m1         ("lcl_m1",          ncol),
    lel_m1         ("lel_m1",          ncol),
    msemax_klev_m1 ("msemax_klev_m1",  ncol),
    jlcl           ("jlcl",            ncol),
    j0             ("j0",              ncol),
    jd             ("jd",              ncol);
  view_1d<bool> active("active", ncol);

  // Workspace for sub-functions (max 20 arrays for zm_cloud_properties + entrainment)
  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(ncol, pver);
  WorkspaceManager wsm(pverp, 20, policy);

  //============================================================================
  // Kernel 1: Initialize outputs/scalars, convert pressures, find PBL top,
  //           set up s_mid / q_mid and initial interface values.
  //============================================================================
  Kokkos::parallel_for("zm_conv_main_init", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Initialize per-column output scalars
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      z_srf(i)               = geos(i) / gravit;
      pbl_top(i)             = pver - 1;
      lcl(i)                 = 0;
      lel(i)                 = pver - 1;
      lcl_m1(i)              = 0;
      lel_m1(i)              = pver - 1;
      msemax_klev_m1(i)      = 0;
      cape_m1(i)             = 0;
      t_pcl_lcl(i)           = 400;
      t_pcl_lcl_m1(i)        = 400;
      cld_base_mass_flux(i)  = 0;
      prec(i)                = 0;
      rliq(i)                = 0;
      dsubcld(i)             = 0;
      jctop(i)               = pver - 1;
      jcbot(i)               = 0;
      cape(i)                = 0;
      dcape(i)               = 0;
      msemax_klev(i)         = 0;
      pflx(i, pver)          = 0;
    });

    // Initialize output arrays
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pver), [&](const Int k) {
      heat(i,k)      = 0;
      qtnd(i,k)      = 0;
      mcon(i,k)      = 0;
      pflx(i,k)      = 0;
      rprd(i,k)      = 0;
      zdu(i,k)       = 0;
      ql(i,k)        = 0;
      dlf(i,k)       = 0;
      p_del(i,k)     = 0;
      mflx_up(i,k)   = 0;
      entr_up(i,k)   = 0;
      detr_up(i,k)   = 0;
      mflx_dn(i,k)   = 0;
      entr_dn(i,k)   = 0;
      t_pcl(i,k)     = 0;
    });
    team.team_barrier();

    // Convert pressures from Pa to mb, heights relative to surface
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pver), [&](const Int k) {
      p_mid(i,k) = p_mid_in(i,k) * ZMC::pa_to_mb;
      z_mid(i,k) = z_mid_in(i,k) + z_srf(i);
      p_int(i,k) = p_int_in(i,k) * ZMC::pa_to_mb;
      z_int(i,k) = z_int_in(i,k) + z_srf(i);
      if (k == pver-1) {
        p_int(i,k+1) = p_int_in(i,k+1) * ZMC::pa_to_mb;
        z_int(i,k+1) = z_int_in(i,k+1) + z_srf(i);
      }
    });
    team.team_barrier();

    // Find PBL top: topmost level where z_mid is within half a layer of pbl_hgt
    Int pbl_top_result = pver - 1;
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, msg + 1, pver - 1),
      [&](const Int k, Int& val) {
        if (std::abs(z_mid(i,k) - pbl_hgt(i)) < (z_int(i,k) - z_int(i,k+1)) * ZMC::half) {
          val = ekat::impl::min(val, k);
        }
      }, Kokkos::Min<Int>(pbl_top_result));
    team.team_barrier();
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      pbl_top(i) = ekat::impl::min(pbl_top_result, pver - 1);
    });
    team.team_barrier();

    // Store input humidity, compute dry static energy, initialize interface values
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pver), [&](const Int k) {
      q_mid(i,k) = q_mid_in(i,k);
      s_mid(i,k) = t_mid(i,k) + (gravit / cpair) * z_mid(i,k);
      s_int(i,k) = s_mid(i,k);
      q_int(i,k) = q_mid(i,k);
    });
  });
  Kokkos::fence();

  //============================================================================
  // Kernel 2: Compute CAPE (standard parcel using current-step state)
  //============================================================================
  Kokkos::parallel_for("zm_conv_main_cape", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    auto ws = wsm.get_workspace(team);
    compute_dilute_cape(team, ws, runtime_opt,
                        pver, pverp,
                        runtime_opt.num_cin, msg,
                        ekat::subview(q_mid, i), ekat::subview(t_mid, i),
                        ekat::subview(z_mid, i), ekat::subview(p_mid, i),
                        ekat::subview(p_int, i),
                        pbl_top(i), tpert(i),
                        true,   // calc_msemax_klev
                        0,      // prev_msemax_klev (unused for first call)
                        false,  // use_input_tq_mx
                        ekat::subview(q_pcl_sat, i),
                        msemax_klev(i), lcl(i), lel(i), cape(i),
                        q_mx(i), t_mx(i),
                        ekat::subview(t_pcl, i), t_pcl_lcl(i));
  });
  Kokkos::fence();

  //============================================================================
  // Kernel 3: Compute DCAPE (using previous-step state) — only when needed
  //============================================================================
  if (!is_first_step && runtime_opt.trig_dcape) {
    Kokkos::parallel_for("zm_conv_main_dcape", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const Int i = team.league_rank();
      const Int prev_msemax_klev_val = msemax_klev(i);
      auto ws = wsm.get_workspace(team);
      compute_dilute_cape(team, ws, runtime_opt,
                          pver, pverp,
                          runtime_opt.num_cin, msg,
                          ekat::subview(q_star, i), ekat::subview(t_star, i),
                          ekat::subview(z_mid, i), ekat::subview(p_mid, i),
                          ekat::subview(p_int, i),
                          pbl_top(i), tpert(i),
                          false, // calc_msemax_klev: use prev level
                          prev_msemax_klev_val,
                          false, // use_input_tq_mx
                          ekat::subview(q_pcl_sat_m1, i),
                          msemax_klev_m1(i), lcl_m1(i), lel_m1(i), cape_m1(i),
                          q_mx_m1(i), t_mx_m1(i),
                          ekat::subview(t_pcl_m1, i), t_pcl_lcl_m1(i));
    });
    Kokkos::fence();

    Kokkos::parallel_for("zm_conv_main_dcape_calc", RangePolicy(0, ncol),
      KOKKOS_LAMBDA(const Int i) {
        dcape(i) = (cape(i) - cape_m1(i)) / time_step;
      });
    Kokkos::fence();
  }

  //============================================================================
  // Host: Determine active columns
  // cape_threshold_loc depends only on runtime_opt / is_first_step (same for all cols)
  //============================================================================
  const Real cape_threshold_loc = (runtime_opt.trig_dcape && !is_first_step)
                                    ? ZMC::cape_threshold_new
                                    : ZMC::cape_threshold_old;
  const bool use_dcape_trigger = runtime_opt.trig_dcape && !is_first_step;
  int inactive_cnt = 0;
  Kokkos::parallel_reduce("zm_conv_main_active", RangePolicy(0, ncol),
                          KOKKOS_LAMBDA(const Int i, Int& local_inactive) {
      active(i) = use_dcape_trigger
        ? (cape(i) > cape_threshold_loc && dcape(i) > ZMC::dcape_threshold)
        : (cape(i) > cape_threshold_loc);
      if (!active(i)) {
        local_inactive+=1;
      }
  }, inactive_cnt);
  // The inactive count may be worth looking at as a high inactive count could indicate
  // an issue with randomly generated input data for testing.
  // std::cout << "INACTIVE COUNT: " << inactive_cnt << std::endl;
  Kokkos::fence();

  //============================================================================
  // Kernel 4: Active columns — convert p_del, compute dsubcld, define s/q interfaces
  //============================================================================
  Kokkos::parallel_for("zm_conv_main_setup_active", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    if (!active(i)) return;

    // Convert layer thickness to mb
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pver), [&](const Int k) {
      p_del(i,k) = p_del_in(i,k) * ZMC::pa_to_mb;
    });
    team.team_barrier();

    // Sub-cloud layer pressure thickness
    Int max_idx = ekat::impl::max(msg, msemax_klev(i));
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, max_idx, pver),
      [&](const Int k, Real& sum) { sum += p_del(i,k); }, dsubcld(i));
    team.team_barrier();

    // Define interfacial s/q values (log-linear interpolation where possible)
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg + 1, pver), [&](const Int k) {
      Real sdifr = 0, qdifr = 0;
      if (s_mid(i,k) > 0 || s_mid(i,k-1) > 0) {
        sdifr = std::abs((s_mid(i,k) - s_mid(i,k-1)) /
                         ekat::impl::max(s_mid(i,k-1), s_mid(i,k)));
      }
      if (q_mid(i,k) > 0 || q_mid(i,k-1) > 0) {
        qdifr = std::abs((q_mid(i,k) - q_mid(i,k-1)) /
                         ekat::impl::max(q_mid(i,k-1), q_mid(i,k)));
      }
      if (sdifr > ZMC::interp_diff_min) {
        s_int(i,k) = std::log(s_mid(i,k-1) / s_mid(i,k)) *
                     s_mid(i,k-1) * s_mid(i,k) / (s_mid(i,k-1) - s_mid(i,k));
      } else {
        s_int(i,k) = ZMC::half * (s_mid(i,k) + s_mid(i,k-1));
      }
      if (qdifr > ZMC::interp_diff_min) {
        q_int(i,k) = std::log(q_mid(i,k-1) / q_mid(i,k)) *
                     q_mid(i,k-1) * q_mid(i,k) / (q_mid(i,k-1) - q_mid(i,k));
      } else {
        q_int(i,k) = ZMC::half * (q_mid(i,k) + q_mid(i,k-1));
      }
    });
  });
  Kokkos::fence();

  //============================================================================
  // Kernel 5: Updraft/downdraft cloud properties
  //============================================================================
  Kokkos::parallel_for("zm_conv_main_cloud_props", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    if (!active(i)) return;
    auto ws = wsm.get_workspace(team);
    zm_cloud_properties(team, ws, runtime_opt,
                        pver, pverp, msg, runtime_opt.limcnv,
                        ekat::subview(p_mid, i), ekat::subview(z_mid, i),
                        ekat::subview(z_int, i),
                        ekat::subview(t_mid, i), ekat::subview(s_mid, i),
                        ekat::subview(s_int, i), ekat::subview(q_mid, i),
                        landfrac(i), tpert(i),
                        msemax_klev(i), lel(i),
                        jt(i), jlcl(i), j0(i), jd(i),
                        ekat::subview(mflx_up, i), ekat::subview(entr_up, i),
                        ekat::subview(detr_up, i), ekat::subview(mflx_dn, i),
                        ekat::subview(entr_dn, i), ekat::subview(mflx_net, i),
                        ekat::subview(s_upd, i), ekat::subview(q_upd, i),
                        ekat::subview(ql, i), ekat::subview(s_dnd, i),
                        ekat::subview(q_dnd, i), ekat::subview(qst, i),
                        ekat::subview(cu, i), ekat::subview(evp, i),
                        ekat::subview(pflx, i), ekat::subview(rprd, i));
  });
  Kokkos::fence();

  //============================================================================
  // Kernel 6: Unit conversion — per-length [1/m] to per-pressure [1/mb]
  //============================================================================
  Kokkos::parallel_for("zm_conv_main_unit_conv", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    if (!active(i)) return;
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg + 1, pver), [&](const Int k) {
      const Real dz = z_int(i,k) - z_int(i,k+1);
      detr_up(i,k) = detr_up(i,k) * dz / p_del(i,k);
      entr_up(i,k) = entr_up(i,k) * dz / p_del(i,k);
      entr_dn(i,k) = entr_dn(i,k) * dz / p_del(i,k);
      cu(i,k)      = cu(i,k)      * dz / p_del(i,k);
      rprd(i,k)    = rprd(i,k)    * dz / p_del(i,k);
      evp(i,k)     = evp(i,k)     * dz / p_del(i,k);
    });
  });
  Kokkos::fence();

  //============================================================================
  // Host: Compute lel_min and msemax_klev_min across active columns.
  // These cross-column minimums are available to zm_closure for consistent
  // loop bounds across the column ensemble.
  //============================================================================
  Int lel_min = pver - 1, msemax_klev_min = pver - 1;
  Kokkos::parallel_reduce("zm_conv_main_lel_min",
    RangePolicy(0, ncol),
    KOKKOS_LAMBDA(const Int i, Int& mn) {
      if (active(i)) mn = ekat::impl::min(mn, lel(i));
    }, Kokkos::Min<Int>(lel_min));
  Kokkos::parallel_reduce("zm_conv_main_msemax_klev_min",
    RangePolicy(0, ncol),
    KOKKOS_LAMBDA(const Int i, Int& mn) {
      if (active(i)) mn = ekat::impl::min(mn, msemax_klev(i));
    }, Kokkos::Min<Int>(msemax_klev_min));
  Kokkos::fence();

  //============================================================================
  // Kernel 7: Closure — cloud base mass flux
  //============================================================================
  Kokkos::parallel_for("zm_conv_main_closure", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    if (!active(i)) return;
    auto ws = wsm.get_workspace(team);
    zm_closure(team, ws, runtime_opt,
               pver, pverp, msg, cape_threshold_loc,
               lcl(i), lel(i), jt(i), msemax_klev(i), lel_min, msemax_klev_min, dsubcld(i),
               ekat::subview(z_mid, i), ekat::subview(z_int, i),
               ekat::subview(p_mid, i), ekat::subview(p_del, i),
               ekat::subview(t_mid, i),
               ekat::subview(s_mid, i), ekat::subview(q_mid, i),
               ekat::subview(qst, i), ekat::subview(ql, i),
               ekat::subview(s_int, i), ekat::subview(q_int, i),
               t_pcl_lcl(i), ekat::subview(t_pcl, i),
               ekat::subview(q_pcl_sat, i),
               ekat::subview(s_upd, i), ekat::subview(q_upd, i),
               ekat::subview(mflx_net, i), ekat::subview(detr_up, i),
               ekat::subview(mflx_up, i), ekat::subview(mflx_dn, i),
               ekat::subview(q_dnd, i), ekat::subview(s_dnd, i),
               cape(i),
               cld_base_mass_flux(i));
  });
  Kokkos::fence();

  //============================================================================
  // Kernel 8: Limit cloud base mass flux and scale all flux arrays
  //============================================================================
  Kokkos::parallel_for("zm_conv_main_scale", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    if (!active(i)) return;

    // Upper bound: mflx_up/p_del must not exceed 1/(dt * max_val)
    Real mflx_up_max_val = 0;
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, msg + 1, pver),
      [&](const Int k, Real& mx) {
        mx = ekat::impl::max(mx, mflx_up(i,k) / p_del(i,k));
      }, Kokkos::Max<Real>(mflx_up_max_val));
    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      if (mflx_up_max_val > 0) {
        cld_base_mass_flux(i) = ekat::impl::min(cld_base_mass_flux(i),
                                                 1 / (time_step * mflx_up_max_val));
      } else {
        cld_base_mass_flux(i) = 0;
      }
      if (runtime_opt.clos_dyn_adj) {
        cld_base_mass_flux(i) = ekat::impl::max(
          cld_base_mass_flux(i) - omega(i, pbl_top(i)) * ZMC::pa_to_mb, Real(0));
      }
      if (runtime_opt.no_deep_pbl && z_mid_in(i, jt(i)) < pbl_hgt(i)) {
        cld_base_mass_flux(i) = 0;
      }
    });
    team.team_barrier();

    // Scale all flux arrays by cloud base mass flux
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg + 1, pver), [&](const Int k) {
      mflx_up(i,k)  *= cld_base_mass_flux(i);
      mflx_dn(i,k)  *= cld_base_mass_flux(i);
      mflx_net(i,k) *= cld_base_mass_flux(i);
      detr_up(i,k)  *= cld_base_mass_flux(i);
      entr_up(i,k)  *= cld_base_mass_flux(i);
      entr_dn(i,k)  *= cld_base_mass_flux(i);
      rprd(i,k)     *= cld_base_mass_flux(i);
      cu(i,k)       *= cld_base_mass_flux(i);
      evp(i,k)      *= cld_base_mass_flux(i);
      pflx(i,k+1)   *= cld_base_mass_flux(i) * ZMC::mb_to_pa / gravit;
    });
  });
  Kokkos::fence();

  //============================================================================
  // Host: Compute ktm = min(jt) and ktb = min(msemax_klev) across active columns.
  // These are used in zm_calc_output_tend to set the loop bounds consistently.
  //============================================================================
  Int ktm_val = pver - 1, ktb_val = pver - 1;
  Kokkos::parallel_reduce("zm_conv_main_ktm",
    RangePolicy(0, ncol),
    KOKKOS_LAMBDA(const Int i, Int& mn) {
      if (active(i)) mn = ekat::impl::min(mn, jt(i));
    }, Kokkos::Min<Int>(ktm_val));
  Kokkos::parallel_reduce("zm_conv_main_ktb",
    RangePolicy(0, ncol),
    KOKKOS_LAMBDA(const Int i, Int& mn) {
      if (active(i)) mn = ekat::impl::min(mn, msemax_klev(i));
    }, Kokkos::Min<Int>(ktb_val));
  Kokkos::fence();

  //============================================================================
  // Kernel 9: Compute temperature and moisture tendencies
  //============================================================================
  Kokkos::parallel_for("zm_conv_main_output_tend", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    if (!active(i)) return;
    zm_calc_output_tend(team,
                        pver, pverp, msg,
                        jt(i), msemax_klev(i), ktm_val, ktb_val, dsubcld(i),
                        ekat::subview(p_del, i),
                        ekat::subview(s_int, i), ekat::subview(q_int, i),
                        ekat::subview(s_upd, i), ekat::subview(q_upd, i),
                        ekat::subview(mflx_up, i), ekat::subview(detr_up, i),
                        ekat::subview(mflx_dn, i), ekat::subview(s_dnd, i),
                        ekat::subview(q_dnd, i),
                        ekat::subview(ql, i), ekat::subview(evp, i),
                        ekat::subview(cu, i),
                        ekat::subview(dsdt, i), ekat::subview(dqdt, i),
                        ekat::subview(dl, i));
  });
  Kokkos::fence();

  //============================================================================
  // Kernel 10: Scatter results, compute precipitation and reserved liquid
  //============================================================================
  Kokkos::parallel_for("zm_conv_main_scatter", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();
    if (!active(i)) return;

    // Scatter tendencies and fluxes
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg + 1, pver), [&](const Int k) {
      q_mid(i,k) = q_mid_in(i,k) + time_step * dqdt(i,k);
      qtnd(i,k)  = dqdt(i,k);
      zdu(i,k)   = detr_up(i,k);
      mcon(i,k)  = mflx_net(i,k);
      heat(i,k)  = dsdt(i,k) * cpair;
      dlf(i,k)   = dl(i,k);
    });
    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      jctop(i) = jt(i);
      jcbot(i) = msemax_klev(i);
    });
    team.team_barrier();

    // Precipitation: integrate change in water vapor minus detrained cloud water
    Real prec_sum = 0;
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, msg + 1, pver),
      [&](const Int k, Real& local_sum) {
        const Real local = p_del_in(i,k) * (q_mid(i,k) - q_mid_in(i,k))
                         - p_del_in(i,k) * dlf(i,k) * time_step;
        local_sum -= local;
      }, prec_sum);
    team.team_barrier();
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      prec(i) = (1 / gravit) * ekat::impl::max(prec_sum, Real(0)) / time_step / 1000;
    });
    team.team_barrier();

    // Reserved liquid: flux out bottom, added back later
    Real rliq_sum = 0;
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, pver),
      [&](const Int k, Real& local_sum) {
        local_sum += dlf(i,k) * p_del_in(i,k) / gravit;
      }, rliq_sum);
    team.team_barrier();
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      rliq(i) = rliq_sum / 1000;
    });
  });
  Kokkos::fence();

  return active;
}

} // namespace zm
} // namespace scream

#endif
