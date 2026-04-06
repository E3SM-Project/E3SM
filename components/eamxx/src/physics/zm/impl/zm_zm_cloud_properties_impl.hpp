#ifndef ZM_ZM_CLOUD_PROPERTIES_IMPL_HPP
#define ZM_ZM_CLOUD_PROPERTIES_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_cloud_properties. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_cloud_properties(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const ZmRuntimeOpt& runtime_opt,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& msg, // number of levels to ignore at model top
  const Int& limcnv, // convection limiting level
  const uview_1d<const Real>& p_mid, // env pressure at mid-point
  const uview_1d<const Real>& z_mid, // env altitude at mid-point
  const uview_1d<const Real>& z_int, // env altitude at interface
  const uview_1d<const Real>& t_mid, // env temperature
  const uview_1d<const Real>& s_mid, // env dry static energy of env [K] (normalized)
  const uview_1d<const Real>& s_int, // interface values of dry stat energy
  const uview_1d<const Real>& q_mid, // env specific humidity
  const Real& landfrac, // Land fraction
  const Real& tpert_g, // PBL temperature perturbation
  const Int& jb, // updraft base level
  const Int& lel, // updraft parcel launch level
  // Outputs
  Int& jt, // updraft plume top
  Int& jlcl, // updraft lifting cond level
  Int& j0, // level where detrainment begins (starting at h_env_min)
  Int& jd, // level of downdraft
  const uview_1d<Real>& mflx_up, // updraft mass flux
  const uview_1d<Real>& entr_up, // entrainment rate of updraft
  const uview_1d<Real>& detr_up, // detrainement rate of updraft
  const uview_1d<Real>& mflx_dn, // downdraft mass flux
  const uview_1d<Real>& entr_dn, // downdraft entrainment rate
  const uview_1d<Real>& mflx_net, // net mass flux
  const uview_1d<Real>& s_upd, // updraft dry static energy [K] (normalized)
  const uview_1d<Real>& q_upd, // updraft specific humidity [kg/kg]
  const uview_1d<Real>& ql, // updraft liq water
  const uview_1d<Real>& s_dnd, // dndraft dry static energy [K] (normalized)
  const uview_1d<Real>& q_dnd, // dndraft specific humidity [kg/kg]
  const uview_1d<Real>& qst, // env saturation mixing ratio
  const uview_1d<Real>& cu, // condensation rate
  const uview_1d<Real>& evp, // evaporation rate
  const uview_1d<Real>& pflx, // precipitation flux thru layer (pverp-sized)
  const uview_1d<Real>& rprd) // rate of production of precip at that layer
{
  // Physics constants (local aliases for readability)
  const Real latvap = PC::LatVap.value;
  const Real cpair  = PC::Cpair.value;
  const Real grav   = PC::gravit.value;
  const Real epsilo = PC::ep_2.value;
  const Real rdair  = PC::Rair.value;

  // Minimum updraft mass flux below which the layer is considered inactive
  constexpr Real mu_min      = 0.02;
  // If updraft MSE undershoots saturated MSE by more than this, cloud top is 1 level up
  constexpr Real hu_diff_min = -2000.0;

  // =========================================================================
  // Allocate temporary workspace arrays (all pver-sized; pflxs used as pverp)
  // =========================================================================
  uview_1d<Real> gamma, dz, h_env, h_env_sat, h_upd, h_dnd, qsthat, hsthat,
                 gamhat, q_dnd_sat, lambda, fice, tug, tmp_frz, pflxs;
  workspace.template take_many_contiguous_unsafe<15>(
    {"gamma", "dz", "h_env", "h_env_sat", "h_upd", "h_dnd", "qsthat", "hsthat",
     "gamhat", "q_dnd_sat", "lambda", "fice", "tug", "tmp_frz", "pflxs"},
    {&gamma, &dz, &h_env, &h_env_sat, &h_upd, &h_dnd, &qsthat, &hsthat,
     &gamhat, &q_dnd_sat, &lambda, &fice, &tug, &tmp_frz, &pflxs});

  // Scalar locals (captured by reference in Kokkos lambdas)
  Real c0mask = 0, totpcp = 0, totevp = 0, h_env_min = 0, lambda_max = 0, tot_frz = 0;
  Int  jto = 0, khighest = 0, klowest = 0, kount = 0;

  // Initialize scalar parameters (serial)
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    totpcp  = 0;
    totevp  = 0;
    // Land/ocean blend of autoconversion coefficient
    c0mask  = runtime_opt.c0_ocn * (1 - landfrac) + runtime_opt.c0_lnd * landfrac;
    h_env_min  = 1.e6;
    lambda_max = 0;
    tot_frz    = 0;
    j0         = msg;
  });

  // =========================================================================
  // 1. Initialize 2D variables (each k independent)
  // =========================================================================
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pver), [&] (const Int& k) {
    // Saturation specific humidity at mid-level
    Real est, qs_val;
    qsat_hPa(t_mid(k), p_mid(k), est, qs_val);
    qst(k) = qs_val;

    // Lapse rate parameter gamma = (Lv/Cp) * d(qsat)/dT  - see eq (4.118)
    gamma(k) = qst(k) * (1 + qst(k)/epsilo) * epsilo *
               latvap / (rdair * t_mid(k)*t_mid(k)) * latvap / cpair;

    // Layer thickness [m]
    dz(k) = z_int(k) - z_int(k+1);

    // Environmental moist static energy
    h_env(k)     = cpair*t_mid(k) + grav*z_mid(k) + latvap*q_mid(k);
    h_env_sat(k) = cpair*t_mid(k) + grav*z_mid(k) + latvap*qst(k);

    // Initialize updraft/downdraft quantities to environment values
    h_upd(k)     = h_env(k);
    h_dnd(k)     = h_env(k);
    s_upd(k)     = s_mid(k);
    s_dnd(k)     = s_mid(k);
    q_upd(k)     = q_mid(k);
    q_dnd(k)     = q_mid(k);
    q_dnd_sat(k) = q_mid(k);

    // Zero mass fluxes and tendency arrays
    mflx_up(k) = 0;  entr_up(k) = 0;  detr_up(k) = 0;
    mflx_dn(k) = 0;  entr_dn(k) = 0;  mflx_net(k) = 0;
    cu(k)      = 0;  evp(k)     = 0;
    ql(k)      = 0;  rprd(k)    = 0;
    pflx(k)    = 0;

    // Zero microphysics and workspace arrays
    lambda(k)  = 0;  fice(k)    = 0;
    tug(k)     = 0;  tmp_frz(k) = 0;
    pflxs(k)   = 0;
  });
  team.team_barrier();

  // =========================================================================
  // 2. Interpolate thermodynamic quantities to interfaces (each k independent)
  //    Geometric mean interpolation - see eq (4.109), (4.118), (4.119)
  // =========================================================================
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg, pver), [&] (const Int& k) {
    // Default: use mid-level values at and above msg
    hsthat(k) = h_env_sat(k);
    qsthat(k) = qst(k);
    gamhat(k) = gamma(k);

    if (k > msg) {
      // Geometric mean interpolation for qst - see eq (4.109)
      if (std::abs(qst(k-1) - qst(k)) > ZMC::interp_diff_min) {
        qsthat(k) = std::log(qst(k-1)/qst(k)) * qst(k-1)*qst(k) / (qst(k-1) - qst(k));
      } else {
        qsthat(k) = qst(k);
      }
      // Interface saturated moist static energy
      hsthat(k) = cpair*s_int(k) + latvap*qsthat(k);

      // Geometric mean interpolation for gamma - see eq (4.118)
      if (std::abs(gamma(k-1) - gamma(k)) > ZMC::interp_diff_min) {
        gamhat(k) = std::log(gamma(k-1)/gamma(k)) * gamma(k-1)*gamma(k) / (gamma(k-1) - gamma(k));
      } else {
        gamhat(k) = gamma(k);
      }
    }
  });
  team.team_barrier();

  // =========================================================================
  // 3. Initialize jt, jd, jlcl; find j0 (level of min saturated MSE)
  //    Running minimum: serial
  // =========================================================================
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    // Initial cloud top estimate: at least one level below the convection limit
    jt   = ekat::impl::max(lel, limcnv + 1);
    jt   = ekat::impl::min(jt,  pver - 1);
    jd   = pver - 1;
    jlcl = lel;
    h_env_min = 1.e6;

    // Find level of minimum h_env_sat between jt and jb (detrainment onset level)
    for (Int k = msg; k < pver; ++k) {
      if (h_env_sat(k) <= h_env_min && k >= jt && k <= jb) {
        h_env_min = h_env_sat(k);
        j0 = k;
      }
    }

    // Constrain j0 to a physically valid detrainment range
    j0 = ekat::impl::min(j0, jb - 2);
    j0 = ekat::impl::max(j0, jt + 2);
    j0 = ekat::impl::min(j0, pver - 1);
  });
  team.team_barrier();

  // =========================================================================
  // 4. Initialize updraft MSE with PBL temperature perturbation (each k independent)
  // =========================================================================
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg, pver), [&] (const Int& k) {
    if (k >= jt && k <= jb) {
      if (runtime_opt.tpert_fix) {
        // No PBL perturbation applied (tpert_fix disables it)
        h_upd(k) = h_env(jb) + cpair * runtime_opt.tiedke_add;
        s_upd(k) = s_mid(jb) + runtime_opt.tiedke_add;
      } else {
        h_upd(k) = h_env(jb) + cpair * (runtime_opt.tiedke_add + runtime_opt.tpert_fac * tpert_g);
        s_upd(k) = s_mid(jb) + runtime_opt.tiedke_add + runtime_opt.tpert_fac * tpert_g;
      }
    }
  });
  team.team_barrier();

  // =========================================================================
  // 5. Compute fractional entrainment rate profile
  // =========================================================================
  zm_calc_fractional_entrainment(team, pver, pverp, msg, jb, jt, j0,
                                  z_mid, z_int, dz, h_env, h_env_sat,
                                  h_env_min, lambda, lambda_max);
  team.team_barrier();

  // =========================================================================
  // ITERATION (itnum = 1 in Fortran: no outer loop needed)
  // =========================================================================

  // 6a. Zero condensation and liquid water arrays (each k independent)
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pver), [&] (const Int& k) {
    cu(k) = 0;
    ql(k) = 0;
  });
  team.team_barrier();

  // =========================================================================
  // 6b. Compute updraft mass flux profile from cloud base upward
  //     (serial: entr_up(k) depends on mflx_up(k+1) from previous iteration)
  // =========================================================================
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    if (lambda_max > 0) {
      // Cloud base: unit mass flux, all entrainment
      mflx_up(jb) = 1;
      entr_up(jb) = mflx_up(jb) / dz(jb);

      // Integrate upward using the fractional entrainment lambda profile
      for (Int k = pver - 1; k >= msg; --k) {
        if (k >= jt && k < jb) {
          // Height above cloud base interface
          const Real zuef = z_int(k) - z_int(jb);
          // Mass flux at top of layer (interface between k and k+1)
          const Real rmue  = (1/lambda_max) * (std::exp(lambda(k+1)*zuef) - 1) / zuef;
          // Mass flux at bottom of layer
          mflx_up(k) = (1/lambda_max) * (std::exp(lambda(k)*zuef) - 1) / zuef;
          entr_up(k) = (rmue - mflx_up(k+1)) / dz(k);
          detr_up(k) = (rmue - mflx_up(k))   / dz(k);
        }
      }
    }
  });
  team.team_barrier();

  // =========================================================================
  // 6c. Update updraft MSE profile from cloud base upward
  //     (serial: h_upd(k) depends on h_upd(k+1) from previous iteration)
  // =========================================================================
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    khighest = lel;
    klowest  = jb;
    for (Int k = klowest - 1; k >= khighest; --k) {
      if (k <= jb - 1 && k >= lel && lambda_max > 0) {
        if (mflx_up(k) < mu_min) {
          // Below threshold: detrain all upward flux, reset to environment
          h_upd(k)   = h_env(k);
          mflx_up(k) = 0;
          entr_up(k) = 0;
          detr_up(k) = mflx_up(k+1) / dz(k);
        } else {
          // MSE budget: updraft MSE is flux-weighted mix of below and environment
          h_upd(k) = mflx_up(k+1)/mflx_up(k) * h_upd(k+1)
                   + dz(k)/mflx_up(k) * (entr_up(k)*h_env(k) - detr_up(k)*h_env_sat(k));
        }
      }
    }
  });
  team.team_barrier();

  // =========================================================================
  // 6d. Accumulate freezing contribution for cloud top determination
  //     (tmp_frz currently zero; parallel reduce for future microphysics)
  // =========================================================================
  Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, pver),
                           [&] (const Int& k, Real& val) {
                             val += tmp_frz(k) * dz(k);
                           }, tot_frz);
  team.team_barrier();

  // =========================================================================
  // 6e. Determine cloud top by comparing updraft and saturated env MSE
  //     (serial: doit flag creates sequential dependency)
  // =========================================================================
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    bool doit = true;
    for (Int k = klowest - 2; k >= khighest - 1; --k) {
      if (doit && k <= jb - 2 && k >= lel - 1) {
        if (h_upd(k) <= hsthat(k) && h_upd(k+1) > hsthat(k+1)
            && mflx_up(k) >= mu_min) {
          // Updraft MSE crosses saturation MSE from above: this is the cloud top
          if (h_upd(k) - hsthat(k) < hu_diff_min) {
            jt = k + 1;  // Large undershoot: cloud top is one level higher
          } else {
            jt = k;
          }
          doit = false;
        } else if ((h_upd(k) > h_upd(jb) && tot_frz <= 0)
                   || mflx_up(k) < mu_min) {
          // Neutral or negative buoyancy, or sub-threshold mass flux
          jt   = k + 1;
          doit = false;
        }
      }
    }
    // Save cloud top from this (only) iteration
    jto = jt;
  });
  team.team_barrier();

  // =========================================================================
  // 6f. Zero mass fluxes at and above updated cloud top (each k independent)
  // =========================================================================
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg, pver), [&] (const Int& k) {
    if (lambda_max > 0) {
      // Zero all updraft quantities at and above cloud top
      if (k >= lel && k <= jt) {
        mflx_up(k) = 0;
        entr_up(k) = 0;
        detr_up(k) = 0;
        h_upd(k)   = h_env(k);
      }
      // At the cloud top level: non-zero detrainment to conserve mass
      if (k == jt) {
        detr_up(k) = mflx_up(k+1) / dz(k);
        entr_up(k) = 0;
        mflx_up(k) = 0;
      }
    }
  });
  team.team_barrier();

  // =========================================================================
  // 6g. LCL determination: compute s_upd, q_upd from cloud base upward and
  //     find the level where the updraft first reaches saturation
  //     (serial: done flag and upward mass-balance chain)
  // =========================================================================
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    bool done = false;
    kount = 0;
    for (Int k = pver - 1; k >= msg + 1; --k) {
      // Initialize s_upd, q_upd at cloud base from h_upd
      if (k == jb && lambda_max > 0) {
        q_upd(k) = q_mid(jb);
        s_upd(k) = (h_upd(k) - latvap * q_upd(k)) / cpair;
      }
      // Compute updraft thermodynamics above cloud base up to cloud top
      if (!done && k > jt && k < jb && lambda_max > 0) {
        s_upd(k) = mflx_up(k+1)/mflx_up(k) * s_upd(k+1)
                 + dz(k)/mflx_up(k) * (entr_up(k) - detr_up(k)) * s_mid(k);
        q_upd(k) = mflx_up(k+1)/mflx_up(k) * q_upd(k+1)
                 + dz(k)/mflx_up(k) * (entr_up(k)*q_mid(k) - detr_up(k)*qst(k));

        // Updraft temperature at the layer interface
        const Real tu    = s_upd(k) - grav/cpair * z_int(k);
        // Saturation check at interface pressure (average of adjacent mid-levels)
        const Real p_ifc = (p_mid(k) + p_mid(k-1)) / 2;
        Real estu, qstu;
        qsat_hPa(tu, p_ifc, estu, qstu);

        // LCL is the highest level where updraft is just saturated on ascent
        if (q_upd(k) >= qstu) {
          jlcl = k;
          ++kount;
          done = true;  // Fortran equivalent: goto 690
        }
      }
    }
  });
  team.team_barrier();

  // =========================================================================
  // 6h. Above LCL: compute s_upd/q_upd from MSE and interface saturation
  //     (each k independent given jt and jlcl are known)
  // =========================================================================
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg + 1, pver), [&] (const Int& k) {
    if (k > jt && k <= jlcl && lambda_max > 0) {
      // Saturated updraft: s and q follow from MSE and saturation constraint
      s_upd(k) = s_int(k)  + (h_upd(k) - hsthat(k)) / (cpair  * (1 + gamhat(k)));
      q_upd(k) = qsthat(k) + gamhat(k) * (h_upd(k) - hsthat(k)) / (latvap * (1 + gamhat(k)));
    }
  });
  team.team_barrier();

  // =========================================================================
  // 6i. Condensation rate from DSE budget (each k independent given mflx_up)
  // =========================================================================
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg + 1, pver), [&] (const Int& k) {
    if (lambda_max > 0 && k >= jt && k < jb) {
      // Condensation = latent heat release implied by DSE flux divergence
      cu(k) = ((mflx_up(k)*s_upd(k) - mflx_up(k+1)*s_upd(k+1))/dz(k)
               - (entr_up(k) - detr_up(k))*s_mid(k))
              / (latvap/cpair);
      if (k == jt) cu(k) = 0;                           // No condensation at cloud top
      cu(k) = ekat::impl::max(Real(0), cu(k));           // Condensation must be non-negative
    }
  });
  team.team_barrier();

  // =========================================================================
  // 6j. Liquid water budget, rain production, total precipitation
  //     (serial: ql(k) depends on ql(k+1) from level below)
  // =========================================================================
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    totpcp = 0;
    for (Int k = pver - 1; k >= msg + 1; --k) {
      rprd(k) = 0;
      if (k >= jt && k < jb && lambda_max > 0 && mflx_up(k) >= 0) {
        if (mflx_up(k) > 0) {
          // Liquid water mass budget (without autoconversion sink)
          const Real ql1 = (mflx_up(k+1)*ql(k+1)
                            - dz(k)*detr_up(k)*ql(k+1)
                            + dz(k)*cu(k)) / mflx_up(k);
          // Apply autoconversion sink (c0mask * dz)
          ql(k) = ql1 / (1 + dz(k)*c0mask);
        } else {
          ql(k) = 0;
        }
        // Net column precipitation = condensation minus detrained liquid
        totpcp += dz(k) * (cu(k) - detr_up(k)*ql(k+1));
        // Rain production rate by autoconversion
        rprd(k) = c0mask * mflx_up(k) * ql(k);
      }
    }
  });
  team.team_barrier();

  // =========================================================================
  // 7. Downdraft properties
  // =========================================================================
  zm_downdraft_properties(team, runtime_opt, pver, pverp, msg, jb, jt, j0, jd,
                           z_int, dz, s_mid, q_mid, h_env,
                           lambda, lambda_max,
                           qsthat, hsthat, gamhat,
                           rprd, mflx_up,
                           mflx_dn, entr_dn, s_dnd, q_dnd, h_dnd,
                           q_dnd_sat, evp, totevp);
  team.team_barrier();

  // =========================================================================
  // 8. Scale downdraft by precipitation availability - see eq (4.106)
  // =========================================================================
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    totpcp = ekat::impl::max(Real(0), totpcp);
    totevp = ekat::impl::max(Real(0), totevp);
  });
  team.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg + 1, pver), [&] (const Int& k) {
    if (totevp > 0 && totpcp > 0) {
      // Scale downdraft so total evaporation does not exceed available precipitation
      const Real scale = ekat::impl::min(Real(1), totpcp / (totevp + totpcp));
      mflx_dn(k) *= scale;
      entr_dn(k) *= scale;
      evp(k)     *= scale;
    } else {
      mflx_dn(k) = 0;
      entr_dn(k) = 0;
      evp(k)     = 0;
    }
    // Net precipitation production after evaporation
    rprd(k) = rprd(k) - evp(k);
  });
  team.team_barrier();

  // =========================================================================
  // 9. Cumulative precipitation flux at interfaces (serial: running sum)
  // =========================================================================
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    pflx(0) = 0;
    for (Int k = 1; k <= pver; ++k) {
      pflx(k) = pflx(k-1) + rprd(k-1) * dz(k-1);
    }
  });
  team.team_barrier();

  // =========================================================================
  // 10. Net mass flux (each k independent)
  // =========================================================================
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg, pver), [&] (const Int& k) {
    mflx_net(k) = mflx_up(k) + mflx_dn(k);
  });
  team.team_barrier();

  // Release workspace
  workspace.template release_many_contiguous<15>(
    {&gamma, &dz, &h_env, &h_env_sat, &h_upd, &h_dnd, &qsthat, &hsthat,
     &gamhat, &q_dnd_sat, &lambda, &fice, &tug, &tmp_frz, &pflxs});
}

} // namespace zm
} // namespace scream

#endif
