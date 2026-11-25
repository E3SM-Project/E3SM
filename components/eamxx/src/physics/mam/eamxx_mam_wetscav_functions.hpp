#ifndef EAMXX_MAM_WETSCAV_FUNCTIONS_HPP
#define EAMXX_MAM_WETSCAV_FUNCTIONS_HPP

// MAMWetscav member helpers for ZM convective tracer transport
// (EAM: convect_deep_tend_2 -> zm_conv_tend_2 -> zm_transport_tracer).
// Only available when the ZM deep convection library is built.
#ifdef EAMXX_HAS_ZM

#include <physics/mam/eamxx_mam_wetscav_process_interface.hpp>

#include <mam4xx/mam4.hpp>

#include <ekat_subview_utils.hpp>

namespace scream {

// ---------------------------------------------------------------------------
// Per-column scratch length (in Reals) needed by zm_tracer_transport_column:
//   [ q (nlev*pcnst) | dqdt (nlev*pcnst) | mu (nlev) | md (nlev) | dp (nlev) ]
// ---------------------------------------------------------------------------
inline int MAMWetscav::zm_tracer_transport_work_len(const int nlev) {
  constexpr int pcnst = mam4::aero_model::pcnst;
  return 2 * nlev * pcnst + 3 * nlev;
}

// ---------------------------------------------------------------------------
// Domain-wide cloud top/base bounds over columns with active deep convection
// (EAM: ktm/kbm computed over the gathered columns inside
// zm_transport_tracer). Indices follow ZM's convention (0 = model top).
// Returns nlev-1 for both when no column is active.
// ---------------------------------------------------------------------------
inline void MAMWetscav::zm_cloud_level_bounds(const int ncol, const int nlev,
                                              const const_view_1d &cldtop_klev,
                                              const const_view_1d &msemax_klev,
                                              const const_view_1d &activity,
                                              int &ktm, int &kbm) {
  int ktm_red = nlev - 1;
  int kbm_red = nlev - 1;
  Kokkos::parallel_reduce("MAMWetscav::zm_ktm", ncol,
      KOKKOS_LAMBDA(const int icol, int &lmin) {
        if(Kokkos::round(activity(icol)) > 0) {
          lmin = Kokkos::min(lmin, static_cast<int>(Kokkos::round(cldtop_klev(icol))));
        }
      }, Kokkos::Min<int>(ktm_red));
  Kokkos::parallel_reduce("MAMWetscav::zm_kbm", ncol,
      KOKKOS_LAMBDA(const int icol, int &lmin) {
        if(Kokkos::round(activity(icol)) > 0) {
          lmin = Kokkos::min(lmin, static_cast<int>(Kokkos::round(msemax_klev(icol))));
        }
      }, Kokkos::Min<int>(kbm_red));
  ktm = Kokkos::min(nlev - 1, ktm_red);
  kbm = Kokkos::min(nlev - 1, kbm_red);
}

// ---------------------------------------------------------------------------
// Gather everything the column kernel needs for ZM tracer transport.
// Returns enabled=false (and empty views) when zm_tracer_transport is off.
// ---------------------------------------------------------------------------
inline MAMWetscav::ZmTransportInputs
MAMWetscav::get_zm_transport_inputs() const {
  ZmTransportInputs zm;
  // Only zm_microp is read by zm_transport_tracer; s_zm_opts is filled by
  // the ZM process (zero-initialized, i.e. zm_microp=false, otherwise).
  zm.runtime_opt = ZMF::s_zm_opts;
  zm.ktm         = nlev_ - 1;
  zm.kbm         = nlev_ - 1;
  zm.enabled     = zm_tracer_transport_;
  if(!zm.enabled) return zm;

  zm.cldtop_klev = get_field_in("zm_cldtop_klev").get_view<const Real *>();
  zm.msemax_klev = get_field_in("zm_msemax_klev").get_view<const Real *>();
  zm.activity    = get_field_in("zm_activity").get_view<const Real *>();
  zm.mflx_up     = get_field_in("zm_mflx_up").get_view<const Real **>();
  zm.mflx_dn     = get_field_in("zm_mflx_dn").get_view<const Real **>();
  zm.entr_up     = get_field_in("zm_entr_up").get_view<const Real **>();
  zm.detr_up     = get_field_in("zm_detr_up").get_view<const Real **>();
  zm.entr_dn     = get_field_in("zm_entr_dn").get_view<const Real **>();
  zm.doconvtran  = zm_doconvtran_;
  zm.fracis      = zm_fracis_;
  zm.work        = zm_transport_work_;
  zm.wsm         = *zm_wsm_;

  // Domain-wide cloud top/base bounds over active columns
  zm_cloud_level_bounds(ncol_, nlev_, zm.cldtop_klev, zm.msemax_klev,
                        zm.activity, zm.ktm, zm.kbm);
  return zm;
}

// ---------------------------------------------------------------------------
// ZM convective tracer transport for one column, applied to the MAM dry
// state after aero_model_wetdep so the transport sees the scavenged state.
//
// Gases and interstitial aerosols are packed into the pcnst-ordered state
// array, transported for every constituent with doconvtran(m) = true, and
// the tendencies are applied in place to progs (negative values clipped to
// zero, as EAM's qneg3 does in physics_update).
//
// Must be called by the whole team, and only for columns with active ZM
// convection (cloud top/base indices are undefined otherwise).
// ---------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
void MAMWetscav::zm_tracer_transport_column(const ThreadTeam &team,
                                            const ZmTransportInputs &zm,
                                            const int icol, const int nlev,
                                            const Real dt,
                                            const mam4::Atmosphere &atm,
                                            const mam4::Prognostics &progs,
                                            const const_view_2d &p_del) {
  constexpr int pcnst    = mam4::aero_model::pcnst;
  constexpr Real mb_to_pa = ZMF::ZMC::mb_to_pa;
  // ZM exports mass fluxes in [kg/m2/s]; zm_transport_tracer expects [mb/s]
  // together with dp in [mb] (inverse of 100/g used in the ZM interface)
  constexpr Real kgm2s_to_mb = ZMF::PC::gravit.value / mb_to_pa;

  // carve per-column scratch: q | dqdt | mu | md | dp
  Real *work = ekat::subview(zm.work, icol).data();
  const ZMF::uview_2d<Real> q(work, nlev, pcnst);
  work += nlev * pcnst;
  const ZMF::uview_2d<Real> dqdt(work, nlev, pcnst);
  work += nlev * pcnst;
  const ZMF::uview_1d<Real> mu(work, nlev);
  work += nlev;
  const ZMF::uview_1d<Real> md(work, nlev);
  work += nlev;
  const ZMF::uview_1d<Real> dp(work, nlev);

  const auto mflx_up = ekat::subview(zm.mflx_up, icol);
  const auto mflx_dn = ekat::subview(zm.mflx_dn, icol);
  const auto p_del_c = ekat::subview(p_del, icol);
  const ZMF::uview_1d<const bool> doconvtran(zm.doconvtran.data(), pcnst);

  // pack gases and interstitial aerosols (dry mmr) and convert ZM inputs
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](int kk) {
    auto q_kk = ekat::subview(q, kk);
    for(int m = 0; m < pcnst; ++m) {
      q_kk(m)     = 0;
      dqdt(kk, m) = 0;
    }
    mam4::utils::extract_stateq_from_prognostics(progs, atm, q_kk, kk);
    mu(kk) = mflx_up(kk) * kgm2s_to_mb;
    md(kk) = mflx_dn(kk) * kgm2s_to_mb;
    dp(kk) = p_del_c(kk) / mb_to_pa;
  });
  team.team_barrier();

  // dpdry is unused by zm_transport_tracer (moist formulation only)
  ZMF::zm_transport_tracer(
      team, zm.wsm.get_workspace(team), zm.runtime_opt, nlev,
      doconvtran, q, pcnst, mu, md,
      ekat::subview(zm.detr_up, icol),
      ekat::subview(zm.entr_up, icol),
      ekat::subview(zm.entr_dn, icol),
      dp,
      static_cast<int>(Kokkos::round(zm.cldtop_klev(icol))),
      static_cast<int>(Kokkos::round(zm.msemax_klev(icol))),
      zm.ktm, zm.kbm,
      ZMF::uview_2d<const Real>(zm.fracis.data(), nlev, pcnst),
      dp, dt, dqdt);
  team.team_barrier();

  // apply tendencies to transported species, using the same gas/aerosol
  // ordering as extract_stateq_from_prognostics
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](int kk) {
    int s_idx = mam4::utils::gasses_start_ind();
    for(int g = 0; g < mam4::AeroConfig::num_gas_ids(); ++g, ++s_idx) {
      if(doconvtran(s_idx) && progs.q_gas[g].data()) {
        progs.q_gas[g](kk) =
            Kokkos::max(Real(0), progs.q_gas[g](kk) + dqdt(kk, s_idx) * dt);
      }
    }
    s_idx = mam4::utils::aero_start_ind();
    for(int m = 0; m < mam4::AeroConfig::num_modes(); ++m) {
      for(int a = 0; a < mam4::num_species_mode(m); ++a, ++s_idx) {
        if(doconvtran(s_idx) && progs.q_aero_i[m][a].data()) {
          progs.q_aero_i[m][a](kk) = Kokkos::max(
              Real(0), progs.q_aero_i[m][a](kk) + dqdt(kk, s_idx) * dt);
        }
      }
      if(doconvtran(s_idx)) {
        progs.n_mode_i[m](kk) =
            Kokkos::max(Real(0), progs.n_mode_i[m](kk) + dqdt(kk, s_idx) * dt);
      }
      ++s_idx;
    }
  });
  team.team_barrier();
}

}  // namespace scream

#endif  // EAMXX_HAS_ZM
#endif  // EAMXX_MAM_WETSCAV_FUNCTIONS_HPP
