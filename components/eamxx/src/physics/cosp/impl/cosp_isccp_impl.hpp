#ifndef SCREAM_COSP_ISCCP_IMPL_HPP
#define SCREAM_COSP_ISCCP_IMPL_HPP

#include "physics/cosp/cosp_cxx_functions.hpp"

namespace scream {
namespace cosp {

/*
 * ICARUS, the ISCCP simulator (icarus.F90), for the configuration EAMxx uses:
 * isccp_top_height = 1 (cloud top from the IR brightness temperature, adjusted for
 * partially transparent clouds) and isccp_top_height_direction = 2 (pick the highest
 * level matching the brightness temperature). Levels go from the top down.
 *
 * The Fortran works on all points at once with (npoints,ncol) and (npoints,nlev) work
 * arrays. Here the per-point part (tropopause, clear-sky flux) is done once per point,
 * and each subcolumn is then independent. Per-level quantities that the Fortran stores
 * (dem_wv, bb) are recomputed where needed, which gives identical values and avoids
 * scratch arrays.
 */

// Water vapor continuum emissivity of one layer at 10.47 microns (the dem_wv loop in
// icarus_subcolumn), following Schwarzkopf and Ramaswamy (JGR 1999)
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::icarus_dem_wv(const S pfull, const S phalf_top, const S phalf_bot,
                                const S qv, const S at)
{
  const S press  = pfull*S(10);
  const S dpress = (phalf_bot - phalf_top)*10;
  const S atmden = dpress/(C::grav*S(100));
  const S rvh20  = qv*C::amd/C::amw;
  const S wk     = rvh20*C::avo*atmden/C::amd;
  const S rhoave = (press/C::pstd)*(C::isccp_t0/at);
  const S rh20s  = rvh20*rhoave;
  const S rfrgn  = rhoave - rh20s;
  const S tmpexp = Kokkos::exp(S(-0.02)*(at - C::isccp_t0));
  // 1.e-20 and 3.41817e-7 are single precision literals in the Fortran
  const S tauwv  = wk*S(1.e-20f)*((S(0.0224697)*rh20s*tmpexp) + (S(3.41817e-7f)*rfrgn))*S(0.98);
  return 1 - Kokkos::exp(-1*tauwv);
}

// Planck-like 10.5 micron emission used throughout ICARUS
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::icarus_bb(const S t)
{
  return 1 / (Kokkos::exp(S(1307.27)/t) - 1);
}

// Brightness temperature from a TOA flux
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::icarus_tb(const S flux)
{
  return S(1307.27)/(Kokkos::log(1 + (1/flux)));
}

// Per-point part of icarus_subcolumn: tropopause, clear-sky TOA flux and brightness
// temperature, and the minimum cloud brightness flux
template <typename S, typename D>
template <typename VQ, typename VA, typename VP, typename VPH>
KOKKOS_INLINE_FUNCTION
typename Functions<S,D>::IcarusPoint
Functions<S,D>::icarus_point(const int nlev, const S skt, const S emsfc_lw, const VQ& qv,
                             const VA& at, const VP& pfull, const VPH& phalf)
{
  IcarusPoint pt;

  // Tropopause: the coldest level between 50 and 400 hPa (bounds are single precision
  // literals in the Fortran, but exact)
  pt.ptrop  = 5000;
  pt.attrop = 120;
  pt.itrop  = 0;
  S attropmin = 400;
  for (int k = 0; k < nlev; ++k) {
    if (pfull(k) < 40000 && pfull(k) > 5000 && at(k) < attropmin) {
      pt.ptrop   = pfull(k);
      attropmin  = at(k);
      pt.attrop  = attropmin;
      pt.itrop   = k;
    }
  }
  pt.atmax = 0;
  for (int k = 0; k < nlev; ++k) {
    if (at(k) > pt.atmax && k >= pt.itrop) pt.atmax = at(k);
  }

  // Clear-sky radiance
  S flux = 0, trans = 1;
  for (int k = 0; k < nlev; ++k) {
    const S bb     = icarus_bb(at(k));
    const S dem_wv = icarus_dem_wv(pfull(k), phalf(k), phalf(k+1), qv(k), at(k));
    flux  = flux + dem_wv*bb*trans;
    trans = trans*(1 - dem_wv);
  }
  flux = flux + emsfc_lw*icarus_bb(skt)*trans;
  pt.fluxtop_clrsky = flux;
  pt.meantbclr      = icarus_tb(flux);

  // Minimum brightness temperature flux
  pt.btcmin = 1 / (Kokkos::exp(S(1307.27)/(pt.attrop - 5)) - 1);

  return pt;
}

// icarus_subcolumn for one subcolumn of one point; dtau and dem are indexed by level.
// Outputs the subcolumn optical depth and cloud-top pressure (Pa) as ISCCP would see
// them (R_UNDEF if clear or not sunlit), and the all-sky brightness temperature.
template <typename S, typename D>
template <typename VT, typename VE, typename VQ, typename VA, typename VP, typename VPH>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::icarus_subcolumn(const int nlev, const int sunlit, const IcarusPoint& pt,
                                      const VT& dtau, const VE& dem, const S skt,
                                      const S emsfc_lw, const VQ& qv, const VA& at,
                                      const VP& pfull, const VPH& phalf,
                                      S& boxtau, S& boxptop, S& boxttop)
{
  constexpr S tauchk = C::isccp_tauchk;

  // Subcolumn optical depth
  S tau = 0;
  for (int k = 0; k < nlev; ++k) tau += dtau(k);

  // All-sky radiance
  S fluxtop = 0, trans = 1;
  for (int k = 0; k < nlev; ++k) {
    const S bb     = icarus_bb(at(k));
    const S dem_wv = icarus_dem_wv(pfull(k), phalf(k), phalf(k+1), qv(k), at(k));
    const S demk   = dem(k) == 0 ? dem_wv : 1 - (1 - dem(k))*(1 - dem_wv);
    fluxtop = fluxtop + demk*bb*trans;
    trans   = trans*(1 - demk);
  }
  fluxtop = fluxtop + emsfc_lw*icarus_bb(skt)*trans;
  boxttop = icarus_tb(fluxtop);

  // Cloud-top temperature, accounting for partially transmitting clouds as ISCCP does.
  // The 0.001 and 0.9999999 in the transmax test are single precision literals.
  const S transmax = (fluxtop - pt.btcmin)/(pt.fluxtop_clrsky - pt.btcmin);
  S tauir = tau/S(2.13);
  const S taumin = -Kokkos::log(Kokkos::max(Kokkos::min(transmax, S(0.9999999)), S(0.001)));
  const bool partial = transmax > S(0.001f) && transmax <= S(0.9999999f);
  S fluxtopinit = 0;
  if (partial) {
    fluxtopinit = fluxtop;
    tauir = tau/S(2.13);
  }
  S tb = 0;
  for (int icycle = 0; icycle < 2; ++icycle) {
    if (tau > tauchk && partial) {
      const S emcld = 1 - Kokkos::exp(-1*tauir);
      fluxtop = fluxtopinit - ((1 - emcld)*pt.fluxtop_clrsky);
      fluxtop = Kokkos::max(S(1.E-06), (fluxtop/emcld));
      tb = icarus_tb(fluxtop);
      if (tb > 260) {
        tauir = tau / S(2.56);
      }
    }
  }
  if (tau > tauchk) {
    tb = icarus_tb(fluxtop);
    if (tauir < taumin) {
      tb  = pt.attrop - 5;
      tau = S(2.13)*taumin;
    }
  }
  if (tau <= tauchk) {
    tb = pt.meantbclr;
  }

  // Cloud-top pressure: find the highest level (searching up from the surface, and not
  // above the tropopause) whose temperature brackets tb, and interpolate in log(p)
  int kmatch = -1;
  for (int k1 = 1; k1 <= nlev-1; ++k1) {
    const int k = nlev-1-k1;
    if (k >= pt.itrop &&
        ((at(k) >= tb && at(k+1) <= tb) || (at(k) <= tb && at(k+1) >= tb))) {
      kmatch = k;
    }
  }
  // Note: if nothing matches and attrop < tb < atmax, the Fortran leaves ptop unset
  // (garbage). That requires tb outside the temperature range below the tropopause,
  // which realistic profiles never produce; here ptop is 0 in that case.
  S ptop = 0;
  if (kmatch >= 0) {
    const int k1 = kmatch, k2 = kmatch+1;
    const S logp1 = Kokkos::log(pfull(k1));
    const S logp2 = Kokkos::log(pfull(k2));
    const S atd   = Kokkos::max(tauchk, Kokkos::abs(at(k2) - at(k1)));
    const S logp  = logp1 + (logp2 - logp1)*Kokkos::abs(tb - at(k1))/atd;
    ptop = Kokkos::exp(logp);
  } else {
    if (tb <= pt.attrop) ptop = pt.ptrop;
    if (tb >= pt.atmax)  ptop = pfull(nlev-1);
  }
  if (tau <= tauchk) ptop = 0;

  // Subcolumn optical depth and cloud-top pressure
  boxtau  = C::R_UNDEF;
  boxptop = C::R_UNDEF;
  if (tau > tauchk && ptop > 0 && sunlit == 1) {
    boxtau  = tau;
    boxptop = ptop;
  }
}

// icarus_column for one point: the tau/CTP joint histogram (percent, indexed (itau, ip)
// with pressure bins in the Fortran order, from low to high pressure) and the total
// cloud area (percent). boxptop is in Pa. The other column means are not computed,
// since EAMxx does not use them.
template <typename S, typename D>
template <typename VT, typename VP, typename H>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::icarus_column(const int ncol, const int sunlit, const VT& boxtau,
                                   const VP& boxptop, const H& fq, S& totalcldarea)
{
  constexpr int ntau  = C::numISCCPTauBins;
  constexpr int npres = C::numISCCPPresBins;

  if (sunlit != 1) {
    for (int i = 0; i < ntau; ++i)
      for (int k = 0; k < npres; ++k)
        fq(i,k) = C::R_UNDEF;
    totalcldarea = C::R_UNDEF;
    return;
  }

  // Joint histogram, with cloud-top pressure in hPa
  for (int i = 0; i < ntau; ++i) {
    for (int k = 0; k < npres; ++k) {
      int count = 0;
      for (int n = 0; n < ncol; ++n) {
        const S p = boxptop(n)/S(100);
        if (boxtau(n) >= C::tau_binBounds(i)  && boxtau(n) < C::tau_binBounds(i+1) &&
            p         >= C::pres_binBounds(k) && p         < C::pres_binBounds(k+1)) {
          ++count;
        }
      }
      fq(i,k) = (S(count)/ncol)*S(100);
    }
  }

  // Cloud area; real(count(...))/ncol is single precision in the Fortran
  int count = 0;
  for (int n = 0; n < ncol; ++n) {
    if (boxtau(n) > C::isccp_tauchk && boxptop(n)/S(100) > 0 && boxtau(n) > C::isccp_taumin) {
      ++count;
    }
  }
  totalcldarea = S(static_cast<float>(count)/static_cast<float>(ncol))*S(100);
}

} // namespace cosp
} // namespace scream

#endif // SCREAM_COSP_ISCCP_IMPL_HPP
