#ifndef P3_TABLE_ICE_IMPL_HPP
#define P3_TABLE_ICE_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

#include <fstream>

namespace scream {
namespace p3 {

/*
 * Implementation of ice table functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
void Functions<S,D>
::get_global_ice_lookup_tables(view_ice_table& ice_table_vals, view_collect_table& collect_table_vals) {
  auto tables = p3_init();
  ice_table_vals = tables.ice_table_vals;
  collect_table_vals = tables.collect_table_vals;
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::lookup_ice (const Spack& qi, const Spack& ni,
              const Spack& qm, const Spack& rhop, TableIce& tab,
              const Smask& context)
{
  // find index for qi (normalized ice mass mixing ratio = qi/ni)
  //   dum1 = (log10(qi)+16.)/0.70757  !orig
  //   dum1 = (log10(qi)+16.)*1.41328
  // we are inverting this equation from the lookup table to solve for i:
  // qi/ni=261.7**((i+10)*0.1)*1.e-18
  // dum1 = (log10(qi/ni)+18.)/(0.1*log10(261.7))-10.

  if (!context.any()) return;

  const auto lookup_table_1a_dum1_c = P3C::lookup_table_1a_dum1_c;
  tab.dum1 = (log10(qi/ni)+18) * lookup_table_1a_dum1_c - 10;
  tab.dumi = IntSmallPack(tab.dum1);

  // set limits (to make sure the calculated index doesn't exceed range of lookup table)
  tab.dum1 = min(tab.dum1, static_cast<Scalar>(P3C::isize));
  tab.dum1 = max(tab.dum1, sp(1.));
  tab.dumi = max(1, tab.dumi);
  tab.dumi = min(P3C::isize-1, tab.dumi);

  // find index for rime mass fraction
  tab.dum4  = (qm/qi)*3 + 1;
  tab.dumii = IntSmallPack(tab.dum4);

  // set limits
  tab.dum4  = min(tab.dum4, static_cast<Scalar>(P3C::rimsize));
  tab.dum4  = max(tab.dum4, sp(1.));
  tab.dumii = max(1, tab.dumii);
  tab.dumii = min(P3C::rimsize-1, tab.dumii);

  // find index for bulk rime density
  // (account for uneven spacing in lookup table for density)
  const auto rhop_le_650 = context && (rhop <= 650);
  const auto rhop_gt_650 = !rhop_le_650 && context;
  tab.dum5.set(rhop_le_650, (rhop-50)*sp(0.005) + 1);
  tab.dum5.set(rhop_gt_650, (rhop-650)*sp(0.004) + 4);

  // set limits
  tab.dumjj = IntSmallPack(tab.dum5);
  tab.dum5  = min(tab.dum5, static_cast<Scalar>(P3C::densize));
  tab.dum5  = max(tab.dum5, sp(1.));
  tab.dumjj = max(1, tab.dumjj);
  tab.dumjj = min(P3C::densize-1, tab.dumjj);

  tab.dum6  = -99;
  tab.dumzz = -99;

  // adjust for 0-based indexing
  tab.dumi  -= 1;
  tab.dumjj -= 1;
  tab.dumii -= 1;
  tab.dumzz -= 1;
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::lookup_rain(const Spack& qr, const Spack& nr, TableRain& tab,
              const Smask& context)
{
  if (!context.any()) return;

  // find index for scaled mean rain size
  // if no rain, then just choose dumj = 1 and do not calculate rain-ice collection processes
  const auto qsmall = C::QSMALL;
  const auto gt_small = qr > qsmall && nr > 0.0 && context;
  const auto lt_small = !gt_small && context;

  // calculate scaled mean size for consistency with ice lookup table
  Spack dumlr(1);
  if (gt_small.any()) {
    dumlr.set(gt_small, cbrt(qr/(C::Pi * C::RHO_H2O * nr)));
  }
  tab.dum3 = (log10(1*dumlr) + 5)*sp(10.70415);
  tab.dumj = IntSmallPack(tab.dum3);

  // set limits
  tab.dum3 = min(tab.dum3, static_cast<Scalar>(P3C::rcollsize));
  tab.dum3 = max(tab.dum3, 1);
  tab.dumj = max(1, tab.dumj);
  tab.dumj = min(P3C::rcollsize-1, tab.dumj);

  tab.dumj.set(lt_small, 1);
  tab.dum3.set(lt_small, 1);

  // adjust for 0-based indexing
  tab.dumj -= 1;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table_ice(const int& idx, const view_ice_table& ice_table_vals, const TableIce& tab,
                  const Smask& context)
{
  using ekat::index;

  Spack proc;
  IntSmallPack idxpk(idx);

  if (!context.any()) return proc;

  // get value at current density index

  // first interpolate for current rimed fraction index
  auto iproc1 = index(ice_table_vals, tab.dumjj, tab.dumii, tab.dumi, idxpk) + (tab.dum1-Spack(tab.dumi)-1) *
    (index(ice_table_vals, tab.dumjj, tab.dumii, tab.dumi+1, idxpk) - index(ice_table_vals, tab.dumjj, tab.dumii, tab.dumi, idxpk));

  // linearly interpolate to get process rates for rimed fraction index + 1
  auto gproc1 = index(ice_table_vals, tab.dumjj, tab.dumii+1, tab.dumi, idxpk) + (tab.dum1-Spack(tab.dumi)-1) *
    (index(ice_table_vals, tab.dumjj, tab.dumii+1, tab.dumi+1, idxpk) - index(ice_table_vals, tab.dumjj, tab.dumii+1, tab.dumi, idxpk));

  const auto tmp1   = iproc1 + (tab.dum4-Spack(tab.dumii)-1) * (gproc1-iproc1);

  // get value at density index + 1

  // first interpolate for current rimed fraction index

  iproc1 = index(ice_table_vals, tab.dumjj+1, tab.dumii, tab.dumi, idxpk) + (tab.dum1-Spack(tab.dumi)-1) *
    (index(ice_table_vals, tab.dumjj+1, tab.dumii, tab.dumi+1, idxpk) - index(ice_table_vals, tab.dumjj+1, tab.dumii, tab.dumi, idxpk));

  // linearly interpolate to get process rates for rimed fraction index + 1

  gproc1 = index(ice_table_vals, tab.dumjj+1, tab.dumii+1, tab.dumi, idxpk) + (tab.dum1-Spack(tab.dumi)-1) *
    (index(ice_table_vals, tab.dumjj+1, tab.dumii+1, tab.dumi+1, idxpk)-index(ice_table_vals, tab.dumjj+1, tab.dumii+1, tab.dumi, idxpk));

  const auto tmp2 = iproc1+(tab.dum4 - Spack(tab.dumii) - 1) * (gproc1-iproc1);

  // get final process rate
  proc = tmp1 + (tab.dum5 - Spack(tab.dumjj) - 1) * (tmp2-tmp1);
  return proc;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table_coll(const int& idx, const view_collect_table& collect_table_vals,
                   const TableIce& ti, const TableRain& tr,
                   const Smask& context)
{
  using ekat::index;

  Spack proc;
  IntSmallPack idxpk(idx);

  if (!context.any()) return proc;

  // current density index

  // current rime fraction index
  auto dproc1  = index(collect_table_vals, ti.dumjj, ti.dumii, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (index(collect_table_vals, ti.dumjj, ti.dumii, ti.dumi+1, tr.dumj, idxpk) -
     index(collect_table_vals, ti.dumjj, ti.dumii, ti.dumi, tr.dumj, idxpk));

  auto dproc2  = index(collect_table_vals, ti.dumjj, ti.dumii, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1)*
    (index(collect_table_vals, ti.dumjj, ti.dumii, ti.dumi+1, tr.dumj+1, idxpk) -
     index(collect_table_vals, ti.dumjj, ti.dumii, ti.dumi, tr.dumj+1, idxpk));

  auto iproc1  = dproc1+(tr.dum3 - Spack(tr.dumj) - 1) * (dproc2 - dproc1);

  // rime fraction index + 1

  dproc1  = index(collect_table_vals, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (index(collect_table_vals, ti.dumjj, ti.dumii+1, ti.dumi+1, tr.dumj, idxpk) -
     index(collect_table_vals, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj, idxpk));

  dproc2  = index(collect_table_vals, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (index(collect_table_vals, ti.dumjj, ti.dumii+1, ti.dumi+1, tr.dumj+1, idxpk) -
     index(collect_table_vals, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk));

  auto gproc1  = dproc1+(tr.dum3 - Spack(tr.dumj) - 1) * (dproc2-dproc1);
  const auto tmp1    = iproc1+(ti.dum4 - Spack(ti.dumii) - 1) * (gproc1-iproc1);

  // density index + 1

  // current rime fraction index

  dproc1  = index(collect_table_vals, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1)*
    (index(collect_table_vals, ti.dumjj+1, ti.dumii, ti.dumi+1, tr.dumj, idxpk) -
     index(collect_table_vals, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj, idxpk));

  dproc2  = index(collect_table_vals, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (index(collect_table_vals, ti.dumjj+1, ti.dumii, ti.dumi+1, tr.dumj+1, idxpk) -
     index(collect_table_vals, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj+1, idxpk));

  iproc1  = dproc1 + (tr.dum3 - Spack(tr.dumj) - 1) * (dproc2-dproc1);

  // rime fraction index + 1

  dproc1  = index(collect_table_vals, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1)*
    (index(collect_table_vals, ti.dumjj+1, ti.dumii+1, ti.dumi+1, tr.dumj, idxpk) -
     index(collect_table_vals, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj, idxpk));

  dproc2  = index(collect_table_vals, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (index(collect_table_vals, ti.dumjj+1, ti.dumii+1, ti.dumi+1, tr.dumj+1, idxpk) -
     index(collect_table_vals, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk));

  gproc1  = dproc1 + (tr.dum3 - Spack(tr.dumj) - 1) * (dproc2-dproc1);
  const auto tmp2    = iproc1 + (ti.dum4 - Spack(ti.dumii) - 1) * (gproc1-iproc1);

  // interpolate over density to get final values
  proc    = tmp1+(ti.dum5 - Spack(ti.dumjj) - 1) * (tmp2-tmp1);
  return proc;
}

} // namespace p3
} // namespace scream

#endif // P3_TABLE_ICE_IMPL_HPP
