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
::init_kokkos_ice_lookup_tables(view_ice_table& ice_table_vals, view_collect_table& collect_table_vals) {

  using DeviceIcetable = typename view_ice_table::non_const_type;
  using DeviceColtable = typename view_collect_table::non_const_type;

  const auto ice_table_vals_d     = DeviceIcetable("ice_table_vals");
  const auto collect_table_vals_d = DeviceColtable("collect_table_vals");

  const auto ice_table_vals_h    = Kokkos::create_mirror_view(ice_table_vals_d);
  const auto collect_table_vals_h = Kokkos::create_mirror_view(collect_table_vals_d);

  //
  // read in ice microphysics table into host views
  //

  std::string filename = std::string(P3C::p3_lookup_base) + std::string(P3C::p3_version);

  std::ifstream in(filename);

  // read header
  std::string version, version_val;
  in >> version >> version_val;
  EKAT_REQUIRE_MSG(version == "VERSION", "Bad " << filename << ", expected VERSION X.Y.Z header");
  EKAT_REQUIRE_MSG(version_val == P3C::p3_version, "Bad " << filename << ", expected version " << P3C::p3_version << ", but got " << version_val);

  // read tables
  double dum_s; int dum_i; // dum_s needs to be double to stream correctly
  for (int jj = 0; jj < P3C::densize; ++jj) {
    for (int ii = 0; ii < P3C::rimsize; ++ii) {
      for (int i = 0; i < P3C::isize; ++i) {
        in >> dum_i >> dum_i;
        int j_idx = 0;
        for (int j = 0; j < 15; ++j) {
          in >> dum_s;
          if (j > 1 && j != 10) {
            ice_table_vals_h(jj, ii, i, j_idx++) = dum_s;
          }
        }
      }

      for (int i = 0; i < P3C::isize; ++i) {
        for (int j = 0; j < P3C::rcollsize; ++j) {
          in >> dum_i >> dum_i;
          int k_idx = 0;
          for (int k = 0; k < 6; ++k) {
            in >> dum_s;
            if (k == 3 || k == 4) {
              collect_table_vals_h(jj, ii, i, j, k_idx++) = std::log10(dum_s);
            }
          }
        }
      }
    }
  }

  // deep copy to device
  Kokkos::deep_copy(ice_table_vals_d, ice_table_vals_h);
  Kokkos::deep_copy(collect_table_vals_d, collect_table_vals_h);
  ice_table_vals    = ice_table_vals_d;
  collect_table_vals = collect_table_vals_d;
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
