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
::init_kokkos_ice_lookup_tables(view_itab_table& itab, view_itabcol_table& itabcol) {

  using DeviceItab    = typename view_itab_table::non_const_type;
  using DeviceItabcol = typename view_itabcol_table::non_const_type;

  const auto itab_d    = DeviceItab("itab");
  const auto itabcol_d = DeviceItabcol("itabcol");

  const auto itab_h    = Kokkos::create_mirror_view(itab_d);
  const auto itabcol_h = Kokkos::create_mirror_view(itabcol_d);

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
            itab_h(jj, ii, i, j_idx++) = dum_s;
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
              itabcol_h(jj, ii, i, j, k_idx++) = std::log10(dum_s);
            }
          }
        }
      }
    }
  }

  // deep copy to device
  Kokkos::deep_copy(itab_d, itab_h);
  Kokkos::deep_copy(itabcol_d, itabcol_h);
  itab    = itab_d;
  itabcol = itabcol_d;
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::lookup_ice (const Spack& qi, const Spack& ni,
              const Spack& qm, const Spack& rhop, TableIce& T_atm,
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
  T_atm.dum1 = (log10(qi/ni)+18) * lookup_table_1a_dum1_c - 10;
  T_atm.dumi = IntSmallPack(T_atm.dum1);

  // set limits (to make sure the calculated index doesn't exceed range of lookup table)
  T_atm.dum1 = min(T_atm.dum1, static_cast<Scalar>(P3C::isize));
  T_atm.dum1 = max(T_atm.dum1, sp(1.));
  T_atm.dumi = max(1, T_atm.dumi);
  T_atm.dumi = min(P3C::isize-1, T_atm.dumi);

  // find index for rime mass fraction
  T_atm.dum4  = (qm/qi)*3 + 1;
  T_atm.dumii = IntSmallPack(T_atm.dum4);

  // set limits
  T_atm.dum4  = min(T_atm.dum4, static_cast<Scalar>(P3C::rimsize));
  T_atm.dum4  = max(T_atm.dum4, sp(1.));
  T_atm.dumii = max(1, T_atm.dumii);
  T_atm.dumii = min(P3C::rimsize-1, T_atm.dumii);

  // find index for bulk rime density
  // (account for uneven spacing in lookup table for density)
  const auto rhop_le_650 = context && (rhop <= 650);
  const auto rhop_gt_650 = !rhop_le_650 && context;
  T_atm.dum5.set(rhop_le_650, (rhop-50)*sp(0.005) + 1);
  T_atm.dum5.set(rhop_gt_650, (rhop-650)*sp(0.004) + 4);

  // set limits
  T_atm.dumjj = IntSmallPack(T_atm.dum5);
  T_atm.dum5  = min(T_atm.dum5, static_cast<Scalar>(P3C::densize));
  T_atm.dum5  = max(T_atm.dum5, sp(1.));
  T_atm.dumjj = max(1, T_atm.dumjj);
  T_atm.dumjj = min(P3C::densize-1, T_atm.dumjj);

  T_atm.dum6  = -99;
  T_atm.dumzz = -99;

  // adjust for 0-based indexing
  T_atm.dumi  -= 1;
  T_atm.dumjj -= 1;
  T_atm.dumii -= 1;
  T_atm.dumzz -= 1;
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::lookup_rain(const Spack& qr, const Spack& nr, TableRain& T_atm,
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
  T_atm.dum3 = (log10(1*dumlr) + 5)*sp(10.70415);
  T_atm.dumj = IntSmallPack(T_atm.dum3);

  // set limits
  T_atm.dum3 = min(T_atm.dum3, static_cast<Scalar>(P3C::rcollsize));
  T_atm.dum3 = max(T_atm.dum3, 1);
  T_atm.dumj = max(1, T_atm.dumj);
  T_atm.dumj = min(P3C::rcollsize-1, T_atm.dumj);

  T_atm.dumj.set(lt_small, 1);
  T_atm.dum3.set(lt_small, 1);

  // adjust for 0-based indexing
  T_atm.dumj -= 1;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table_ice(const int& idx, const view_itab_table& itab, const TableIce& T_atm,
                  const Smask& context)
{
  using ekat::pack::index;

  Spack proc;
  IntSmallPack idxpk(idx);

  if (!context.any()) return proc;

  // get value at current density index

  // first interpolate for current rimed fraction index
  auto iproc1 = index(itab, T_atm.dumjj, T_atm.dumii, T_atm.dumi, idxpk) + (T_atm.dum1-Spack(T_atm.dumi)-1) *
    (index(itab, T_atm.dumjj, T_atm.dumii, T_atm.dumi+1, idxpk) - index(itab, T_atm.dumjj, T_atm.dumii, T_atm.dumi, idxpk));

  // linearly interpolate to get process rates for rimed fraction index + 1
  auto gproc1 = index(itab, T_atm.dumjj, T_atm.dumii+1, T_atm.dumi, idxpk) + (T_atm.dum1-Spack(T_atm.dumi)-1) *
    (index(itab, T_atm.dumjj, T_atm.dumii+1, T_atm.dumi+1, idxpk) - index(itab, T_atm.dumjj, T_atm.dumii+1, T_atm.dumi, idxpk));

  const auto tmp1   = iproc1 + (T_atm.dum4-Spack(T_atm.dumii)-1) * (gproc1-iproc1);

  // get value at density index + 1

  // first interpolate for current rimed fraction index

  iproc1 = index(itab, T_atm.dumjj+1, T_atm.dumii, T_atm.dumi, idxpk) + (T_atm.dum1-Spack(T_atm.dumi)-1) *
    (index(itab, T_atm.dumjj+1, T_atm.dumii, T_atm.dumi+1, idxpk) - index(itab, T_atm.dumjj+1, T_atm.dumii, T_atm.dumi, idxpk));

  // linearly interpolate to get process rates for rimed fraction index + 1

  gproc1 = index(itab, T_atm.dumjj+1, T_atm.dumii+1, T_atm.dumi, idxpk) + (T_atm.dum1-Spack(T_atm.dumi)-1) *
    (index(itab, T_atm.dumjj+1, T_atm.dumii+1, T_atm.dumi+1, idxpk)-index(itab, T_atm.dumjj+1, T_atm.dumii+1, T_atm.dumi, idxpk));

  const auto tmp2 = iproc1+(T_atm.dum4 - Spack(T_atm.dumii) - 1) * (gproc1-iproc1);

  // get final process rate
  proc = tmp1 + (T_atm.dum5 - Spack(T_atm.dumjj) - 1) * (tmp2-tmp1);
  return proc;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table_coll(const int& idx, const view_itabcol_table& itabcoll,
                   const TableIce& ti, const TableRain& tr,
                   const Smask& context)
{
  using ekat::pack::index;

  Spack proc;
  IntSmallPack idxpk(idx);

  if (!context.any()) return proc;

  // current density index

  // current rime fraction index
  auto dproc1  = index(itabcoll, ti.dumjj, ti.dumii, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (index(itabcoll, ti.dumjj, ti.dumii, ti.dumi+1, tr.dumj, idxpk) -
     index(itabcoll, ti.dumjj, ti.dumii, ti.dumi, tr.dumj, idxpk));

  auto dproc2  = index(itabcoll, ti.dumjj, ti.dumii, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1)*
    (index(itabcoll, ti.dumjj, ti.dumii, ti.dumi+1, tr.dumj+1, idxpk) -
     index(itabcoll, ti.dumjj, ti.dumii, ti.dumi, tr.dumj+1, idxpk));

  auto iproc1  = dproc1+(tr.dum3 - Spack(tr.dumj) - 1) * (dproc2 - dproc1);

  // rime fraction index + 1

  dproc1  = index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi+1, tr.dumj, idxpk) -
     index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj, idxpk));

  dproc2  = index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi+1, tr.dumj+1, idxpk) -
     index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk));

  auto gproc1  = dproc1+(tr.dum3 - Spack(tr.dumj) - 1) * (dproc2-dproc1);
  const auto tmp1    = iproc1+(ti.dum4 - Spack(ti.dumii) - 1) * (gproc1-iproc1);

  // density index + 1

  // current rime fraction index

  dproc1  = index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1)*
    (index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi+1, tr.dumj, idxpk) -
     index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj, idxpk));

  dproc2  = index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi+1, tr.dumj+1, idxpk) -
     index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj+1, idxpk));

  iproc1  = dproc1 + (tr.dum3 - Spack(tr.dumj) - 1) * (dproc2-dproc1);

  // rime fraction index + 1

  dproc1  = index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1)*
    (index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi+1, tr.dumj, idxpk) -
     index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj, idxpk));

  dproc2  = index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi+1, tr.dumj+1, idxpk) -
     index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk));

  gproc1  = dproc1 + (tr.dum3 - Spack(tr.dumj) - 1) * (dproc2-dproc1);
  const auto tmp2    = iproc1 + (ti.dum4 - Spack(ti.dumii) - 1) * (gproc1-iproc1);

  // interpolate over density to get final values
  proc    = tmp1+(ti.dum5 - Spack(ti.dumjj) - 1) * (tmp2-tmp1);
  return proc;
}

} // namespace p3
} // namespace scream

#endif // P3_TABLE_ICE_IMPL_HPP
