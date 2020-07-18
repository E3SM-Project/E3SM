#ifndef P3_TABLE_ICE_IMPL_HPP
#define P3_TABLE_ICE_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

#include <fstream>

namespace scream {
namespace p3 {

/*
 * Implementation of p3 ice table functions. Clients should NOT #include
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
  scream_require_msg(version == "VERSION", "Bad " << filename << ", expected VERSION X.Y.Z header");
  scream_require_msg(version_val == P3C::p3_version, "Bad " << filename << ", expected version " << P3C::p3_version << ", but got " << version_val);

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
::lookup_ice (const Spack& qitot, const Spack& nitot,
              const Spack& qirim, const Spack& rhop, TableIce& t,
              const Smask& context)
{
  // find index for qi (normalized ice mass mixing ratio = qitot/nitot)
  //   dum1 = (log10(qitot)+16.)/0.70757  !orig
  //   dum1 = (log10(qitot)+16.)*1.41328
  // we are inverting this equation from the lookup table to solve for i:
  // qitot/nitot=261.7**((i+10)*0.1)*1.e-18
  // dum1 = (log10(qitot/nitot)+18.)/(0.1*log10(261.7))-10.

  if (!context.any()) return;

  const auto lookup_table_1a_dum1_c = P3C::lookup_table_1a_dum1_c;
  t.dum1 = (pack::log10(qitot/nitot)+18) * lookup_table_1a_dum1_c - 10;
  t.dumi = IntSmallPack(t.dum1);

  // set limits (to make sure the calculated index doesn't exceed range of lookup table)
  t.dum1 = pack::min(t.dum1, static_cast<Scalar>(P3C::isize));
  t.dum1 = pack::max(t.dum1, sp(1.));
  t.dumi = pack::max(1, t.dumi);
  t.dumi = pack::min(P3C::isize-1, t.dumi);

  // find index for rime mass fraction
  t.dum4  = (qirim/qitot)*3 + 1;
  t.dumii = IntSmallPack(t.dum4);

  // set limits
  t.dum4  = pack::min(t.dum4, static_cast<Scalar>(P3C::rimsize));
  t.dum4  = pack::max(t.dum4, sp(1.));
  t.dumii = pack::max(1, t.dumii);
  t.dumii = pack::min(P3C::rimsize-1, t.dumii);

  // find index for bulk rime density
  // (account for uneven spacing in lookup table for density)
  const auto rhop_le_650 = context && (rhop <= 650);
  const auto rhop_gt_650 = !rhop_le_650 && context;
  t.dum5.set(rhop_le_650, (rhop-50)*sp(0.005) + 1);
  t.dum5.set(rhop_gt_650, (rhop-650)*sp(0.004) + 4);

  // set limits
  t.dumjj = IntSmallPack(t.dum5);
  t.dum5  = pack::min(t.dum5, static_cast<Scalar>(P3C::densize));
  t.dum5  = pack::max(t.dum5, sp(1.));
  t.dumjj = pack::max(1, t.dumjj);
  t.dumjj = pack::min(P3C::densize-1, t.dumjj);

  t.dum6  = -99;
  t.dumzz = -99;

  // adjust for 0-based indexing
  t.dumi  -= 1;
  t.dumjj -= 1;
  t.dumii -= 1;
  t.dumzz -= 1;
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::lookup_rain(const Spack& qr, const Spack& nr, TableRain& t,
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
    dumlr.set(gt_small, pack::cbrt(qr/(C::Pi * C::RHOW * nr)));
  }
  t.dum3 = (pack::log10(1*dumlr) + 5)*sp(10.70415);
  t.dumj = IntSmallPack(t.dum3);

  // set limits
  t.dum3 = pack::min(t.dum3, static_cast<Scalar>(P3C::rcollsize));
  t.dum3 = pack::max(t.dum3, 1);
  t.dumj = pack::max(1, t.dumj);
  t.dumj = pack::min(P3C::rcollsize-1, t.dumj);

  t.dumj.set(lt_small, 1);
  t.dum3.set(lt_small, 1);

  // adjust for 0-based indexing
  t.dumj -= 1;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table_ice(const int& index, const view_itab_table& itab, const TableIce& t,
                  const Smask& context)
{
  Spack proc;
  IntSmallPack idxpk(index);

  if (!context.any()) return proc;

  // get value at current density index

  // first interpolate for current rimed fraction index
  auto iproc1 = pack::index(itab, t.dumjj, t.dumii, t.dumi, idxpk) + (t.dum1-Spack(t.dumi)-1) *
    (pack::index(itab, t.dumjj, t.dumii, t.dumi+1, idxpk) - pack::index(itab, t.dumjj, t.dumii, t.dumi, idxpk));

  // linearly interpolate to get process rates for rimed fraction index + 1
  auto gproc1 = pack::index(itab, t.dumjj, t.dumii+1, t.dumi, idxpk) + (t.dum1-Spack(t.dumi)-1) *
    (pack::index(itab, t.dumjj, t.dumii+1, t.dumi+1, idxpk) - pack::index(itab, t.dumjj, t.dumii+1, t.dumi, idxpk));

  const auto tmp1   = iproc1 + (t.dum4-Spack(t.dumii)-1) * (gproc1-iproc1);

  // get value at density index + 1

  // first interpolate for current rimed fraction index

  iproc1 = pack::index(itab, t.dumjj+1, t.dumii, t.dumi, idxpk) + (t.dum1-Spack(t.dumi)-1) *
    (pack::index(itab, t.dumjj+1, t.dumii, t.dumi+1, idxpk) - pack::index(itab, t.dumjj+1, t.dumii, t.dumi, idxpk));

  // linearly interpolate to get process rates for rimed fraction index + 1

  gproc1 = pack::index(itab, t.dumjj+1, t.dumii+1, t.dumi, idxpk) + (t.dum1-Spack(t.dumi)-1) *
    (pack::index(itab, t.dumjj+1, t.dumii+1, t.dumi+1, idxpk)-pack::index(itab, t.dumjj+1, t.dumii+1, t.dumi, idxpk));

  const auto tmp2 = iproc1+(t.dum4 - Spack(t.dumii) - 1) * (gproc1-iproc1);

  // get final process rate
  proc = tmp1 + (t.dum5 - Spack(t.dumjj) - 1) * (tmp2-tmp1);
  return proc;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table_coll(const int& index, const view_itabcol_table& itabcoll,
                   const TableIce& ti, const TableRain& tr,
                   const Smask& context)
{
  Spack proc;
  IntSmallPack idxpk(index);

  if (!context.any()) return proc;

  // current density index

  // current rime fraction index
  auto dproc1  = pack::index(itabcoll, ti.dumjj, ti.dumii, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (pack::index(itabcoll, ti.dumjj, ti.dumii, ti.dumi+1, tr.dumj, idxpk) -
     pack::index(itabcoll, ti.dumjj, ti.dumii, ti.dumi, tr.dumj, idxpk));

  auto dproc2  = pack::index(itabcoll, ti.dumjj, ti.dumii, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1)*
    (pack::index(itabcoll, ti.dumjj, ti.dumii, ti.dumi+1, tr.dumj+1, idxpk) -
     pack::index(itabcoll, ti.dumjj, ti.dumii, ti.dumi, tr.dumj+1, idxpk));

  auto iproc1  = dproc1+(tr.dum3 - Spack(tr.dumj) - 1) * (dproc2 - dproc1);

  // rime fraction index + 1

  dproc1  = pack::index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (pack::index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi+1, tr.dumj, idxpk) -
     pack::index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj, idxpk));

  dproc2  = pack::index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (pack::index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi+1, tr.dumj+1, idxpk) -
     pack::index(itabcoll, ti.dumjj, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk));

  auto gproc1  = dproc1+(tr.dum3 - Spack(tr.dumj) - 1) * (dproc2-dproc1);
  const auto tmp1    = iproc1+(ti.dum4 - Spack(ti.dumii) - 1) * (gproc1-iproc1);

  // density index + 1

  // current rime fraction index

  dproc1  = pack::index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1)*
    (pack::index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi+1, tr.dumj, idxpk) -
     pack::index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj, idxpk));

  dproc2  = pack::index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (pack::index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi+1, tr.dumj+1, idxpk) -
     pack::index(itabcoll, ti.dumjj+1, ti.dumii, ti.dumi, tr.dumj+1, idxpk));

  iproc1  = dproc1 + (tr.dum3 - Spack(tr.dumj) - 1) * (dproc2-dproc1);

  // rime fraction index + 1

  dproc1  = pack::index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1)*
    (pack::index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi+1, tr.dumj, idxpk) -
     pack::index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj, idxpk));

  dproc2  = pack::index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk) +
    (ti.dum1 - Spack(ti.dumi) - 1) *
    (pack::index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi+1, tr.dumj+1, idxpk) -
     pack::index(itabcoll, ti.dumjj+1, ti.dumii+1, ti.dumi, tr.dumj+1, idxpk));

  gproc1  = dproc1 + (tr.dum3 - Spack(tr.dumj) - 1) * (dproc2-dproc1);
  const auto tmp2    = iproc1 + (ti.dum4 - Spack(ti.dumii) - 1) * (gproc1-iproc1);

  // interpolate over density to get final values
  proc    = tmp1+(ti.dum5 - Spack(ti.dumjj) - 1) * (tmp2-tmp1);
  return proc;
}

} // namespace p3
} // namespace scream

#endif
