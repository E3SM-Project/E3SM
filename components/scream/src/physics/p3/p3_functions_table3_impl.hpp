#ifndef P3_FUNCTIONS_TABLE3_IMPL_HPP
#define P3_FUNCTIONS_TABLE3_IMPL_HPP

#include <fstream>

#include "p3_functions.hpp"
#include "p3_constants.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 table functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::lookup (const Smask& qr_gt_small, const Spack& mu_r,
          const Spack& lamr, Table3& t) {
  // find location in scaled mean size space
  const auto dum1 = (mu_r+1.) / lamr;
  const auto dum1_lt = qr_gt_small && (dum1 <= 195.e-6);
  t.dumii = 1;
  if (dum1_lt.any()) {
    scream_masked_loop(dum1_lt, s) {
      const auto inv_dum3 = 0.1;
      auto rdumii = (dum1[s]*1.e6+5.)*inv_dum3;
      rdumii = util::max<Scalar>(rdumii,  1.);
      rdumii = util::min<Scalar>(rdumii, 20.);
      Int dumii = rdumii;
      dumii = util::max(dumii,  1);
      dumii = util::min(dumii, 20);
      t.rdumii[s] = rdumii;
      t.dumii[s] = dumii;
    }
  }
  const auto dum1_gte = qr_gt_small && ! dum1_lt;
  if (dum1_gte.any()) {
    scream_masked_loop(dum1_gte, s) {
      const auto inv_dum3 = C::THIRD*0.1;
      auto rdumii = (dum1[s]*1.e+6-195.)*inv_dum3 + 20.;
      rdumii = util::max<Scalar>(rdumii, 20.);
      rdumii = util::min<Scalar>(rdumii,300.);
      Int dumii = rdumii;
      dumii = util::max(dumii, 20);
      dumii = util::min(dumii,299);
      t.rdumii[s] = rdumii;
      t.dumii[s] = dumii;
    }
  }

  // find location in mu_r space
  {
    auto rdumjj = mu_r+1.;
    rdumjj = max(rdumjj,1.);
    rdumjj = min(rdumjj,10.);
    IntSmallPack dumjj(rdumjj);
    dumjj  = max(dumjj,1);
    dumjj  = min(dumjj,9);
    t.rdumjj.set(qr_gt_small, rdumjj);
    t.dumjj = 1;
    t.dumjj.set(qr_gt_small, dumjj);
  }
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table (const Smask& qr_gt_small, const view_2d_table& table,
               const Table3& t) {
  const auto rdumii_m_dumii = t.rdumii - Spack(t.dumii);
  const auto t_im1_jm1 = index(table, t.dumii-1, t.dumjj-1);
  // Linear interpolant.
  const auto dum1 = (t_im1_jm1 + rdumii_m_dumii *
                     (index(table, t.dumii, t.dumjj-1) - t_im1_jm1));
  const auto t_im1_j = index(table, t.dumii-1, t.dumjj);
  // Linear interpolant.
  const auto dum2 = (t_im1_j + rdumii_m_dumii *
                     (index(table, t.dumii, t.dumjj) - t_im1_j));
  // Linear interpolation in other direction to complete bilinear interpolant.
  return dum1 + (t.rdumjj - Spack(t.dumjj)) * (dum2 - dum1);
}

template <typename S, typename D>
void Functions<S,D>
::init_kokkos_tables (view_2d_table& vn_table, view_2d_table& vm_table, view_1d_table& mu_r_table) {
  // initialize on host

  using DeviceTable1 = typename view_1d_table::non_const_type;
  using DeviceTable2 = typename view_2d_table::non_const_type;

  const auto vn_table_d = DeviceTable2("vn_table");
  const auto vm_table_d = DeviceTable2("vm_table");
  const auto mu_r_table_d = DeviceTable1("mu_r_table");
  const auto vn_table_h = Kokkos::create_mirror_view(vn_table_d);
  const auto vm_table_h = Kokkos::create_mirror_view(vm_table_d);
  const auto mu_table_h = Kokkos::create_mirror_view(mu_r_table_d);

  scream_require(G::VN_TABLE.size()    == vn_table_h.extent(0) && G::VN_TABLE.size() > 0);
  scream_require(G::VN_TABLE[0].size() == vn_table_h.extent(1));
  scream_require(G::VM_TABLE.size()    == vm_table_h.extent(0) && G::VM_TABLE.size() > 0);
  scream_require(G::VM_TABLE[0].size() == vm_table_h.extent(1));
  scream_require(G::MU_R_TABLE.size()  == mu_table_h.extent(0));

  for (size_t i = 0; i < vn_table_h.extent(0); ++i) {
    for (size_t k = 0; k < vn_table_h.extent(1); ++k) {
      vn_table_h(i, k) = G::VN_TABLE[i][k];
      vm_table_h(i, k) = G::VM_TABLE[i][k];
    }
  }

  for (size_t i = 0; i < mu_table_h.extent(0); ++i) {
    mu_table_h(i) = G::MU_R_TABLE[i];
  }

  // deep copy to device
  Kokkos::deep_copy(vn_table_d, vn_table_h);
  Kokkos::deep_copy(vm_table_d, vm_table_h);
  Kokkos::deep_copy(mu_r_table_d, mu_table_h);
  vn_table = vn_table_d;
  vm_table = vm_table_d;
  mu_r_table = mu_r_table_d;
}

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
::lookup_ice (const Smask& qiti_gt_small, const Spack& qitot, const Spack& nitot,
              const Spack& qirim, const Spack& rhop, TableIce& t)
{
  // find index for qi (normalized ice mass mixing ratio = qitot/nitot)
  //   dum1 = (log10(qitot)+16.)/0.70757  !orig
  //   dum1 = (log10(qitot)+16.)*1.41328
  // we are inverting this equation from the lookup table to solve for i:
  // qitot/nitot=261.7**((i+10)*0.1)*1.e-18
  // dum1 = (log10(qitot/nitot)+18.)/(0.1*log10(261.7))-10.

  if (qiti_gt_small.any()) {
    scream_masked_loop(qiti_gt_small, s) {
      t.dum1[s] = (std::log10(qitot[s]/nitot[s])+18) * P3C::lookup_table_1a_dum1_c - 10; // For computational efficiency
      t.dumi[s] = static_cast<int>(t.dum1[s]);

      // set limits (to make sure the calculated index doesn't exceed range of lookup table)
      t.dum1[s] = util::min<Scalar>(t.dum1[s], P3C::isize);
      t.dum1[s] = util::max(t.dum1[s], static_cast<Scalar>(1.));
      t.dumi[s] = util::max(1, t.dumi[s]);
      t.dumi[s] = util::min(P3C::isize, t.dumi[s]);

      // find index for rime mass fraction
      t.dum4[s]  = (qirim[s]/qitot[s])*3 + 1;
      t.dumii[s] = static_cast<int>(t.dum4[s]);

      // set limits
      t.dum4[s]  = util::min<Scalar>(t.dum4[s], P3C::rimsize);
      t.dum4[s]  = util::max<Scalar>(t.dum4[s], static_cast<Scalar>(1.));
      t.dumii[s] = util::max(1, t.dumii[s]);
      t.dumii[s] = util::min(P3C::rimsize-1, t.dumii[s]);
    }

    // find index for bulk rime density
    // (account for uneven spacing in lookup table for density)
    const auto rhop_le_650 = qiti_gt_small && (rhop <= 650);
    const auto rhop_gt_650 = qiti_gt_small && (rhop > 650);
    scream_masked_loop(rhop_le_650, s) {
      t.dum5[s] = (rhop[s]-50)*0.005 + 1;
    }
    scream_masked_loop(rhop_gt_650, s) {
      t.dum5[s] = (rhop[s]-650)*0.004 + 4;
    }

    // set limits
    scream_masked_loop(qiti_gt_small, s) {
      t.dumjj[s] = static_cast<int>(t.dum5[s]);
      t.dum5[s]  = util::min<Scalar>(t.dum5[s], P3C::densize);
      t.dum5[s]  = util::max<Scalar>(t.dum5[s], 1);
      t.dumjj[s] = util::max(1, t.dumjj[s]);
      t.dumjj[s] = util::min(P3C::densize-1, t.dumjj[s]);

      t.dum6[s]  = -99;
      t.dumzz[s] = -99;
    }

    // adjust for 0-based indexing
    t.dumi  -= 1;
    t.dumjj -= 1;
    t.dumii -= 1;
    t.dumzz -= 1;
  }
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::lookup_rain(const Smask& qiti_gt_small, const Spack& qr, const Spack& nr, TableRain& t)
{
  if (qiti_gt_small.any()) {
    // find index for scaled mean rain size
    // if no rain, then just choose dumj = 1 and do not calculate rain-ice collection processes
    Smask gt_small = qr > C::QSMALL && nr > 0.0;
    Smask gt_both = gt_small && qiti_gt_small;
    scream_masked_loop(gt_both, s) {
      // calculate scaled mean size for consistency with ice lookup table
      Scalar dumlr = std::pow(qr[s]/(C::Pi * C::RHOW * nr[s]), C::THIRD);
      t.dum3[s] = (std::log10(1.0*dumlr) + 5.0)*10.70415;
      t.dumj[s] = static_cast<int>(t.dum3[s]);

      // set limits
      t.dum3[s] = util::min<Scalar>(t.dum3[s], P3C::rcollsize);
      t.dum3[s] = util::max(t.dum3[s], static_cast<Scalar>(1.));
      t.dumj[s] = util::max(1, t.dumj[s]);
      t.dumj[s] = util::min(P3C::rcollsize-1, t.dumj[s]);
    }

    t.dumj.set(qiti_gt_small && !gt_small, 1);
    t.dum3.set(qiti_gt_small && !gt_small, 1.);

    // adjust for 0-based indexing
    t.dumj -= 1;
  }
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table_ice(const Smask& qiti_gt_small, const int& index, const view_itab_table& itab,
                  const TableIce& t)
{
  Spack proc;
  scream_masked_loop(qiti_gt_small, s) {
    // get value at current density index

    // first interpolate for current rimed fraction index
    Scalar iproc1 = itab(t.dumjj[s], t.dumii[s], t.dumi[s], index) + (t.dum1[s]-t.dumi[s]-1) *
      (itab(t.dumjj[s], t.dumii[s], t.dumi[s]+1, index) - itab(t.dumjj[s], t.dumii[s], t.dumi[s], index));

    // linearly interpolate to get process rates for rimed fraction index + 1
    Scalar gproc1 = itab(t.dumjj[s], t.dumii[s]+1, t.dumi[s], index) + (t.dum1[s]-t.dumi[s]-1) *
      (itab(t.dumjj[s], t.dumii[s]+1, t.dumi[s]+1, index) - itab(t.dumjj[s], t.dumii[s]+1, t.dumi[s], index));

    Scalar tmp1   = iproc1 + (t.dum4[s]-t.dumii[s]-1) * (gproc1-iproc1);

    // get value at density index + 1

    // first interpolate for current rimed fraction index

    iproc1 = itab(t.dumjj[s]+1, t.dumii[s], t.dumi[s], index) + (t.dum1[s]-t.dumi[s]-1) *
      (itab(t.dumjj[s]+1, t.dumii[s], t.dumi[s]+1, index) - itab(t.dumjj[s]+1, t.dumii[s], t.dumi[s], index));

    // linearly interpolate to get process rates for rimed fraction index + 1

    gproc1 = itab(t.dumjj[s]+1, t.dumii[s]+1, t.dumi[s], index) + (t.dum1[s]-t.dumi[s]-1) *
      (itab(t.dumjj[s]+1, t.dumii[s]+1, t.dumi[s]+1, index)-itab(t.dumjj[s]+1, t.dumii[s]+1, t.dumi[s], index));

    Scalar tmp2 = iproc1+(t.dum4[s] - t.dumii[s] - 1) * (gproc1-iproc1);

    // get final process rate
    proc[s] = tmp1 + (t.dum5[s] - t.dumjj[s] - 1) * (tmp2-tmp1);
  }
  return proc;
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table_coll(const Smask& qiti_gt_small, const int& index, const view_itabcol_table& itabcoll,
                   const TableIce& ti, const TableRain& tr)
{
  Spack proc;
  scream_masked_loop(qiti_gt_small, s) {

    // current density index

    // current rime fraction index
    Scalar dproc1  = itabcoll(ti.dumjj[s], ti.dumii[s], ti.dumi[s], tr.dumj[s], index) +
      (ti.dum1[s] - ti.dumi[s] - 1) *
      (itabcoll(ti.dumjj[s], ti.dumii[s], ti.dumi[s]+1, tr.dumj[s], index) -
       itabcoll(ti.dumjj[s], ti.dumii[s], ti.dumi[s], tr.dumj[s], index));

    Scalar dproc2  = itabcoll(ti.dumjj[s], ti.dumii[s], ti.dumi[s], tr.dumj[s]+1, index) +
      (ti.dum1[s] - ti.dumi[s] - 1)*
      (itabcoll(ti.dumjj[s], ti.dumii[s], ti.dumi[s]+1, tr.dumj[s]+1, index) -
       itabcoll(ti.dumjj[s], ti.dumii[s], ti.dumi[s], tr.dumj[s]+1, index));

    Scalar iproc1  = dproc1+(tr.dum3[s] - tr.dumj[s] - 1) * (dproc2 - dproc1);

    // rime fraction index + 1

    dproc1  = itabcoll(ti.dumjj[s], ti.dumii[s]+1, ti.dumi[s], tr.dumj[s], index) +
      (ti.dum1[s] - ti.dumi[s] - 1) *
      (itabcoll(ti.dumjj[s], ti.dumii[s]+1, ti.dumi[s]+1, tr.dumj[s], index) -
       itabcoll(ti.dumjj[s], ti.dumii[s]+1, ti.dumi[s], tr.dumj[s], index));

    dproc2  = itabcoll(ti.dumjj[s], ti.dumii[s]+1, ti.dumi[s], tr.dumj[s]+1, index) +
      (ti.dum1[s] - ti.dumi[s] - 1) *
      (itabcoll(ti.dumjj[s], ti.dumii[s]+1, ti.dumi[s]+1, tr.dumj[s]+1, index) -
       itabcoll(ti.dumjj[s], ti.dumii[s]+1, ti.dumi[s], tr.dumj[s]+1, index));

    Scalar gproc1  = dproc1+(tr.dum3[s] - tr.dumj[s] - 1) * (dproc2-dproc1);
    Scalar tmp1    = iproc1+(ti.dum4[s] - ti.dumii[s] - 1) * (gproc1-iproc1);

    // density index + 1

    // current rime fraction index

    dproc1  = itabcoll(ti.dumjj[s]+1, ti.dumii[s], ti.dumi[s], tr.dumj[s], index) +
      (ti.dum1[s] - ti.dumi[s] - 1)*
      (itabcoll(ti.dumjj[s]+1, ti.dumii[s], ti.dumi[s]+1, tr.dumj[s], index) -
       itabcoll(ti.dumjj[s]+1, ti.dumii[s], ti.dumi[s], tr.dumj[s], index));

    dproc2  = itabcoll(ti.dumjj[s]+1, ti.dumii[s], ti.dumi[s], tr.dumj[s]+1, index) +
      (ti.dum1[s] - ti.dumi[s] - 1) *
      (itabcoll(ti.dumjj[s]+1, ti.dumii[s], ti.dumi[s]+1, tr.dumj[s]+1, index) -
       itabcoll(ti.dumjj[s]+1, ti.dumii[s], ti.dumi[s], tr.dumj[s]+1,index));

    iproc1  = dproc1 + (tr.dum3[s] - tr.dumj[s] - 1) * (dproc2-dproc1);

    // rime fraction index + 1

    dproc1  = itabcoll(ti.dumjj[s]+1, ti.dumii[s]+1, ti.dumi[s], tr.dumj[s], index) +
      (ti.dum1[s] - ti.dumi[s] - 1)*
      (itabcoll(ti.dumjj[s]+1, ti.dumii[s]+1, ti.dumi[s]+1, tr.dumj[s], index) -
       itabcoll(ti.dumjj[s]+1, ti.dumii[s]+1, ti.dumi[s], tr.dumj[s], index));

    dproc2  = itabcoll(ti.dumjj[s]+1, ti.dumii[s]+1, ti.dumi[s], tr.dumj[s]+1, index) +
      (ti.dum1[s] - ti.dumi[s] - 1) *
      (itabcoll(ti.dumjj[s]+1, ti.dumii[s]+1, ti.dumi[s]+1, tr.dumj[s]+1, index) -
       itabcoll(ti.dumjj[s]+1, ti.dumii[s]+1, ti.dumi[s], tr.dumj[s]+1, index));

    gproc1  = dproc1 + (tr.dum3[s] - tr.dumj[s] - 1) * (dproc2-dproc1);
    Scalar tmp2    = iproc1 + (ti.dum4[s] - ti.dumii[s] - 1) * (gproc1-iproc1);

    // interpolate over density to get final values
    proc[s]    = tmp1+(ti.dum5[s] - ti.dumjj[s] - 1) * (tmp2-tmp1);
  }
  return proc;
}

} // namespace p3
} // namespace scream

#endif
