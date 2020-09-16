#ifndef P3_TABLE3_IMPL_HPP
#define P3_TABLE3_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of table functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::lookup (const Spack& mu_r,
          const Spack& lamr, Table3& t,
          const Smask& context)
{
  // find location in scaled mean size space
  const auto dum1 = (mu_r+1) / lamr;
  const auto dum1_lt = context && (dum1 <= sp(195.e-6));
  t.dumii = 1;
  if (dum1_lt.any()) {
    ekat_masked_loop(dum1_lt, s) {
      const auto inv_dum3 = sp(0.1);
      auto rdumii = (dum1[s] * sp(1.e6) + 5) * inv_dum3;
      rdumii = ekat::impl::max(rdumii, sp( 1.0));
      rdumii = ekat::impl::min(rdumii, sp(20.0));
      Int dumii = rdumii;
      dumii = ekat::impl::max(dumii,  1);
      dumii = ekat::impl::min(dumii, 20);
      t.rdumii[s] = rdumii;
      t.dumii[s] = dumii;
    }
  }
  const auto dum1_gte = context && !dum1_lt;
  if (dum1_gte.any()) {
    ekat_masked_loop(dum1_gte, s) {
      const auto inv_dum3 = C::THIRD * sp(0.1);
      auto rdumii = (dum1[s] * sp(1.e+6) - 195) * inv_dum3 + 20;
      rdumii = ekat::impl::max(rdumii,sp( 20.0));
      rdumii = ekat::impl::min(rdumii,sp(300.0));
      Int dumii = rdumii;
      dumii = ekat::impl::max(dumii, 20);
      dumii = ekat::impl::min(dumii,299);
      t.rdumii[s] = rdumii;
      t.dumii[s] = dumii;
    }
  }

  // find location in mu_r space
  {
    auto rdumjj = mu_r + 1;
    rdumjj = max(rdumjj, 1);
    rdumjj = min(rdumjj, 10);
    IntSmallPack dumjj(rdumjj);
    dumjj  = max(dumjj, 1);
    dumjj  = min(dumjj, 9);
    t.rdumjj.set(context, rdumjj);
    t.dumjj = 1;
    t.dumjj.set(context, dumjj);
  }
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table (const view_2d_table& table,
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
::init_kokkos_tables (view_2d_table& vn_table, view_2d_table& vm_table,
                      view_2d_table& revap_table, view_1d_table& mu_r_table,
                      view_dnu_table& dnu) {
  // initialize on host

  using DeviceTable1   = typename view_1d_table::non_const_type;
  using DeviceTable2   = typename view_2d_table::non_const_type;
  using DeviceDnuTable = typename view_dnu_table::non_const_type;

  const auto vn_table_d    = DeviceTable2("vn_table");
  const auto vm_table_d    = DeviceTable2("vm_table");
  const auto revap_table_d = DeviceTable2("revap_table");
  const auto mu_r_table_d  = DeviceTable1("mu_r_table");
  const auto dnu_table_d   = DeviceDnuTable("dnu");
  const auto vn_table_h    = Kokkos::create_mirror_view(vn_table_d);
  const auto vm_table_h    = Kokkos::create_mirror_view(vm_table_d);
  const auto revap_table_h = Kokkos::create_mirror_view(revap_table_d);
  const auto mu_table_h    = Kokkos::create_mirror_view(mu_r_table_d);
  const auto dnu_table_h   = Kokkos::create_mirror_view(dnu_table_d);

  // Need 2d-tables with fortran-style layout
  using P3F         = Functions<Real, HostDevice>;
  using LHostTable2 = typename P3F::KT::template lview<Real[C::VTABLE_DIM0][C::VTABLE_DIM1]>;
  LHostTable2 vn_table_lh("vn_table_lh"), vm_table_lh("vm_table_lh"), revap_table_lh("revap_table_lh");
  init_tables_from_f90_c(vn_table_lh.data(), vm_table_lh.data(), revap_table_lh.data(), mu_table_h.data());
  for (int i = 0; i < C::VTABLE_DIM0; ++i) {
    for (int j = 0; j < C::VTABLE_DIM1; ++j) {
      vn_table_h(i, j) = vn_table_lh(i, j);
      vm_table_h(i, j) = vm_table_lh(i, j);
      revap_table_h(i, j) = revap_table_lh(i, j);
    }
  }

  dnu_table_h(0)  =  0.000;
  dnu_table_h(1)  = -0.557;
  dnu_table_h(2)  = -0.430;
  dnu_table_h(3)  = -0.307;
  dnu_table_h(4)  = -0.186;
  dnu_table_h(5)  = -0.067;
  dnu_table_h(6)  = -0.050;
  dnu_table_h(7)  = -0.167;
  dnu_table_h(8)  = -0.282;
  dnu_table_h(9)  = -0.397;
  dnu_table_h(10) = -0.512;
  dnu_table_h(11) = -0.626;
  dnu_table_h(12) = -0.739;
  dnu_table_h(13) = -0.853;
  dnu_table_h(14) = -0.966;
  dnu_table_h(15) = -0.966;

  // deep copy to device
  Kokkos::deep_copy(vn_table_d, vn_table_h);
  Kokkos::deep_copy(vm_table_d, vm_table_h);
  Kokkos::deep_copy(revap_table_d, revap_table_h);
  Kokkos::deep_copy(mu_r_table_d, mu_table_h);
  Kokkos::deep_copy(dnu_table_d, dnu_table_h);
  vn_table   = vn_table_d;
  vm_table   = vm_table_d;
  revap_table   = revap_table_d;
  mu_r_table = mu_r_table_d;
  dnu        = dnu_table_d;
}

} // namespace p3
} // namespace scream

#endif
