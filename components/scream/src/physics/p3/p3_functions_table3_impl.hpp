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

} // namespace p3
} // namespace scream

#endif
