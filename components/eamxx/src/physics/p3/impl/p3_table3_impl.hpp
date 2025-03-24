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
          const Spack& lamr, Table3& tab,
          const Smask& context)
{
  // find location in scaled mean size space
  const auto dum1 = (mu_r+1) / lamr;
  const auto dum1_lt = context && (dum1 <= sp(195.e-6));
  tab.dumii = 1;
  if (dum1_lt.any()) {
    ekat_masked_loop(dum1_lt, s) {
      const auto inv_dum3 = sp(0.1);
      auto rdumii = (dum1[s] * sp(1.e6) + 5) * inv_dum3;
      rdumii = ekat::impl::max(rdumii, sp( 1.0));
      rdumii = ekat::impl::min(rdumii, sp(20.0));
      Int dumii = rdumii;
      dumii = ekat::impl::max(dumii,  1);
      dumii = ekat::impl::min(dumii, 20);
      tab.rdumii[s] = rdumii;
      tab.dumii[s] = dumii;
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
      tab.rdumii[s] = rdumii;
      tab.dumii[s] = dumii;
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
    tab.rdumjj.set(context, rdumjj);
    tab.dumjj = 1;
    tab.dumjj.set(context, dumjj);
  }
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table (const view_2d_table& table,
               const Table3& tab3) {
  const auto rdumii_m_dumii = tab3.rdumii - Spack(tab3.dumii);
  const auto t_im1_jm1 = index(table, tab3.dumii-1, tab3.dumjj-1);
  // Linear interpolant.
  const auto dum1 = (t_im1_jm1 + rdumii_m_dumii *
                     (index(table, tab3.dumii, tab3.dumjj-1) - t_im1_jm1));
  const auto t_im1_j = index(table, tab3.dumii-1, tab3.dumjj);
  // Linear interpolant.
  const auto dum2 = (t_im1_j + rdumii_m_dumii *
                     (index(table, tab3.dumii, tab3.dumjj) - t_im1_j));
  // Linear interpolation in other direction to complete bilinear interpolant.
  return dum1 + (tab3.rdumjj - Spack(tab3.dumjj)) * (dum2 - dum1);
}

template <typename S, typename D>
void Functions<S,D>
::get_global_tables (view_2d_table& vn_table_vals, view_2d_table& vm_table_vals,
                     view_2d_table& revap_table_vals, view_1d_table& mu_r_table_vals,
                     view_dnu_table& dnu) {
  auto tables = p3_init();
  vn_table_vals = tables.vn_table_vals;
  vm_table_vals = tables.vm_table_vals;
  revap_table_vals = tables.revap_table_vals;
  mu_r_table_vals = tables.mu_r_table_vals;
  dnu = tables.dnu_table_vals;
}

} // namespace p3
} // namespace scream

#endif
