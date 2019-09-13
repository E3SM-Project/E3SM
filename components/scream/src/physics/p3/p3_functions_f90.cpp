#include "p3_functions_f90.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"
#include "p3_f90.hpp"

using scream::Real;
using scream::Int;

//
// A C++ interface to micro_p3 fortran calls and vice versa
//


extern "C" {

void p3_init_a_c(Real* itab, Real* itabcol);

void find_lookuptable_indices_1a_c(Int* dumi, Int* dumjj, Int* dumii, Int* dumzz,
                                   Real* dum1, Real* dum4, Real* dum5, Real* dum6,
                                   Real qitot, Real nitot, Real qirim, Real rhop);

void find_lookuptable_indices_1b_c(Int* dumj, Real* dum3, Real qr, Real nr);

void access_lookup_table_c(Int dumjj, Int dumii, Int dumi, Int index,
                           Real dum1, Real dum4, Real dum5, Real* proc);

void access_lookup_table_coll_c(Int dumjj, Int dumii, Int dumj, Int dumi, Int index,
                                Real dum1, Real dum3, Real dum4, Real dum5, Real* proc);

void get_cloud_dsd2_c(Real qc, Real* nc, Real* mu_c, Real rho, Real* nu, Real* lamc,
                      Real* cdist, Real* cdist1, Real lcldm);

void get_rain_dsd2_c(Real qr, Real* nr, Real* mu_r, Real* lamr, Real* cdistr, Real* logn0r, Real rcldm);

}

namespace scream {
namespace p3 {

//
// In all C++ -> Fortran bridge functions you should see p3_init(true). P3 needs
// to be initialized since most of its function depend on global tables to be
// populated. The 'true' argument is to set p3 to use its fortran implementations
// instead of calling back to C++. We want this behavior since it doesn't make much
// sense for C++ to bridge over to fortran only to have fortran bridge back to C++.
// If the client wanted the C++ implementation, they should just call it directly.
//

void p3_init_a(P3InitAFortranData& d)
{
  p3_init(true); // need to initialize p3 first so that tables are loaded
  p3_init_a_c(d.itab.data(), d.itabcol.data());
}

void find_lookuptable_indices_1a(LookupIceData& d)
{
  p3_init(true); // need to initialize p3 first so that tables are loaded
  find_lookuptable_indices_1a_c(&d.dumi, &d.dumjj, &d.dumii, &d.dumzz,
                                &d.dum1, &d.dum4, &d.dum5, &d.dum6,
                                d.qitot, d.nitot, d.qirim, d.rhop);
}

void find_lookuptable_indices_1b(LookupIceDataB& d)
{
  p3_init(true);
  find_lookuptable_indices_1b_c(&d.dumj, &d.dum3, d.qr, d.nr);
}

void access_lookup_table(AccessLookupTableData& d)
{
  p3_init(true); // need to initialize p3 first so that tables are loaded
  access_lookup_table_c(d.lid.dumjj, d.lid.dumii, d.lid.dumi, d.index,
                        d.lid.dum1, d.lid.dum4, d.lid.dum5, &d.proc);
}

void access_lookup_table_coll(AccessLookupTableCollData& d)
{
  p3_init(true);
  access_lookup_table_coll_c(d.lid.dumjj, d.lid.dumii, d.lidb.dumj, d.lid.dumi, d.index,
                             d.lid.dum1, d.lidb.dum3, d.lid.dum4, d.lid.dum5, &d.proc);
}

void get_cloud_dsd2(GetCloudDsd2Data& d)
{
  p3_init(true);
  Real nc_in = d.nc_in;
  get_cloud_dsd2_c(d.qc, &nc_in, &d.mu_c, d.rho, &d.nu, &d.lamc, &d.cdist, &d.cdist1, d.lcldm);
  d.nc_out = nc_in;
}

void get_rain_dsd2(GetRainDsd2Data& d)
{
  p3_init(true);
  Real nr_in = d.nr_in;
  get_rain_dsd2_c(d.qr, &nr_in, &d.mu_r, &d.lamr, &d.cdistr, &d.logn0r, d.rcldm);
  d.nr_out = nr_in;
}

std::shared_ptr<P3GlobalForFortran::Views> P3GlobalForFortran::s_views;

const P3GlobalForFortran::Views& P3GlobalForFortran::get()
{
  if (!P3GlobalForFortran::s_views) {
    P3GlobalForFortran::s_views = std::make_shared<Views>();
    P3F::init_kokkos_ice_lookup_tables(s_views->m_itab, s_views->m_itabcol);
    P3F::init_kokkos_tables(s_views->m_vn_table, s_views->m_vm_table, s_views->m_mu_r_table, s_views->m_dnu);
  }
  return *P3GlobalForFortran::s_views;
}

void P3GlobalForFortran::deinit()
{
  P3GlobalForFortran::s_views = nullptr;
}

void find_lookuptable_indices_1a_f(Int* dumi, Int* dumjj, Int* dumii, Int* dumzz,
                                   Real* dum1, Real* dum4, Real* dum5, Real* dum6,
                                   Real qitot_, Real nitot_, Real qirim_, Real rhop_)
{
  using P3F = Functions<Real, HostDevice>;

  typename P3F::Smask qiti_gt_small(qitot_ > P3F::C::QSMALL);
  typename P3F::Spack qitot(qitot_), nitot(nitot_), qirim(qirim_), rhop(rhop_);
  typename P3F::TableIce t;
  P3F::lookup_ice(qiti_gt_small, qitot, nitot, qirim, rhop, t);

  // adjust for 1-based indexing
  *dumi  = t.dumi[0]  + 1;
  *dumjj = t.dumjj[0] + 1;
  *dumii = t.dumii[0] + 1;
  *dumzz = t.dumzz[0] + 1;

  *dum1 = t.dum1[0];
  *dum4 = t.dum4[0];
  *dum5 = t.dum5[0];
  *dum6 = t.dum6[0];
}

void find_lookuptable_indices_1b_f(Int* dumj, Real* dum3, Real qr_, Real nr_)
{
  using P3F = Functions<Real, HostDevice>;

  // we can assume fortran would not be calling this routine if qiti_gt_small was not true
  typename P3F::Smask qiti_gt_small(true);

  typename P3F::Spack qr(qr_), nr(nr_);
  typename P3F::TableRain t;
  P3F::lookup_rain(qiti_gt_small, qr, nr, t);

  // adjust for 1-based indexing
  *dumj = t.dumj[0] + 1;

  *dum3 = t.dum3[0];
}

void access_lookup_table_f(Int dumjj, Int dumii, Int dumi, Int index,
                           Real dum1, Real dum4, Real dum5, Real* proc)
{
  using P3F = Functions<Real, HostDevice>;

  // we can assume fortran would not be calling this routine if qiti_gt_small was not true
  typename P3F::Smask qiti_gt_small(true);
  typename P3F::TableIce t;

  // Adjust for 0-based indexing
  t.dumi  = dumi  - 1;
  t.dumjj = dumjj - 1;
  t.dumii = dumii - 1;

  int adjusted_index = index - 1;

  t.dum1 = dum1;
  t.dum4 = dum4;
  t.dum5 = dum5;

  *proc = P3F::apply_table_ice(qiti_gt_small, adjusted_index, P3GlobalForFortran::itab(), t)[0];
}

void access_lookup_table_coll_f(Int dumjj, Int dumii, Int dumj, Int dumi, Int index,
                                Real dum1, Real dum3, Real dum4, Real dum5, Real* proc)
{
  using P3F = Functions<Real, HostDevice>;

  // we can assume fortran would not be calling this routine if qiti_gt_small was not true
  typename P3F::Smask qiti_gt_small(true);

  typename P3F::TableIce ti;
  typename P3F::TableRain tr;

  // Adjust for 0-based indexing
  ti.dumi  = dumi  - 1;
  ti.dumjj = dumjj - 1;
  ti.dumii = dumii - 1;
  tr.dumj  = dumj  - 1;

  int adjusted_index = index - 1;

  ti.dum1 = dum1;
  ti.dum4 = dum4;
  ti.dum5 = dum5;
  tr.dum3 = dum3;

  *proc = P3F::apply_table_coll(qiti_gt_small, adjusted_index, P3GlobalForFortran::itabcol(), ti, tr)[0];
}

void get_cloud_dsd2_f(Real qc_, Real* nc_, Real* mu_c_, Real rho_, Real* nu_, Real* lamc_,
                      Real* cdist_, Real* cdist1_, Real lcldm_)
{
  using P3F = Functions<Real, HostDevice>;

  typename P3F::Smask qc_gt_small(qc_ > P3F::C::QSMALL);
  typename P3F::Spack qc(qc_), nc(*nc_), rho(rho_), lcldm(lcldm_);

  typename P3F::Spack mu_c, nu, lamc, cdist, cdist1;

  P3F::get_cloud_dsd2(qc_gt_small, qc, nc, mu_c, rho, nu, P3GlobalForFortran::dnu(), lamc, cdist, cdist1, lcldm);

  *nc_     = nc[0];
  *mu_c_   = mu_c[0];
  *nu_     = nu[0];
  *lamc_   = lamc[0];
  *cdist_  = cdist[0];
  *cdist1_ = cdist1[0];
}

void get_rain_dsd2_f(Real qr_, Real* nr_, Real* mu_r_, Real* lamr_, Real* cdistr_, Real* logn0r_, Real rcldm_)
{
  using P3F = Functions<Real, HostDevice>;

  typename P3F::Smask qr_gt_small(qr_ > P3F::C::QSMALL);
  typename P3F::Spack qr(qr_), rcldm(rcldm_), nr(*nr_);
  typename P3F::Spack lamr, mu_r, cdistr, logn0r;

  P3F::get_rain_dsd2(qr_gt_small, qr, nr, mu_r, lamr, cdistr, logn0r, rcldm);

  *nr_     = nr[0];
  *mu_r_   = mu_r[0];
  *lamr_   = lamr[0];
  *cdistr_ = cdistr[0];
  *logn0r_ = logn0r[0];
}

Real cxx_pow(Real base, Real exp)
{
  return std::pow(base, exp);
}

Real cxx_gamma(Real input)
{
  return std::tgamma(input);
}


} // namespace p3
} // namespace scream
