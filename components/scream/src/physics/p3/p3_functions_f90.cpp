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

void cloud_water_conservation_c(Real qc, Real qcnuc, Real dt, Real* qcaut, Real* qcacc, Real* qccol,
  Real* qcheti, Real* qcshd, Real* qiberg, Real* qisub, Real* qidep);

void rain_water_conservation_c(Real qr, Real qcaut, Real qcacc, Real qimlt, Real qcshd,
  Real dt, Real* qrevp, Real* qrcol, Real* qrheti);

void ice_water_conservation_c(Real qitot, Real qidep, Real qinuc, Real qiberg, Real qrcol, Real qccol,
  Real qrheti, Real qcheti, Real dt, Real* qisub, Real* qimlt);
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

std::shared_ptr<P3GlobalForFortran::Views> P3GlobalForFortran::s_views;

const P3GlobalForFortran::Views& P3GlobalForFortran::get()
{
  if (!P3GlobalForFortran::s_views) {
    P3GlobalForFortran::s_views = std::make_shared<Views>();
    P3F::init_kokkos_ice_lookup_tables(s_views->m_itab, s_views->m_itabcol);
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

  typename P3F::Smask qiti_gt_small(true);
  typename P3F::Spack qitot(qitot_);
  typename P3F::Spack nitot(nitot_);
  typename P3F::Spack qirim(qirim_);
  typename P3F::Spack rhop(rhop_);
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

  typename P3F::Smask qiti_gt_small(true);
  typename P3F::Spack qr(qr_);
  typename P3F::Spack nr(nr_);
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

void cloud_water_conservation_f(Real qc_, Real qcnuc_, Real dt, Real* qcaut_, Real* qcacc_, Real* qccol_,
  Real* qcheti_, Real* qcshd_, Real* qiberg_, Real* qisub_, Real* qidep_)
  {
    using P3F = Functions<Real, HostDevice>;

    typename P3F::Spack qc(qc_);
    typename P3F::Spack qcnuc(qcnuc_);
    typename P3F::Spack qcaut(*qcaut_);
    typename P3F::Spack qcacc(*qcacc_);
    typename P3F::Spack qccol(*qccol_);
    typename P3F::Spack qcheti(*qcheti_);
    typename P3F::Spack qcshd(*qcshd_);
    typename P3F::Spack qiberg(*qiberg_);
    typename P3F::Spack qisub(*qisub_);
    typename P3F::Spack qidep(*qidep_);

    P3F::cloud_water_conservation(qc, qcnuc, dt, qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep);
  }

void rain_water_conservation_f(Real qr_, Real qcaut_, Real qcacc_, Real qimlt_, Real qcshd_,
  Real dt, Real* qrevp_, Real* qrcol_, Real* qrheti_)
  {
    using P3F = Functions<Real, HostDevice>;

    typename P3F::Spack qr(qr_);
    typename P3F::Spack qcaut(qcaut_);
    typename P3F::Spack qcacc(qcacc_);
    typename P3F::Spack qimlt(qimlt_);
    typename P3F::Spack qcshd(qcshd_);
    typename P3F::Spack qrevp(*qrevp_);
    typename P3F::Spack qrcol(*qrcol_);
    typename P3F::Spack qrheti(*qrheti_);

    P3F::rain_water_conservation(qr, qcaut, qcacc, qimlt, qcshd, dt, qrevp, qrcol, qrheti);
    *qrevp_ = qrevp[0];
    *qrcol_ = qrcol[0];
    *qrheti_ = qrheti[0];
  }

void ice_water_conservation_f(Real qitot_, Real qidep_, Real qinuc_, Real qiberg_, Real qrcol_, Real qccol_,
  Real qrheti_, Real qcheti_, Real dt, Real* qisub_, Real* qimlt_)
{
    using P3F = Functions<Real, HostDevice>;
    typename P3F::Spack qitot(qitot_);
    typename P3F::Spack qidep(qidep_);
    typename P3F::Spack qinuc(qinuc_);
    typename P3F::Spack qiberg(qiberg_);
    typename P3F::Spack qrcol(qrcol_);
    typename P3F::Spack qccol(qccol_);
    typename P3F::Spack qrheti(qrheti_);
    typename P3F::Spack qcheti(qcheti_);
    typename P3F::Spack qisub(*qisub_);
    typename P3F::Spack qimlt(*qimlt_);

    P3F::ice_water_conservation(qitot, qidep, qinuc, qiberg, qrcol, qccol, qrheti, qcheti, dt, qisub, qimlt);
    *qisub_ = qisub[0]; 
    *qimlt_ = qimlt[0]; 

}

} // namespace p3
} // namespace scream
