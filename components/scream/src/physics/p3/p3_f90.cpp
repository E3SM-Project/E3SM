#include "p3_f90.hpp"
#include "p3_ic_cases.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"
#include "share/util/math_utils.hpp"


using scream::Real;
using scream::Int;
extern "C" {
  void p3_init_c(const char** lookup_file_dir, int* info);
  void p3_main_c(Real* qc, Real* nc, Real* qr, Real* nr, Real* th_old, Real* th,
                 Real* qv_old, Real* qv, Real dt, Real* qitot, Real* qirim,
                 Real* nitot, Real* birim, Real* ssat, Real* pres,
                 Real* dzq, Int it, Real* prt_liq, Real* prt_sol, Int its,
                 Int ite, Int kts, Int kte, Real* diag_ze,
                 Real* diag_effc, Real* diag_effi, Real* diag_vmi,
                 Real* diag_di, Real* diag_rhoi, 
                 bool log_predictNc);
}

namespace scream {
namespace p3 {

FortranData::FortranData (Int ncol_, Int nlev_)
  : ncol(ncol_), nlev(nlev_)
{

  dt = -1; // model time step, s; set to invalid -1
  it = 1;  // seems essentially unused
  // In/out
  qc = Array2("cloud liquid water mixing ratio, kg/kg", ncol, nlev);
  nc = Array2("cloud liquid drop number, #/kg", ncol, nlev);
  qr = Array2("rain water mixing ratio, kg/kg", ncol, nlev);
  nr = Array2("rain drop number, #/kg", ncol, nlev);
  qitot = Array2("total ice mass mixing ratio, kg/kg", ncol, nlev);
  nitot = Array2("total ice number, #/kg", ncol, nlev);
  qirim = Array2("rime ice mass mixing ratio, kg/kg", ncol, nlev);
  birim = Array2("rime ice volume mixing ratio, m3/kg", ncol, nlev);
  ssat = Array2("supersaturation (qv - qs), kg/kg", ncol, nlev);
  qv = Array2("water vapor mixing ratio, kg/kg", ncol, nlev);
  th = Array2("potential temperature, K", ncol, nlev);
  qv_old = Array2("qv at beginning of timestep, kg/kg", ncol, nlev);
  th_old = Array2("theta at beginning of timestep, K", ncol, nlev);
  pres = Array2("pressure, Pa", ncol, nlev);
  dzq = Array2("vertical grid spacing, m", ncol, nlev);
  // Out
  prt_liq = Array1("precipitation rate, liquid  m/s", ncol);
  prt_sol = Array1("precipitation rate, solid   m/s", ncol);
  diag_ze = Array2("equivalent reflectivity, dBZ", ncol, nlev);
  diag_effc = Array2("effective radius, cloud, m", ncol, nlev);
  diag_effi = Array2("effective radius, ice, m", ncol, nlev);
  diag_vmi = Array2("mass-weighted fall speed of ice, m/s", ncol, nlev);
  diag_di = Array2("mean diameter of ice, m", ncol, nlev);
  diag_rhoi = Array2("bulk density of ice, kg/m", ncol, nlev);
}

FortranDataIterator::FortranDataIterator (const FortranData::Ptr& d) {
  init(d);
}

void FortranDataIterator::init (const FortranData::Ptr& dp) {
  d_ = dp;
#define fdipb(name)                                                     \
  fields_.push_back({#name,                                             \
        2,                                                              \
        {d_->name.extent_int(0), d_->name.extent_int(1), d_->name.extent_int(2)}, \
        d_->name.data(),                                                \
        d_->name.size()})
  fdipb(qv); fdipb(th); fdipb(qv_old); fdipb(th_old); fdipb(pres);
  fdipb(dzq); fdipb(qc); fdipb(nc); fdipb(qr); fdipb(nr);
  fdipb(ssat); fdipb(qitot); fdipb(nitot);
  fdipb(qirim); fdipb(birim); fdipb(prt_liq); fdipb(prt_sol);
  fdipb(diag_ze); fdipb(diag_effc); fdipb(diag_effi);
  fdipb(diag_vmi); fdipb(diag_di); fdipb(diag_rhoi);
#undef fdipb
}

const FortranDataIterator::RawArray&
FortranDataIterator::getfield (Int i) const {
  scream_assert(i >= 0 || i < nfield());
  return fields_[i];
}

void p3_init () {
  static const char* dir = ".";
  Int info;
  p3_init_c(&dir, &info);
  scream_require_msg(info == 0, "p3_init_c returned info " << info);
}

void p3_main (const FortranData& d) {
  p3_main_c(d.qc.data(), d.nc.data(), d.qr.data(), d.nr.data(), d.th_old.data(),
            d.th.data(), d.qv_old.data(), d.qv.data(), d.dt, d.qitot.data(),
            d.qirim.data(), d.nitot.data(), d.birim.data(), d.ssat.data(),
            d.pres.data(), d.dzq.data(), d.it, d.prt_liq.data(),
            d.prt_sol.data(), 1, d.ncol, 1, d.nlev, d.diag_ze.data(),
            d.diag_effc.data(), d.diag_effi.data(), d.diag_vmi.data(),
            d.diag_di.data(), d.diag_rhoi.data(),
            d.log_predictnc);
}

Int check_against_python (const FortranData& d) {
  Int nerr = 0;
  if (util::is_single_precision<Real>::value) {
    const double tol = 0;
    if (util::reldif<double>(d.birim(0,d.nlev-1), 7.237245824853744e-08) > tol)
      ++nerr;
    if (util::reldif<double>(d.qirim(0,d.nlev-1), 9.047746971191373e-06) > tol)
      ++nerr;
    if (util::reldif<double>(d.nr(0,d.nlev-1), 3.177030468750000e+04) > tol)
      ++nerr;
  }
  return nerr;
}

int test_FortranData () {
  FortranData d(11, 72);
  return 0;
}

int test_p3_init () {
  p3_init();
  return 0;
}

int test_p3_main () {
  FortranData d(11, 72);
  d.dt = 0;
  p3_main(d);
  return 0;
}

int test_p3_ic () {
  const auto d = ic::Factory::create(ic::Factory::mixed);
  p3_init();
  p3_main(*d);
  return check_against_python(*d);
}

} // namespace p3
} // namespace scream
