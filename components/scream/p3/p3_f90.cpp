#include "p3_f90.hpp"
#include "p3_ic_cases.hpp"
#include "scream_util.hpp"

using scream::Real;
using scream::Int;
extern "C" {
  void p3_init_c(const char** lookup_file_dir, int ncat, int* info);
  void p3_main_c(Real* qc, Real* nc, Real* qr, Real* nr, Real* th_old, Real* th,
                 Real* qv_old, Real* qv, Real dt, Real* qitot, Real* qirim,
                 Real* nitot, Real* birim, Real* ssat, Real* uzpl, Real* pres,
                 Real* dzq, Int it, Real* prt_liq, Real* prt_sol, Int its,
                 Int ite, Int kts, Int kte, Int nCat, Real* diag_ze,
                 Real* diag_effc, Real* diag_effi, Real* diag_vmi,
                 Real* diag_di, Real* diag_rhoi, Int n_diag_2d, Real* diag_2d,
                 Int n_diag_3d, Real* diag_3d, bool log_predictNc,
                 bool typeDiags_ON, const char** model, Real* prt_drzl,
                 Real* prt_rain, Real* prt_crys, Real* prt_snow, Real* prt_grpl,
                 Real* prt_pell, Real* prt_hail, Real* prt_sndp);
}

namespace scream {
namespace p3 {

FortranData::FortranData (Int ncol_, Int nlev_)
  : ncol(ncol_), nlev(nlev_)
{
  const Int n_diag_2d = 1, n_diag_3d = 1;

  dt = -1; // model time step, s; set to invalid -1
  it = 1;  // seems essentially unused
  // In/out
  qc = Array2("cloud liquid water mixing ratio, kg/kg", ncol, nlev);
  nc = Array2("cloud liquid drop number, #/kg", ncol, nlev);
  qr = Array2("rain water mixing ratio, kg/kg", ncol, nlev);
  nr = Array2("rain drop number, #/kg", ncol, nlev);
  qitot = Array3("total ice mass mixing ratio, kg/kg", ncol, nlev, ncat);
  nitot = Array3("total ice number, #/kg", ncol, nlev, ncat);
  qirim = Array3("rime ice mass mixing ratio, kg/kg", ncol, nlev, ncat);
  birim = Array3("rime ice volume mixing ratio, m3/kg", ncol, nlev, ncat);
  ssat = Array2("supersaturation (qv - qs), kg/kg", ncol, nlev);
  qv = Array2("water vapor mixing ratio, kg/kg", ncol, nlev);
  th = Array2("potential temperature, K", ncol, nlev);
  qv_old = Array2("qv at beginning of timestep, kg/kg", ncol, nlev);
  th_old = Array2("theta at beginning of timestep, K", ncol, nlev);
  uzpl = Array2("vertical air velocity, m/s", ncol, nlev);
  pres = Array2("pressure, Pa", ncol, nlev);
  dzq = Array2("vertical grid spacing, m", ncol, nlev);
  // Out
  prt_liq = Array1("precipitation rate, liquid  m/s", ncol);
  prt_sol = Array1("precipitation rate, solid   m/s", ncol);
  prt_drzl = Array1("precip rate, drizzle       m/s", ncol);
  prt_rain = Array1("precip rate, rain          m/s", ncol);
  prt_crys = Array1("precip rate, ice cystals   m/s", ncol);
  prt_snow = Array1("precip rate, snow          m/s", ncol);
  prt_grpl = Array1("precip rate, graupel       m/s", ncol);
  prt_pell = Array1("precip rate, ice pellets   m/s", ncol);
  prt_hail = Array1("precip rate, hail          m/s", ncol);
  prt_sndp = Array1("precip rate, unmelted snow m/s", ncol);
  diag_ze = Array2("equivalent reflectivity, dBZ", ncol, nlev);
  diag_effc = Array2("effective radius, cloud, m", ncol, nlev);
  diag_2d = Array2("user-defined 2D diagnostic fields", ncol, n_diag_2d);
  diag_effi = Array3("effective radius, ice, m", ncol, nlev, ncat);
  diag_vmi = Array3("mass-weighted fall speed of ice, m/s", ncol, nlev, ncat);
  diag_di = Array3("mean diameter of ice, m", ncol, nlev, ncat);
  diag_rhoi = Array3("bulk density of ice, kg/m", ncol, nlev, ncat);
  diag_3d = Array3("user-defined 3D diagnostic fields", ncol, nlev, n_diag_3d);
}

FortranDataIterator::FortranDataIterator (const FortranData::Ptr& d) {
  init(d);
}

void FortranDataIterator::init (const FortranData::Ptr& dp) {
  d_ = dp;
  fields_.reserve(34);
#define fdipb(name)                                                     \
  fields_.push_back({#name,                                             \
        2,                                                              \
        {d_->name.extent_int(0), d_->name.extent_int(1), d_->name.extent_int(2)}, \
        d_->name.data(),                                                \
        d_->name.size()})
  fdipb(qv); fdipb(th); fdipb(qv_old); fdipb(th_old); fdipb(pres);
  fdipb(dzq); fdipb(qc); fdipb(nc); fdipb(qr); fdipb(nr);
  fdipb(ssat); fdipb(uzpl); fdipb(qitot); fdipb(nitot);
  fdipb(qirim); fdipb(birim); fdipb(prt_liq); fdipb(prt_sol);
  fdipb(prt_drzl); fdipb(prt_rain); fdipb(prt_crys); fdipb(prt_snow);
  fdipb(prt_grpl); fdipb(prt_pell); fdipb(prt_hail); fdipb(prt_sndp);
  fdipb(diag_ze); fdipb(diag_effc); fdipb(diag_2d); fdipb(diag_effi);
  fdipb(diag_vmi); fdipb(diag_di); fdipb(diag_rhoi); fdipb(diag_3d);
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
  p3_init_c(&dir, 1, &info);
  scream_throw_if(info != 0, "p3_init_c returned info " << info);
}

void p3_main (const FortranData& d) {
  static const char* model = "GEM";
  p3_main_c(d.qc.data(), d.nc.data(), d.qr.data(), d.nr.data(), d.th_old.data(),
            d.th.data(), d.qv_old.data(), d.qv.data(), d.dt, d.qitot.data(),
            d.qirim.data(), d.nitot.data(), d.birim.data(), d.ssat.data(),
            d.uzpl.data(), d.pres.data(), d.dzq.data(), d.it, d.prt_liq.data(),
            d.prt_sol.data(), 1, d.ncol, 1, d.nlev, d.ncat, d.diag_ze.data(),
            d.diag_effc.data(), d.diag_effi.data(), d.diag_vmi.data(),
            d.diag_di.data(), d.diag_rhoi.data(), d.diag_2d.extent_int(1),
            d.diag_2d.data(), d.diag_3d.extent_int(2), d.diag_3d.data(),
            d.log_predictnc, d.typediags_on, &model, d.prt_drzl.data(),
            d.prt_rain.data(), d.prt_crys.data(), d. prt_snow.data(),
            d.prt_grpl.data(), d.prt_pell.data(), d.prt_hail.data(),
            d.prt_sndp.data());
}

Int check_against_python (const FortranData& d) {
  Int nerr = 0;
  if (util::is_single_precision<Real>::value) {
    const double tol = 0;
    if (util::reldif<double>(d.birim(0,d.nlev-1,0), 7.237245824853744e-08) > tol)
      ++nerr;
    if (util::reldif<double>(d.qirim(0,d.nlev-1,0), 9.047746971191373e-06) > tol)
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
