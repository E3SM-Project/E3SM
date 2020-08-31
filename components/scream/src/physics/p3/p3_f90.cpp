#include "p3_f90.hpp"
#include "p3_functions_f90.hpp"
#include "physics_constants.hpp"
#include "p3_ic_cases.hpp"

#include "ekat/ekat_assert.hpp"

using scream::Real;
using scream::Int;
extern "C" {
  void micro_p3_utils_init_c(Real Cpair, Real Rair, Real RH2O, Real RHO_H2O,
                 Real MWH2O, Real MWdry, Real gravit, Real LatVap, Real LatIce,
                 Real CpLiq, Real Tmelt, Real Pi, Int iulog, bool masterproc);
  void p3_init_c(const char** lookup_file_dir, int* info);
  void p3_main_c(Real* qc, Real* nc, Real* qr, Real* nr, Real* th,
                 Real* qv, Real dt, Real* qi, Real* qm,
                 Real* ni, Real* bm, Real* pres,
                 Real* dz, Real* nc_nuceat_tend, Real* ni_activated, Real* inv_qc_relvar,
                 Int it, Real* precip_liq_surf, Real* precip_ice_surf, Int its,
                 Int ite, Int kts, Int kte, Real* diag_effc, Real* diag_effi,
                 Real* rho_qi, bool do_predict_nc, Real* dpres, Real* exner,
                 Real* cmeiout, Real* precip_total_tend, Real* nevapr, Real* qr_evap_tend,
                 Real* precip_liq_flux, Real* precip_ice_flux, // 1 extra column size
                 Real* cld_frac_r, Real* cld_frac_l, Real* cld_frac_i, Real* mu_c, Real* lamc,
                 Real* liq_ice_exchange, Real* vap_liq_exchange,
                 Real* vap_ice_exchange);
}

namespace scream {
namespace p3 {

FortranData::FortranData (Int ncol_, Int nlev_)
  : ncol(ncol_), nlev(nlev_)
{
  do_predict_nc = true;
  dt = -1; // model time step, s; set to invalid -1
  it = 1;  // seems essentially unused
  // In/out
  qc              = Array2("cloud liquid water mixing ratio, kg/kg", ncol, nlev);
  nc             = Array2("cloud liquid drop number, #/kg", ncol, nlev);
  qr             = Array2("rain water mixing ratio, kg/kg", ncol, nlev);
  nr             = Array2("rain drop number, #/kg", ncol, nlev);
  qi             = Array2("total ice mass mixing ratio, kg/kg", ncol, nlev);
  ni             = Array2("total ice number, #/kg", ncol, nlev);
  qm             = Array2("rime ice mass mixing ratio, kg/kg", ncol, nlev);
  bm             = Array2("rime ice volume mixing ratio, m3/kg", ncol, nlev);
  qv             = Array2("water vapor mixing ratio, kg/kg", ncol, nlev);
  th             = Array2("potential temperature, K", ncol, nlev);
  pres           = Array2("pressure, Pa", ncol, nlev);
  dz             = Array2("vertical grid spacing, m", ncol, nlev);
  nc_nuceat_tend = Array2("ccn activated number tendency, kg-1 s-1", ncol, nlev);
  ni_activated   = Array2("activated nuclei concentration, kg-1", ncol, nlev);
  inv_qc_relvar  = Array2("Assumed SGS 1/(var(qc)/mean(qc)), kg2/kg2", ncol, nlev);
  dpres          = Array2("pressure thickness, Pa", ncol, nlev);
  exner          = Array2("Exner expression", ncol, nlev);
  // Out
  precip_liq_surf   = Array1("precipitation rate, liquid  m/s", ncol);
  precip_ice_surf   = Array1("precipitation rate, solid   m/s", ncol);
  diag_effc         = Array2("effective radius, cloud, m", ncol, nlev);
  diag_effi         = Array2("effective radius, ice, m", ncol, nlev);
  rho_qi            = Array2("bulk density of ice, kg/m", ncol, nlev);
  cmeiout           = Array2("qitend due to deposition/sublimation ", ncol, nlev);
  precip_total_tend = Array2("Total precipitation (rain + snow)", ncol, nlev);
  nevapr            = Array2("evaporation of total precipitation (rain + snow)", ncol, nlev);
  qr_evap_tend      = Array2("evaporation of rain", ncol, nlev);
  precip_liq_flux   = Array2("grid-box average rain flux (kg m^-2 s^-1), pverp", ncol, nlev+1);
  precip_ice_flux   = Array2("grid-box average ice/snow flux (kg m^-2 s^-1), pverp", ncol, nlev+1);
  cld_frac_r        = Array2("Rain cloud fraction", ncol, nlev);
  cld_frac_l        = Array2("Liquid cloud fraction", ncol, nlev);
  cld_frac_i        = Array2("Ice cloud fraction", ncol, nlev);
  mu_c              = Array2("Size distribution shape paramter", ncol, nlev);
  lamc              = Array2("Size distribution slope paramter", ncol, nlev);
  liq_ice_exchange  = Array2("sum of liq-ice phase change tendenices", ncol, nlev);
  vap_liq_exchange  = Array2("sum of vap-liq phase change tendenices", ncol, nlev);
  vap_ice_exchange  = Array2("sum of vap-ice phase change tendenices", ncol, nlev);
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
  fdipb(qv); fdipb(th); fdipb(pres);
  fdipb(dz); fdipb(nc_nuceat_tend); fdipb(ni_activated); fdipb(inv_qc_relvar); fdipb(qc);
  fdipb(nc); fdipb(qr); fdipb(nr); fdipb(qi); fdipb(ni);
  fdipb(qm); fdipb(bm); fdipb(precip_liq_surf); fdipb(precip_ice_surf);
  fdipb(diag_effc); fdipb(diag_effi); fdipb(rho_qi);
  fdipb(dpres); fdipb(exner); fdipb(cmeiout); fdipb(precip_total_tend);
  fdipb(nevapr); fdipb(qr_evap_tend); fdipb(precip_liq_flux); fdipb(precip_ice_flux);
  fdipb(cld_frac_r); fdipb(cld_frac_l); fdipb(cld_frac_i);
  fdipb(mu_c); fdipb(lamc), fdipb(liq_ice_exchange); fdipb(vap_liq_exchange);
  fdipb(vap_ice_exchange);;
#undef fdipb
}

const FortranDataIterator::RawArray&
FortranDataIterator::getfield (Int i) const {
  EKAT_ASSERT(i >= 0 || i < nfield());
  return fields_[i];
}

void micro_p3_utils_init () {
  using c = scream::physics::Constants<Real>;
  micro_p3_utils_init_c(c::Cpair, c::Rair, c::RH2O, c::RHO_H2O,
                 c::MWH2O, c::MWdry, c::gravit, c::LatVap, c::LatIce,
                 c::CpLiq, c::Tmelt, c::Pi, c::iulog, c::masterproc);
}

void p3_init () {
  static bool is_init = false;
  if (!is_init) {
    micro_p3_utils_init();
    static const char* dir = "./data";
    Int info;
    p3_init_c(&dir, &info);
    EKAT_REQUIRE_MSG(info == 0, "p3_init_c returned info " << info);
    is_init = true;
  }
}

void p3_main (const FortranData& d, bool use_fortran) {
  if (use_fortran) {
    p3_main_c(d.qc.data(), d.nc.data(), d.qr.data(), d.nr.data(),
              d.th.data(), d.qv.data(), d.dt, d.qi.data(),
              d.qm.data(), d.ni.data(), d.bm.data(),
              d.pres.data(), d.dz.data(), d.nc_nuceat_tend.data(), d.ni_activated.data(), d.inv_qc_relvar.data(),
              d.it, d.precip_liq_surf.data(), d.precip_ice_surf.data(), 1, d.ncol, 1, d.nlev,
              d.diag_effc.data(), d.diag_effi.data(), d.rho_qi.data(),
              d.do_predict_nc, d.dpres.data(), d.exner.data(), d.cmeiout.data(),
              d.precip_total_tend.data(), d.nevapr.data(), d.qr_evap_tend.data(),
              d.precip_liq_flux.data(), d.precip_ice_flux.data(), d.cld_frac_r.data(), d.cld_frac_l.data(),
              d.cld_frac_i.data(), d.mu_c.data(), d.lamc.data(),
              d.liq_ice_exchange.data(),
              d.vap_liq_exchange.data(),d.vap_ice_exchange.data());
  }
  else {
    p3_main_f(d.qc.data(), d.nc.data(), d.qr.data(), d.nr.data(), d.th.data(),
              d.qv.data(), d.dt, d.qi.data(), d.qm.data(), d.ni.data(),
              d.bm.data(), d.pres.data(), d.dz.data(), d.nc_nuceat_tend.data(),
              d.ni_activated.data(), d.inv_qc_relvar.data(), d.it, d.precip_liq_surf.data(),
              d.precip_ice_surf.data(), 1, d.ncol, 1, d.nlev, d.diag_effc.data(),
              d.diag_effi.data(), d.rho_qi.data(), d.do_predict_nc,
              d.dpres.data(), d.exner.data(), d.cmeiout.data(), d.precip_total_tend.data(),
              d.nevapr.data(), d.qr_evap_tend.data(), d.precip_liq_flux.data(), d.precip_ice_flux.data(),
              d.cld_frac_r.data(), d.cld_frac_l.data(), d.cld_frac_i.data(), d.mu_c.data(),
              d.lamc.data(), d.liq_ice_exchange.data(), d.vap_liq_exchange.data(),
              d.vap_ice_exchange.data());
  }
}

int test_FortranData () {
  FortranData d(11, 72);
  return 0;
}

int test_p3_init () {
  p3_init();
  P3GlobalForFortran::deinit();
  return 0;
}

int test_p3_ic (bool use_fortran) {
  const auto d = ic::Factory::create(ic::Factory::mixed);
  p3_init();
  p3_main(*d, use_fortran);
  P3GlobalForFortran::deinit();
  return 0;
}

} // namespace p3
} // namespace scream
