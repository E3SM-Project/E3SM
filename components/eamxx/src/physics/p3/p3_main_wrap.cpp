#include "p3_main_wrap.hpp"
#include "p3_f90.hpp"
#include "p3_functions_f90.hpp"
#include "physics_constants.hpp"
#include "p3_ic_cases.hpp"

#include "ekat/ekat_assert.hpp"

using scream::Real;
using scream::Int;
extern "C" {
  void p3_main_c(Real* qc, Real* nc, Real* qr, Real* nr, Real* th_atm,
                 Real* qv, Real dt, Real* qi, Real* qm,
                 Real* ni, Real* bm, Real* pres,
                 Real* dz, Real* nc_nuceat_tend, Real* nccn_prescribed, Real* ni_activated, Real* inv_qc_relvar,
                 Int it, Real* precip_liq_surf, Real* precip_ice_surf, Int its,
                 Int ite, Int kts, Int kte, Real* diag_eff_radius_qc, Real* diag_eff_radius_qi,
                 Real* rho_qi, bool do_predict_nc, bool do_prescribed_CCN, Real* dpres, Real* inv_exner,
                 Real* qv2qi_depos_tend,
                 Real* precip_liq_flux, Real* precip_ice_flux, // 1 extra column size
                 Real* cld_frac_r, Real* cld_frac_l, Real* cld_frac_i,
                 Real* liq_ice_exchange, Real* vap_liq_exchange,
                 Real* vap_ice_exchange, Real* qv_prev, Real* t_prev, Real* elapsed_s);
}

namespace scream {
namespace p3 {

Int p3_main_wrap(const FortranData& d, bool use_fortran) {
  EKAT_REQUIRE_MSG(d.dt > 0, "invalid dt");
  if (use_fortran) {
    Real elapsed_s;
    p3_main_c(d.qc.data(), d.nc.data(), d.qr.data(), d.nr.data(),
              d.th_atm.data(), d.qv.data(), d.dt, d.qi.data(),
              d.qm.data(), d.ni.data(), d.bm.data(),
              d.pres.data(), d.dz.data(), d.nc_nuceat_tend.data(), d.nccn_prescribed.data(), d.ni_activated.data(), d.inv_qc_relvar.data(),
              d.it, d.precip_liq_surf.data(), d.precip_ice_surf.data(), 1, d.ncol, 1, d.nlev,
              d.diag_eff_radius_qc.data(), d.diag_eff_radius_qi.data(), d.rho_qi.data(),
              d.do_predict_nc, d.do_prescribed_CCN, d.dpres.data(), d.inv_exner.data(), d.qv2qi_depos_tend.data(),
              d.precip_liq_flux.data(), d.precip_ice_flux.data(), d.cld_frac_r.data(), d.cld_frac_l.data(), d.cld_frac_i.data(),
              d.liq_ice_exchange.data(), d.vap_liq_exchange.data(),d.vap_ice_exchange.data(),d.qv_prev.data(),d.t_prev.data(), &elapsed_s);
    return static_cast<Int>(elapsed_s * 1000000);
  }
  else {
    return p3_main_f(d.qc.data(), d.nc.data(), d.qr.data(), d.nr.data(), d.th_atm.data(),
                     d.qv.data(), d.dt, d.qi.data(), d.qm.data(), d.ni.data(),
                     d.bm.data(), d.pres.data(), d.dz.data(), d.nc_nuceat_tend.data(), d.nccn_prescribed.data(),
                     d.ni_activated.data(), d.inv_qc_relvar.data(), d.it, d.precip_liq_surf.data(),
                     d.precip_ice_surf.data(), 1, d.ncol, 1, d.nlev, d.diag_eff_radius_qc.data(),
                     d.diag_eff_radius_qi.data(), d.rho_qi.data(), d.do_predict_nc, d.do_prescribed_CCN,
                     d.dpres.data(), d.inv_exner.data(), d.qv2qi_depos_tend.data(),
                     d.precip_liq_flux.data(), d.precip_ice_flux.data(),
                     d.cld_frac_r.data(), d.cld_frac_l.data(), d.cld_frac_i.data(),
                     d.liq_ice_exchange.data(), d.vap_liq_exchange.data(),
                     d.vap_ice_exchange.data(),d.qv_prev.data(),d.t_prev.data() );

  }
}

int test_p3_init () {
  p3_init();
  P3GlobalForFortran::deinit();
  return 0;
}

int test_p3_ic (bool use_fortran) {
  const auto d = ic::Factory::create(ic::Factory::mixed);
  d->dt = 300.0;
  p3_init();
  p3_main_wrap(*d, use_fortran);
  P3GlobalForFortran::deinit();
  return 0;
}

} // namespace p3
} // namespace scream
