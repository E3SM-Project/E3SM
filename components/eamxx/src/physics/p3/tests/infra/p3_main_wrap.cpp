#include "p3_main_wrap.hpp"
#include "p3_data.hpp"
#include "p3_test_data.hpp"
#include "physics_constants.hpp"
#include "p3_ic_cases.hpp"

#include "ekat/ekat_assert.hpp"

using scream::Real;
using scream::Int;

namespace scream {
namespace p3 {

Int p3_main_wrap(const P3Data& d) {
  EKAT_REQUIRE_MSG(d.dt > 0, "invalid dt");
  return p3_main_host(d.qc.data(), d.nc.data(), d.qr.data(), d.nr.data(), d.th_atm.data(),
                      d.qv.data(), d.dt, d.qi.data(), d.qm.data(), d.ni.data(),
                      d.bm.data(), d.pres.data(), d.dz.data(), d.nc_nuceat_tend.data(), d.nccn_prescribed.data(),
                      d.ni_activated.data(), d.inv_qc_relvar.data(), d.it, d.precip_liq_surf.data(),
                      d.precip_ice_surf.data(), 1, d.ncol, 1, d.nlev, d.diag_eff_radius_qc.data(),
                      d.diag_eff_radius_qi.data(), d.diag_eff_radius_qr.data(), d.rho_qi.data(), d.do_predict_nc, d.do_prescribed_CCN,
                      d.dpres.data(), d.inv_exner.data(), d.qv2qi_depos_tend.data(),
                      d.precip_liq_flux.data(), d.precip_ice_flux.data(),
                      d.cld_frac_r.data(), d.cld_frac_l.data(), d.cld_frac_i.data(),
                      d.liq_ice_exchange.data(), d.vap_liq_exchange.data(),
                      d.vap_ice_exchange.data(),d.qv_prev.data(),d.t_prev.data() );
}

int test_p3_init () {
  using P3F = Functions<Real, DefaultDevice>;

  P3F::p3_init();
  P3GlobalForFortran::deinit();
  return 0;
}

int test_p3_ic () {
  using P3F = Functions<Real, DefaultDevice>;

  const auto d = ic::Factory::create(ic::Factory::mixed);

  d->dt = 300.0;
  P3F::p3_init();
  p3_main_wrap(*d);
  P3GlobalForFortran::deinit();
  return 0;
}

} // namespace p3
} // namespace scream
