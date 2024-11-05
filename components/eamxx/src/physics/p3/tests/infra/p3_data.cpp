#include "p3_data.hpp"
#include "physics_constants.hpp"
#include "p3_ic_cases.hpp"

#include "ekat/ekat_assert.hpp"

using scream::Real;
using scream::Int;

namespace scream {
namespace p3 {

P3Data::P3Data (Int ncol_, Int nlev_)
  : ncol(ncol_), nlev(nlev_)
{
  do_predict_nc = true;
  do_prescribed_CCN = true;
  dt = -1; // model time step, s; set to invalid -1
  it = 1;
  // In/out
  qc              = Array2("cloud liquid water mixing ratio, kg/kg", ncol, nlev);
  nc              = Array2("cloud liquid drop number, #/kg", ncol, nlev);
  qr              = Array2("rain water mixing ratio, kg/kg", ncol, nlev);
  nr              = Array2("rain drop number, #/kg", ncol, nlev);
  qi              = Array2("total ice mass mixing ratio, kg/kg", ncol, nlev);
  ni              = Array2("total ice number, #/kg", ncol, nlev);
  qm              = Array2("rime ice mass mixing ratio, kg/kg", ncol, nlev);
  bm              = Array2("rime ice volume mixing ratio, m3/kg", ncol, nlev);
  qv              = Array2("water vapor mixing ratio, kg/kg", ncol, nlev);
  th_atm          = Array2("potential temperature, K", ncol, nlev);
  qv_prev         = Array2("prev-step water vapor mixing ratio, kg/kg", ncol, nlev);
  t_prev          = Array2("prev-step temperature, K", ncol, nlev);
  pres            = Array2("pressure, Pa", ncol, nlev);
  dz              = Array2("vertical grid spacing, m", ncol, nlev);
  nc_nuceat_tend  = Array2("ccn activated number tendency, kg-1 s-1", ncol, nlev);
  nccn_prescribed = Array2("CCN concentration, kg-1", ncol, nlev);
  ni_activated    = Array2("activated nuclei concentration, kg-1", ncol, nlev);
  inv_qc_relvar   = Array2("Assumed SGS 1/(var(qc)/mean(qc)), kg2/kg2", ncol, nlev);
  dpres           = Array2("pressure thickness, Pa", ncol, nlev);
  inv_exner       = Array2("Exner expression", ncol, nlev);
  // Out
  precip_liq_surf    = Array1("precipitation rate, liquid  m/s", ncol);
  precip_ice_surf    = Array1("precipitation rate, solid   m/s", ncol);
  diag_eff_radius_qc = Array2("effective radius, cloud, m", ncol, nlev);
  diag_eff_radius_qi = Array2("effective radius, ice, m", ncol, nlev);
  diag_eff_radius_qr = Array2("effective radius, rain, m", ncol, nlev);
  rho_qi             = Array2("bulk density of ice, kg/m", ncol, nlev);
  qv2qi_depos_tend   = Array2("qitend due to deposition/sublimation ", ncol, nlev);
  precip_liq_flux    = Array2("grid-box average rain flux (kg m^-2 s^-1), pverp", ncol, nlev+1);
  precip_ice_flux    = Array2("grid-box average ice/snow flux (kg m^-2 s^-1), pverp", ncol, nlev+1);
  cld_frac_r         = Array2("Rain cloud fraction", ncol, nlev);
  cld_frac_l         = Array2("Liquid cloud fraction", ncol, nlev);
  cld_frac_i         = Array2("Ice cloud fraction", ncol, nlev);
  liq_ice_exchange   = Array2("sum of liq-ice phase change tendenices", ncol, nlev);
  vap_liq_exchange   = Array2("sum of vap-liq phase change tendenices", ncol, nlev);
  vap_ice_exchange   = Array2("sum of vap-ice phase change tendenices", ncol, nlev);
}

P3DataIterator::P3DataIterator (const P3Data::Ptr& d) {
  init(d);
}

void P3DataIterator::init (const P3Data::Ptr& dp) {
  d_ = dp;
#define fdipb(name)                                                     \
  fields_.push_back({#name,                                             \
        2,                                                              \
        {d_->name.extent_int(0), d_->name.extent_int(1), d_->name.extent_int(2)}, \
        d_->name.data(),                                                \
        d_->name.size()})
  fdipb(qv); fdipb(th_atm); fdipb(pres);
  fdipb(dz); fdipb(nc_nuceat_tend); fdipb(nccn_prescribed); fdipb(ni_activated); fdipb(inv_qc_relvar); fdipb(qc);
  fdipb(nc); fdipb(qr); fdipb(nr); fdipb(qi); fdipb(ni);
  fdipb(qm); fdipb(bm); fdipb(precip_liq_surf); fdipb(precip_ice_surf);
  fdipb(diag_eff_radius_qc); fdipb(diag_eff_radius_qi); fdipb(diag_eff_radius_qr); fdipb(rho_qi);
  fdipb(dpres); fdipb(inv_exner); fdipb(qv2qi_depos_tend);
  fdipb(precip_liq_flux); fdipb(precip_ice_flux);
  fdipb(cld_frac_r); fdipb(cld_frac_l); fdipb(cld_frac_i);
  fdipb(liq_ice_exchange); fdipb(vap_liq_exchange);
  fdipb(vap_ice_exchange); fdipb(qv_prev); fdipb(t_prev);;
#undef fdipb
}

const P3DataIterator::RawArray&
P3DataIterator::getfield (Int i) const {
  EKAT_ASSERT(i >= 0 || i < nfield());
  return fields_[i];
}

int test_P3Data () {
  P3Data d(11, 72);
  return 0;
}

} // namespace p3
} // namespace scream
