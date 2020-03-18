#include "shoc_f90.hpp"
#include "shoc_functions_f90.hpp"
#include "shoc_constants.hpp"
#include "shoc_ic_cases.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"

using scream::Real;
using scream::Int;
extern "C" {
  void shoc_init_c(int nlev, Real gravit, Real rair, Real rh2o, Real cpair,
                   Real zvir, Real latvap, Real latice, Real karman,
                   Real* pref_mid, int nbot_shoc, int ntop_shoc);
  void shoc_use_cxx_c(bool use_cxx);
  void shoc_main_c(int shcol, int nlev, int nlevi, Real dtime, int nadv,
                   Real* host_dx, Real* host_dy, Real* thv, Real* zt_grid,
                   Real* zi_grid, Real* pres, Real* presi, Real* pdel,
                   Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc,
                   Real* wtracer_sfc, int num_qtracers, Real* w_field,
                   Real* exner, Real* phis, Real* host_dse, Real* tke,
                   Real* thetal, Real* qw, Real* u_wind, Real* v_wind,
                   Real* qtracers, Real* wthv_sec, Real* tkh, Real* tk,
                   Real* shoc_ql, Real* shoc_cldfrac, Real* pblh,
                   Real* shoc_mix, Real* isotropy, Real* w_sec, Real* thl_sec,
                   Real* qw_sec, Real* qwthl_sec, Real* wthl_sec, Real* wqw_sec,
                   Real* wtke_sec, Real* uw_sec, Real* vw_sec, Real* w3,
                   Real* wqls_sec, Real* brunt);
}

namespace scream {
namespace shoc {

FortranData::FortranData(Int shcol_, Int nlev_, Int nlevi_,
                         Int num_qtracers_, Int nadv_)
  : shcol(shcol_), nlev(nlev_), nlevi(nlevi_), num_qtracers(num_qtracers_),
    nadv(nadv_)
{
  dtime = -1; // model time step [s]; set to invalid -1

  // In variables
  host_dx = Array1("grid spacing of host model in x direction [m]", shcol);
  host_dy = Array1("grid spacing of host model in y direction [m]", shcol);
  zt_grid = Array2("heights for thermo grid [m]", shcol, nlev);
  zi_grid = Array2("heights for interface grid [m]", shcol, nlevi);
  pres = Array2("Pressure levels on thermo grid [Pa]", shcol, nlev);
  presi = Array2("Pressure levels on interface grid [Pa]", shcol, nlevi);
  pdel = Array2("Differences in pressure levels [Pa]", shcol, nlev);
  thv = Array2("Virtual potential temperature [K]", shcol, nlev);
  w_field = Array2("Large-scale vertical velocity [m/s]", shcol, nlev);
  wthl_sfc = Array1("Surface sensible heat flux [K m/s]", shcol);
  wqw_sfc = Array1("Surface latent heat flux [kg/kg m/s]", shcol);
  uw_sfc = Array1("Surface momentum flux (u direction) [m2/s2]", shcol);
  vw_sfc = Array1("Surface momentum flux (v direction) [m2/s2]", shcol);
  wtracer_sfc = Array2("Surface tracer flux [various]", shcol, num_qtracers);
  exner = Array2("Exner function [-]", shcol, nlev);
  phis = Array1("Host model surface geopotential height [m?]", shcol);

  // In/out variables
  host_dse = Array2("Dry static energy [J/kg]", shcol, nlev);
  tke = Array2("Turbulent kinetic energy [m2/s2]", shcol, nlev);
  thetal = Array2("Liquid water potential temperature [K]", shcol, nlev);
  qw = Array2("Total water mixing ratio [kg/kg]", shcol, nlev);
  u_wind = Array2("U wind component [m/s]", shcol, nlev);
  v_wind = Array2("V wind component [m/s]", shcol, nlev);
  wthv_sec = Array2("Buoyancy flux [K m/s]", shcol, nlev);
  qtracers = Array3("Tracers [varies]", shcol, nlev, num_qtracers);
  tk = Array2("Eddy coefficient for momentum [m2/s]", shcol, nlev);
  tkh = Array2("Eddy coefficent for heat [m2/s]", shcol, nlev);

  // Out variables (including diagnostics)
  shoc_cldfrac = Array2("Cloud fraction [-]", shcol, nlev);
  shoc_ql = Array2("Cloud liquid mixing ratio [kg/kg]", shcol, nlev);
  pblh = Array1("Planetary boundary layer depth [m]", shcol);
  shoc_mix = Array2("Turbulent length scale [m]", shcol, nlev);
  w_sec = Array2("Vertical velocity variance [m2/s2]", shcol, nlev);
  thl_sec = Array2("Temperature variance [K^2]", shcol, nlevi);
  qw_sec = Array2("Moisture variance [kg2/kg2]", shcol, nlevi);
  qwthl_sec = Array2("Temp moisture covariance [K kg/kg]", shcol, nlevi);
  wthl_sec = Array2("Vertical heat flux [K m/s]", shcol, nlevi);
  wqw_sec = Array2("Vertical moisture flux [K m/s]", shcol, nlevi);
  wtke_sec = Array2("Vertical tke flux [m3/s3]", shcol, nlevi);
  uw_sec = Array2("Vertical zonal momentum flux [m2/s2]", shcol, nlevi);
  vw_sec = Array2("Vertical meridional momentum flux [m2/s2]", shcol, nlevi);
  w3 = Array2("Third moment vertical velocity [m3/s3]", shcol, nlevi);
  wqls_sec = Array2("Liquid water flux [kg/kg m/s]", shcol, nlev);
  brunt = Array2("Brunt vaisala frequency [s-1]", shcol, nlev);
  isotropy = Array2("Return to isotropic timescale [s]", shcol, nlev);
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
  fdipb(host_dx); fdipb(host_dy);
  fdipb(zt_grid); fdipb(zi_grid); fdipb(pres); fdipb(presi); fdipb(pdel);
  fdipb(thv); fdipb(w_field);
  fdipb(wthl_sfc); fdipb(wqw_sfc); fdipb(uw_sfc); fdipb(vw_sfc);
  fdipb(wtracer_sfc); fdipb(exner);
  fdipb(phis);

  fdipb(host_dse); fdipb(tke); fdipb(thetal); fdipb(qw);
  fdipb(u_wind); fdipb(v_wind); fdipb(wthv_sec);
  fdipb(qtracers);
  fdipb(tk); fdipb(tkh);

  fdipb(shoc_cldfrac); fdipb(shoc_ql);
  fdipb(pblh);
  fdipb(shoc_mix); fdipb(w_sec); fdipb(thl_sec); fdipb(qw_sec);
  fdipb(qwthl_sec); fdipb(wthl_sec); fdipb(wqw_sec); fdipb(wtke_sec);
  fdipb(uw_sec); fdipb(vw_sec); fdipb(w3); fdipb(wqls_sec); fdipb(isotropy);
  fdipb(brunt);
#undef fdipb
}

const FortranDataIterator::RawArray&
FortranDataIterator::getfield (Int i) const {
  scream_assert(i >= 0 || i < nfield());
  return fields_[i];
}

void shoc_init(Int nlev, bool use_fortran) {
  static bool is_init = false;
  if (!is_init) {
    using KT     = KokkosTypes<HostDevice>;
    using Scalar = Real;
    using Array1 = typename KT::template lview<Scalar*>;
    using C = Constants<Scalar>;

    auto pref_mid = Array1("reference pressures (unused)", nlev);
    shoc_init_c((int)nlev, C::GGr, C::RAir, C::RH2O, C::CpAir, C::ZVir, C::LatVap,
                C::LatIce, C::Karman, pref_mid.data(), (int)nlev, 1);
    is_init = true;
  }
  shoc_use_cxx_c(!use_fortran);
}

void shoc_main(FortranData& d) {
  shoc_main_c((int)d.shcol, (int)d.nlev, (int)d.nlevi, d.dtime, (int)d.nadv,
              d.host_dx.data(), d.host_dy.data(), d.thv.data(),
              d.zt_grid.data(), d.zi_grid.data(), d.pres.data(), d.presi.data(),
              d.pdel.data(), d.wthl_sfc.data(), d.wqw_sfc.data(), d.uw_sfc.data(),
              d.vw_sfc.data(), d.wtracer_sfc.data(), (int)d.num_qtracers,
              d.w_field.data(), d.exner.data(), d.phis.data(), d.host_dse.data(),
              d.tke.data(), d.thetal.data(), d.qw.data(), d.u_wind.data(),
              d.v_wind.data(), d.qtracers.data(), d.wthv_sec.data(), d.tkh.data(),
              d.tk.data(), d.shoc_ql.data(), d.shoc_cldfrac.data(), d.pblh.data(),
              d.shoc_mix.data(), d.isotropy.data(), d.w_sec.data(),
              d.thl_sec.data(), d.qw_sec.data(), d.qwthl_sec.data(),
              d.wthl_sec.data(), d.wqw_sec.data(), d.wtke_sec.data(),
              d.uw_sec.data(), d.vw_sec.data(), d.w3.data(), d.wqls_sec.data(),
              d.brunt.data());
}

int test_FortranData () {
  Int shcol = 1;
  Int nlev = 128, num_tracers = 1;
  Int nadv = 1;
  FortranData d(shcol, nlev, nlev+1, num_tracers, nadv);
  return 0;
}

int test_shoc_init (bool use_fortran) {
  shoc_init(use_fortran);
  return 0;
}

int test_shoc_ic (bool use_fortran) {
  const auto d = ic::Factory::create();
  shoc_init(use_fortran);
  shoc_main(*d);
  return 0;
}

} // namespace shoc
} // namespace scream
