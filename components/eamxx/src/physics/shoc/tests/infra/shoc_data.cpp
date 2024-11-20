#include "shoc_data.hpp"
#include "physics_constants.hpp"
#include "shoc_ic_cases.hpp"

#include "ekat/ekat_assert.hpp"

using scream::Real;
using scream::Int;

namespace scream {
namespace shoc {

FortranData::FortranData(Int shcol_, Int nlev_, Int nlevi_,
                         Int num_qtracers_)
  : shcol(shcol_), nlev(nlev_), nlevi(nlevi_), num_qtracers(num_qtracers_)
{
  dtime = -1; // model time step [s]; set to invalid -1
  nadv = -1;  // depends on timestep, so also invalid

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
  inv_exner = Array2("Inverse of the Exner function [-]", shcol, nlev);
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
  brunt = Array2("Brunt-Vaisala frequency [s-1]", shcol, nlev);
  isotropy = Array2("Return to isotropic timescale [s]", shcol, nlev);
  shoc_ql2 = Array2("Variance in liquid water [kg2/kg2]", shcol,nlev);
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
  fdipb(wtracer_sfc); fdipb(inv_exner);
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
  fdipb(brunt); fdipb(shoc_ql2);
#undef fdipb
}

const FortranDataIterator::RawArray&
FortranDataIterator::getfield (Int i) const {
  EKAT_ASSERT(i >= 0 || i < nfield());
  return fields_[i];
}

int test_FortranData () {
  Int shcol = 1;
  Int nlev = 128, num_tracers = 1;
  FortranData d(shcol, nlev, nlev+1, num_tracers);
  return 0;
}

} // namespace shoc
} // namespace scream
