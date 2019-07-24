#include "p3_ic_cases.hpp"
#include "p3_constants.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_assert.hpp"

namespace scream {
namespace p3 {
namespace ic {

// From mixed_case_data.py in scream-docs at commit 4bbea4.
FortranData::Ptr make_mixed () {
  const Int nk = 72;
  Int k;
  const auto dp = std::make_shared<FortranData>(1, nk);
  auto& d = *dp;

  // We use double throughout b/c the python code uses double in the r.h.s.; the
  // result of the r.h.s. expression then is cast to Real.

  // max cld at ~700mb, decreasing to 0 at 900mb.
  for (k = 0; k < 15; ++k) d.qc(0,nk-20+k) = 1e-4*(1 - double(k)/14);
  for (k = 0; k < nk; ++k) d.nc(0,k) = 1e6;
  // max rain at 700mb, decreasing to zero at surf.
  for (k = 0; k < 20; ++k) d.qr(0,nk-20+k) = 1e-5*(1 - double(k)/19);
  for (k = 0; k < nk; ++k) d.nr(0,k) = 1e6;

  //                                                      v (in the python)
  for (k = 0; k < 15; ++k) d.qitot(0,nk-20+k) = 1e-4; //*(1 - double(k)/14)
  for (k = 0; k < nk; ++k) d.nitot(0,k) = 1e6;
  for (k = 0; k < 15; ++k) d.qirim(0,nk-20+k) = 1e-4*(1 - double(k)/14);
  // guess at reasonable value based on: m3/kg is 1/density and liquid water has
  // a density of 1000 kg/m3
  for (k = 0; k < 15; ++k) d.birim(0,nk-20+k) = 1e-2;

  // qv goes to zero halfway through profile (to avoid condensate near model
  // top)
  for (k = 0; k < nk; ++k) {
    const auto tmp = -5e-4 + 1e-3/double(nk)*k;
    d.qv(0,k) = tmp > 0 ? tmp : 0;
  }
  // make layer with qc saturated.
  for (k = 0; k < 15; ++k) d.qv(0,nk-20+k) = 5e-3;

  // pres is actually an input variable, but needed here to compute theta.
  for (k = 0; k < nk; ++k) d.pres(0,k) = 100 + 1e5/double(nk)*k;
  // pdel is actually an input variable, but needed here to compute theta.
  for (k = 0; k < nk; ++k) d.pdel(0,k) = 1e5/double(nk);
  // exner is actually an input variable, but needed here to compute theta.
  for (k = 0; k < nk; ++k) d.exner(0,k) = std::pow((1e5/d.pres(0,k)), (287.15/1005.0));
  // cloud fraction is an input variable, just set to 1 everywhere
  for (k = 0; k < nk; ++k) d.icldm(0,k) = 1.0;
  for (k = 0; k < nk; ++k) d.lcldm(0,k) = 1.0;
  for (k = 0; k < nk; ++k) d.rcldm(0,k) = 1.0;

  // To get potential temperature, start by making absolute temperature vary
  // between 150K at top of atmos and 300k at surface, then convert to potential
  // temp.
  FortranData::Array1 T("T", nk);
  for (k = 0; k < nk; ++k) {
    T(k) = 150 + 150/double(nk)*k;
    d.th(0,k) = T(k)*std::pow(Real(Constants<double>::P0/d.pres(0,k)), Real(Constants<double>::RD/Constants<double>::CP));
  }

  // The next section modifies inout variables to satisfy weird conditions
  // needed for code coverage.
  d.qitot(0,nk-1) = 1e-9;
  d.qv(0,nk-1) = 5e-2; // also needs to be supersaturated to avoid getting set
                       // to 0 earlier.

  // make lowest-level qc and qr>0 to trigger surface rain and drizzle
  // calculation.
  d.qr(0,nk-1) = 1e-6;
  d.qc(0,nk-1) = 1e-6;

  // make qitot>1e-8 where qr=0 to test rain collection conditional.
  d.qitot(0,nk-25) = 5e-8;

  // make qc>0 and qr>0 where T<233.15 to test homogeneous freezing.
  d.qc(0,35) = 1e-7;
  d.qv(0,35) = 1e-6;

  // deposition/condensation-freezing needs t<258.15 and >5% supersat.
  d.qv(0,33) = 1e-4;

  // input variables.
  d.dt = 1800;

  // compute vertical grid spacing dzq (in m) from pres and theta.
  static constexpr double
    g = 9.8; // gravity, m/s^2
  for (k = 0; k < nk; ++k) {
    double plo, phi; // pressure at cell edges, Pa
    plo = (k == 0   ) ?
      std::max<double>(0, d.pres(0,0) - 0.5*(d.pres(0,1) - d.pres(0,0))/(1 - 0)) :
      0.5*(d.pres(0,k-1) + d.pres(0,k));
    phi = (k == nk-1) ?
      d.pres(0,nk-1) + 0.5*(d.pres(0,nk-1) - d.pres(0,nk-2))/(1 - 0) :
      0.5*(d.pres(0,k) + d.pres(0,k+1));
    const auto dpres = phi - plo;
    d.dzq(0,k) = Constants<double>::RD*T(k)/(g*d.pres(0,k))*dpres;
  }

  return dp;
}

FortranData::Ptr Factory::create (IC ic) {
 switch (ic) {
   case mixed: return make_mixed();
 default:
   scream_require_msg(false, "Not an IC: " << ic);
 }
}

} // namespace ic
} // namespace p3
} // namespace scream
