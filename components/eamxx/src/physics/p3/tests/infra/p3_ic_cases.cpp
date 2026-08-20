#include "p3_ic_cases.hpp"
#include "physics_constants.hpp"

#include <algorithm>
#include <cmath>
#include <ekat_assert.hpp>
#include <limits>

namespace scream {
namespace p3 {
namespace ic {

// From mixed_case_data.py in scream-docs at commit 4bbea4.
P3Data::Ptr make_mixed (const Int ncol, const Int nlev) {
  using consts = scream::physics::Constants<Real>;

  const Int nk = nlev;
  EKAT_REQUIRE_MSG(nk >= 2, "The mixed P3 IC requires at least two levels");
  Int k;
  const auto dp = std::make_shared<P3Data>(ncol, nk);
  auto& d = *dp;

  for (Int i = 0; i < ncol; ++i) {
    // For column i = 0, use the ICs as originally coded in python and
    // subsequently modified here. For columns i > 0, introduce some small
    // variations.

    // nccn_presribed for the cases where do_prescribed_CCN=true
    for (k=0; k < nk; ++k) d.nccn_prescribed(i,k) = (100.0 + double(i)/ncol - double(k)/nk) * 1e6;

    const Int k_profile_start = std::max(0, nk-20);
    const Int n_cloud_levels = std::min(15, nk-1);
    // Leave at least one level without rain so the collection trigger can
    // be placed after the profile has been initialized.
    const Int n_rain_levels = std::min(20, std::max(0, nk-2));

    // max cld at ~700mb, decreasing to 0 at 900mb.
    for (k = 0; k < n_cloud_levels; ++k) {
      const auto fraction = n_cloud_levels > 1 ? double(k)/(n_cloud_levels-1) : 0;
      d.qc(i,k_profile_start+k) = 1e-4*(1 - fraction);
    }
    for (k = 0; k < nk; ++k) d.nc(i,k) = 1e6;
    // max rain at 700mb, decreasing to zero at surf.
    for (k = 0; k < n_rain_levels; ++k) {
      const auto fraction = n_rain_levels > 1 ? double(k)/(n_rain_levels-1) : 0;
      d.qr(i,k_profile_start+k) = 1e-5*(1 - fraction);
    }
    for (k = 0; k < nk; ++k) d.nr(i,k) = 1e6;

    //                                                      v (in the python)
    for (k = 0; k < n_cloud_levels; ++k) d.qi(i,k_profile_start+k) = 1e-4;
    for (k = 0; k < nk; ++k) d.ni(i,k) = 1e6;
    for (k = 0; k < n_cloud_levels; ++k) {
      const auto fraction = n_cloud_levels > 1 ? double(k)/(n_cloud_levels-1) : 0;
      d.qm(i,k_profile_start+k) = 1e-4*(1 - fraction);
    }
    // guess at reasonable value based on: m3/kg is 1/density and liquid water has
    // a density of 1000 kg/m3
    for (k = 0; k < n_cloud_levels; ++k) d.bm(i,k_profile_start+k) = 1e-2;

    // qv goes to zero halfway through profile (to avoid condensate near model
    // top)
    for (k = 0; k < nk; ++k) {
      const auto tmp = -5e-4 + 1e-3/double(nk)*k;
      d.qv(i,k) = tmp > 0 ? tmp : 0;
    }
    // make layer with qc saturated.
    for (k = 0; k < n_cloud_levels; ++k) d.qv(i,k_profile_start+k) = 5e-3;

    // pres is actually an input variable, but needed here to compute theta.
    for (k = 0; k < nk; ++k) d.pres(i,k) = 100 + 1e5/double(nk)*k;
    // dpres is actually an input variable, but needed here to compute theta.
    for (k = 0; k < nk; ++k) d.dpres(i,k) = 1e5/double(nk);
    // inv_exner is actually an input variable, but needed here to compute theta.
    for (k = 0; k < nk; ++k) d.inv_exner(i,k) = std::pow((1e5/d.pres(i,k)), (287.15/1005.0));
    // cloud fraction is an input variable, just set to 1 everywhere
    for (k = 0; k < nk; ++k) d.cld_frac_i(i,k) = 1.0;
    for (k = 0; k < nk; ++k) d.cld_frac_l(i,k) = 1.0;
    for (k = 0; k < nk; ++k) d.cld_frac_r(i,k) = 1.0;
    // inv_qc_relvar=mean(qc)/var(qc) measures subgrid qc variability. It is computed in SHOC
    // and used by P3. It can range between 0.1 and 10.0. Setting to a typical value of 1.0
    // here.
    for (k = 0; k < nk; ++k) d.inv_qc_relvar(i,k) = 1.0;

    // To get potential temperature, start by making absolute temperature vary
    // between 150K at top of atmos and 300k at surface, then convert to potential
    // temp.
    P3Data::Array1 T_atm("T", nk);
    for (k = 0; k < nk; ++k) {
      T_atm(k) = 150 + 150/double(nk)*k;
      if (i > 0) T_atm(k) += ((i % 3) - 0.5)/double(nk)*k;
      d.th_atm(i,k) = T_atm(k)*std::pow(Real(consts::P0.value/d.pres(i,k)), Real(consts::RD.value/consts::CP.value));
    }

    // The next section modifies inout variables to satisfy weird conditions
    // needed for code coverage.
    d.qi(i,nk-1) = 1e-9;
    d.qv(i,nk-1) = 5e-2; // also needs to be supersaturated to avoid getting set
    // to 0 earlier.

    // make lowest-level qc and qr>0 to trigger surface rain and drizzle
    // calculation.
    d.qr(i,nk-1) = 1e-6;
    d.qc(i,nk-1) = 1e-6;

    // Apply coverage triggers only after the common profiles and temperature
    // have been initialized.  Search the profile rather than relying on a
    // particular vertical resolution.
    Int k_rain_collection = -1;
    for (k = nk-1; k >= 0; --k) {
      if (k_rain_collection < 0 && d.qr(i,k) == 0) {
        k_rain_collection = k;
      }
    }
    EKAT_REQUIRE_MSG(k_rain_collection >= 0, "Could not place the rain-collection trigger");

    // Select the levels nearest the target temperatures.  Using the nearest
    // level rather than a fixed index keeps the IC usable at any resolution.
    auto nearest_temperature_level = [&](const Real target, const Int exclude) {
      Int best = -1;
      Real best_error = std::numeric_limits<Real>::max();
      for (Int j = 0; j < nk; ++j) {
        if (j == exclude) continue;
        const Real error = std::abs(T_atm(j) - target);
        if (error < best_error) {
          best = j;
          best_error = error;
        }
      }
      return best;
    };
    const Int k_homog_freeze = nearest_temperature_level(233.15, -1);
    const Int k_dep_cond_freeze = nearest_temperature_level(258.15, k_homog_freeze);

    // qr is intentionally zero at the collection level.
    d.qi(i,k_rain_collection) = 5e-8;
    d.qc(i,k_homog_freeze) = 1e-7;
    d.qr(i,k_homog_freeze) = 1e-7;
    d.qv(i,k_homog_freeze) = 1e-6;
    // Deposition/condensation-freezing needs T < 258.15 and supersaturation.
    d.qv(i,k_dep_cond_freeze) = 1e-4;

    // set qv_prev and t_prev to qv and T vals
    for (k = 0; k < nk; ++k){
      d.qv_prev(i,k) = d.qv(i,k);
      d.t_prev(i,k) = T_atm(k);
    }

    // compute vertical grid spacing dz (in m) from pres and theta.
    static constexpr double
      g = 9.8; // gravity, m/s^2
    for (k = 0; k < nk; ++k) {
      double plo, phi; // pressure at cell edges, Pa
      plo = (k == 0   ) ?
        std::max<double>(0, d.pres(i,0) - 0.5*(d.pres(i,1) - d.pres(i,0))/(1 - 0)) :
        0.5*(d.pres(i,k-1) + d.pres(i,k));
      phi = (k == nk-1) ?
        d.pres(i,nk-1) + 0.5*(d.pres(i,nk-1) - d.pres(i,nk-2))/(1 - 0) :
        0.5*(d.pres(i,k) + d.pres(i,k+1));
      const auto dpres = phi - plo;
      d.dz(i,k) = consts::RD.value*T_atm(k)/(g*d.pres(i,k))*dpres;
    }
    for (k = 0; k < nk; ++k) {
      d.hetfrz_immersion_nucleation_tend(i,k) = 0.01;
      d.hetfrz_contact_nucleation_tend(i,k) = 0.02;
      d.hetfrz_deposition_nucleation_tend(i,k) = 0.03;
    }
  }

  return dp;
}

P3Data::Ptr Factory::create (IC ic, Int ncol, Int nlev)
{
  P3Data::Ptr ret;
  switch (ic) {
    case mixed: ret = make_mixed(ncol, nlev); break;
    default:
      EKAT_REQUIRE_MSG(false, "Not an IC: " << ic);
  }
  return ret;
}

} // namespace ic
} // namespace p3
} // namespace scream
