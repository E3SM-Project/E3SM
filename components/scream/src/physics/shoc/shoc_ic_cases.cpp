#include "shoc_ic_cases.hpp"
#include "shoc_constants.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_assert.hpp"

namespace scream {
namespace shoc {
namespace ic {

namespace {

Real interpolate_theta(Real z) {
}

Real interpolate_qw(Real z) {
}

Real interpolate_ql(Real z) {
}

Real calc_hydrostatic_p(Real p0, Real theta0, Real z, theta_z) {
}

}

// From scream-docs/shoc-port/shocintr.py.
FortranData::Ptr Factory::create (Int shcol, Int nlev) {
  using consts = Constants<Real>;

  Int num_qtracers = 1, nadv = 1;
  const auto dp = std::make_shared<FortranData>(shcol, nlev, nlev+1,
                                                num_qtracers, nadv);
  auto& d = *dp;

  d.dtime = 12;
  const Real ztop = 2400.0;
  const Real dz = ztop/nlev;
  const Real ps = 1015e2;
  const Real theta0 = interpolate_theta(0.0);
  for (Int i = 0; i < shcol; ++i) {

    for (Int k = 0; k < nlev; ++k) {
      // Set up the grid.
      Real zi = k * dz;
      Real zt = (k+0.5) * dz;
      d.host_dx(i, k) = 5300.0;
      d.host_dy(i, k) = 5300.0;
      d.zi_grid(i, k) = zi;
      d.zt_grid(i, k) = zt;

      // Set the potential temperature and pressure.
      const Real theta_zi = interpolate_theta(zi);
      const Real theta_zt = interpolate_theta(zt);
      const Real qw = interpolate_qw(zt);
      const Real ql = interpolate_ql(zt);
      d.qw(i, k) = qw;
      d.ql(i, k) = ql;
      d.thv(i, k) = theta_to_thetav(theta_zt, qw);
      d.thetal(i, k) = theta_zt;
      d.pres(i, k) = calc_hydrostatic_p(ps, theta0, zi, theta_zt);
      d.presi(i, k) = calc_hydrostatic_p(ps, theta0, zi, theta_zi);
      const Real pres1 = calc_hydrostatic_p(ps, theta0, zi+dz, interpolate_theta(zi+dz));
      d.pdel(i, k) = pres1 - d.pres(i, k);

      // Set wind speeds.
      d.u_wind(i, k) = interpolate_u(zt);
      d.v_wind(i, k) = interpolate_v(zt);
      d.w_field(i, k) = interpolate_w(zt);

      // Set turbulent kinetic energy.
      // TODO: Not done yet!

      // Set tracers.
      for (Int q = 0; q < num_qtracers; ++q)
      {
        d.qtracers(i, k, q) = sin(q/3. + 0.1 * zt * 0.2 * (k+1));
        d.wtracer_sfc(i, q) = 0.1 * k;
      }

      // Set surface fluxes and forcings.
      // (Not clear whether these can be set independent of mesh)
      d.wthl_sfc[i] = 1e-4;
      d.wqw_sfc[i] = 1e-6;
      d.uw_sfc[i] = 1e-2;
      d.vw_sfc[i] = 1e-4;
    }
  }

  return dp;
}

} // namespace ic
} // namespace shoc
} // namespace scream
