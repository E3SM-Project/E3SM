#include "shoc_ic_cases.hpp"
#include "shoc_constants.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_assert.hpp"

namespace scream {
namespace shoc {
namespace ic {

namespace {

// Reference elevations for data interpolation.
const std::array<Real, 5> z_ref = {0.0, 520.0, 1480.0, 2000.0, 3000.0};
const std::array<Real, 5> qw_ref = {17.0, 16.3, 10.7, 4.2, 3.0};
const std::array<Real, 5> ql_ref = {0.0, 5.0, 7.0, 6.0, 0.0};
const std::array<Real, 5> theta_ref = {299.7, 298.7, 302.4, 308.2, 312.85};

// Wind speed interpolation data.
const std::array<Real, 3> wind_z_ref = {0.0, 700.0, 3000.0};
const std::array<Real, 3> u_ref = {-7.75, -8.75, -4.61};
const std::array<Real, 3> v_ref = {0.0, 0.1, 0.0};
const std::array<Real, 3> w_ref = {0.1, 0.1, 0.0};

// Linearly interpolates data at a specific elevation using elevation and
// data arrays.
template <size_t N>
Real interpolate_data(const std::array<Real, N>& ref_elevations,
                      const std::array<Real, N>& ref_data,
                      Real z)
{
  auto pos = std::lower_bound(ref_elevations.begin(), ref_elevations.end(), z);
  Int index = pos - ref_elevations.begin();
  if (index < N) {
    const Real f1 = ref_data[index];
    const Real f2 = ref_data[index+1];
    const Real z1 = ref_elevations[index];
    const Real z2 = ref_elevations[index+1];
    const Real dfdz = (f2 - f1) / (z2 - z1);
    const Real dz = z - z1;
    return f1 + dfdz * dz;
  }
  else {
    // Don't extrapolate off the end of the table.
    return ref_data[N-1];
  }
}

using KT     = KokkosTypes<HostDevice>;
using Scalar = Real;
using Array1 = typename KT::template lview<Scalar*>;
using Array2 = typename KT::template lview<Scalar**>;

// Calculates hydrostatic pressure for a specific column given elevation
// data.
void compute_column_pressure(Int col, Int nlev, const Array2& z,
                             Array2& pres) {
  using consts = Constants<Real>;
  const Int i = col;
  const Real k = consts::RAir / consts::CpAir;
  const Real c = -consts::GGr * pow(consts::P0, k) / consts::RAir;
  const Real p_s = 1015e2; // surface pressure
  const Real theta_s = theta_ref[0]; // surface temperature

  // Compute the cell and interface pressures at the bottom of the column.
  {
    Real dz = z(i, 1);
    Real theta1 = interpolate_data(z_ref, theta_ref, dz);
    if (std::abs(theta_s - theta1) < 1e-14 * theta_s) {
      pres(i, 0) = pow(pow(p_s, k) + k*c*dz/theta_s, 1.0/k);
    }
    else {
      Real ra = dz/(theta1 - theta_s);
      pres(i, 0) = pow(pow(p_s, k) + k*c*ra*log(theta1/theta_s), 1.0/k);
    }
  }

  // Move up the column, computing the pressures at cell centers and
  // interfaces.
  for (Int k = 1; k < nlev; ++k) {
    Real z0 = z(i, k-1);
    Real z1 = z(i, k);
    Real th0 = interpolate_data(z_ref, theta_ref, z0);
    Real th1 = interpolate_data(z_ref, theta_ref, z1);
    Real p0 = pres(i, k-1);
    if (std::abs(th0 - th1) < 1e-14 * th0) {
      pres(i, k) = pow(pow(p0, k) + k*c*(z1 - z0)/th0, 1.0/k);
    }
    else {
      Real ra = (z1 - z0)/(th1 - th0);
      pres(i, k) = pow(pow(p0, k) + k*c*ra*log(th1/th0), 1.0/k);
    }
  }
}

void interpolate_column_data(Real ztop, Int col, FortranData& d) {
  using consts = Constants<Real>;
  const Int i = col;
  const Int nlev = d.nlev;
  const Real dz = ztop/nlev;

  for (Int k = 0; k < nlev; ++k) {

    // Set up the grid.
    Real zi = k * dz;
    Real zt = (k+0.5) * dz;
    d.host_dx(i, k) = 5300.0;
    d.host_dy(i, k) = 5300.0;
    d.zi_grid(i, k) = zi;
    d.zt_grid(i, k) = zt;

    // Interpolate the potential temperature.
    const Real theta_zt = interpolate_data(z_ref, theta_ref, zt);
    const Real qw = interpolate_data(z_ref, qw_ref, zt);
    const Real ql = interpolate_data(z_ref, ql_ref, zt);
    d.qw(i, k) = qw;
    d.shoc_ql(i, k) = ql;
    Real zvir = (consts::RH2O / consts::RAir) - 1.0;
    d.thv(i, k) = theta_zt * (1.0 + zvir * qw);
    d.thetal(i, k) = theta_zt;

    // Set cell-centered wind speeds.
    d.u_wind(i, k) = interpolate_data(wind_z_ref, u_ref, zt);
    d.v_wind(i, k) = interpolate_data(wind_z_ref, v_ref, zt);
    d.w_field(i, k) = interpolate_data(wind_z_ref, w_ref, zt);

    // Set turbulent kinetic energy.
    // TODO: Not done yet!

    // Set tracers.
    for (Int q = 0; q < d.num_qtracers; ++q)
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

  // Compute hydrostatic pressure at cell centers and interfaces.
  compute_column_pressure(i, d.nlev, d.zt_grid, d.pres);
  compute_column_pressure(i, d.nlevi, d.zi_grid, d.presi);

  // Compute pressure differences.
  for (Int k = 0; k < nlev; ++k) {
    d.pdel(i, k) = d.presi(i, k+1) - d.presi(i, k);
  }
}

} // end anonymous namespace

// From scream-docs/shoc-port/shocintr.py.
FortranData::Ptr Factory::create (Int shcol, Int nlev) {

  Int num_qtracers = 1, nadv = 1;
  const auto dp = std::make_shared<FortranData>(shcol, nlev, nlev+1,
                                                num_qtracers, nadv);
  auto& d = *dp;

  d.dtime = 12;
  const Real ztop = 2400.0;
  for (Int col = 0; col < shcol; ++col) {
    interpolate_column_data(ztop, col, d);
  }

  return dp;
}

} // namespace ic
} // namespace shoc
} // namespace scream
