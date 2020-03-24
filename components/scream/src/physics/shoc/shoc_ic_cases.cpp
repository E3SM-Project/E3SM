#include "shoc_ic_cases.hpp"
#include "physics_constants.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_assert.hpp"

namespace scream {
namespace shoc {
namespace ic {

namespace {

// Reference elevations for data interpolation.
const std::array<Real, 5> z_ref = {0.0, 520.0, 1480.0, 2000.0, 3000.0};
const std::array<Real, 5> qw_ref = {0.017, 0.0163, 0.0107, 0.0042, 0.003};
const std::array<Real, 5> ql_ref = {0.0, 0.005, 0.007, 0.006, 0.0};
const std::array<Real, 5> theta_ref = {299.7, 298.7, 302.4, 308.2, 312.85};

// Wind speed interpolation data.
const std::array<Real, 3> wind_z_ref = {0.0, 700.0, 3000.0};
const std::array<Real, 3> u_ref = {-7.75, -8.75, -4.61};
const std::array<Real, 3> v_ref = {0.0, 0.1, 0.0};
const std::array<Real, 3> w_ref = {0.1, 0.1, 0.0};

const Real surface_pressure = 1015e2;

// Linearly interpolates data at a specific elevation using elevation and
// data arrays.
template <size_t N>
Real interpolate_data(const std::array<Real, N>& ref_elevations,
                      const std::array<Real, N>& ref_data,
                      Real z)
{
  auto pos = std::lower_bound(ref_elevations.begin(), ref_elevations.end(), z);
  Int index = pos - ref_elevations.begin();
  if (index == 0)
    return ref_data[0];
  else if (index < (Int)N) {
    const Real a = (z - ref_elevations[index-1]) /
                   (ref_elevations[index] - ref_elevations[index-1]);
    return (1.0 - a) * ref_data[index-1] + a * ref_data[index];
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
  using consts = scream::physics::Constants<Real>;
  const Int i = col;
  const Real k = consts::Rair / consts::Cpair;
  const Real c = -consts::gravit * pow(consts::P0, k) / consts::Rair;
  const Real p_s = surface_pressure;

  // Move up the column, computing the pressures at each elevation.
  for (Int j = 0; j < nlev; ++j) {
    Real z0 = (j == 0) ? 0.0 : z(i, j-1);
    Real z1 = z(i, j);
    Real th0 = interpolate_data(z_ref, theta_ref, z0);
    Real th1 = interpolate_data(z_ref, theta_ref, z1);
    Real p0 = (j == 0) ? p_s : pres(i, j-1);
    if (std::abs(th0 - th1) < 1e-14 * th0) {
      pres(i, j) = pow(pow(p0, k) + k*c*(z1 - z0)/th0, 1.0/k);
    }
    else {
      Real ra = (z1 - z0)/(th1 - th0);
      pres(i, j) = pow(pow(p0, k) + k*c*ra*log(th1/th0), 1.0/k);
    }
  }
}

FortranData::Ptr make_standard(const Int shcol, Int nlev, Int num_qtracers) {
  using consts = scream::physics::Constants<Real>;

  const auto dp = std::make_shared<FortranData>(shcol, nlev, nlev+1,
                                                num_qtracers);

  auto& d = *dp;

  const Real ztop = 2400.0;
  const Real dz = ztop/nlev;
  for (Int i = 0; i < shcol; ++i) {
    // Set the horizontal grid spacing.
    d.host_dx(i) = 5300.0;
    d.host_dy(i) = 5300.0;

    d.zt_grid(i, 0) = 0;
    for (Int k = 0; k < nlev; ++k) {

      // Set up the vertical grid.
      Real zi1 = (k+1) * dz;
      d.zi_grid(i, k+1) = zi1;
      Real zt = (k+0.5) * dz;
      d.zt_grid(i, k) = zt;

      // Interpolate the potential temperature, introducing small variations
      // between columns.
      Real theta_zt = interpolate_data(z_ref, theta_ref, zt);
      if (i > 0)
        theta_zt += ((i % 3) - 0.5)/double(nlev)*k;
      const Real qw = interpolate_data(z_ref, qw_ref, zt);
      const Real ql = interpolate_data(z_ref, ql_ref, zt);
      d.qw(i, k) = qw;
      d.shoc_ql(i, k) = ql;
      Real zvir = (consts::RH2O / consts::Rair) - 1.0;
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
        d.wtracer_sfc(i, q) = 0.1 * i;
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
      d.pdel(i, k) = std::abs(d.presi(i, k+1) - d.presi(i, k));
    }

    // Initialize host_dse and exner, which are assumed to be valid inputs to
    // shoc_main.
    for (Int k = 0; k < nlev; ++k) {
      d.exner(i, k) = pow(d.pres(i, k)/consts::P0, consts::Rair/consts::Cpair);
      d.host_dse(i, k) = consts::Cpair * d.exner(i, k) * d.thv(i, k) +
        consts::gravit * d.zt_grid(i, k);
    }
  }

  return dp;
}

} // end anonymous namespace

// From scream-docs/shoc-port/shocintr.py.
FortranData::Ptr Factory::create (IC ic, Int shcol, Int nlev, Int num_qtracers) {
  switch (ic) {
    case standard: return make_standard(shcol, nlev, num_qtracers);
    default: scream_require_msg(false, "Not an IC: " << ic);
  }
}

} // namespace ic
} // namespace shoc
} // namespace scream
