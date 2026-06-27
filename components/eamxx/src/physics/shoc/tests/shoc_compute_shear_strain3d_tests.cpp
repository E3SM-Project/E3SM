#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "share/core/eamxx_types.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <thread>
#include <vector>

namespace scream {
namespace shoc {
namespace unit_test {

namespace {

inline Real interp_interface_to_midpoint(const std::vector<Real>& x1,
                                         const std::vector<Real>& y1,
                                         const std::vector<Real>& x2,
                                         const Int k2)
{
  const Int km1 = x1.size();
  Int idx = k2 + 1;
  if (idx >= km1) {
    idx = km1 - 1;
  }

  const Real x = x1[idx];
  const Real xs = x1[idx-1];
  const Real y = y1[idx];
  const Real ys = y1[idx-1];

  return ys + (y - ys) * (x2[k2] - xs) / (x - xs);
}

inline std::vector<Real> reference_vertical_shear_component(const std::vector<Real>& dz_zi,
                                                            const std::vector<Real>& field,
                                                            const std::vector<Real>& zt_grid,
                                                            const std::vector<Real>& zi_grid)
{
  const Int nlev = field.size();
  const Int nlevi = dz_zi.size();

  std::vector<Real> grad_i(nlevi, 0);
  std::vector<Real> grad_m(nlev, 0);

  for (Int k = 1; k < nlev; ++k) {
    grad_i[k] = (field[k-1] - field[k]) / dz_zi[k];
  }

  grad_i[0] = 0;
  grad_i[nlevi-1] = 0;

  for (Int k = 0; k < nlev; ++k) {
    grad_m[k] = interp_interface_to_midpoint(zi_grid, grad_i, zt_grid, k);
    grad_m[k] = std::max(grad_m[k], Real(0));
  }

  return grad_m;
}

inline void build_midpoint_grid_from_interfaces(const std::vector<Real>& zi_grid,
                                                std::vector<Real>& zt_grid)
{
  const Int nlev = zt_grid.size();
  for (Int k = 0; k < nlev; ++k) {
    zt_grid[k] = 0.5 * (zi_grid[k] + zi_grid[k+1]);
  }
}

inline void build_dz_zi_from_midpoints(const std::vector<Real>& zt_grid,
                                       std::vector<Real>& dz_zi)
{
  const Int nlev = zt_grid.size();
  dz_zi[0] = 0;
  for (Int k = 1; k < nlev; ++k) {
    dz_zi[k] = zt_grid[k-1] - zt_grid[k];
  }
  dz_zi[nlev] = 0;
}

} // namespace

template <typename D>
struct UnitWrap::UnitTest<D>::TestComputeVerticalShearTerms : public UnitWrap::UnitTest<D>::Base {
  static void require_close(const Real a, const Real b, const Real scale = 1)
  {
    const Real tol = 100 * std::numeric_limits<Real>::epsilon() * scale;
    REQUIRE(a == Approx(b).margin(tol));
  }

  void run_property()
  {
    static constexpr Int shcol  = 2;
    static constexpr Int nlev   = 5;
    static constexpr Int nlevi  = nlev + 1;

    ComputeVerticalShearTermsData d(shcol, nlev, nlevi);

    const std::vector<Real> dz_zi = {0, 500, 200, 100, 50, 0};
    const std::vector<Real> zi_grid = {860, 360, 160, 60, 10, 0};
    std::vector<Real> zt_grid(nlev);
    build_midpoint_grid_from_interfaces(zi_grid, zt_grid);

    const std::vector<Real> u_wind_shr = {2, 1, 0, -1, -2};
    const std::vector<Real> v_wind_shr = {1, 2, 3, 4, 5};
    const std::vector<Real> w_field_shr = {5, 4, 3, 2, 1};

    const auto exp_du = reference_vertical_shear_component(dz_zi, u_wind_shr, zt_grid, zi_grid);
    const auto exp_dv = reference_vertical_shear_component(dz_zi, v_wind_shr, zt_grid, zi_grid);
    const auto exp_dw = reference_vertical_shear_component(dz_zi, w_field_shr, zt_grid, zi_grid);

    REQUIRE(d.shcol == shcol);
    REQUIRE(d.nlev == nlev);
    REQUIRE(d.nlevi == nlevi);

    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlevi; ++k) {
        const auto offset = k + s * nlevi;
        d.dz_zi[offset] = dz_zi[k];
        d.zi_grid[offset] = zi_grid[k];
      }

      for (Int k = 0; k < nlev; ++k) {
        const auto offset = k + s * nlev;
        d.u_wind[offset] = u_wind_shr[k];
        d.v_wind[offset] = v_wind_shr[k];
        d.w_field[offset] = w_field_shr[k];
        d.zt_grid[offset] = zt_grid[k];
      }
    }

    compute_vertical_shear_terms(d);

    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlev; ++k) {
        const auto offset = k + s * nlev;
        require_close(d.du_dz_m[offset], exp_du[k], std::abs(exp_du[k]) + 1);
        require_close(d.dv_dz_m[offset], exp_dv[k], std::abs(exp_dv[k]) + 1);
        require_close(d.dw_dz_m[offset], exp_dw[k], std::abs(exp_dw[k]) + 1);

        REQUIRE(d.du_dz_m[offset] > 0);
        REQUIRE(d.dv_dz_m[offset] == 0);
        REQUIRE(d.dw_dz_m[offset] > 0);
      }

      for (Int k = 0; k < nlev - 2; ++k) {
        const auto offset = k + s * nlev;
        REQUIRE(std::abs(d.du_dz_m[offset]) < std::abs(d.du_dz_m[offset+1]));
        if (d.dv_dz_m[offset+1] > 0) {
          REQUIRE(std::abs(d.dv_dz_m[offset]) < std::abs(d.dv_dz_m[offset+1]));
        }
        REQUIRE(std::abs(d.dw_dz_m[offset]) < std::abs(d.dw_dz_m[offset+1]));
      }

      REQUIRE(std::abs(d.du_dz_m[s*nlev + nlev-1]) < std::abs(d.du_dz_m[s*nlev + nlev-2]));
      if (d.dv_dz_m[s*nlev + nlev-2] > 0) {
        REQUIRE(std::abs(d.dv_dz_m[s*nlev + nlev-1]) < std::abs(d.dv_dz_m[s*nlev + nlev-2]));
      }
      REQUIRE(std::abs(d.dw_dz_m[s*nlev + nlev-1]) < std::abs(d.dw_dz_m[s*nlev + nlev-2]));
    }

    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlev; ++k) {
        const auto offset = k + s * nlev;
        d.u_wind[offset] = 10;
        d.v_wind[offset] = -5;
        d.w_field[offset] = 3;
      }
    }

    compute_vertical_shear_terms(d);

    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlev; ++k) {
        const auto offset = k + s * nlev;
        REQUIRE(d.du_dz_m[offset] == 0);
        REQUIRE(d.dv_dz_m[offset] == 0);
        REQUIRE(d.dw_dz_m[offset] == 0);
      }
    }

    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlev; ++k) {
        const auto offset = k + s * nlev;
        d.u_wind[offset] = 7;
        d.v_wind[offset] = -2;
        d.w_field[offset] = Real(k);
      }
    }

    compute_vertical_shear_terms(d);

    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlev; ++k) {
        const auto offset = k + s * nlev;
        REQUIRE(d.du_dz_m[offset] == 0);
        REQUIRE(d.dv_dz_m[offset] == 0);
        REQUIRE(d.dw_dz_m[offset] == 0);
      }
    }
  }

  void run_property_random()
  {
    auto engine = Base::get_engine();

    std::uniform_int_distribution<Int> shcol_dist(1, 6);
    std::uniform_int_distribution<Int> nlev_dist(7, 19);
    std::uniform_real_distribution<Real> dz_dist(10, 300);
    std::uniform_real_distribution<Real> frac_dist(0.2, 0.8);
    std::uniform_real_distribution<Real> field_dist(-25, 25);

    for (Int trial = 0; trial < 20; ++trial) {
      const Int shcol = shcol_dist(engine);
      const Int nlev = nlev_dist(engine);
      const Int nlevi = nlev + 1;

      ComputeVerticalShearTermsData d(shcol, nlev, nlevi);
      std::vector<std::vector<Real>> exp_du(shcol);
      std::vector<std::vector<Real>> exp_dv(shcol);
      std::vector<std::vector<Real>> exp_dw(shcol);

      for (Int s = 0; s < shcol; ++s) {
        std::vector<Real> zi_grid(nlevi);
        std::vector<Real> zt_grid(nlev);
        std::vector<Real> dz_zi(nlevi);
        std::vector<Real> u_wind(nlev);
        std::vector<Real> v_wind(nlev);
        std::vector<Real> w_field(nlev);

        zi_grid[nlevi-1] = 0;
        for (Int k = nlevi - 2; k >= 0; --k) {
          zi_grid[k] = zi_grid[k+1] + dz_dist(engine);
        }

        for (Int k = 0; k < nlev; ++k) {
          const Real upper = zi_grid[k];
          const Real lower = zi_grid[k+1];
          zt_grid[k] = lower + frac_dist(engine) * (upper - lower);
        }

        build_dz_zi_from_midpoints(zt_grid, dz_zi);

        for (Int k = 0; k < nlev; ++k) {
          u_wind[k] = field_dist(engine);
          v_wind[k] = field_dist(engine);
          w_field[k] = field_dist(engine);
        }

        exp_du[s] = reference_vertical_shear_component(dz_zi, u_wind, zt_grid, zi_grid);
        exp_dv[s] = reference_vertical_shear_component(dz_zi, v_wind, zt_grid, zi_grid);
        exp_dw[s] = reference_vertical_shear_component(dz_zi, w_field, zt_grid, zi_grid);

        for (Int k = 0; k < nlevi; ++k) {
          const auto offset = k + s * nlevi;
          d.dz_zi[offset] = dz_zi[k];
          d.zi_grid[offset] = zi_grid[k];
        }

        for (Int k = 0; k < nlev; ++k) {
          const auto offset = k + s * nlev;
          d.u_wind[offset] = u_wind[k];
          d.v_wind[offset] = v_wind[k];
          d.w_field[offset] = w_field[k];
          d.zt_grid[offset] = zt_grid[k];
        }
      }

      compute_vertical_shear_terms(d);

      for (Int s = 0; s < shcol; ++s) {
        for (Int k = 0; k < nlev; ++k) {
          const auto offset = k + s * nlev;
          require_close(d.du_dz_m[offset], exp_du[s][k], std::abs(exp_du[s][k]) + 1);
          require_close(d.dv_dz_m[offset], exp_dv[s][k], std::abs(exp_dv[s][k]) + 1);
          require_close(d.dw_dz_m[offset], exp_dw[s][k], std::abs(exp_dw[s][k]) + 1);
        }
      }
    }
  }
};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("compute_vertical_shear_terms_property", "shoc")
{
  using TestStruct =
    scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeVerticalShearTerms;

  TestStruct().run_property();
  TestStruct().run_property_random();
}

} // namespace
