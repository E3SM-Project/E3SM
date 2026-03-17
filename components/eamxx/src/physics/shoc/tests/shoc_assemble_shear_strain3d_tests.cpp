#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_test_data.hpp"
#include "share/core/eamxx_types.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>

namespace scream {
namespace shoc {
namespace unit_test {

namespace {

inline Real reference_shear_strain3d(const Real a00, const Real a01, const Real a10,
                                     const Real a11, const Real a20, const Real a21,
                                     const Real du_dz_m, const Real dv_dz_m, const Real dw_dz_m)
{
  const Real s00 = a00;
  const Real s11 = a11;
  const Real s22 = dw_dz_m;
  const Real s01 = 0.5 * (a01 + a10);
  const Real s02 = 0.5 * (du_dz_m + a20);
  const Real s12 = 0.5 * (dv_dz_m + a21);

  return 2.0 * (s00*s00 + s11*s11 + s22*s22
              + 2.0*s01*s01 + 2.0*s02*s02 + 2.0*s12*s12);
}

} // namespace

template <typename D>
struct UnitWrap::UnitTest<D>::TestAssembleShocShearStrain3d : public UnitWrap::UnitTest<D>::Base {

  static void require_close(const Real a, const Real b, const Real scale = 1)
  {
    const Real tol = 200 * std::numeric_limits<Real>::epsilon() * scale;
    REQUIRE(a == Approx(b).margin(tol));
  }

  void run_property()
  {
    static constexpr Int shcol = 2;
    static constexpr Int nlev = 5;

    AssembleShocShearStrain3dData d(shcol, nlev);

    REQUIRE(d.shcol == shcol);
    REQUIRE(d.nlev == nlev);

    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlev; ++k) {
        const Int off2 = k + s*nlev;
        d.du_dz_m[off2] = 0;
        d.dv_dz_m[off2] = 0;
        d.dw_dz_m[off2] = 0;
        d.shear_strain3d[off2] = -1;
        for (Int j = 0; j < 6; ++j) {
          const Int off3 = k + nlev*(j + 6*s);
          d.shear_strain3d_components[off3] = 0;
        }
      }
    }

    assemble_shoc_shear_strain3d(d);

    for (Int i = 0; i < shcol*nlev; ++i) {
      REQUIRE(d.shear_strain3d[i] == 0);
    }

    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlev; ++k) {
        const Int off2 = k + s*nlev;
        d.du_dz_m[off2] = 1 + k;
        d.dv_dz_m[off2] = 2 + k;
        d.dw_dz_m[off2] = 3 + k;
        const Real vals[6] = {1, -2, 4, -3, 5, -6};
        for (Int j = 0; j < 6; ++j) {
          const Int off3 = k + nlev*(j + 6*s);
          d.shear_strain3d_components[off3] = vals[j];
        }
      }
    }

    assemble_shoc_shear_strain3d(d);

    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlev; ++k) {
        const Int off2 = k + s*nlev;
        const Real expected = reference_shear_strain3d(1, -2, 4, -3, 5, -6,
                                                       1 + k, 2 + k, 3 + k);
        require_close(d.shear_strain3d[off2], expected, expected + 1);
        REQUIRE(d.shear_strain3d[off2] >= 0);
      }
    }

    // Antisymmetric horizontal off-diagonal part should cancel.
    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlev; ++k) {
        const Int off2 = k + s*nlev;
        d.du_dz_m[off2] = 0;
        d.dv_dz_m[off2] = 0;
        d.dw_dz_m[off2] = 0;
        const Real vals[6] = {0, 7, -7, 0, 0, 0};
        for (Int j = 0; j < 6; ++j) {
          const Int off3 = k + nlev*(j + 6*s);
          d.shear_strain3d_components[off3] = vals[j];
        }
      }
    }

    assemble_shoc_shear_strain3d(d);

    for (Int i = 0; i < shcol*nlev; ++i) {
      REQUIRE(d.shear_strain3d[i] == 0);
    }

    // The vertical cross term should depend only on the symmetric average 0.5*(A20 + du_dz_m).
    for (Int s = 0; s < shcol; ++s) {
      for (Int k = 0; k < nlev; ++k) {
        const Int off2 = k + s*nlev;
        d.du_dz_m[off2] = 3;
        d.dv_dz_m[off2] = 0;
        d.dw_dz_m[off2] = 0;
        const Real vals[6] = {0, 0, 0, 0, -3, 0};
        for (Int j = 0; j < 6; ++j) {
          const Int off3 = k + nlev*(j + 6*s);
          d.shear_strain3d_components[off3] = vals[j];
        }
      }
    }

    assemble_shoc_shear_strain3d(d);

    for (Int i = 0; i < shcol*nlev; ++i) {
      REQUIRE(d.shear_strain3d[i] == 0);
    }
  }

  void run_property_random()
  {
    auto engine = Base::get_engine();

    std::uniform_int_distribution<Int> shcol_dist(1, 5);
    std::uniform_int_distribution<Int> nlev_dist(7, 19);
    std::uniform_real_distribution<Real> val_dist(-10, 10);

    for (Int trial = 0; trial < 20; ++trial) {
      const Int shcol = shcol_dist(engine);
      const Int nlev = nlev_dist(engine);

      AssembleShocShearStrain3dData d(shcol, nlev);

      for (Int s = 0; s < shcol; ++s) {
        for (Int k = 0; k < nlev; ++k) {
          const Int off2 = k + s*nlev;
          d.du_dz_m[off2] = val_dist(engine);
          d.dv_dz_m[off2] = val_dist(engine);
          d.dw_dz_m[off2] = val_dist(engine);
          for (Int j = 0; j < 6; ++j) {
            const Int off3 = k + nlev*(j + 6*s);
            d.shear_strain3d_components[off3] = val_dist(engine);
          }
        }
      }

      assemble_shoc_shear_strain3d(d);

      for (Int s = 0; s < shcol; ++s) {
        for (Int k = 0; k < nlev; ++k) {
          const Int off2 = k + s*nlev;
          const auto comp = [&](const Int j) {
            return d.shear_strain3d_components[k + nlev*(j + 6*s)];
          };
          const Real expected = reference_shear_strain3d(comp(0), comp(1), comp(2), comp(3), comp(4), comp(5),
                                                         d.du_dz_m[off2], d.dv_dz_m[off2], d.dw_dz_m[off2]);
          require_close(d.shear_strain3d[off2], expected, std::abs(expected) + 1);
          REQUIRE(d.shear_strain3d[off2] >= 0);
        }
      }
    }
  }
};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("assemble_shoc_shear_strain3d_property", "shoc")
{
  using TestStruct =
    scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestAssembleShocShearStrain3d;

  TestStruct().run_property();
  TestStruct().run_property_random();
}

} // namespace
