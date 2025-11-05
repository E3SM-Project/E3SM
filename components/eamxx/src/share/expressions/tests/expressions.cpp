#include <catch2/catch.hpp>
#include <numeric>

#include "share/expressions/base.hpp"
#include "share/expressions/binary_op.hpp"
#include "share/expressions/compare.hpp"
#include "share/expressions/conditional.hpp"
#include "share/expressions/field.hpp"

#include "share/expressions/evaluate.hpp"
#include "share/expressions/helpers.hpp"

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include "share/core/eamxx_setup_random_test.hpp"

#include <chrono>

namespace scream {

TEST_CASE("core_expressions", "") {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  std::uniform_real_distribution<Real> pdf(0, 1);
  auto engine = setup_random_test();

  FieldLayout fl ({COL},{10});

  FieldIdentifier fid ("", fl, kg, "some_grid");
  Field f1 (fid.alias("f1"),true);
  Field f2 (fid.alias("f2"),true);
  randomize(f1, engine, pdf);
  randomize(f2, engine, pdf);

  SECTION ("field") {
    auto expression = RealFieldExpression(f1);

    // Evaluate
    Field f3(fid.alias("f3"),true);
    f3.deep_copy(-1);
    expression.evaluate(f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      REQUIRE (v3h(i)==v1h(i));
    }
  }
  SECTION ("f_plus_f") {
    auto sum = f1 + f2;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    sum.evaluate(f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      REQUIRE (v3h(i)==(v1h(i)+v2h(i)));
    }
  }
  SECTION ("f_plus_prod") {
    auto expression = f1 * f2 + f1;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    expression.evaluate(f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      auto tgt = v1h(i)*v2h(i) + v1h(i);
      REQUIRE (v3h(i)==tgt);
    }
  }
  SECTION ("prod_plus_div") {
    auto expression = f1 * f2 + f1 / f2;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    expression.evaluate(f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      auto tgt = v1h(i)*v2h(i) + v1h(i)/v2h(i);
      REQUIRE (v3h(i)==tgt);
    }
  }
  SECTION ("prod_plus_div_minus_f") {
    auto expression = (f1 * f2 + f1 / f2) - f1;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    expression.evaluate(f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      auto tgt = v1h(i)*v2h(i) + v1h(i)/v2h(i) - v1h(i);
      REQUIRE (v3h(i)==tgt);
    }
  }
  SECTION ("cmp") {
    auto f1e = RealFieldExpression(f1);
    auto f2e = RealFieldExpression(f2);
    auto cmp = f1e>=f2e;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    cmp.evaluate(f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      REQUIRE (v3h(i)==static_cast<Real>(v1h(i)>=v2h(i)));
    }
  }
  SECTION ("conditional") {
    auto f1e = RealFieldExpression(f1);
    auto f2e = RealFieldExpression(f2);
    auto expression = conditional(f1e>=f2e,f1+f2,f1-f2);

    // Evaluate
    Field f3(fid.alias("f3"),true);
    expression.evaluate(f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      auto tgt = v1h(i)>=v2h(i) ? v1h(i)+v2h(i) : v1h(i)-v2h(i);
      REQUIRE (v3h(i)==tgt);
    }
  }
}

TEST_CASE("math_expressions", "") {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  std::uniform_real_distribution<Real> pdf(0.1, 1);
  auto engine = setup_random_test();

  FieldLayout fl ({COL},{10000});

  FieldIdentifier fid ("", fl, kg, "some_grid");
  Field f1 (fid.alias("f1"),true);
  Field f2 (fid.alias("f2"),true);
  randomize(f1, engine, pdf);
  randomize(f2, engine, pdf);

  SECTION ("exp") {
    // auto expression = exp(f1) + log(f2)*sin(f1);
    auto expression = f1*f2 - f1/f2;
    // Evaluate
    Field f3(fid.alias("f3"),true);
    f3.deep_copy(-1);
    double tot_exp = 0;
    double tot_man = 0;
    double tot_fop = 0;

    Field tmp = f3.clone();

    for (int i=0; i<1000; ++i) {
      auto beg = std::chrono::high_resolution_clock::now();
      expression.evaluate(f3);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> t_exp = end - beg;
      tot_exp += t_exp.count();

      auto v1d = f1.get_view<const Real*>();
      auto v2d = f2.get_view<const Real*>();
      auto v3d = f3.get_view<Real*>();
      auto manual = KOKKOS_LAMBDA(int i) {
        // v3d(i) = Kokkos::exp(v1d(i)) + Kokkos::log(v2d(i))*Kokkos::sin(v1d(i));
        v3d(i) = v1d(i)*v2d(i) - v1d(i)/v2d(i);
      };
      Kokkos::RangePolicy<> p(0,v1d.size());
      beg = std::chrono::high_resolution_clock::now();
      Kokkos::parallel_for(p,manual);
      end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> t_man = end - beg;
      tot_man += t_man.count();

      beg = std::chrono::high_resolution_clock::now();
      f3.deep_copy(f1);
      f3.update<CombineMode::Multiply>(f2,1,1);
      tmp.deep_copy(f1);
      tmp.update<CombineMode::Divide>(f2,1,1);
      f3.update(tmp,-1,1);
      end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> t_fop = end - beg;
      tot_fop += t_fop.count();
    }

    std::cout << "expressions: " << tot_exp << "\n";
    std::cout << "manual loop: " << tot_man << "\n";
    std::cout << "field ops  : " << tot_fop << "\n";
    f1.sync_to_host();
    f2.sync_to_host();
    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      auto x = v1h(i);
      auto y = v2h(i);
      auto z = v3h(i);
      REQUIRE (z==(std::exp(x) + std::log(y)*std::sin(x)));
    }
  }
}

} // namespace scream
