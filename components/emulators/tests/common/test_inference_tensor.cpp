// Catch2 v2 single header
#include "field_types.hpp"
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "inference_error.hpp"
#include "tensor.hpp"

#include <algorithm>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

namespace e3sm::emulator::inference::test {

TEST_CASE("Owning tensor allocates and zero-fills", "[tensor]") {
  Tensor t("dT", {4, 3});

  REQUIRE(t.name() == "dT");
  REQUIRE(std::ranges::equal(t.dims(), std::vector<std::int64_t>{4, 3}));
  REQUIRE(t.size() == 12);
  REQUIRE(t.writable());
  REQUIRE(t.as<double>()[11] == 0.0);
  REQUIRE(t.to_string() == "dT[4,3]");
}

TEST_CASE("Owning tensor survives a move", "[tensor]") {
  Tensor a("dT", {2});
  a.as<double>()[1] = 7.0;

  Tensor b = std::move(a);
  REQUIRE(b.as<double>()[1] == 7.0);
}

TEST_CASE("View shares the caller's memory", "[tensor]") {
  std::vector<double> field(6, 1.5);
  Tensor t = Tensor::from_view("T", field, {2, 3});

  REQUIRE(t.data() == static_cast<void*>(field.data()));

  t.as<double>()[0] = 42.0;
  REQUIRE(field[0] == 42.0);
}

TEST_CASE("Const view cannot be written through", "[tensor]") {
  const std::vector<double> field(6, 1.5);
  Tensor t = Tensor::from_const_view("T", field, {2, 3});
  REQUIRE_FALSE(t.writable());
  REQUIRE_THROWS_AS(t.as<double>()[3] == 1.5, InferenceError);
}

TEST_CASE("Invalid views are refused", "[tensor]") {
  std::vector<double> field(6, 0.0);

  REQUIRE_THROWS_AS(Tensor::from_view("T", FieldView{}, {2, 3}),
                    InferenceError);
  REQUIRE_THROWS_AS(Tensor::from_view("T", field, {-2, 3}), InferenceError);
}

TEST_CASE("Empty tensors are legal", "[tensor]") {
  // e.g. an MPI rank that owns no columns
  Tensor owned("dT", {0, 3});
  REQUIRE(owned.size() == 0);
  REQUIRE(owned.writable());

  Tensor v = Tensor::from_view("dT", FieldView{}, {0, 3});
  REQUIRE(v.writable());
  REQUIRE(v.size() == 0);

  REQUIRE_FALSE(
      Tensor::from_const_view("T", ConstFieldView{}, {0, 3}).writable());
}

TEST_CASE("TensorMap preserves order and rejects duplicates", "[tensor]") {
  const std::vector<double> a(4, 1.0);
  std::vector<double> b(4, 2.0);

  TensorMap map;
  map.wrap("T", a, {4});
  map.wrap("q", FieldView{b}, {4});
  map.add(Tensor("p", {4}));

  std::vector<std::string> names;
  for (const auto& t : map) {
    names.push_back(t.name());
  }
  REQUIRE(names == std::vector<std::string>{"T", "q", "p"});

  // const and non-const pointers select the matching view
  REQUIRE_FALSE(map.begin()->writable());
  REQUIRE((map.begin() + 1)->writable());

  REQUIRE_THROWS_AS(map.wrap("T", FieldView{b}, {4}), InferenceError);
  REQUIRE(map.size() == 3);
}

} // namespace e3sm::emulator::inference::test
