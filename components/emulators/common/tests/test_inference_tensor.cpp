// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "inference_error.hpp"
#include "tensor.hpp"

#include <vector>

namespace emulator {
namespace inference {
namespace test {

TEST_CASE("An owning tensor allocates and zero-fills", "[tensor]") {
  Tensor t("dT", {4, 3});

  REQUIRE(t.name() == "dT");
  REQUIRE(t.rank() == 2);
  REQUIRE(t.dim(0) == 4);
  REQUIRE(t.size() == 12);
  REQUIRE(t.nbytes() == 12 * sizeof(double));
  REQUIRE(t.owns_data());
  REQUIRE(t.writable());
  REQUIRE(t.cdata()[7] == 0.0);
  REQUIRE(t.to_string() == "dT[4,3]");
}

TEST_CASE("A view shares the caller's memory", "[tensor]") {
  std::vector<double> field(6, 1.5);
  Tensor t = Tensor::view("T", field.data(), {2, 3});

  REQUIRE_FALSE(t.owns_data());
  REQUIRE(t.writable());
  REQUIRE(t.cdata() == field.data());

  t.data()[0] = 42.0;
  REQUIRE(field[0] == 42.0); // no copy happened anywhere
}

TEST_CASE("A const view cannot be written through", "[tensor]") {
  const std::vector<double> field(6, 1.5);
  Tensor t = Tensor::const_view("T", field.data(), {2, 3});

  REQUIRE_FALSE(t.writable());
  REQUIRE(t.cdata()[3] == 1.5);
  REQUIRE_THROWS_AS(t.data(), InferenceError);
}

TEST_CASE("clone() is the only way to deep copy", "[tensor]") {
  std::vector<double> field{1.0, 2.0, 3.0};
  Tensor view = Tensor::view("T", field.data(), {3});
  Tensor copy = view.clone();

  REQUIRE(copy.owns_data());
  copy.data()[0] = 99.0;
  REQUIRE(field[0] == 1.0); // the original is untouched
}

TEST_CASE("TensorMap keeps its order and its names distinct", "[tensor]") {
  std::vector<double> a(4, 1.0);
  std::vector<double> b(4, 2.0);

  TensorMap map;
  map.wrap("T", static_cast<const double *>(a.data()), {4});
  map.wrap("q", b.data(), {4});

  REQUIRE(map.size() == 2);

  // Order is preserved, because a positional model signature depends on it.
  std::vector<std::string> names;
  for (const auto &t : map) {
    names.push_back(t.name());
  }
  REQUIRE(names == std::vector<std::string>{"T", "q"});

  // A duplicate would replace the first in the dict handed to the model, so
  // that field would silently never arrive.
  REQUIRE_THROWS_AS(map.wrap("T", b.data(), {4}), InferenceError);
}

TEST_CASE("An empty tensor is legal, and still writable", "[tensor]") {
  // A rank owning no columns is normal on a large layout, and inference is
  // collective, so such a rank still takes part with zero-length fields.
  // std::vector::data() may be null there, so writability cannot be inferred
  // from the pointer: an empty output would arrive read-only and throw.
  Tensor t("dT", {0, 3});
  REQUIRE(t.size() == 0);
  REQUIRE(t.nbytes() == 0);
  REQUIRE(t.writable());
  // Owning nothing is still owning: allocation, not the size of the result,
  // is what decides who frees the memory.
  REQUIRE(t.owns_data());

  Tensor v = Tensor::view("dT", nullptr, {0, 3});
  REQUIRE(v.writable());
  REQUIRE(v.data() == nullptr);
  REQUIRE_FALSE(v.owns_data());

  REQUIRE_FALSE(Tensor::const_view("T", nullptr, {0, 3}).writable());

  // Copying nothing is legal too, and yields an owning tensor.
  Tensor copy = v.clone();
  REQUIRE(copy.size() == 0);
  REQUIRE(copy.owns_data());
}

} // namespace test
} // namespace inference
} // namespace emulator
