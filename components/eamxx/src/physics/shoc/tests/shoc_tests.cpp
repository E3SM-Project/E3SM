#include "catch2/catch.hpp"

#include "shoc_main_wrap.hpp"
#include "shoc_ic_cases.hpp"

#include "ekat/util/ekat_test_utils.hpp"

namespace {

TEST_CASE("FortranData", "shoc") {
  int val = scream::shoc::test_FortranData();
  REQUIRE(val == 0);
}

TEST_CASE("FortranDataIterator", "shoc") {
  using scream::shoc::ic::Factory;
  const auto d = Factory::create(Factory::standard);
  scream::shoc::FortranDataIterator fdi(d);
  REQUIRE(fdi.nfield() == 44);
  const auto& f = fdi.getfield(0);
  REQUIRE(f.dim == 2);
  REQUIRE(f.extent[0] == d->shcol);
  REQUIRE(f.extent[1] == 1);
  REQUIRE(f.extent[2] == 1);
  REQUIRE(f.data == d->host_dx.data());
  REQUIRE(static_cast<int>(f.size) == d->shcol);
}

TEST_CASE("shoc_init_f", "shoc") {
  int nerr = scream::shoc::test_shoc_init(true);
  REQUIRE(nerr == 0);
}

// This helper returns true if we've been asked to generate Python
// plotting scripts, false otherwise.
bool generating_plot_scripts() {
  bool gen_plot_scripts = false;
  auto& ts = ekat::TestSession::get();
  auto iter = ts.params.find("gen_plot_scripts");
  if (iter != ts.params.end()) {
    // Here's val, passed as gen_plot_scripts=val
    std::string val = iter->second;
    // Low-case the thing. Isn't C++ a friendly language??
    std::transform(val.begin(), val.end(), val.begin(),
      [](unsigned char c){ return std::tolower(c); });

    // Now decide if the value is true or not. Use CMake sensibilities.
    gen_plot_scripts = ((val == "1") or (val == "true") or
                        (val == "yes") or (val == "on"));
  }
  return gen_plot_scripts;
}

TEST_CASE("shoc_ic_f", "shoc") {
  int nerr = scream::shoc::test_shoc_ic(true, generating_plot_scripts());
  REQUIRE(nerr == 0);
}

TEST_CASE("shoc_init_c", "shoc") {
  int nerr = scream::shoc::test_shoc_init(false);
  REQUIRE(nerr == 0);
}

TEST_CASE("shoc_ic_c", "shoc") {
  int nerr = scream::shoc::test_shoc_ic(false);
  REQUIRE(nerr == 0);
}

} // anonymous namespace

