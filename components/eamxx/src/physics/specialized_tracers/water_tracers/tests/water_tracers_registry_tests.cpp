#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_registry.hpp"

#include <fstream>
#include <string>

namespace scream {
namespace {

namespace wt = water_tracers;
using Catch::Matchers::Contains;

// Writes a fixture into the working directory and returns its name
std::string write_file (const std::string& name, const std::string& body) {
  std::ofstream(name) << body;
  return name;
}

std::string tracers (const std::string& entries) {
  return "%YAML 1.1\n---\ntracers:\n" + entries + "...\n";
}

const wt::RegistryOptions no_iso{};
const wt::RegistryOptions iso{true, {"H216O", "HDO", "H218O", "H217O", "HTO"}};

TEST_CASE("registry_valid") {
  const auto flow = write_file("wtrc_flow.yaml", tracers(
    "  wtrc_bp:  {type: generic, ocean_flux: bulk_proportional, r_ocn: 1}\n"
    "  wtrc_bf:  {type: generic, ocean_flux: bulk_formula, r_ocn: 0.5, init: ratio}\n"
    "  wtrc_hdo: {type: isotope, species: HDO}\n"));
  const auto block = write_file("wtrc_block.yaml", tracers(
    "  wtrc_bp:\n    type: generic\n"));

  const auto specs = wt::parse_tracers_file(flow, iso);
  REQUIRE(specs.size() == 3);
  // Sorted by name
  REQUIRE(specs[0].name == "wtrc_bf");
  REQUIRE(specs[1].name == "wtrc_bp");
  REQUIRE(specs[2].name == "wtrc_hdo");

  REQUIRE(specs[0].ocean_flux == wt::OceanFlux::BulkFormula);
  REQUIRE(specs[0].r_ocn == Real(0.5));
  REQUIRE(specs[0].r_init == Real(0.5));  // default r_init = r_ocn

  REQUIRE(specs[1].r_ocn == 1);            // int in the file
  REQUIRE(specs[1].ocean_flux == wt::OceanFlux::BulkProportional);

  REQUIRE(specs[2].type == wt::TracerType::Isotope);
  REQUIRE(specs[2].species == "HDO");
  REQUIRE(specs[2].init == wt::InitType::Ratio);
  REQUIRE(not specs[2].r_init_given);

  const auto bspecs = wt::parse_tracers_file(block, no_iso);
  REQUIRE(bspecs.size() == 1);
  REQUIRE(bspecs[0].name == "wtrc_bp");
  REQUIRE(bspecs[0].r_ocn == 1);           // default
  REQUIRE(bspecs[0].init == wt::InitType::Ratio);
}

TEST_CASE("registry_errors") {
  auto fails = [](const std::string& tag, const std::string& entries, const std::string& msg,
                  const wt::RegistryOptions& opts = no_iso) {
    const auto f = write_file("wtrc_bad_" + tag + ".yaml", tracers(entries));
    REQUIRE_THROWS_WITH(wt::parse_tracers_file(f, opts), Contains(msg));
  };

  fails("type",    "  wtrc_a: {type: bogus}\n", "unknown type");
  fails("notype",  "  wtrc_a: {r_ocn: 1}\n", "missing required key 'type'");
  fails("flux",    "  wtrc_a: {type: generic, ocean_flux: bogus}\n", "unknown ocean_flux");
  fails("key",     "  wtrc_a: {type: generic, rocn: 1}\n", "unknown key 'rocn'");
  fails("init",    "  wtrc_a: {type: generic, init: bogus}\n", "unknown init");
  fails("iso",     "  wtrc_a: {type: isotope, species: HDO}\n", "EAMXX_ENABLE_WATER_ISOTOPES");
  fails("iso2",    "  wtrc_a: {type: isotope, species: HDO}\n", "Rebuild with isotope support");
  fails("species", "  wtrc_a: {type: isotope, species: D2O}\n", "unknown isotope species", iso);
  fails("nospec",  "  wtrc_a: {type: isotope}\n", "missing required key 'species'", iso);
  fails("pfx_qv",  "  qv: {type: generic}\n", "must start with 'wtrc_'");
  fails("pfx_sst", "  sst: {type: generic}\n", "must start with 'wtrc_'");
  fails("pfx_iso", "  hdo: {type: isotope, species: HDO}\n", "must start with 'wtrc_'", iso);
  fails("pfx_bare", "  wtrc_: {type: generic}\n", "must start with 'wtrc_'");
  fails("sfx",     "  wtrc_a_surf_flux: {type: generic}\n", "cannot end in '_surf_flux'");

  const auto notracers = write_file("wtrc_bad_top.yaml", "%YAML 1.1\n---\nfoo: 1\n...\n");
  REQUIRE_THROWS_WITH(wt::parse_tracers_file(notracers, no_iso), Contains("top-level 'tracers'"));
  const auto notmap = write_file("wtrc_bad_entry.yaml", tracers("  wtrc_a: 1\n"));
  REQUIRE_THROWS_WITH(wt::parse_tracers_file(notmap, no_iso), Contains("must be a map"));
}

} // anonymous namespace
} // namespace scream
