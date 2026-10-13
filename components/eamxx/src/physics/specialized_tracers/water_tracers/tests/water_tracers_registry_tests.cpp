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
    "  wtrc_bp: {type: generic, ocean_flux: bulk_proportional, r_ocn: 1}\n"
    "  wtrc_bf: {type: generic, ocean_flux: bulk_formula, r_ocn: 0.5, init: ratio}\n"
    "  hdo:     {type: isotope, species: HDO}\n"
    "  p_sst:   {type: property, carrier: wtrc_bf, source: sst, scale: 1.0, offset: 0.0}\n"
    "  p_lat:   {type: property, carrier: hdo, source: lat, offset: 90}\n"));
  const auto block = write_file("wtrc_block.yaml", tracers(
    "  wtrc_bp:\n    type: generic\n"
    "  p_fld:\n    type: property\n    carrier: wtrc_bp\n    source: field\n"
    "    field_name: T_2m\n    field_units: K\n    init: ratio\n    r_init: 2\n"));

  const auto specs = wt::parse_tracers_file(flow, iso);
  REQUIRE(specs.size() == 5);
  // Sorted by name
  REQUIRE(specs[0].name == "hdo");
  REQUIRE(specs[1].name == "p_lat");
  REQUIRE(specs[2].name == "p_sst");
  REQUIRE(specs[3].name == "wtrc_bf");
  REQUIRE(specs[4].name == "wtrc_bp");

  REQUIRE(specs[0].type == wt::TracerType::Isotope);
  REQUIRE(specs[0].species == "HDO");
  REQUIRE(specs[0].init == wt::InitType::Ratio);
  REQUIRE(not specs[0].r_init_given);

  REQUIRE(specs[1].type == wt::TracerType::Property);
  REQUIRE(specs[1].offset == 90);
  REQUIRE(specs[1].scale == 1);
  REQUIRE(specs[1].init == wt::InitType::Zero);

  REQUIRE(specs[3].ocean_flux == wt::OceanFlux::BulkFormula);
  REQUIRE(specs[3].r_ocn == Real(0.5));
  REQUIRE(specs[3].r_init == Real(0.5));  // default r_init = r_ocn

  REQUIRE(specs[4].r_ocn == 1);            // int in the file
  REQUIRE(specs[4].ocean_flux == wt::OceanFlux::BulkProportional);

  const auto bspecs = wt::parse_tracers_file(block, no_iso);
  REQUIRE(bspecs.size() == 2);
  REQUIRE(bspecs[0].name == "p_fld");
  REQUIRE(bspecs[0].source == "field");
  REQUIRE(bspecs[0].field_name == "T_2m");
  REQUIRE(bspecs[0].field_units == "K");
  REQUIRE(bspecs[0].init == wt::InitType::Ratio);
  REQUIRE(bspecs[0].r_init == 2);
  REQUIRE(bspecs[1].r_ocn == 1);           // default
  REQUIRE(bspecs[1].init == wt::InitType::Ratio);
}

TEST_CASE("registry_errors") {
  auto fails = [](const std::string& tag, const std::string& entries, const std::string& msg,
                  const wt::RegistryOptions& opts = no_iso) {
    const auto f = write_file("wtrc_bad_" + tag + ".yaml", tracers(entries));
    REQUIRE_THROWS_WITH(wt::parse_tracers_file(f, opts), Contains(msg));
  };

  fails("type",    "  a: {type: bogus}\n", "unknown type");
  fails("notype",  "  a: {r_ocn: 1}\n", "missing required key 'type'");
  fails("flux",    "  a: {type: generic, ocean_flux: bogus}\n", "unknown ocean_flux");
  fails("key",     "  a: {type: generic, rocn: 1}\n", "unknown key 'rocn'");
  fails("init",    "  a: {type: generic, init: bogus}\n", "unknown init");
  fails("iso",     "  a: {type: isotope, species: HDO}\n", "water_isotopes");
  fails("iso2",    "  a: {type: isotope, species: HDO}\n", "atm_procs_list");
  fails("species", "  a: {type: isotope, species: D2O}\n", "unknown isotope species", iso);
  fails("nospec",  "  a: {type: isotope}\n", "missing required key 'species'", iso);
  fails("rsv_qv",  "  qv: {type: generic}\n", "reserved");
  fails("rsv_tke", "  tke: {type: generic}\n", "reserved");
  fails("rsv_sfx", "  a_surf_flux: {type: generic}\n", "reserved");
  fails("source",  "  a: {type: generic}\n  p: {type: property, carrier: a, source: bogus}\n", "unknown source");
  fails("nosrc",   "  a: {type: generic}\n  p: {type: property, carrier: a}\n", "missing required key 'source'");
  fails("nocar",   "  p: {type: property, source: sst}\n", "missing required key 'carrier'");
  fails("carrier", "  p: {type: property, carrier: zz, source: sst}\n", "not defined");
  fails("pofp",    "  a: {type: generic}\n  p: {type: property, carrier: a, source: sst}\n"
                   "  q: {type: property, carrier: p, source: sst}\n", "is a property tracer");
  fails("range",   "  a: {type: generic}\n  p: {type: property, carrier: a, source: lat}\n", "can be negative");
  fails("range2",  "  a: {type: generic}\n  p: {type: property, carrier: a, source: z, scale: -2, offset: 1}\n",
                   "can be negative");
  fails("units",   "  a: {type: generic}\n  p: {type: property, carrier: a, source: field,"
                   " field_name: T_2m, field_units: degC}\n", "Unknown field_units");
  fails("fname",   "  a: {type: generic}\n  p: {type: property, carrier: a, source: field, field_units: K}\n",
                   "missing required key 'field_name'");
  fails("funits",  "  a: {type: generic}\n  p: {type: property, carrier: a, source: field, field_name: T_2m}\n",
                   "missing required key 'field_units'");

  const auto notracers = write_file("wtrc_bad_top.yaml", "%YAML 1.1\n---\nfoo: 1\n...\n");
  REQUIRE_THROWS_WITH(wt::parse_tracers_file(notracers, no_iso), Contains("top-level 'tracers'"));
  const auto notmap = write_file("wtrc_bad_entry.yaml", tracers("  a: 1\n"));
  REQUIRE_THROWS_WITH(wt::parse_tracers_file(notmap, no_iso), Contains("must be a map"));

  // A bounded source shifted into range is fine; so is lon (always >= 0)
  const auto ok = write_file("wtrc_ok_range.yaml", tracers(
    "  a: {type: generic}\n"
    "  p: {type: property, carrier: a, source: z, scale: -2, offset: 2}\n"
    "  q: {type: property, carrier: a, source: lon}\n"
    "  u: {type: property, carrier: a, source: field, field_name: ocnfrac, field_units: 1}\n"));
  REQUIRE(wt::parse_tracers_file(ok, no_iso).size() == 4);
}

TEST_CASE("registry_field_units") {
  using namespace ekat::units;
  REQUIRE(wt::property_field_units("K") == K);
  REQUIRE(wt::property_field_units("kg/m2/s") == kg/(m*m*s));
  REQUIRE(wt::property_field_units("1") == none);
  REQUIRE_THROWS(wt::property_field_units("degC"));
}

} // anonymous namespace
} // namespace scream
