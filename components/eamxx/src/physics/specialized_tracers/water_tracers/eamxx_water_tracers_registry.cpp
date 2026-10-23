#include "eamxx_water_tracers_registry.hpp"

#include <ekat_assert.hpp>
#include <ekat_std_utils.hpp>
#include <ekat_string_utils.hpp>
#include <ekat_yaml.hpp>

namespace scream {
namespace water_tracers {

namespace {

// Every water tracer name must start with this prefix, so that it cannot
// alias a field registered by this or any other process (qv, sst, o3, ...).
const std::string name_prefix = "wtrc_";

// Each tracer also registers <name>_surf_flux, so no tracer name may end in
// it (wtrc_a_surf_flux would collide with the flux field of wtrc_a).
const std::string flux_suffix = "_surf_flux";

const std::vector<std::string> generic_keys  = {"type", "ocean_flux", "r_ocn", "init", "r_init"};
const std::vector<std::string> isotope_keys  = {"type", "species", "init", "r_init"};

bool ends_with (const std::string& s, const std::string& suffix) {
  return s.size()>=suffix.size() and
         s.compare(s.size()-suffix.size(), suffix.size(), suffix)==0;
}

// Context prefix for error messages
std::string where (const std::string& path, const std::string& tracer) {
  return "Error! In water_tracers_file '" + path + "', tracer '" + tracer + "': ";
}

// ekat's YAML parser stores a scalar as the first of bool, int, double,
// string that fits it, so `r_ocn: 1` is an int and `r_ocn: 1.0` a double.
// Returns false if the key is absent.
bool get_real (const ekat::ParameterList& entry, const std::string& key,
               Real& value, const std::string& ctx)
{
  if (not entry.isParameter(key)) return false;
  if (entry.isType<int>(key)) {
    value = static_cast<Real>(entry.get<int>(key));
  } else if (entry.isType<double>(key)) {
    value = static_cast<Real>(entry.get<double>(key));
  } else {
    EKAT_ERROR_MSG(ctx + "'" + key + "' must be a number.\n");
  }
  return true;
}

// String-valued keys. An int is accepted and converted, so that the units
// token `1` (parsed as an int) can be used. Returns "" if the key is absent.
std::string get_token (const ekat::ParameterList& entry, const std::string& key,
                       const std::string& ctx)
{
  if (not entry.isParameter(key)) return "";
  if (entry.isType<std::string>(key)) return entry.get<std::string>(key);
  if (entry.isType<int>(key)) return std::to_string(entry.get<int>(key));
  EKAT_ERROR_MSG(ctx + "'" + key + "' must be a string.\n");
  return "";
}

std::string require_token (const ekat::ParameterList& entry, const std::string& key,
                           const std::string& ctx, const std::string& why)
{
  const auto value = get_token(entry, key, ctx);
  EKAT_REQUIRE_MSG(value!="", ctx + "missing required key '" + key + "' (" + why + ").\n");
  return value;
}

TracerSpec parse_entry (const std::string& path, const std::string& name,
                        const ekat::ParameterList& entry, const RegistryOptions& opts)
{
  const auto ctx = where(path, name);

  EKAT_REQUIRE_MSG(name.size()>name_prefix.size() and name.compare(0, name_prefix.size(), name_prefix)==0,
      ctx + "tracer names must start with '" + name_prefix + "' (e.g. " + name_prefix + "hdo).\n");
  EKAT_REQUIRE_MSG(not ends_with(name, flux_suffix),
      ctx + "tracer names cannot end in '" + flux_suffix + "'.\n");

  TracerSpec spec;
  spec.name  = name;
  spec.extra = entry;

  const auto type = require_token(entry, "type", ctx, "generic or isotope");
  const std::vector<std::string>* valid_keys = nullptr;
  if (type=="generic") {
    spec.type = TracerType::Generic;
    valid_keys = &generic_keys;
  } else if (type=="isotope") {
    EKAT_REQUIRE_MSG(opts.accept_isotopes,
        ctx + "type: isotope entries require EAMXX_ENABLE_WATER_ISOTOPES=ON at build\n"
        "  time; this build was compiled with it OFF. Rebuild with isotope support\n"
        "  enabled, or remove the isotope entries from the tracers file.\n");
    spec.type = TracerType::Isotope;
    valid_keys = &isotope_keys;
  } else {
    EKAT_ERROR_MSG(ctx + "unknown type '" + type + "'. Valid types: generic, isotope.\n");
  }

  // Reject misspelled or misplaced keys rather than silently using defaults.
  for (const auto& key : entry.param_names()) {
    EKAT_REQUIRE_MSG(ekat::contains(*valid_keys, key),
        ctx + "unknown key '" + key + "' for a " + type + " tracer. Valid keys: " +
        ekat::join(*valid_keys, ", ") + ".\n");
  }
  EKAT_REQUIRE_MSG(entry.sublist_names().empty(),
      ctx + "unexpected nested map '" + entry.sublist_names().front() + "'.\n");

  // Initialization (default: ratio)
  spec.init = InitType::Ratio;
  const auto init = get_token(entry, "init", ctx);
  if (init=="ratio") {
    spec.init = InitType::Ratio;
  } else if (init=="zero") {
    spec.init = InitType::Zero;
  } else if (init!="") {
    EKAT_ERROR_MSG(ctx + "unknown init '" + init + "'. Valid values: ratio, zero.\n");
  }
  spec.r_init_given = get_real(entry, "r_init", spec.r_init, ctx);
  EKAT_REQUIRE_MSG(spec.r_init>=0, ctx + "r_init must be >= 0.\n");

  switch (spec.type) {
    case TracerType::Generic:
    {
      get_real(entry, "r_ocn", spec.r_ocn, ctx);
      EKAT_REQUIRE_MSG(spec.r_ocn>=0, ctx + "r_ocn must be >= 0.\n");
      if (not spec.r_init_given) {
        spec.r_init = spec.r_ocn;
      }

      const auto flux = get_token(entry, "ocean_flux", ctx);
      if (flux=="" or flux=="bulk_proportional") {
        spec.ocean_flux = OceanFlux::BulkProportional;
      } else if (flux=="bulk_formula") {
        spec.ocean_flux = OceanFlux::BulkFormula;
      } else {
        EKAT_ERROR_MSG(ctx + "unknown ocean_flux '" + flux + "'. "
                       "Valid values: bulk_proportional, bulk_formula.\n");
      }
      break;
    }
    case TracerType::Isotope:
    {
      // r_ocn (and the default r_init) come from the isotope parameter
      // tables; the water_isotopes process sets them after parsing.
      spec.species = require_token(entry, "species", ctx, "the isotopologue");
      EKAT_REQUIRE_MSG(ekat::contains(opts.isotope_species, spec.species),
          ctx + "unknown isotope species '" + spec.species + "'. Valid species: " +
          ekat::join(opts.isotope_species, ", ") + ".\n");
      break;
    }
  }

  return spec;
}

} // anonymous namespace

std::vector<TracerSpec> parse_tracers_file (const std::string& path,
                                            const RegistryOptions& opts)
{
  const auto file = ekat::parse_yaml_file(path);
  EKAT_REQUIRE_MSG(file.isSublist("tracers"),
      "Error! water_tracers_file '" + path + "' needs a top-level 'tracers' map, keyed by tracer name.\n");
  const auto& tracers = file.sublist("tracers");

  const auto scalars = tracers.param_names();
  EKAT_REQUIRE_MSG(scalars.empty(),
      where(path, scalars.empty() ? "" : scalars.front()) +
      "each entry of 'tracers' must be a map, e.g. {type: generic}.\n");

  std::vector<TracerSpec> specs;
  for (const auto& name : tracers.sublist_names()) {
    specs.push_back(parse_entry(path, name, tracers.sublist(name), opts));
  }

  return specs;
}

} // namespace water_tracers
} // namespace scream
