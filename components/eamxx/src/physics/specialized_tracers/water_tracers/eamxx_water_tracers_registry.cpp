#include "eamxx_water_tracers_registry.hpp"

#include <ekat_assert.hpp>
#include <ekat_std_utils.hpp>
#include <ekat_string_utils.hpp>
#include <ekat_yaml.hpp>

#include <algorithm>
#include <map>

namespace scream {
namespace water_tracers {

namespace {

// Fields that other processes already register as tracers. A water tracer
// with one of these names would silently alias them.
const std::vector<std::string> reserved_names =
  {"qv", "qc", "qr", "qi", "qm", "nc", "nr", "ni", "bm", "tke"};

const std::string flux_suffix = "_surf_flux";

const std::vector<std::string> generic_keys  = {"type", "ocean_flux", "r_ocn", "init", "r_init"};
const std::vector<std::string> isotope_keys  = {"type", "species", "init", "r_init"};
const std::vector<std::string> property_keys = {"type", "carrier", "source", "field_name",
                                                "field_units", "scale", "offset", "init", "r_init"};

const std::vector<std::string> property_sources =
  {"sst", "t_bot", "rh_sfc", "lat", "lon", "x", "y", "z", "field"};

bool ends_with (const std::string& s, const std::string& suffix) {
  return s.size()>=suffix.size() and
         s.compare(s.size()-suffix.size(), suffix.size(), suffix)==0;
}

// Context prefix for error messages
std::string where (const std::string& path, const std::string& tracer) {
  return "Error! In tracers_file '" + path + "', tracer '" + tracer + "': ";
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

  EKAT_REQUIRE_MSG(not ekat::contains(reserved_names, name) and not ends_with(name, flux_suffix),
      ctx + "the name is reserved. Tracer names cannot be one of {" +
      ekat::join(reserved_names, ", ") + "} or end in '" + flux_suffix + "'.\n");

  TracerSpec spec;
  spec.name  = name;
  spec.extra = entry;

  const auto type = require_token(entry, "type", ctx, "generic, isotope or property");
  const std::vector<std::string>* valid_keys = nullptr;
  if (type=="generic") {
    spec.type = TracerType::Generic;
    valid_keys = &generic_keys;
  } else if (type=="isotope") {
    EKAT_REQUIRE_MSG(opts.accept_isotopes,
        ctx + "type: isotope entries need the water_isotopes process, but this file\n"
        "  is being read by water_tracers. Put water_isotopes (not water_tracers) in\n"
        "  atm_procs_list, or remove the isotope entries from the file.\n");
    spec.type = TracerType::Isotope;
    valid_keys = &isotope_keys;
  } else if (type=="property") {
    spec.type = TracerType::Property;
    valid_keys = &property_keys;
  } else {
    EKAT_ERROR_MSG(ctx + "unknown type '" + type + "'. Valid types: generic, isotope, property.\n");
  }

  // Reject misspelled or misplaced keys rather than silently using defaults.
  for (const auto& key : entry.param_names()) {
    EKAT_REQUIRE_MSG(ekat::contains(*valid_keys, key),
        ctx + "unknown key '" + key + "' for a " + type + " tracer. Valid keys: " +
        ekat::join(*valid_keys, ", ") + ".\n");
  }
  EKAT_REQUIRE_MSG(entry.sublist_names().empty(),
      ctx + "unexpected nested map '" + entry.sublist_names().front() + "'.\n");

  // Initialization (default: ratio for mass tracers, zero for properties)
  spec.init = spec.is_mass_tracer() ? InitType::Ratio : InitType::Zero;
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
    case TracerType::Property:
    {
      spec.carrier = require_token(entry, "carrier", ctx, "the mass tracer that carries this property");
      spec.source  = require_token(entry, "source",  ctx, ekat::join(property_sources, ", "));
      EKAT_REQUIRE_MSG(ekat::contains(property_sources, spec.source),
          ctx + "unknown source '" + spec.source + "'. Valid sources: " +
          ekat::join(property_sources, ", ") + ".\n");

      if (spec.source=="field") {
        spec.field_name  = require_token(entry, "field_name",  ctx, "needed by source: field");
        spec.field_units = require_token(entry, "field_units", ctx, "needed by source: field");
        property_field_units(spec.field_units);  // validates the token
      } else {
        EKAT_REQUIRE_MSG(not entry.isParameter("field_name") and not entry.isParameter("field_units"),
            ctx + "field_name and field_units are only valid with source: field.\n");
      }

      get_real(entry, "scale",  spec.scale,  ctx);
      get_real(entry, "offset", spec.offset, ctx);

      // P = scale*source + offset must be >= 0. For bounded sources this can
      // be checked now; the others are checked on device in debug builds.
      Real lo, hi;
      if (bounded_source_range(spec.source, lo, hi)) {
        const Real p_min = std::min(spec.scale*lo, spec.scale*hi) + spec.offset;
        EKAT_REQUIRE_MSG(p_min>=0,
            ctx + "scale*" + spec.source + " + offset can be negative (minimum " +
            std::to_string(p_min) + " for " + spec.source + " in [" + std::to_string(lo) +
            ", " + std::to_string(hi) + "]). Property values must be >= 0.\n");
      }
      // With init: ratio, r_init is the initial property value (default 1).
      break;
    }
  }

  return spec;
}

} // anonymous namespace

ekat::units::Units property_field_units (const std::string& token)
{
  using namespace ekat::units;
  static const std::map<std::string, Units> table = {
    {"1",       none},
    {"K",       K},
    {"kg/kg",   kg/kg},
    {"m/s",     m/s},
    {"Pa",      Pa},
    {"W/m2",    W/(m*m)},
    {"kg/m2/s", kg/(m*m*s)},
    {"m",       m},
  };
  const auto it = table.find(token);
  if (it!=table.end()) {
    return it->second;
  }

  std::vector<std::string> tokens;
  for (const auto& entry : table) tokens.push_back(entry.first);
  EKAT_ERROR_MSG("Error! Unknown field_units '" + token + "'. Valid units: " +
                 ekat::join(tokens, ", ") + ".\n");
  return none;  // unreachable, unless EKAT_ERROR_MSG is hit while unwinding
}

bool bounded_source_range (const std::string& source, Real& lo, Real& hi)
{
  if (source=="lat") {
    lo = -90; hi = 90;
  } else if (source=="lon") {
    lo = 0; hi = 360;
  } else if (source=="x" or source=="y" or source=="z") {
    lo = -1; hi = 1;
  } else {
    return false;
  }
  return true;
}

std::vector<TracerSpec> parse_tracers_file (const std::string& path,
                                            const RegistryOptions& opts)
{
  const auto file = ekat::parse_yaml_file(path);
  EKAT_REQUIRE_MSG(file.isSublist("tracers"),
      "Error! tracers_file '" + path + "' needs a top-level 'tracers' map, keyed by tracer name.\n");
  const auto& tracers = file.sublist("tracers");

  const auto scalars = tracers.param_names();
  EKAT_REQUIRE_MSG(scalars.empty(),
      where(path, scalars.empty() ? "" : scalars.front()) +
      "each entry of 'tracers' must be a map, e.g. {type: generic}.\n");

  std::vector<TracerSpec> specs;
  for (const auto& name : tracers.sublist_names()) {
    specs.push_back(parse_entry(path, name, tracers.sublist(name), opts));
  }

  // A property is carried by a mass tracer defined in the same file.
  for (const auto& spec : specs) {
    if (spec.type!=TracerType::Property) continue;
    const auto carrier = std::find_if(specs.begin(), specs.end(),
        [&](const TracerSpec& s) { return s.name==spec.carrier; });
    EKAT_REQUIRE_MSG(carrier!=specs.end(),
        where(path, spec.name) + "carrier '" + spec.carrier + "' is not defined in this file.\n");
    EKAT_REQUIRE_MSG(carrier->is_mass_tracer(),
        where(path, spec.name) + "carrier '" + spec.carrier + "' is a property tracer; "
        "the carrier must be a generic or isotope tracer.\n");
  }

  return specs;
}

} // namespace water_tracers
} // namespace scream
