#ifndef EAMXX_WATER_TRACERS_REGISTRY_HPP
#define EAMXX_WATER_TRACERS_REGISTRY_HPP

#include "eamxx_water_tracers_functions.hpp"  // OceanFlux

#include "share/core/eamxx_types.hpp"

#include <ekat_parameter_list.hpp>
#include <ekat_units.hpp>

#include <string>
#include <vector>

namespace scream {
namespace water_tracers {

/*
 * Host-only parser for the water tracer YAML file (the tracers_file parameter,
 * shared by the water_tracers and water_isotopes processes).
 *
 * The file has one top-level key, `tracers`, a map keyed by tracer name.
 * Flow style and block style are both accepted:
 *
 *   tracers:
 *     wtrc_bp: {type: generic, ocean_flux: bulk_proportional, r_ocn: 1.0}
 *     wtrc_bf: {type: generic, ocean_flux: bulk_formula, r_ocn: 1.0, init: ratio}
 *     hdo:     {type: isotope, species: HDO}
 *     p_sst:   {type: property, carrier: wtrc_bf, source: sst, scale: 1.0, offset: 0.0}
 *     p_t2m:
 *       type: property
 *       carrier: wtrc_bp
 *       source: field
 *       field_name: T_2m
 *       field_units: K
 *
 * Keys by type (anything else is an error, so a misspelled key cannot be
 * silently ignored):
 *   generic:  type, ocean_flux (bulk_proportional [default] | bulk_formula),
 *             r_ocn (default 1), init (ratio [default] | zero),
 *             r_init (default r_ocn)
 *   isotope:  type, species (required), init (default ratio), r_init (default:
 *             the species' ocean ratio, set by the water_isotopes process)
 *   property: type, carrier (required; a generic or isotope tracer in this
 *             file), source (required; sst, t_bot, rh_sfc, lat, lon, x, y, z,
 *             or field), field_name and field_units (required iff source is
 *             field), scale (default 1), offset (default 0), init (zero
 *             [default] | ratio), r_init (the initial property value for
 *             init: ratio, default 1)
 *
 * The returned specs are sorted by tracer name: ekat::ParameterList stores
 * sublists in a std::map. Nothing downstream depends on the order.
 */

// What the calling process can handle.
struct RegistryOptions {
  bool accept_isotopes = false;
  // Valid `species:` values for isotope entries (used only if accept_isotopes)
  std::vector<std::string> isotope_species;
};

enum class TracerType { Generic, Isotope, Property };

// How a tracer is initialized on an Initial run:
//  - Ratio: q_i = r_init*qv for mass tracers; q_p = r_init*q_c for property
//    tracers (q_c = the carrier after its own initialization).
//  - Zero:  the tracer keeps the initial condition the driver gave it (buildnml
//    sets initial_conditions::<name> to 0.0 unless the user sets it).
enum class InitType { Ratio, Zero };

struct TracerSpec {
  std::string name;
  TracerType  type = TracerType::Generic;

  // Isotope tracers only: one of RegistryOptions::isotope_species
  std::string species;

  // Ocean ratio R_ocn. Generic: from the file (default 1). Isotope: set by the
  // water_isotopes process from its parameter tables. Unused for properties.
  Real r_ocn = 1;

  // Initialization (see InitType). r_init_given records whether the file set
  // r_init, so a subclass that resolves r_ocn later knows whether to default
  // r_init to it.
  InitType init = InitType::Ratio;
  bool r_init_given = false;
  Real r_init = 1;

  // Generic tracers only
  OceanFlux ocean_flux = OceanFlux::BulkProportional;

  // Property tracers only: P = scale*source + offset, carried by `carrier`.
  // field_name/field_units are set only when source == "field".
  std::string carrier, source, field_name, field_units;
  Real scale = 1;
  Real offset = 0;

  // The tracer's raw entry from the file. Extension point for keys added by
  // later changes; nothing reads it yet.
  ekat::ParameterList extra;

  bool is_mass_tracer () const { return type != TracerType::Property; }
};

// Parse and validate a tracer file. Errors (std::runtime_error, via
// EKAT_REQUIRE_MSG) name the file, the tracer and the offending key.
std::vector<TracerSpec> parse_tracers_file (const std::string& path,
                                            const RegistryOptions& opts);

// Units of a `field_units` token: 1, K, kg/kg, m/s, Pa, W/m2, kg/m2/s, m.
// Errors on any other token.
ekat::units::Units property_field_units (const std::string& token);

// Range of a bounded built-in property source (lat [-90,90] degN, lon [0,360]
// degE, x/y/z [-1,1]). Returns false if the source has no fixed bounds.
bool bounded_source_range (const std::string& source, Real& lo, Real& hi);

} // namespace water_tracers
} // namespace scream

#endif // EAMXX_WATER_TRACERS_REGISTRY_HPP
