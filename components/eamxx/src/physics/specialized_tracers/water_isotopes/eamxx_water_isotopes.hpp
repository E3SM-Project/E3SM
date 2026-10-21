#ifndef EAMXX_WATER_ISOTOPES_HPP
#define EAMXX_WATER_ISOTOPES_HPP

#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_registry.hpp"
#include "eamxx_water_isotopes_parameters.hpp"

#include <ekat_parameter_list.hpp>

#include <string>
#include <vector>

namespace scream
{
namespace wiso
{

// One resolved isotope tracer: its tracers_file name, species, and ocean
// ratio R_ocn (resolved from the parameter tables by WaterIsotopes::resolve).
struct IsotopeTracer {
  std::string        name;
  WaterIsotopologues species;
  Real               r_ocn;
};

/*
 * Isotope-specific state and configuration composed into the water_tracers
 * process: which species names are valid in the tracers file, the runtime
 * formulation options (liquid/vapor and ice/vapor fractionation, standard
 * ratio, ocean enrichment), and the isotope tracers resolved from a parsed
 * tracers_file.
 *
 * This is a plain helper class, not an AtmosphereProcess: it owns no fields
 * and makes no field-manager calls. WaterTracers owns one unconditionally
 * (empty/inert unless EAMXX_HAS_WATER_ISOTOPES is defined and the tracers
 * file actually defines isotope entries) and does its own field lookups and
 * kernel launches around the state and math this class provides.
 */
class WaterIsotopes
{
public:
  // Allows `type: isotope` entries in the tracers file and lists every
  // species name parse_tracers_file should accept. Called only when
  // EAMXX_HAS_WATER_ISOTOPES is defined; otherwise opts.accept_isotopes stays
  // false and the registry rejects any isotope entry outright.
  static void enable (water_tracers::RegistryOptions& opts);

  // Reads the runtime formulation options from the process's parameters
  // (liquid_vapor_formulation, ice_vapor_formulation, standard_ratio_formulation,
  // ocean_enrichment_formulation; see eamxx_water_tracers_registry.hpp for
  // defaults). Safe to call unconditionally: with no isotope tracers these
  // options are simply never read.
  void configure (ekat::ParameterList& params);

  // Resolves r_ocn (and the default r_init) for every Isotope spec from the
  // parameter tables, and records {name, species, r_ocn} for each. No-op if
  // specs contains no isotope entries (which the registry guarantees unless
  // enable() was called first).
  void resolve (std::vector<water_tracers::TracerSpec>& specs);

  const std::vector<IsotopeTracer>& tracers () const { return m_isotopes; }

  const WaterIsotopeParameters<Real>& parameters () const { return m_iso_params; }

private:
  WaterIsotopeRuntimeOptions   m_runtime_options;
  WaterIsotopeParameters<Real> m_iso_params;
  std::vector<IsotopeTracer>   m_isotopes;

}; // class WaterIsotopes

} // namespace wiso
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_HPP
