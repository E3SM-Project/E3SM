#ifndef EAMXX_WATER_ISOTOPES_HPP
#define EAMXX_WATER_ISOTOPES_HPP

#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_registry.hpp"
#include "eamxx_water_isotopes_parameters.hpp"

#include <ekat_parameter_list.hpp>

#include <string>
#include <vector>

namespace scream
{
namespace water_isotopes
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
 */
class WaterIsotopes
{
public:
  // Allows `type: isotope` entries in the tracers file and lists every
  // species name parse_tracers_file should accept. Called only when
  // EAMXX_HAS_WATER_ISOTOPES is true.
  static void enable (water_tracers::RegistryOptions& opts);

  // Reads the runtime formulation options from the process's parameters.
  // with no isotope tracers these options are just never read.
  void configure (ekat::ParameterList& params);

  // Resolves r_ocn (and the default r_init) for every Isotope spec from the
  // parameter tables, and records {name, species, r_ocn} for each.
  void resolve (std::vector<water_tracers::TracerSpec>& specs);

  const std::vector<IsotopeTracer>& tracers () const { return m_isotopes; }

  const WaterIsotopeParameters<Real>& parameters () const { return m_iso_params; }

private:
  WaterIsotopeRuntimeOptions   m_runtime_options;
  WaterIsotopeParameters<Real> m_iso_params;
  std::vector<IsotopeTracer>   m_isotopes;

}; // class WaterIsotopes

} // namespace water_isotopes
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_HPP
