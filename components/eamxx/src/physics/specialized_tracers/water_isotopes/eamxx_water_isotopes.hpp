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

/*
 * Isotope-specific state and configuration composed into the water_tracers
 * process: which species names are valid in the tracers file, and the runtime
 * formulation options (liquid/vapor and ice/vapor fractionation, standard
 * ratio, ocean enrichment). Per-tracer data (name, species, R_ocn) lives only
 * in the water_tracers::TracerSpec entries, which resolve() completes.
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

  // Completes every Isotope spec: r_ocn (and the default r_init) from the
  // parameter tables, and ocean_flux = BulkFormula (the only form that takes
  // the fractionation factor; see water_tracers::ocean_flux).
  void resolve (std::vector<water_tracers::TracerSpec>& specs) const;

  // The species enum of an Isotope spec, from its `species:` name. Host-only.
  static WaterIsotopologues species (const water_tracers::TracerSpec& spec);

  const WaterIsotopeParameters<Real>& parameters () const { return m_iso_params; }

private:
  WaterIsotopeRuntimeOptions   m_runtime_options;
  WaterIsotopeParameters<Real> m_iso_params;

}; // class WaterIsotopes

} // namespace water_isotopes
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_HPP
