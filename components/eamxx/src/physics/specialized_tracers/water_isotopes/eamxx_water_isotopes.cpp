#include "eamxx_water_isotopes.hpp"

#include <ekat_assert.hpp>
#include <ekat_string_utils.hpp>

#include <iterator>
#include <utility>

namespace scream {
namespace water_isotopes {
namespace {

// generic parse_option utility function.
template <typename E>
E parse_option (ekat::ParameterList& p,
                const std::string& key,
                const std::string& defaultValue,
                std::initializer_list<std::pair<const char*, E>> choices)
{
  const std::string value = p.get<std::string>(key, defaultValue);

  for (const auto& [str, enum_val] : choices) {
    if (value == str) return enum_val;
  }

  // Build valid options list
  std::vector<std::string> valid_opts;
  for (const auto& [str, _] : choices) {
    valid_opts.push_back(str);
  }

  EKAT_ERROR_MSG("Invalid " + key + ": '" + value +
                 "'. Valid options: " + ekat::join(valid_opts, ", "));
  return choices.begin()->second;  // unreachable; silences -Wreturn-type warning
}

// Species names accepted in `species:` of a water_tracers_file isotope entry.
// Host-only: device code receives the resolved enum, never this table.
constexpr std::pair<const char*, WaterIsotopologues> species_names[] = {
  {"H216O", WaterIsotopologues::H216O},
  {"HDO",   WaterIsotopologues::HDO},
  {"H218O", WaterIsotopologues::H218O},
  {"H217O", WaterIsotopologues::H217O},
  {"HTO",   WaterIsotopologues::HTO},
};
static_assert(std::size(species_names) == etoi(WaterIsotopologues::Count),
              "Every WaterIsotopologues species needs a water_tracers_file name.");

WaterIsotopologues species_from_name (const std::string& name)
{
  for (const auto& [str, species] : species_names) {
    if (name == str) return species;
  }
  // The registry already rejects names missing from isotope_species.
  EKAT_ERROR_MSG("Error! Unknown water isotope species '" + name + "'.\n");
  return WaterIsotopologues::H216O;  // unreachable; see parse_option
}

} // anonymous namespace

// =========================================================================================
void WaterIsotopes::enable (water_tracers::RegistryOptions& opts)
{
  opts.accept_isotopes = true;
  for (const auto& [str, _] : species_names) {
    opts.isotope_species.push_back(str);
  }
}

// =========================================================================================
void WaterIsotopes::configure (ekat::ParameterList& params)
{
  // Liquid/vapor fractionation formulation
  m_runtime_options.liquid_vapor = parse_option<LiquidVaporFractionation>(
      params, "liquid_vapor_formulation", "horita_wesolowski_1994", {
        {"horita_wesolowski_1994", LiquidVaporFractionation::HoritaWesolowski1994},
        {"majoube_1971",           LiquidVaporFractionation::Majoube1971}
      });

  // Standard ratio formulation
  m_runtime_options.standard_ratio = parse_option<StandardRatioFormulation>(
      params, "standard_ratio_formulation", "normalized", {
        {"normalized",         StandardRatioFormulation::Normalized},
        {"natural_abundance",  StandardRatioFormulation::NaturalAbundance}
      });

  // Ocean enrichment formulation
  m_runtime_options.ocean_enrichment = parse_option<OceanEnrichmentFormulation>(
      params, "ocean_enrichment_formulation", "modern", {
        {"modern", OceanEnrichmentFormulation::Modern},
        {"LGM",    OceanEnrichmentFormulation::LGM}
      });

  // Ice/vapor fractionation formulation
  m_runtime_options.ice_vapor = parse_option<IceVaporFractionation>(
      params, "ice_vapor_formulation", "merlivat_nief_1967", {
        {"merlivat_nief_1967", IceVaporFractionation::MerlivatNief1967},
        {"isocam3",            IceVaporFractionation::IsoCAM3}
      });

  m_iso_params = WaterIsotopeParameters<Real>(m_runtime_options);
}

// =========================================================================================
void WaterIsotopes::resolve (std::vector<water_tracers::TracerSpec>& specs)
{
  // Resolve each isotope's ocean ratio here, on the host, so kernels get a
  // plain Real. The tracer's own `init: ratio` is applied later, in
  // WaterTracers::initialize_impl, so setting r_init here is still in time.
  for (auto& spec : specs) {
    if (spec.type != water_tracers::TracerType::Isotope) continue;

    const auto species = species_from_name(spec.species);
    spec.r_ocn = m_iso_params.std_iso_ratio(species) *
                 m_iso_params.mean_ocean_enrichment(species);
    if (!spec.r_init_given) {
      spec.r_init = spec.r_ocn;
    }
    m_isotopes.push_back({spec.name, species, spec.r_ocn});
  }
}

} // namespace water_isotopes
} // namespace scream
