#include "eamxx_water_isotopes_process_interface.hpp"

namespace scream {
namespace {
  // generic parse_option utility function.
  template <typename E>
  E parse_option(ekat::ParameterList& p, 
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
  }
}
// =========================================================================================
WaterIsotopes::WaterIsotopes(const ekat::Comm& comm, const ekat::ParameterList& params)
  : WaterTracers(comm, params)
{
  // Water isotopes will inherit all tracer handling from WaterTracers
  // Read runtime formulation choices from parameters with sensible defaults

  // Liquid/vapor fractionation formulation
  m_runtime_options.liquid_vapor = parse_option<wiso::LiquidVaporFractionation>(
      m_params, "liquid_vapor_formulation", "horita_wesolowski_1994", {
        {"horita_wesolowski_1994",
  wiso::LiquidVaporFractionation::HoritaWesolowski1994},
        {"majoube_1971", wiso::LiquidVaporFractionation::Majoube1971}
      });

  // Standard ratio formulation
  m_runtime_options.standard_ratio = parse_option<wiso::StandardRatioFormulation>(
    m_params,"standard_ratio_formulation", "normalized", {
      {"normalized",wiso::StandardRatioFormulation::Normalized},
      {"natural_abundance",wiso::StandardRatioFormulation::NaturalAbundance}
    });

  // Ocean enrichment formulation
  m_runtime_options.ocean_enrichment = parse_option<wiso::OceanEnrichmentFormulation>(
    m_params,"ocean_enrichment_formulation","modern", {
      {"modern",wiso::OceanEnrichmentFormulation::Modern},
      {"LGM",wiso::OceanEnrichmentFormulation::LGM}
    });

  // Ice/vapor fractionation formulation
  m_runtime_options.ice_vapor = parse_option<wiso::IceVaporFractionation>(
    m_params,"ice_vapor_formulation", "merlivat_nief_1967", {
      {"merlivat_nief_1967",wiso::IceVaporFractionation::MerlivatNief1967},
      {"isocam3",wiso::IceVaporFractionation::IsoCAM3}
    });

}

// =========================================================================================
void WaterIsotopes::run_impl(const double dt)
{
  // Call base class tracer physics (currently a no-op)
  WaterTracers::run_impl(dt);

  // TODO: Add fractionation physics here
}

} // namespace scream
