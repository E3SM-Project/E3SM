#include "eamxx_water_isotopes_process_interface.hpp"

namespace scream
{

// =========================================================================================
WaterIsotopes::WaterIsotopes(const ekat::Comm& comm, const ekat::ParameterList& params)
  : WaterTracers(comm, params)
{
  // Water isotopes will inherit all tracer handling from WaterTracers

  // Read runtime formulation choices from parameters with sensible defaults

  // Liquid/vapor fractionation formulation
  const std::string lv_form = m_params.get<std::string>(
    "liquid_vapor_formulation", "horita_wesolowski_1994");
  if (lv_form == "majoube_1971a") {
    m_runtime_options.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971a;
  } else if (lv_form == "horita_wesolowski_1994") {
    m_runtime_options.liquid_vapor = wiso::LiquidVaporFractionation::HoritaWesolowski1994;
  } else {
    EKAT_ERROR_MSG("Invalid liquid_vapor_formulation: " + lv_form +
                   ". Valid options: horita_wesolowski_1994, majoube_1971a");
  }

  // Diffusivity formulation
  const std::string diff_form = m_params.get<std::string>(
    "diffusivity_formulation", "merlivat_1978");
  if (diff_form == "cappa_2003") {
    m_runtime_options.diffusivity = wiso::DiffusivityFormulation::Cappa2003;
  } else if (diff_form == "merlivat_1978") {
    m_runtime_options.diffusivity = wiso::DiffusivityFormulation::Merlivat1978;
  } else {
    EKAT_ERROR_MSG("Invalid diffusivity_formulation: " + diff_form +
                   ". Valid options: merlivat_1978, cappa_2003");
  }

  // Standard ratio formulation
  const std::string ratio_form = m_params.get<std::string>(
    "standard_ratio_formulation", "normalized");
  if (ratio_form == "natural_abundance") {
    m_runtime_options.standard_ratio = wiso::StandardRatioFormulation::NaturalAbundance;
  } else if (ratio_form == "normalized") {
    m_runtime_options.standard_ratio = wiso::StandardRatioFormulation::Normalized;
  } else {
    EKAT_ERROR_MSG("Invalid standard_ratio_formulation: " + ratio_form +
                   ". Valid options: normalized, natural_abundance");
  }

  // Ocean enrichment formulation
  const std::string ocean_form = m_params.get<std::string>(
    "ocean_enrichment_formulation", "none");
  if (ocean_form == "lgm") {
    m_runtime_options.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
  } else if (ocean_form == "none") {
    m_runtime_options.ocean_enrichment = wiso::OceanEnrichmentFormulation::None;
  } else {
    EKAT_ERROR_MSG("Invalid ocean_enrichment_formulation: " + ocean_form +
                   ". Valid options: none, lgm");
  }

  // Ice/vapor fractionation formulation
  const std::string iv_form = m_params.get<std::string>(
    "ice_vapor_formulation", "merlivat_nief_1967");
  if (iv_form == "isocam3") {
    m_runtime_options.ice_vapor = wiso::IceVaporFractionation::IsoCAM3;
  } else if (iv_form == "merlivat_nief_1967") {
    m_runtime_options.ice_vapor = wiso::IceVaporFractionation::MerlivatNief1967;
  } else {
    EKAT_ERROR_MSG("Invalid ice_vapor_formulation: " + iv_form +
                   ". Valid options: merlivat_nief_1967, isocam3");
  }
}

// =========================================================================================
void WaterIsotopes::run_impl(const double dt)
{
  // Call base class tracer physics (currently a no-op)
  WaterTracers::run_impl(dt);

  // TODO: Add fractionation physics here
}

} // namespace scream
