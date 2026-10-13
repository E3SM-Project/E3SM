#include "eamxx_water_isotopes_process_interface.hpp"

#include <ekat_assert.hpp>
#include <ekat_pack.hpp>

#include <iterator>
#include <utility>

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
    return choices.begin()->second;  // unreachable; silences -Wreturn-type warning
                                   
  }

  // Species names accepted in `species:` of a tracers_file isotope entry.
  // Host-only: device code receives the resolved enum, never this table.
  constexpr std::pair<const char*, wiso::WaterIsotopologues> species_names[] = {
    {"H216O", wiso::WaterIsotopologues::H216O},
    {"HDO",   wiso::WaterIsotopologues::HDO},
    {"H218O", wiso::WaterIsotopologues::H218O},
    {"H217O", wiso::WaterIsotopologues::H217O},
    {"HTO",   wiso::WaterIsotopologues::HTO},
  };
  static_assert(std::size(species_names) == etoi(wiso::WaterIsotopologues::Count),
                "Every WaterIsotopologues species needs a tracers_file name.");

  water_tracers::RegistryOptions isotope_registry_options()
  {
    water_tracers::RegistryOptions opts;
    opts.accept_isotopes = true;
    for (const auto& [str, _] : species_names) {
      opts.isotope_species.push_back(str);
    }
    return opts;
  }

  wiso::WaterIsotopologues species_from_name(const std::string& name)
  {
    for (const auto& [str, species] : species_names) {
      if (name == str) return species;
    }
    // The registry already rejects names missing from isotope_species.
    EKAT_ERROR_MSG("Error! Unknown water isotope species '" + name + "'.\n");
    return wiso::WaterIsotopologues::H216O;  // unreachable; see parse_option
  }
}
// =========================================================================================
WaterIsotopes::WaterIsotopes(const ekat::Comm& comm, const ekat::ParameterList& params)
  : WaterTracers(comm, params, isotope_registry_options())
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

  m_iso_params = wiso::WaterIsotopeParameters<Real>(m_runtime_options);

  // Resolve each isotope's ocean ratio here, on the host, so kernels get a
  // plain Real. The base class applies `init: ratio` later, in
  // initialize_impl, so setting r_init here is still in time.
  for (auto& spec : m_specs) {
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

// =========================================================================================
void WaterIsotopes::compute_isotope_fluxes()
{
  using Pack1 = ekat::Pack<Real, 1>;

  // Local copies, so the lambda captures views and values rather than `this`.
  const auto sfc_f     = m_sfc_f;
  const auto sfc_rho   = m_sfc_rho;
  const auto sfc_ustar = m_sfc_ustar;
  const auto sfc_re    = m_sfc_re;
  const auto sfc_ssq   = m_sfc_ssq;
  const auto sfc_qv    = m_sfc_qv;
  const auto sfc_evap  = m_sfc_evap;
  const auto sfc_sst   = m_sfc_sst;
  const auto iso_params = m_iso_params;
  const int  ilev_bot  = m_num_levs - 1;

  // Tracers may also live on another grid (e.g. physics_gll for initial
  // conditions), so look them up on this process' physics grid explicitly.
  const auto& grid_name = m_grid->name();

  // One kernel per isotope: there are at most a handful, and each reads a
  // different tracer field.
  for (const auto& iso : m_isotopes) {
    const auto q_i  = get_field_in(iso.name, grid_name).get_view<const Real**>();
    const auto flux = get_field_out(iso.name + "_surf_flux", grid_name).get_view<Real*>();
    const auto species = iso.species;
    const Real r_ocn   = iso.r_ocn;

    Kokkos::parallel_for("water_isotopes_ocean_flux_" + iso.name,
                         KT::RangePolicy(0, m_num_cols),
                         KOKKOS_LAMBDA(const int icol) {
      flux(icol) = wiso::isotope_ocean_flux(
          Pack1(sfc_f(icol)), Pack1(sfc_rho(icol)), Pack1(sfc_ustar(icol)),
          Pack1(sfc_re(icol)), Pack1(sfc_ssq(icol)), Pack1(sfc_qv(icol)),
          Pack1(q_i(icol, ilev_bot)), Pack1(sfc_evap(icol)), Pack1(sfc_sst(icol)),
          r_ocn, species, iso_params)[0];
    });
  }
}

} // namespace scream
