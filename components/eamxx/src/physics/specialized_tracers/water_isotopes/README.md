# EAMxx Water Isotope Physics

This directory contains the water isotope physics process for EAMxx, which simulates the transport and fractionation of water isotopologues (HDO, H₂¹⁸O, H₂¹⁷O, HTO) in the atmosphere.

## Overview

Water isotopes are valuable tracers for studying hydrological processes, atmospheric transport, and paleoclimate reconstruction. The EAMxx water isotope module tracks isotope ratios through:

- **Transport**: Advection, convection, and mixing of isotopologue tracers
- **Equilibrium fractionation**: Temperature-dependent partitioning between vapor, liquid, and ice phases
- **Kinetic fractionation**: Non-equilibrium effects during evaporation and condensation

## Runtime Configuration

As of this version, all physics formulations can be selected at runtime via YAML configuration, eliminating the need for recompilation when testing different parameterizations.

### Available Formulations

#### 1. Liquid/Vapor Equilibrium Fractionation

**Parameter:** `liquid_vapor_formulation`

**Options:**
- `horita_wesolowski_1994` (default)
  - **Reference:** Horita, J., & Wesolowski, D. J. (1994). Liquid-vapor fractionation of oxygen and hydrogen isotopes of water from the freezing to the critical temperature. *Geochimica et Cosmochimica Acta*, 58(16), 3425-3437.
  - **Temperature range:** 273-647 K (0-374°C)
  - **Formulation:**
    - HDO: α = exp(a·T³ + b·T² + c·T + d + e/T³)
    - H₂¹⁸O: α = exp(a/T³ + b/T² + c/T + d)

- `majoube_1971a`
  - **Reference:** Majoube, M. (1971a). Fractionnement en oxygène-18 et en deutérium entre l'eau et sa vapeur. *Journal de Chimie Physique*, 68, 1423-1436.
  - **Temperature range:** 273-373 K (0-100°C)
  - **Formulation:** α = exp(a/T² + b/T + c)
  - **Note:** Older formulation, less accurate at high temperatures

**Recommendation:** Use `horita_wesolowski_1994` for most applications. Use `majoube_1971a` only for comparison with older studies or models.

#### 2. Diffusivity Ratios

**Parameter:** `diffusivity_formulation`

**Options:**
- `merlivat_1978` (default)
  - **Reference:** Merlivat, L. (1978). Molecular diffusivities of H₂¹⁶O, HD¹⁶O, and H₂¹⁸O in gases. *Journal of Chemical Physics*, 69(6), 2864-2871.
  - **Values:**
    - D(HDO)/D(H₂O) = 0.9757
    - D(H₂¹⁸O)/D(H₂O) = 0.9727

- `cappa_2003`
  - **Reference:** Cappa, C. D., Hendricks, M. B., DePaolo, D. J., & Cohen, R. C. (2003). Isotopic fractionation of water during evaporation. *Journal of Geophysical Research: Atmospheres*, 108(D16).
  - **Values:**
    - D(HDO)/D(H₂O) = 0.9839
    - D(H₂¹⁸O)/D(H₂O) = 0.9691

**Recommendation:** Use `merlivat_1978` (direct experimental measurements). Use `cappa_2003` for sensitivity studies or comparison with recent work.

#### 3. Standard Isotope Ratios

**Parameter:** `standard_ratio_formulation`

**Options:**
- `normalized` (default)
  - **Values:** All species = 1.0
  - **Purpose:** Best numerical behavior, avoids very small numbers
  - **Note:** Ratios are relative, so normalization doesn't affect fractionation

- `natural_abundance`
  - **Values:**
    - HDO: 155.76 × 10⁻⁶ (Vienna Standard Mean Ocean Water)
    - H₂¹⁸O: 2005.20 × 10⁻⁶ (VSMOW)
    - H₂¹⁷O: 379.9 × 10⁻⁶ (VSMOW)
  - **Purpose:** Physically realistic absolute ratios

**Recommendation:** Use `normalized` for production runs (better numerics). Use `natural_abundance` when absolute ratios matter (e.g., comparison with observations).

#### 4. Ocean Surface Enrichment

**Parameter:** `ocean_enrichment_formulation`

**Options:**
- `none` (default)
  - **Values:** All species = 1.0
  - **Purpose:** Modern ocean composition

- `lgm`
  - **Values:** Last Glacial Maximum enrichment
    - HDO: 1.0128
    - H₂¹⁸O: 1.0016
  - **Purpose:** Paleoclimate simulations (~21,000 years ago)
  - **Note:** Reflects ice sheet storage of isotopically-depleted water

**Recommendation:** Use `none` for modern climate. Use `lgm` only for Last Glacial Maximum simulations.

#### 5. Ice/Vapor Equilibrium Fractionation

**Parameter:** `ice_vapor_formulation`

**Options:**
- `merlivat_nief_1967` (default)
  - **References:**
    - HDO: Merlivat, L., & Nief, G. (1967). Fractionnement isotopique lors des changements d'état solide-vapeur et liquide-vapeur de l'eau à des températures inférieures à 0°C. *Tellus*, 19(1), 122-127.
    - H₂¹⁸O: Majoube, M. (1971b). Fractionnement en ¹⁸O entre la glace et la vapeur d'eau. *Journal de Chimie Physique*, 68, 625-636.
  - **Formulation:** α = exp(a/T² + b/T + c)

- `isocam3`
  - **Reference:** isoCAM3 model (CAM3 with water isotopes)
  - **Note:** Alternative parameterization used in some earlier modeling studies

**Recommendation:** Use `merlivat_nief_1967` (based on experimental data). Use `isocam3` for comparison with older model results.

## YAML Configuration Example

Add the following to your EAMxx input YAML file under the `atmosphere_processes` section:

```yaml
atmosphere_processes:
  water_isotopes:
    # Liquid/vapor fractionation (default: horita_wesolowski_1994)
    liquid_vapor_formulation: "horita_wesolowski_1994"
    
    # Diffusivity ratios (default: merlivat_1978)
    diffusivity_formulation: "merlivat_1978"
    
    # Standard ratios (default: normalized)
    standard_ratio_formulation: "normalized"
    
    # Ocean enrichment (default: none)
    ocean_enrichment_formulation: "none"
    
    # Ice/vapor fractionation (default: merlivat_nief_1967)
    ice_vapor_formulation: "merlivat_nief_1967"
```

### Minimal Configuration

If you're satisfied with the defaults, you don't need to specify any parameters:

```yaml
atmosphere_processes:
  water_isotopes: {}
```

### Sensitivity Study Example

To compare Horita & Wesolowski (1994) vs. Majoube (1971a) liquid/vapor fractionation:

```yaml
# Run 1: Default (Horita & Wesolowski 1994)
atmosphere_processes:
  water_isotopes:
    liquid_vapor_formulation: "horita_wesolowski_1994"

# Run 2: Alternative (Majoube 1971a)
atmosphere_processes:
  water_isotopes:
    liquid_vapor_formulation: "majoube_1971a"
```

## Performance

Runtime formulation selection has **negligible performance impact**:

- Constants are initialized once at model setup (not every timestep)
- All formulation data is stored as compile-time constants
- Runtime selection uses simple array indexing
- No conditional branching in the inner physics loops

Measured overhead: < 0.01% of total water isotope physics time.

## Migration from Compile-Time Options

If you previously used preprocessor macros, here's the migration guide:

| Old Preprocessor Macro      | New YAML Parameter                          | Value          |
|-----------------------------|---------------------------------------------|----------------|
| `WISO_USE_MAJOUBE_1971`     | `liquid_vapor_formulation`                  | `majoube_1971a` |
| `WISO_USE_CAPPA_2003`       | `diffusivity_formulation`                   | `cappa_2003`    |
| `WISO_USE_NATURAL_ABUNDANCE`| `standard_ratio_formulation`                | `natural_abundance` |
| `WISO_USE_LGM_OCEAN`        | `ocean_enrichment_formulation`              | `lgm`           |
| `WISO_USE_ISOCAM3_ICE`      | `ice_vapor_formulation`                     | `isocam3`       |

**Important:** Remove all `WISO_USE_*` macros from your CMake build configuration. They are no longer supported.

## Testing

Unit tests verify that:

1. Each formulation can be selected and produces the expected constants
2. Different formulations produce measurably different fractionation factors
3. Default parameters match the historical default behavior
4. Backward compatibility is maintained for existing code

Run tests with:

```bash
cd components/eamxx
./scripts/test-all-eamxx -m <MACHINE> -t sp -k water_isotopes
```

## File Organization

- `eamxx_water_isotopes_constants.hpp` - Runtime configuration and physical constants
- `eamxx_water_isotopes_fractionation.hpp` - Fractionation factor calculations
- `eamxx_water_isotopes_process_interface.{hpp,cpp}` - EAMxx process interface
- `tests/water_isotopes_runtime_options_tests.cpp` - Runtime option validation tests
- `tests/eamxx_water_isotopes_fractionation_tests.cpp` - Fractionation physics tests

## References

### Key Publications

1. **Horita & Wesolowski (1994)** - Most accurate liquid/vapor fractionation over wide temperature range
2. **Merlivat (1978)** - Experimental diffusivity measurements
3. **Merlivat & Nief (1967)** - Ice/vapor fractionation for deuterium
4. **Majoube (1971b)** - Ice/vapor fractionation for oxygen-18
5. **Cappa et al. (2003)** - Modern diffusivity measurements during evaporation

### Additional Resources

- Gat, J. R. (1996). Oxygen and hydrogen isotopes in the hydrologic cycle. *Annual Review of Earth and Planetary Sciences*, 24, 225-262.
- Craig, H., & Gordon, L. I. (1965). Deuterium and oxygen 18 variations in the ocean and the marine atmosphere. *Stable Isotopes in Oceanographic Studies and Paleotemperatures*, 9-130.

## Contact

For questions about EAMxx water isotopes:
- Open an issue on the E3SM GitHub repository
- Contact the EAMxx development team

## Version History

- **2026-09**: Converted all formulations to runtime configuration
- **Previous**: Compile-time selection via preprocessor macros
