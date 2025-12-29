/**
 * @file atm_coupling.hpp
 * @brief Atmosphere coupling field indices and transfer functions.
 *
 * Defines the mapping between MCT attribute vector field names and
 * internal storage, as well as functions to transfer data between
 * the MCT buffers and the atmosphere field manager.
 */

#ifndef ATM_COUPLING_HPP
#define ATM_COUPLING_HPP

#include "../../../common/src/coupling_fields.hpp"
#include "atm_field_manager.hpp"

namespace emulator {
namespace impl {

/**
 * @brief Coupling field indices for atmosphere component.
 *
 * Stores the index positions of all export (a2x) and import (x2a)
 * fields within the MCT attribute vector buffers. An index of -1
 * indicates the field is not present in the coupling.
 *
 * ## Export Fields (a2x - atmosphere to coupler)
 * 
 * - Sa_*   : State variables at lowest atmosphere level
 * - Faxa_* : Atmosphere fluxes
 *
 * ## Import Fields (x2a - coupler to atmosphere)
 * 
 * - Sx_*   : Surface state variables (merged)
 * - So_*   : Ocean state variables
 * - Faxx_* : Surface fluxes from other components
 * - Sf_*   : Surface fractions
 */
struct AtmCouplingIndices {
  // =========================================================================
  // Export indices (a2x)
  // =========================================================================
  int Sa_z = -1;       ///< Atmospheric height at bottom level [m]
  int Sa_u = -1;       ///< Zonal wind at bottom level [m/s]
  int Sa_v = -1;       ///< Meridional wind at bottom level [m/s]
  int Sa_tbot = -1;    ///< Temperature at bottom level [K]
  int Sa_ptem = -1;    ///< Potential temperature [K]
  int Sa_shum = -1;    ///< Specific humidity [kg/kg]
  int Sa_dens = -1;    ///< Air density [kg/m³]
  int Sa_pbot = -1;    ///< Pressure at bottom level [Pa]
  int Sa_pslv = -1;    ///< Sea level pressure [Pa]
  int Faxa_lwdn = -1;  ///< Downward longwave radiation [W/m²]
  int Faxa_rainc = -1; ///< Convective precipitation [kg/m²/s]
  int Faxa_rainl = -1; ///< Large-scale precipitation [kg/m²/s]
  int Faxa_snowc = -1; ///< Convective snowfall [kg/m²/s]
  int Faxa_snowl = -1; ///< Large-scale snowfall [kg/m²/s]
  int Faxa_swndr = -1; ///< NIR direct shortwave [W/m²]
  int Faxa_swvdr = -1; ///< Visible direct shortwave [W/m²]
  int Faxa_swndf = -1; ///< NIR diffuse shortwave [W/m²]
  int Faxa_swvdf = -1; ///< Visible diffuse shortwave [W/m²]
  int Faxa_swnet = -1; ///< Net shortwave radiation [W/m²]

  // =========================================================================
  // Import indices (x2a)
  // =========================================================================
  int Sx_t = -1;      ///< Merged surface temperature [K]
  int So_t = -1;      ///< Ocean surface temperature [K]
  int Faxx_sen = -1;  ///< Sensible heat flux [W/m²]
  int Faxx_evap = -1; ///< Evaporative flux [kg/m²/s]
  int Faxx_lat = -1;  ///< Latent heat flux [W/m²]
  int Faxx_taux = -1; ///< Zonal wind stress [N/m²]
  int Faxx_tauy = -1; ///< Meridional wind stress [N/m²]
  int Faxx_lwup = -1; ///< Upward longwave radiation [W/m²]
  int Sx_avsdr = -1;  ///< Visible direct albedo [-]
  int Sx_anidr = -1;  ///< NIR direct albedo [-]
  int Sx_avsdf = -1;  ///< Visible diffuse albedo [-]
  int Sx_anidf = -1;  ///< NIR diffuse albedo [-]
  int Sl_snowh = -1;  ///< Snow height over land [m]
  int Si_snowh = -1;  ///< Snow height over ice [m]
  int Sx_tref = -1;   ///< Reference temperature [K]
  int Sx_qref = -1;   ///< Reference specific humidity [kg/kg]
  int Sx_u10 = -1;    ///< 10m wind speed [m/s]
  int Sf_ifrac = -1;  ///< Sea ice fraction [-]
  int Sf_ofrac = -1;  ///< Ocean fraction [-]
  int Sf_lfrac = -1;  ///< Land fraction [-]

  /**
   * @brief Initialize indices from coupling field base.
   *
   * Looks up each field name in the CouplingFieldsBase maps
   * and stores the corresponding index.
   *
   * @param fields Initialized CouplingFieldsBase with parsed field lists
   */
  void initialize(CouplingFieldsBase &fields);
};

/**
 * @brief Import fields from MCT buffer to local field manager.
 *
 * Copies data from the MCT import buffer (x2a fields) to the
 * corresponding vectors in the AtmFieldManager.
 *
 * @param import_data Pointer to MCT import buffer (column-major layout)
 * @param ncols Number of local columns
 * @param nfields Number of import fields
 * @param idx Coupling indices
 * @param fields Field manager to populate
 *
 * @note MCT uses column-major (Fortran) layout: data[col * nfields + field_idx]
 */
void import_atm_fields(const double *import_data, int ncols, int nfields,
                       const AtmCouplingIndices &idx, AtmFieldManager &fields);

/**
 * @brief Export fields from local field manager to MCT buffer.
 *
 * Copies data from the AtmFieldManager vectors to the MCT export
 * buffer (a2x fields).
 *
 * @param export_data Pointer to MCT export buffer (column-major layout)
 * @param ncols Number of local columns
 * @param nfields Number of export fields
 * @param idx Coupling indices
 * @param fields Field manager with data to export
 */
void export_atm_fields(double *export_data, int ncols, int nfields,
                       const AtmCouplingIndices &idx,
                       const AtmFieldManager &fields);

} // namespace impl
} // namespace emulator

#endif // ATM_COUPLING_HPP
