/**
 * @file atm_field_manager.hpp
 * @brief Field storage container for atmosphere emulator.
 *
 * Manages allocation and access to all atmosphere fields, including
 * import fields (from coupler), export fields (to coupler), and
 * AI model input/output tensors.
 */

#ifndef ATM_FIELD_MANAGER_HPP
#define ATM_FIELD_MANAGER_HPP

#include <map>
#include <string>
#include <vector>

namespace emulator {
namespace impl {

/**
 * @brief Field storage container for atmosphere emulator.
 *
 * Manages all field vectors for the atmosphere component, including:
 * 
 * - Import fields (x2a) received from other components via coupler
 * - Export fields (a2x) sent to other components via coupler
 * - AI model tensor buffers (net_inputs, net_outputs)
 * - Dynamic fields for configurable I/O variables
 *
 * ## Allocation
 * Call `allocate(ncols)` to allocate all field vectors for the given
 * number of local columns. The `net_inputs` and `net_outputs` vectors
 * are allocated dynamically by `EmulatorAtm::prepare_inputs()`.
 *
 * ## Field Access
 * Fields can be accessed directly or via `get_field_ptr()` for
 * string-based lookup. Dynamic fields can be registered at runtime
 * via `register_dynamic_field()`.
 */
class AtmFieldManager {
public:
  // =========================================================================
  // Constants
  // =========================================================================

  /** @brief Default number of input channels (legacy, use config instead). */
  static constexpr int N_INPUT_CHANNELS = 39;

  /** @brief Default number of output channels (legacy, use config instead). */
  static constexpr int N_OUTPUT_CHANNELS = 44;

  AtmFieldManager() = default;
  ~AtmFieldManager() = default;

  // =========================================================================
  // Allocation
  // =========================================================================

  /**
   * @brief Allocate all field vectors for given number of columns.
   *
   * Allocates import and export field vectors. The `net_inputs` and
   * `net_outputs` vectors are cleared but not pre-allocated (they are
   * sized dynamically in `prepare_inputs()`).
   *
   * @param ncols Number of local columns
   */
  void allocate(int ncols);

  /**
   * @brief Deallocate all field vectors.
   */
  void deallocate();

  /**
   * @brief Set default climatological values for export fields.
   *
   * Initializes export fields with reasonable default values for testing.
   *
   * @param ncols Number of columns
   */
  void set_defaults(int ncols);

  /**
   * @brief Check if fields are allocated.
   * @return true if allocate() has been called
   */
  bool is_allocated() const { return m_allocated; }

  // =========================================================================
  // Generic Field Access
  // =========================================================================

  /**
   * @brief Get pointer to a field vector by name.
   *
   * Looks up the field in the hardcoded field map, then in dynamic fields.
   *
   * @param name Field name (e.g., "ts", "tbot", "pbot")
   * @return Pointer to the field vector, or nullptr if not found
   */
  std::vector<double> *get_field_ptr(const std::string &name);

  /**
   * @brief Register or create a dynamic field.
   *
   * If the field doesn't exist (not in hardcoded or dynamic fields),
   * creates a new entry in dynamic_fields. If already allocated,
   * resizes the new field to match ncols.
   *
   * @param name Field name to register
   */
  void register_dynamic_field(const std::string &name);

  // =========================================================================
  // AI Model Tensors
  // =========================================================================

  std::vector<double> net_inputs;  ///< Packed input tensor for inference
  std::vector<double> net_outputs; ///< Packed output tensor from inference

  // =========================================================================
  // Imported Fields (x2a - from coupler)
  // =========================================================================

  std::vector<double> shf;          ///< Sensible heat flux [W/m²]
  std::vector<double> cflx;         ///< CO2 flux [kg/m²/s]
  std::vector<double> lhf;          ///< Latent heat flux [W/m²]
  std::vector<double> wsx;          ///< Zonal wind stress [N/m²]
  std::vector<double> wsy;          ///< Meridional wind stress [N/m²]
  std::vector<double> lwup;         ///< Upward longwave radiation [W/m²]
  std::vector<double> asdir;        ///< Visible direct albedo [-]
  std::vector<double> aldir;        ///< NIR direct albedo [-]
  std::vector<double> asdif;        ///< Visible diffuse albedo [-]
  std::vector<double> aldif;        ///< NIR diffuse albedo [-]
  std::vector<double> ts;           ///< Surface temperature [K]
  std::vector<double> sst;          ///< Sea surface temperature [K]
  std::vector<double> snowhland;    ///< Snow height over land [m]
  std::vector<double> snowhice;     ///< Snow height over ice [m]
  std::vector<double> tref;         ///< Reference temperature [K]
  std::vector<double> qref;         ///< Reference specific humidity [kg/kg]
  std::vector<double> u10;          ///< 10m wind speed [m/s]
  std::vector<double> u10withgusts; ///< 10m wind with gusts [m/s]
  std::vector<double> icefrac;      ///< Sea ice fraction [-]
  std::vector<double> ocnfrac;      ///< Ocean fraction [-]
  std::vector<double> lndfrac;      ///< Land fraction [-]

  // =========================================================================
  // Exported Fields (a2x - to coupler)
  // =========================================================================

  std::vector<double> zbot;  ///< Height at bottom level [m]
  std::vector<double> ubot;  ///< Zonal wind at bottom level [m/s]
  std::vector<double> vbot;  ///< Meridional wind at bottom level [m/s]
  std::vector<double> tbot;  ///< Temperature at bottom level [K]
  std::vector<double> ptem;  ///< Potential temperature [K]
  std::vector<double> shum;  ///< Specific humidity [kg/kg]
  std::vector<double> dens;  ///< Air density [kg/m³]
  std::vector<double> pbot;  ///< Pressure at bottom level [Pa]
  std::vector<double> pslv;  ///< Sea level pressure [Pa]
  std::vector<double> lwdn;  ///< Downward longwave radiation [W/m²]
  std::vector<double> rainc; ///< Convective precipitation [kg/m²/s]
  std::vector<double> rainl; ///< Large-scale precipitation [kg/m²/s]
  std::vector<double> snowc; ///< Convective snowfall [kg/m²/s]
  std::vector<double> snowl; ///< Large-scale snowfall [kg/m²/s]
  std::vector<double> swndr; ///< NIR direct shortwave [W/m²]
  std::vector<double> swvdr; ///< Visible direct shortwave [W/m²]
  std::vector<double> swndf; ///< NIR diffuse shortwave [W/m²]
  std::vector<double> swvdf; ///< Visible diffuse shortwave [W/m²]
  std::vector<double> swnet; ///< Net shortwave radiation [W/m²]

  // =========================================================================
  // Dynamic Fields
  // =========================================================================

  /** @brief Runtime-registered fields for configurable I/O. */
  std::map<std::string, std::vector<double>> dynamic_fields;

private:
  bool m_allocated = false; ///< Allocation flag
  int m_ncols = 0;          ///< Number of columns

  void init_field_map();
  std::map<std::string, std::vector<double> *> m_field_map;
};

} // namespace impl
} // namespace emulator

#endif // ATM_FIELD_MANAGER_HPP
