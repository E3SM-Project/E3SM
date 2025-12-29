/**
 * @file emulator_atm.hpp
 * @brief Atmosphere emulator component declaration.
 *
 * Defines the EmulatorAtm class, which implements an AI-based atmosphere
 * component for E3SM. Supports multiple inference backends (STUB, LibTorch)
 * for executing neural network models like ACE2.
 */

#ifndef EMULATOR_ATM_HPP
#define EMULATOR_ATM_HPP

#include "../../common/src/coupling_fields.hpp"
#include "../../common/src/emulator_comp.hpp"
#include "../../common/src/emulator_config.hpp"
#include "../../common/src/emulator_output_manager.hpp"
#include "../../common/src/inference/inference_backend.hpp"
#include "impl/atm_coupling.hpp"
#include "impl/atm_field_data_provider.hpp"
#include "impl/atm_field_manager.hpp"
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace emulator {

/**
 * @brief Atmosphere emulator component.
 *
 * Derived from EmulatorComp, provides atmosphere-specific functionality:
 *
 * - Coupling field mappings (x2a → atm inputs, a2x → atm outputs)
 * - AI model integration via configurable inference backends
 * - Tensor reshaping for CNN vs MLP models (spatial_mode)
 *
 * ## Supported Backends
 *
 * - **STUB**: No-op for testing (returns zeros)
 * - **LIBTORCH**: Native C++ TorchScript inference
 *
 * ## Lifecycle
 *
 * 1. Constructor creates EmulatorAtm with CompType::ATM
 * 2. `create_instance()` initializes MPI and reads grid
 * 3. `init_coupling_indices()` parses MCT field lists
 * 4. `set_inference_config()` configures the backend
 * 5. `initialize()` loads model and reads initial conditions
 * 6. `run()` executes time steps (import → inference → export)
 * 7. `finalize()` cleans up resources
 *
 * ## Spatial Mode
 * When `config.model_io.spatial_mode = true` (for CNN models like ACE2):
 *
 * - Input data is packed as [C, H, W] (channel-first)
 * - Inference runs with batch_size=1
 * - Output is unpacked from [C, H, W] format
 *
 * When `spatial_mode = false` (for pointwise MLP models):
 *
 * - Data is packed as [H*W, C] (sample-first)
 * - Inference runs with batch_size=H*W
 *
 * @see EmulatorComp for base class functionality
 * @see InferenceBackend for backend interface
 */
class EmulatorAtm : public EmulatorComp {
public:
  EmulatorAtm();
  ~EmulatorAtm() override = default;

  /**
   * @brief Initialize coupling field indices from MCT field lists.
   *
   * Parses the colon-separated field lists provided by the MCT layer
   * and sets up internal index mappings for import/export operations.
   *
   * @param export_fields Colon-separated list of a2x field names
   * @param import_fields Colon-separated list of x2a field names
   */
  void init_coupling_indices(const std::string &export_fields,
                             const std::string &import_fields);

  /**
   * @brief Configure inference backend.
   *
   * Must be called before initialize() to select the inference backend
   * and set model parameters.
   *
   * @param config Inference configuration (backend type, model path, etc.)
   */
  void set_inference_config(const inference::InferenceConfig &config);

  /**
   * @brief Get current inference configuration.
   * @return Reference to the current InferenceConfig
   */
  const inference::InferenceConfig &get_inference_config() const {
    return m_inference_config;
  }

  /**
   * @brief Set logging file.
   *
   * Redirects the component logger to the specified file.
   * If empty, logging goes to stdout.
   *
   * @param filename Path to log file
   */
  void set_log_file(const std::string &filename);

protected:
  /** @brief Component-specific initialization (load model, read ICs). */
  void init_impl() override;

  /** @brief Execute one time step of inference. */
  void run_impl(int dt) override;

  /** @brief Component-specific finalization. */
  void final_impl() override;

  /** @brief Run AI inference via the configured backend. */
  void run_inference(const std::vector<double> &inputs,
                     std::vector<double> &outputs) override;

private:
  // =========================================================================
  // Coupling and Fields
  // =========================================================================
  CouplingFieldsBase m_coupling_fields;    ///< Base coupling field parser
  impl::AtmCouplingIndices m_coupling_idx; ///< Atmosphere coupling indices
  impl::AtmFieldManager m_fields;          ///< Field storage manager

  // =========================================================================
  // Configuration and Inference
  // =========================================================================
  EmulatorConfig m_config;                       ///< Component config
  inference::InferenceConfig m_inference_config; ///< Backend config
  std::unique_ptr<inference::InferenceBackend>
      m_inference;                        ///< Backend instance
  EmulatorOutputManager m_output_manager; ///< Diagnostic output manager
  std::unique_ptr<impl::AtmFieldDataProvider>
      m_field_provider; ///< Field data adapter for output

  // =========================================================================
  // Helper Methods
  // =========================================================================

  /** @brief Copy fields from coupler import buffer to internal storage. */
  void import_coupling_fields();

  /** @brief Copy fields from internal storage to coupler export buffer. */
  void export_coupling_fields();

  /**
   * @brief Pack input fields into net_inputs tensor.
   *
   * Handles spatial reshaping based on config.model_io.spatial_mode.
   */
  void prepare_inputs();

  /**
   * @brief Unpack net_outputs tensor into output fields.
   *
   * Handles inverse spatial reshaping based on config.model_io.spatial_mode.
   */
  void process_outputs();

  /**
   * @brief Read initial condition fields from a NetCDF file.
   * @param filename Path to IC file
   * @return true if successful
   */
  bool read_initial_conditions(const std::string &filename);
};

} // namespace emulator

#endif // EMULATOR_ATM_HPP
