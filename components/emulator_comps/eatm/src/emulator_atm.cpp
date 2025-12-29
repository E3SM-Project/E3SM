/**
 * @file emulator_atm.cpp
 * @brief Atmosphere emulator component implementation.
 *
 * Implements the atmospheric component using AI/ML-based inference
 * for atmospheric physics and dynamics emulation.
 */

#include "emulator_atm.hpp"
#include "emulator_config.hpp"
#include "impl/atm_io.hpp"
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>

namespace emulator {

EmulatorAtm::EmulatorAtm() : EmulatorComp(CompType::ATM) {
  m_inference_config.backend = inference::BackendType::STUB;
  m_inference_config.verbose = false;
}

void EmulatorAtm::set_log_file(const std::string &filename) {
  m_logger.set_file(filename);
}

void EmulatorAtm::set_inference_config(
    const inference::InferenceConfig &config) {
  m_inference_config = config;

  if (is_root()) {
    m_logger.info("[EmulatorAtm] " + get_comp_name(m_type) +
                  " inference config set:\n          backend = " +
                  inference::backend_type_to_string(config.backend) +
                  "\n          model_path = " + config.model_path);
  }
}

void EmulatorAtm::init_coupling_indices(const std::string &export_fields,
                                        const std::string &import_fields) {
  m_coupling_fields.initialize(export_fields, import_fields);
  m_coupling_idx.initialize(m_coupling_fields);
}

void EmulatorAtm::init_impl() {
  // Load YAML configuration
  m_config =
      parse_emulator_config_with_defaults(m_input_file, "eatm", is_root());

  // Configure inference backend from build config
  if (m_config.build.inference_backend == "libtorch") {
    m_inference_config.backend = inference::BackendType::LIBTORCH;
  } else {
    m_inference_config.backend = inference::BackendType::STUB;
  }

  // Model path from runtime config
  m_inference_config.model_path = m_config.runtime.model_path;

  if (is_root()) {
    m_logger.info("[EmulatorAtm] inference config:");
    m_logger.info("  backend=" + inference::backend_type_to_string(
                                     m_inference_config.backend));
    m_logger.info("  model_path=" + m_inference_config.model_path);
    m_logger.info("  spatial_mode=" + std::string(m_config.model_io.spatial_mode
                                                      ? "true"
                                                      : "false"));
    m_logger.info("  input_channels=" +
                  std::to_string(m_config.model_io.input_variables.size()) +
                  " output_channels=" +
                  std::to_string(m_config.model_io.output_variables.size()));
  }

  // Allocate all field storage
  m_fields.allocate(m_num_local_cols);

  // Read initial conditions
  if (!m_config.runtime.ic_file.empty()) {
    if (!read_initial_conditions(m_config.runtime.ic_file)) {
      throw std::runtime_error("[EmulatorAtm] Failed to read IC file: " +
                               m_config.runtime.ic_file);
    }
  } else if (m_inference_config.backend == inference::BackendType::STUB) {
    // Test mode: STUB backend without IC file uses defaults
    m_fields.set_defaults(m_num_local_cols);
  } else {
    throw std::runtime_error(
        "[EmulatorAtm] runtime.ic_file is required for non-stub backends");
  }

  // Set up inference channel counts
  // For spatial_mode: input_channels = C (the backend receives [1, C*H*W])
  // For pointwise: input_channels = C (backend receives [H*W, C])
  m_inference_config.input_channels = m_config.model_io.input_variables.size();
  m_inference_config.output_channels =
      m_config.model_io.output_variables.size();

  // Set spatial mode and grid dimensions for proper tensor reshaping
  m_inference_config.spatial_mode = m_config.model_io.spatial_mode;
  m_inference_config.grid_height = m_ny;
  m_inference_config.grid_width = m_nx;

  if (is_root()) {
    m_logger.info("[EmulatorAtm] Creating inference backend...");
  }

  m_inference = inference::create_backend(m_inference_config);

  if (!m_inference->initialize(m_inference_config)) {
    if (is_root()) {
      m_logger.warn(
          "[EmulatorAtm] Failed to initialize inference backend, using stub");
    }
    // Fallback to stub
    m_inference_config.backend = inference::BackendType::STUB;
    m_inference = inference::create_backend(m_inference_config);
    m_inference->initialize(m_inference_config);
  }

  if (is_root()) {
    m_logger.info("[EmulatorAtm] Inference backend initialized.");
  }

  // Create field data provider for output manager
  m_field_provider =
      std::make_unique<impl::AtmFieldDataProvider>(m_fields, m_num_local_cols);

  // Initialize diagnostic output manager
  std::string case_name = "emulator"; // TODO: Get from runtime config
  m_output_manager.initialize(m_config.diagnostics, m_comm, m_col_gids, m_ny,
                              m_nx, case_name, ".", m_logger);
  m_output_manager.setup(*m_field_provider);

  if (is_root()) {
    m_logger.info("[EmulatorAtm] Diagnostic output manager initialized "
                  "(" +
                  std::to_string(m_config.diagnostics.history_streams.size()) +
                  " history stream(s))");
  }

  // Export initial values to coupler
  export_coupling_fields();
}

void EmulatorAtm::run_impl(int dt) {
  // 1. Import fields from coupler
  import_coupling_fields();

  // 2. Prepare AI model inputs (with optional spatial reshape)
  prepare_inputs();

  // 3. Run AI inference
  run_inference(m_fields.net_inputs, m_fields.net_outputs);

  // 4. Process AI outputs (with optional spatial reshape)
  process_outputs();

  // 5. Diagnostic output step
  // Detect any new stacked fields from AI output (e.g., wind_0, wind_1)
  m_field_provider->detect_stacked_fields();

  // Run output manager with current field state
  m_output_manager.init_timestep(m_step_count, dt);
  m_output_manager.run(m_step_count, *m_field_provider);

  // 6. Restart file writing (if this is a restart step)
  if (m_output_manager.is_restart_step(m_step_count)) {
    m_output_manager.write_restart(*m_field_provider, m_step_count);
    m_output_manager.write_history_restart(m_step_count);
  }

  // 7. Export fields to coupler
  export_coupling_fields();
}

void EmulatorAtm::import_coupling_fields() {
  if (m_import_data == nullptr) {
    if (is_root()) {
      m_logger.error(
          "[EmulatorAtm] import_coupling_fields: No import data pointer!");
    }
    return;
  }

  impl::import_atm_fields(m_import_data, m_num_local_cols, m_num_imports,
                          m_coupling_idx, m_fields);
}

void EmulatorAtm::export_coupling_fields() {
  if (m_export_data == nullptr) {
    if (is_root()) {
      m_logger.error("[EmulatorAtm] export_coupling_fields: No export data!");
    }
    return;
  }

  impl::export_atm_fields(m_export_data, m_num_local_cols, m_num_exports,
                          m_coupling_idx, m_fields);
}

/**
 * @brief Pack input fields into flattened tensor for inference.
 *
 * Field data is stored as separate vectors per variable.
 * This function packs them into a format suitable for the inference backend:
 *
 * - For spatial_mode=true (CNN models like ACE2):
 *   Input: separate fields [ncols] each where ncols = H*W
 *   Output: [1, C, H, W] flattened as [C*H*W] in memory
 *   The backend receives this directly as [1, C, H, W]
 *
 *   CRITICAL: The spatial ordering of columns must match the expected (H,W)
 * layout. E3SM's domain decomposition provides columns in a specific order
 * determined by the grid file. For global uniform grids (ne4, ne30, etc.),
 * columns are typically ordered as (lat, lon) pairs. The user must ensure:
 *   1. m_ny, m_nx match the grid dimensions (e.g., 180x360 for ne4)
 *   2. Column ordering in the grid file matches PyTorch's [H, W] row-major
 * convention
 *   3. For non-structured grids, spatial_mode should NOT be used
 *
 *   Memory layout after packing (C++ row-major, PyTorch NCHW):
 *   net_inputs[c*H*W + h*W + w] corresponds to channel c, height h, width w
 *
 * - For spatial_mode=false (pointwise MLP models):
 *   Input: separate fields [ncols] each
 *   Output: [H*W, C] in row-major order
 *   The backend receives this as [batch_size=H*W, channels=C]
 */
void EmulatorAtm::prepare_inputs() {
  const auto &input_vars = m_config.model_io.input_variables;

  if (input_vars.empty()) {
    if (is_root()) {
      m_logger.warn("[EmulatorAtm] No input variables configured!");
    }
    return;
  }

  const int C_in = static_cast<int>(input_vars.size());
  const int H = m_ny;
  const int W = m_nx;
  const int HW = H * W; // Should equal m_num_local_cols

  // Allocate net_inputs and net_outputs
  if (m_config.model_io.spatial_mode) {
    // Spatial mode: flatten to [C, H, W] = C*H*W total elements
    m_fields.net_inputs.resize(static_cast<size_t>(C_in * HW));
    m_fields.net_outputs.resize(m_config.model_io.output_variables.size() *
                                static_cast<size_t>(HW));
  } else {
    // Pointwise mode: [H*W, C]
    m_fields.net_inputs.resize(static_cast<size_t>(HW * C_in));
    m_fields.net_outputs.resize(static_cast<size_t>(HW) *
                                m_config.model_io.output_variables.size());
  }

  if (m_config.model_io.spatial_mode) {
    // SPATIAL MODE: Pack as [C, H, W] (channel-major)
    // For each channel c, copy all H*W spatial values contiguously
    // This matches PyTorch NCHW layout (without the N dimension)
    for (int c = 0; c < C_in; ++c) {
      std::vector<double> *field_ptr = m_fields.get_field_ptr(input_vars[c]);

      if (!field_ptr) {
        if (is_root()) {
          m_logger.error("[EmulatorAtm] Missing input field: " + input_vars[c]);
        }
        continue;
      }

      // Copy entire field as one contiguous block for this channel
      // net_inputs[c * HW ... (c+1)*HW - 1] = field[0 ... HW-1]
      std::memcpy(&m_fields.net_inputs[c * HW], field_ptr->data(),
                  HW * sizeof(double));
    }
  } else {
    // POINTWISE MODE: Pack as [batch_size, channels] = [HW, C]
    // For each grid point, interleave all channels
    for (int c = 0; c < C_in; ++c) {
      std::vector<double> *field_ptr = m_fields.get_field_ptr(input_vars[c]);

      if (!field_ptr) {
        if (is_root()) {
          m_logger.error("[EmulatorAtm] Missing input field: " + input_vars[c]);
        }
        continue;
      }

      for (int col = 0; col < HW; ++col) {
        m_fields.net_inputs[col * C_in + c] = (*field_ptr)[col];
      }
    }
  }
}

/**
 * @brief Unpack output tensor from inference into separate field vectors.
 *
 * The inverse of prepare_inputs():
 *
 * - For spatial_mode=true (CNN models):
 *   Input: [1, C, H, W] flattened (directly from backend)
 *   Memory layout: net_outputs[c*H*W + h*W + w]
 *   Output: separate fields [ncols] each
 *
 *   The unpacking preserves the spatial structure: each output channel is
 *   extracted as a contiguous block of H*W values, which should match the
 *   original column ordering from the grid decomposition.
 *
 * - For spatial_mode=false (pointwise MLP models):
 *   Input: [H*W, C] in row-major order
 *   Output: separate fields [ncols] each
 */
void EmulatorAtm::process_outputs() {
  const auto &output_vars = m_config.model_io.output_variables;

  if (output_vars.empty()) {
    return;
  }

  const int C_out = static_cast<int>(output_vars.size());
  const int H = m_ny;
  const int W = m_nx;
  const int HW = H * W;

  // Verify output size
  size_t expected_size = static_cast<size_t>(C_out * HW);
  if (m_fields.net_outputs.size() < expected_size) {
    if (is_root()) {
      m_logger.error("[EmulatorAtm] net_outputs size mismatch. Expected " +
                     std::to_string(expected_size) + ", got " +
                     std::to_string(m_fields.net_outputs.size()));
    }
    return;
  }

  if (m_config.model_io.spatial_mode) {
    // SPATIAL MODE: Unpack from [C, H, W]
    for (int c = 0; c < C_out; ++c) {
      m_fields.register_dynamic_field(output_vars[c]);
      std::vector<double> *field_ptr = m_fields.get_field_ptr(output_vars[c]);

      if (!field_ptr) {
        continue;
      }

      // Ensure field is correctly sized
      if (field_ptr->size() != static_cast<size_t>(HW)) {
        field_ptr->resize(HW);
      }

      // Copy contiguous channel data back to field
      std::memcpy(field_ptr->data(), &m_fields.net_outputs[c * HW],
                  HW * sizeof(double));
    }
  } else {
    // POINTWISE MODE: Unpack from [HW, C]
    for (int c = 0; c < C_out; ++c) {
      m_fields.register_dynamic_field(output_vars[c]);
      std::vector<double> *field_ptr = m_fields.get_field_ptr(output_vars[c]);

      if (!field_ptr) {
        continue;
      }

      if (field_ptr->size() != static_cast<size_t>(HW)) {
        field_ptr->resize(HW);
      }

      for (int col = 0; col < HW; ++col) {
        (*field_ptr)[col] = m_fields.net_outputs[col * C_out + c];
      }
    }
  }
}

/**
 * @brief Run inference using the configured backend.
 *
 * For spatial_mode, the backend is called with batch_size=1 and the
 * full [C*H*W] flattened tensor. For pointwise mode, batch_size=H*W.
 *
 * @note If inference fails, this function will abort the simulation
 * with a clear error message.
 */
void EmulatorAtm::run_inference(const std::vector<double> &inputs,
                                std::vector<double> &outputs) {
  if (!m_inference || !m_inference->is_initialized()) {
    m_logger.error("[EmulatorAtm] FATAL: run_inference() called but no "
                   "backend is initialized!");
    std::cerr << "\n*** EMULATOR ABORT: Inference backend not initialized ***\n"
              << std::endl;
    MPI_Abort(m_comm, 1);
  }

  const int C_in = static_cast<int>(m_config.model_io.input_variables.size());
  const int C_out = static_cast<int>(m_config.model_io.output_variables.size());
  const int HW = m_ny * m_nx;

  // Ensure output buffer is sized correctly
  size_t required_size = static_cast<size_t>(C_out * HW);
  if (outputs.size() != required_size) {
    outputs.resize(required_size);
  }

  bool success = false;

  if (m_config.model_io.spatial_mode) {
    // SPATIAL MODE for CNN models:
    // We pass batch_size=1 and input_channels = C*H*W
    // The data is in [C, H, W] flattened format
    // The backend (LibTorchBackend) will reshape to [1, C, H, W] before
    // inference
    success = m_inference->infer(inputs.data(), outputs.data(), 1);
  } else {
    // POINTWISE MODE for MLP models:
    // Each grid cell is a separate sample
    // Backend receives [H*W, C]
    success = m_inference->infer(inputs.data(), outputs.data(), HW);
  }

  if (!success) {
    m_logger.error("[EmulatorAtm] FATAL: Inference failed!");
    m_logger.error("[EmulatorAtm] Input shape: [" +
                   std::to_string(m_config.model_io.spatial_mode ? 1 : HW) +
                   ", " + std::to_string(C_in) + "]");
    m_logger.error("[EmulatorAtm] Expected output: [" +
                   std::to_string(m_config.model_io.spatial_mode ? 1 : HW) +
                   ", " + std::to_string(C_out) + "]");
    m_logger.error("[EmulatorAtm] Check model path and input/output "
                   "configuration in atm_in");
    std::cerr << "\n*** EMULATOR ABORT: Inference failed! ***\n"
              << "Check log for details.\n"
              << std::endl;
    MPI_Abort(m_comm, 1);
  }
}

void EmulatorAtm::final_impl() {
  if (is_root()) {
    m_logger.info("[EmulatorAtm] Finalizing...");
  }

  // Write final restart files if enabled
  if (m_field_provider) {
    m_output_manager.write_restart(*m_field_provider, m_step_count);
    m_output_manager.write_history_restart(m_step_count);
  }

  // Finalize output manager
  m_output_manager.finalize();

  if (m_inference) {
    m_inference->finalize();
    m_inference.reset();
  }

  m_fields.deallocate();
}

bool EmulatorAtm::read_initial_conditions(const std::string &filename) {
  return impl::read_atm_initial_conditions(
      filename, m_num_global_cols, m_num_local_cols, m_col_gids, m_lat,
      m_fields, m_config.model_io.input_variables, m_logger, is_root());
}

} // namespace emulator
