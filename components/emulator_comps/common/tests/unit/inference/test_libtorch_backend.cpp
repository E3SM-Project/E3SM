/**
 * @file test_libtorch_backend.cpp
 * @brief Unit tests for the LibTorch inference backend.
 *
 * Tests model loading, inference execution, and basic validation
 * of the LibTorch backend implementation.
 *
 * @note Requires a TorchScript model file (dummy_model.pt) to be present
 *       in the test data directory. Run create_dummy_model.py to generate it.
 */

#include "../../src/inference/inference_backend.hpp"
#include "../../src/inference/libtorch_backend.hpp"
#include "../test_support/test_config.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace emulator::inference;

/**
 * @brief Main test entry point for LibTorch backend.
 *
 * Tests:
 * - Backend type identification
 * - Model loading from .pt file
 * - Single-sample inference
 * - Multi-sample batch inference
 * - Output value validation
 *
 * @return 0 on success, 77 if test should be skipped, 1 on failure
 */
int main(int argc, char **argv) {
  (void)argc;
  (void)argv;

#ifdef EMULATOR_HAS_LIBTORCH
  std::cout << "=== Testing LibTorch Backend ===" << std::endl;

  // Get model path from test configuration
  std::string model_path =
      emulator::testing::get_test_data_file("dummy_model.pt");

  // Check if model exists
  std::ifstream model_file(model_path);
  if (!model_file.good()) {
    std::cerr << "Dummy model not found: " << model_path << std::endl;
    std::cerr << "  Run: python3 tests/test_support/data/create_dummy_model.py"
              << std::endl;
    return 77; // SKIP exit code
  }
  model_file.close();

  // Create and verify backend type
  auto backend = std::make_unique<LibTorchBackend>();
  assert(backend->type() == BackendType::LIBTORCH);
  assert(backend->name() == "LibTorch");

  // Configure backend
  InferenceConfig config;
  config.backend = BackendType::LIBTORCH;
  config.model_path = model_path;
  config.input_channels = 10;
  config.output_channels = 5;
  config.device_id = -1; // CPU
  config.use_fp16 = false;

  // Initialize
  bool init_success = backend->initialize(config);
  assert(init_success);
  assert(backend->is_initialized());
  std::cout << "  Initialization: PASSED" << std::endl;

  // Test single batch inference
  {
    const int batch_size = 1;
    std::vector<double> inputs(batch_size * config.input_channels, 1.0);
    std::vector<double> outputs(batch_size * config.output_channels, -999.0);

    bool infer_success =
        backend->infer(inputs.data(), outputs.data(), batch_size);
    assert(infer_success);
    std::cout << "  Single batch inference: PASSED" << std::endl;

    // Validate output values
    // Model is linear: weight=0.1, bias=0, input=1.0
    // Expected: output = 0.1 * sum(inputs) = 0.1 * 10 = 1.0
    for (int i = 0; i < config.output_channels; ++i) {
      double expected = 1.0;
      double diff = std::abs(outputs[i] - expected);
      if (diff > 0.01) {
        std::cerr << "  Output[" << i << "] = " << outputs[i] << ", expected ~"
                  << expected << std::endl;
        assert(false);
      }
    }
    std::cout << "  Output validation: PASSED" << std::endl;
  }

  // Test multi-batch inference
  {
    const int batch_size = 4;
    std::vector<double> inputs(batch_size * config.input_channels);
    std::vector<double> outputs(batch_size * config.output_channels, -999.0);

    // Fill with varying values per batch item
    for (int b = 0; b < batch_size; ++b) {
      for (int i = 0; i < config.input_channels; ++i) {
        inputs[b * config.input_channels + i] = (b + 1) * 0.5;
      }
    }

    bool infer_success =
        backend->infer(inputs.data(), outputs.data(), batch_size);
    assert(infer_success);
    std::cout << "  Multi-batch inference (batch_size=" << batch_size
              << "): PASSED" << std::endl;
  }

  // Test finalization
  backend->finalize();
  assert(!backend->is_initialized());
  std::cout << "  Finalization: PASSED" << std::endl;

  std::cout << "\n=== All LibTorch Backend Tests Passed ===" << std::endl;
  return 0;

#else
  std::cout << "=== LibTorch Backend Not Available ===" << std::endl;
  std::cout << "Skipping test (EMULATOR_HAS_LIBTORCH not defined)" << std::endl;
  return 77; // SKIP exit code
#endif
}
