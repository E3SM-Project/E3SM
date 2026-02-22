/**
 * @file create_backend.hpp
 * @brief Factory for creating inference backends.
 */

#ifndef E3SM_EMULATOR_CREATE_INFERENCE_BACKEND_HPP
#define E3SM_EMULATOR_CREATE_INFERENCE_BACKEND_HPP

#include <memory>

#include "inference_backend.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Create an inference backend by type.
 * @param type Backend type to create
 * @param config Configuration for the backend
 * @return Shared pointer to new backend instance
 */
std::shared_ptr<InferenceBackend>
create_backend(BackendType type, const InferenceConfig &config);

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_CREATE_INFERENCE_BACKEND_HPP
