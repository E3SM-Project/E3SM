/**
 * @file inference_factory.hpp
 * @brief Factory for creating inference backends.
 */

#ifndef INFERENCE_FACTORY_HPP
#define INFERENCE_FACTORY_HPP

#include <memory>

#include "inference_backend.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Create an inference backend by type.
 * @param type Backend type to create
 * @return Unique pointer to new backend instance
 */
std::unique_ptr<InferenceBackend> create_backend(BackendType type);

} // namespace inference
} // namespace emulator

#endif // INFERENCE_FACTORY_HPP
