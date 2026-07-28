/**
 * @file create_inference_backend.hpp
 * @brief Factory for creating inference backends.
 */

#ifndef E3SM_EMULATOR_CREATE_INFERENCE_BACKEND_HPP
#define E3SM_EMULATOR_CREATE_INFERENCE_BACKEND_HPP

#include <memory>
#include <string>
#include <vector>

#include "inference_backend.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Create and initialize the backend named by `config.backend`.
 *
 * The returned backend is ready to infer.  For a distributed model this is
 * a collective operation -- every rank of the component communicator must
 * call it -- because that is when the model builds its process group.
 *
 * @param config  Backend settings
 * @param context Ranks and decomposition from the coupler
 * @throws InferenceError if the name is unknown or that backend was not built
 */
std::shared_ptr<InferenceBackend>
create_backend(const InferenceConfig &config, const InferenceContext &context);

/// @brief Backend names this build can create, for error messages and tools.
std::vector<std::string> available_backends();

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_CREATE_INFERENCE_BACKEND_HPP
