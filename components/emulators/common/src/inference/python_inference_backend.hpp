/**
 * @file python_inference_backend.hpp
 * @brief Inference backend that runs a Python model in-process.
 */

#ifndef E3SM_EMULATOR_PYTHON_INFERENCE_BACKEND_HPP
#define E3SM_EMULATOR_PYTHON_INFERENCE_BACKEND_HPP

#include <memory>

#include "inference_backend.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Runs a model written in Python through an embedded interpreter.
 *
 * On construction the backend imports a module and calls its factory with
 * a dict holding the configuration (model_path, input_channels,
 * output_channels, verbose, and every option as a string). The factory
 * returns an object with an `infer(inputs, outputs)` method and,
 * optionally, a `finalize()` method:
 *
 * @code{.py}
 *     def create_emulator(config):
 *         return MyModel(config["model_path"])
 *
 *     class MyModel:
 *         def infer(self, inputs, outputs):
 *             outputs["dT"][:] = self.net(inputs["T"])
 * @endcode
 *
 * `inputs` and `outputs` are dicts of float64 numpy arrays, keyed by tensor
 * name, that view the tensors' memory: nothing is copied, inputs are
 * read-only, and the model writes its results in place.
 *
 * Options read from `config.options`:
 * - `python_module`  Module to import (required).
 * - `python_factory` Factory function in it (default `create_emulator`).
 * - `python_path`    Colon-separated directories to prepend to sys.path.
 *
 * Only the CPython C API is used, so numpy is needed at run time but not
 * to build. If the process already runs an interpreter (e.g. EAMxx's), it
 * is shared. The interpreter is never shut down, as extension modules such
 * as numpy cannot be loaded twice in one process.
 *
 * @see InferenceBackend for the base interface
 */
class PythonBackend : public InferenceBackend {
public:
  /// @throws InferenceError if the module or its factory fails
  explicit PythonBackend(const InferenceConfig &config);
  ~PythonBackend() override;

  using InferenceBackend::infer;

  /// @copydoc InferenceBackend::infer
  bool infer(const TensorMap &inputs, TensorMap &outputs) override;

  /// @copydoc InferenceBackend::finalize
  void finalize() override;

  /// @copydoc InferenceBackend::name
  std::string name() const override { return "Python"; }

private:
  struct Impl; ///< Keeps Python.h out of this header
  std::unique_ptr<Impl> m_impl;
};

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_PYTHON_INFERENCE_BACKEND_HPP
