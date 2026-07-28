/**
 * @file python_inference_backend.hpp
 * @brief Inference backend that keeps the model in Python.
 */

#ifndef E3SM_EMULATOR_PYTHON_INFERENCE_BACKEND_HPP
#define E3SM_EMULATOR_PYTHON_INFERENCE_BACKEND_HPP

#include <memory>
#include <string>

#include "inference_backend.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Runs a Python emulator in this process.
 *
 * At initialization the backend imports a module (`python_module`, default
 * `e3sm_emulator.bridge`) and calls a factory (`python_factory`, default
 * `create_emulator`) with a dict holding every configuration option plus a
 * `context` entry describing this rank and the columns it owns.  The factory
 * returns an object with an `infer(inputs, outputs)` method, called once per
 * step:
 *
 * @code{.py}
 *     def infer(self, inputs, outputs):
 *         outputs["dT"][:] = self.model(inputs["T"])
 * @endcode
 *
 * `inputs` and `outputs` are dicts of numpy arrays that *view E3SM memory*
 * directly -- nothing is copied in either direction, and the model writes
 * its results into the component's own arrays.  Input arrays are read-only,
 * so a bug in a model cannot corrupt component state.  Zero-copy ends at the
 * numpy view: a model that moves the array to a GPU or converts it to
 * float32 pays for that copy like anybody else.
 *
 * The backend hands the model the *component's* rank and size (see
 * InferenceContext) rather than the job's, so a model is at least told the
 * truth about which processes it shares work with.  That is as far as this
 * goes today: nothing here establishes a rendezvous or builds a process
 * group, and no distributed model has been run through it.  Python is
 * nonetheless the backend worth having first, because it is the only one in
 * which a model *could* contain a collective -- neither ONNX nor TorchScript
 * can express one.
 */
class PythonBackend : public InferenceBackend {
public:
  PythonBackend(const InferenceConfig &config, const InferenceContext &context);
  ~PythonBackend() override;

  /// @copydoc InferenceBackend::name
  std::string name() const override { return "Python"; }

protected:
  void init_impl() override;
  bool infer_impl(const TensorMap &inputs, TensorMap &outputs) override;
  void final_impl() override;

private:
  /// @brief Import numpy and the module, and build the emulator.
  ///
  /// Separate from init_impl() so that everything it acquires can be rolled
  /// back in one place when any step of it throws.
  void load_model();

  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_PYTHON_INFERENCE_BACKEND_HPP
