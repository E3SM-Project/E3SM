/**
 * @file python_inference_backend.cpp
 * @brief Embedded-Python inference backend implementation.
 *
 * The numpy interaction goes through the ordinary Python API
 * (`numpy.frombuffer` over a memoryview of our own storage) rather than the
 * numpy C API.  That keeps numpy a *runtime* dependency only -- building
 * E3SM never needs numpy headers -- while still handing Python zero-copy
 * views of E3SM memory.
 */

#include "python_inference_backend.hpp"

#include "python_interpreter.hpp"

#include <iostream>
#include <sstream>
#include <vector>

#ifndef EMULATOR_PYTHON_PACKAGE_DIR
#define EMULATOR_PYTHON_PACKAGE_DIR ""
#endif

namespace emulator {
namespace inference {

namespace {

/// Split "a:b:c" into {"a","b","c"}, dropping empties.
std::vector<std::string> split_paths(const std::string &joined) {
  std::vector<std::string> out;
  std::istringstream iss(joined);
  std::string item;
  while (std::getline(iss, item, ':')) {
    if (!item.empty()) {
      out.push_back(item);
    }
  }
  return out;
}

/// A Python tuple of ints, for a shape.
PyRef dims_tuple(const std::vector<std::int64_t> &dims) {
  PyRef tuple(PyTuple_New(static_cast<Py_ssize_t>(dims.size())));
  EMULATOR_INFER_REQUIRE(static_cast<bool>(tuple),
                         "Out of memory building a shape tuple.");
  for (std::size_t i = 0; i < dims.size(); ++i) {
    // PyTuple_SetItem steals the reference it is given.
    PyTuple_SetItem(tuple.get(), static_cast<Py_ssize_t>(i),
                    PyLong_FromLongLong(static_cast<long long>(dims[i])));
  }
  return tuple;
}

void dict_set(PyObject *dict, const char *key, PyRef value) {
  EMULATOR_INFER_REQUIRE(static_cast<bool>(value),
                         "Out of memory building the value for '" << key
                                                                  << "'.");
  PyDict_SetItemString(dict, key, value.get());
}

void dict_set_string(PyObject *dict, const char *key, const std::string &v) {
  dict_set(dict, key, py_string(v));
}

void dict_set_int(PyObject *dict, const char *key, long long v) {
  dict_set(dict, key, PyRef(PyLong_FromLongLong(v)));
}

} // namespace

// ===========================================================================
// Impl
// ===========================================================================

struct PythonBackend::Impl {
  bool interpreter_started = false;

  PyRef numpy;
  PyRef np_frombuffer;
  PyRef emulator;       ///< The object the factory returned
  PyRef infer_callable; ///< Its infer method

  std::string module_name = "e3sm_emulator.bridge";
  std::string infer_method = "infer";
  std::string finalize_method = "finalize";

  /**
   * @brief A numpy array sharing memory with the given buffer.
   *
   * `writable == false` yields a read-only array, which is what a component
   * wants for its input fields: a bug in the model then cannot corrupt
   * state.
   */
  PyRef wrap_buffer(const void *base, std::size_t nbytes,
                    const std::vector<std::int64_t> &dims, const char *dtype,
                    bool writable, const std::string &what) {
    if (nbytes == 0) {
      // frombuffer rejects empty buffers; hand back a correctly shaped empty.
      PyRef args(PyTuple_Pack(1, dims_tuple(dims).get()));
      PyRef kwargs(PyDict_New());
      dict_set_string(kwargs.get(), "dtype", dtype);
      PyRef arr(
          PyObject_Call(numpy.attr("empty").get(), args.get(), kwargs.get()));
      if (!arr) {
        py_throw("allocating an empty array for '" + what + "'");
      }
      // np.empty is writable; an input must not be, even when there is
      // nothing in it, so that the contract does not change on a rank that
      // happens to own no columns.
      if (!writable) {
        PyRef no_args(PyTuple_New(0));
        PyRef flags(PyDict_New());
        PyDict_SetItemString(flags.get(), "write", Py_False);
        PyRef done(PyObject_Call(arr.attr("setflags").get(), no_args.get(),
                                 flags.get()));
        if (!done) {
          py_throw("marking the empty array for '" + what + "' read-only");
        }
      }
      return arr;
    }

    // memoryview over our storage, then frombuffer + reshape.  Both are
    // views: nothing is copied, and writes land in E3SM memory.  The cast
    // away from const is required by the C API even for PyBUF_READ, which is
    // what actually makes the resulting array read-only.
    PyRef view(PyMemoryView_FromMemory(
        const_cast<char *>(static_cast<const char *>(base)),
        static_cast<Py_ssize_t>(nbytes), writable ? PyBUF_WRITE : PyBUF_READ));
    if (!view) {
      py_throw("creating a memoryview for '" + what + "'");
    }

    PyRef args(PyTuple_Pack(1, view.get()));
    PyRef kwargs(PyDict_New());
    dict_set_string(kwargs.get(), "dtype", dtype);
    PyRef flat(PyObject_Call(np_frombuffer.get(), args.get(), kwargs.get()));
    if (!flat) {
      py_throw("wrapping '" + what + "' as a numpy array");
    }
    if (dims.size() <= 1) {
      return flat;
    }

    PyRef reshaped(PyObject_CallMethod(flat.get(), "reshape", "O",
                                       dims_tuple(dims).get()));
    if (!reshaped) {
      py_throw("reshaping '" + what + "'");
    }
    return reshaped;
  }

  /// An owning numpy copy, for metadata Python may outlive or mutate.
  template <typename T>
  PyRef copy_vector(const std::vector<T> &values, const char *dtype,
                    const std::string &what) {
    const auto n = static_cast<std::int64_t>(values.size());
    PyRef view = wrap_buffer(values.data(), values.size() * sizeof(T), {n},
                             dtype, /*writable=*/false, what);
    PyRef copy(PyObject_CallMethod(view.get(), "copy", nullptr));
    if (!copy) {
      py_throw("copying '" + what + "'");
    }
    return copy;
  }

  /// Dict of name -> numpy view for a whole TensorMap.
  PyRef tensor_dict(TensorMap &tensors, bool writable) {
    PyRef dict(PyDict_New());
    EMULATOR_INFER_REQUIRE(static_cast<bool>(dict),
                           "Out of memory building the tensor dict.");
    for (auto &tensor : tensors) {
      const void *base =
          writable ? static_cast<const void *>(tensor.data()) : tensor.cdata();
      PyRef arr = wrap_buffer(base, tensor.nbytes(), tensor.dims(), "float64",
                              writable, tensor.name());
      if (PyDict_SetItemString(dict.get(), tensor.name().c_str(), arr.get()) !=
          0) {
        py_throw("adding '" + tensor.name() + "' to the tensor dict");
      }
    }
    return dict;
  }

  /// The dict handed to the factory: every setting, plus the context.
  PyRef config_dict(const InferenceConfig &config,
                    const InferenceContext &context) {
    PyRef cfg(PyDict_New());
    EMULATOR_INFER_REQUIRE(static_cast<bool>(cfg),
                           "Out of memory building the config dict.");
    dict_set_string(cfg.get(), "backend", config.backend);
    dict_set_string(cfg.get(), "model_path", config.model_path);
    dict_set_int(cfg.get(), "input_channels", config.input_channels);
    dict_set_int(cfg.get(), "output_channels", config.output_channels);
    PyDict_SetItemString(cfg.get(), "verbose",
                         config.verbose ? Py_True : Py_False);
    const auto name_list = [](const std::vector<std::string> &names) {
      PyRef list(PyList_New(0));
      for (const auto &n : names) {
        PyList_Append(list.get(), py_string(n).get());
      }
      return list;
    };
    dict_set(cfg.get(), "inputs", name_list(config.inputs));
    dict_set(cfg.get(), "outputs", name_list(config.outputs));
    for (const auto &kv : config.options) {
      dict_set_string(cfg.get(), kv.first.c_str(), kv.second);
    }

    PyRef ctx(PyDict_New());
    dict_set_int(ctx.get(), "rank", context.rank);
    dict_set_int(ctx.get(), "world_size", context.size);
    dict_set_string(ctx.get(), "node_name", context.node_name);
    dict_set_int(ctx.get(), "nx", context.nx);
    dict_set_int(ctx.get(), "ny", context.ny);
    dict_set_int(ctx.get(), "num_global_cols", context.num_global_cols);
    // Copies, not views: this is small, one-time metadata that the model is
    // free to keep, reorder or sort for as long as it likes.
    dict_set(ctx.get(), "col_gids",
             copy_vector(context.col_gids, "int32", "col_gids"));
    dict_set(ctx.get(), "lat", copy_vector(context.lat, "float64", "lat"));
    dict_set(ctx.get(), "lon", copy_vector(context.lon, "float64", "lon"));
    PyDict_SetItemString(cfg.get(), "context", ctx.get());

    return cfg;
  }
};

// ===========================================================================
// PythonBackend
// ===========================================================================

PythonBackend::PythonBackend(const InferenceConfig &config,
                             const InferenceContext &context)
    : InferenceBackend(config, context), m_impl(new Impl()) {
  m_impl->module_name = config.get("python_module", "e3sm_emulator.bridge");
  m_impl->infer_method = config.get("python_infer_method", "infer");
  m_impl->finalize_method = config.get("python_finalize_method", "finalize");
}

PythonBackend::~PythonBackend() {
  try {
    PythonBackend::finalize();
  } catch (const std::exception &e) {
    // A destructor must not throw; a failed Python teardown is worth a note
    // but not worth aborting a run that is already shutting down.
    std::cerr << "[emulator::inference] warning: Python backend teardown "
              << "failed: " << e.what() << "\n";
  }
}

void PythonBackend::init_impl() {
  PyInterpreter::instance().initialize();
  m_impl->interpreter_started = true;

  // Everything past this point can throw: numpy may be absent, the module may
  // not import, the factory may raise.  InferenceBackend::initialize() only
  // sets m_initialized once init_impl() *returns*, so on that path finalize()
  // would decline to run -- leaking our interpreter customer, and leaving the
  // Python references gathered so far to be released by ~Impl with no GIL
  // held.  That last part is undefined behaviour in a process where somebody
  // else (EAMxx, or a Python host) owns the interpreter, which is exactly the
  // process this backend is meant to survive.  So roll back here.
  try {
    load_model();
  } catch (...) {
    final_impl(); // idempotent, takes the GIL, drops our customer
    throw;
  }
}

void PythonBackend::load_model() {
  PyGilGuard gil;
  auto &interpreter = PyInterpreter::instance();

  // The package shipped alongside this source goes on first so that it ends
  // up *last* in sys.path: a site can override it wholesale.  Within
  // `python_path` the user's first entry should win, so walk it backwards.
  interpreter.add_sys_path(EMULATOR_PYTHON_PACKAGE_DIR);
  if (m_config.get_bool("python_add_cwd", true)) {
    interpreter.add_sys_path(".");
  }
  const auto paths = split_paths(m_config.get("python_path"));
  for (auto it = paths.rbegin(); it != paths.rend(); ++it) {
    interpreter.add_sys_path(*it);
  }

  {
    FpeGuard fpe_guard; // importing numpy raises benign FPEs
    m_impl->numpy = PyRef(PyImport_ImportModule("numpy"));
  }
  if (!m_impl->numpy) {
    throw InferenceError(
        "The python inference backend needs numpy in the interpreter E3SM is "
        "linked against, but importing it failed:\n" +
        py_take_error());
  }
  m_impl->np_frombuffer = m_impl->numpy.attr("frombuffer");
  EMULATOR_INFER_REQUIRE(static_cast<bool>(m_impl->np_frombuffer),
                         "The imported 'numpy' has no frombuffer; is it "
                         "really numpy?");

  PyRef module;
  {
    FpeGuard fpe_guard; // the module imports torch, and torch traps too
    module = PyRef(PyImport_ImportModule(m_impl->module_name.c_str()));
  }
  if (!module) {
    throw InferenceError("Could not import the python emulator module '" +
                         m_impl->module_name +
                         "'. Is its directory on sys.path (setting "
                         "`inference.python_path`)?\n" +
                         py_take_error());
  }

  const std::string factory_name =
      m_config.get("python_factory", "create_emulator");
  PyRef factory = module.attr(factory_name.c_str());
  EMULATOR_INFER_REQUIRE(static_cast<bool>(factory) &&
                             PyCallable_Check(factory.get()) != 0,
                         "Module '" << m_impl->module_name
                                    << "' has no callable '" << factory_name
                                    << "' (setting "
                                       "`inference.python_factory`).");
  {
    PyRef cfg = m_impl->config_dict(m_config, m_context);
    FpeGuard fpe_guard; // the factory loads the model
    m_impl->emulator =
        PyRef(PyObject_CallFunctionObjArgs(factory.get(), cfg.get(), nullptr));
  }
  if (!m_impl->emulator) {
    py_throw("calling " + m_impl->module_name + "." + factory_name + "(config)");
  }

  m_impl->infer_callable = m_impl->emulator.attr(m_impl->infer_method.c_str());
  EMULATOR_INFER_REQUIRE(
      static_cast<bool>(m_impl->infer_callable) &&
          PyCallable_Check(m_impl->infer_callable.get()) != 0,
      "The python emulator from '"
          << m_impl->module_name << "' has no callable '"
          << m_impl->infer_method
          << "(inputs, outputs)' method (setting "
             "`inference.python_infer_method`).");
}

bool PythonBackend::infer_impl(const TensorMap &inputs, TensorMap &outputs) {
  PyGilGuard gil;

  // Read-only views for the inputs, writable views for the outputs.  The
  // const_cast is confined here: tensor_dict never writes through an input,
  // it only needs a non-const TensorMap& to iterate.
  PyRef in_dict =
      m_impl->tensor_dict(const_cast<TensorMap &>(inputs), /*writable=*/false);
  PyRef out_dict = m_impl->tensor_dict(outputs, /*writable=*/true);

  FpeGuard fpe_guard; // model kernels are entitled to raise benign FPEs
  PyRef result(PyObject_CallFunctionObjArgs(
      m_impl->infer_callable.get(), in_dict.get(), out_dict.get(), nullptr));
  if (!result) {
    py_throw("calling " + m_impl->module_name + "." + m_impl->infer_method +
             "(inputs, outputs)");
  }
  return true;
}

void PythonBackend::final_impl() {
  if (!m_impl->interpreter_started) {
    return;
  }

  {
    PyGilGuard gil;
    if (m_impl->emulator && !m_impl->finalize_method.empty()) {
      PyRef fin = m_impl->emulator.attr(m_impl->finalize_method.c_str());
      if (fin && PyCallable_Check(fin.get()) != 0) {
        PyRef result(PyObject_CallObject(fin.get(), nullptr));
        if (!result) {
          // Report, but do not throw: teardown must not break a run that has
          // already produced its answers.
          std::cerr << "[emulator::inference] warning: " << m_impl->module_name
                    << "." << m_impl->finalize_method << "() failed:\n"
                    << py_take_error() << "\n";
        }
      }
    }

    // Drop our references while the interpreter is still alive.  This is what
    // actually releases the model and its weights.
    m_impl->infer_callable = PyRef();
    m_impl->emulator = PyRef();
    m_impl->np_frombuffer = PyRef();
    m_impl->numpy = PyRef();
  }

  PyInterpreter::instance().finalize();
  m_impl->interpreter_started = false;
}

} // namespace inference
} // namespace emulator
