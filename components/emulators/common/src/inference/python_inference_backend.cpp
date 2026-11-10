/**
 * @file python_inference_backend.cpp
 * @brief Python inference backend implementation.
 */

// Python.h must come before any standard header
#include <Python.h>

#include "python_inference_backend.hpp"

#include "fpe_guard.hpp"
#include "inference_error.hpp"

#include <sstream>
#include <vector>

namespace emulator {
namespace inference {

namespace {

/// Holds the GIL for the lifetime of the object; safe to nest.
class GilGuard {
public:
  GilGuard() : m_state(PyGILState_Ensure()) {}
  ~GilGuard() { PyGILState_Release(m_state); }
  GilGuard(const GilGuard &) = delete;
  GilGuard &operator=(const GilGuard &) = delete;

private:
  PyGILState_STATE m_state;
};

/// Owning, move-only handle to a PyObject; takes over a new reference.
class PyRef {
public:
  PyRef() = default;
  explicit PyRef(PyObject *obj) : m_obj(obj) {}
  ~PyRef() { reset(); }

  PyRef(const PyRef &) = delete;
  PyRef &operator=(const PyRef &) = delete;
  PyRef(PyRef &&other) noexcept : m_obj(other.m_obj) { other.m_obj = nullptr; }
  PyRef &operator=(PyRef &&other) noexcept {
    if (this != &other) {
      reset();
      m_obj = other.m_obj;
      other.m_obj = nullptr;
    }
    return *this;
  }

  PyObject *get() const { return m_obj; }
  explicit operator bool() const { return m_obj != nullptr; }

private:
  void reset() {
    if (m_obj != nullptr && Py_IsInitialized() != 0) {
      GilGuard gil;
      Py_DECREF(m_obj);
    }
    m_obj = nullptr;
  }

  PyObject *m_obj = nullptr;
};

std::string to_string(PyObject *obj) {
  PyRef str(obj != nullptr ? PyObject_Str(obj) : nullptr);
  const char *utf8 = str ? PyUnicode_AsUTF8(str.get()) : nullptr;
  if (utf8 == nullptr) {
    PyErr_Clear();
    return "";
  }
  return utf8;
}

/// Clear the pending Python exception and return it with its traceback.
std::string take_error() {
#if PY_VERSION_HEX >= 0x030C0000
  PyRef value(PyErr_GetRaisedException());
#else
  PyObject *t = nullptr, *v = nullptr, *tb = nullptr;
  PyErr_Fetch(&t, &v, &tb);
  PyErr_NormalizeException(&t, &v, &tb);
  if (v != nullptr && tb != nullptr) {
    PyException_SetTraceback(v, tb);
  }
  PyRef type(t), value(v), trace(tb);
#endif
  if (!value) {
    return "";
  }

  // traceback.format_exception(type, value, tb), joined
  PyRef module(PyImport_ImportModule("traceback"));
  PyRef trace_obj(PyException_GetTraceback(value.get()));
  PyRef lines(module ? PyObject_CallMethod(
                           module.get(), "format_exception", "OOO",
                           reinterpret_cast<PyObject *>(Py_TYPE(value.get())),
                           value.get(), trace_obj ? trace_obj.get() : Py_None)
                     : nullptr);
  PyRef empty(PyUnicode_FromString(""));
  PyRef joined(lines && empty ? PyUnicode_Join(empty.get(), lines.get())
                              : nullptr);
  const std::string message = to_string(joined.get());
  PyErr_Clear();
  return message.empty() ? to_string(value.get()) : message;
}

/// Own a new reference, or throw the pending Python exception.
PyRef checked(PyObject *obj, const std::string &doing) {
  if (obj == nullptr) {
    throw InferenceError("Python error while " + doing + ":\n" + take_error());
  }
  return PyRef(obj);
}

void set_item(const PyRef &dict, const std::string &key, const PyRef &value) {
  if (PyDict_SetItemString(dict.get(), key.c_str(), value.get()) != 0) {
    throw InferenceError("Python error while setting '" + key + "':\n" +
                         take_error());
  }
}

PyRef py_string(const std::string &s) {
  return checked(PyUnicode_FromString(s.c_str()), "converting a string");
}

PyRef py_int(long long value) {
  return checked(PyLong_FromLongLong(value), "converting an integer");
}

/// Prepend colon-separated directories to sys.path; the first one wins.
void prepend_sys_path(const std::string &paths) {
  std::vector<std::string> entries;
  std::istringstream stream(paths);
  for (std::string entry; std::getline(stream, entry, ':');) {
    if (!entry.empty()) {
      entries.push_back(entry);
    }
  }
  PyObject *sys_path = PySys_GetObject("path"); // borrowed
  EMULATOR_INFER_REQUIRE(sys_path != nullptr, "Python has no sys.path.");
  for (auto it = entries.rbegin(); it != entries.rend(); ++it) {
    PyRef entry = py_string(*it);
    if (PySequence_Contains(sys_path, entry.get()) == 0) {
      PyList_Insert(sys_path, 0, entry.get());
    }
  }
}

/// The dict handed to the factory.
PyRef config_dict(const InferenceConfig &config) {
  PyRef dict = checked(PyDict_New(), "building the config");
  for (const auto &option : config.options) {
    set_item(dict, option.first, py_string(option.second));
  }
  set_item(dict, "model_path", py_string(config.model_path));
  set_item(dict, "input_channels", py_int(config.input_channels));
  set_item(dict, "output_channels", py_int(config.output_channels));
  set_item(dict, "verbose", PyRef(PyBool_FromLong(config.verbose ? 1 : 0)));
  return dict;
}

/// A numpy array that views (does not copy) a tensor's memory.
PyRef as_numpy(const PyRef &numpy, const Tensor &tensor, const double *data,
               bool writable) {
  const std::string doing = "wrapping tensor '" + tensor.name() + "'";
  const auto &dims = tensor.dims();

  PyRef shape =
      checked(PyTuple_New(static_cast<Py_ssize_t>(dims.size())), doing);
  for (std::size_t i = 0; i < dims.size(); ++i) {
    // PyTuple_SetItem takes over the reference
    PyTuple_SetItem(shape.get(), static_cast<Py_ssize_t>(i),
                    PyLong_FromLongLong(dims[i]));
  }

  if (tensor.size() == 0) {
    // There is no memory to view
    PyRef array = checked(
        PyObject_CallMethod(numpy.get(), "empty", "(O)", shape.get()), doing);
    if (!writable) {
      PyRef flags =
          checked(PyObject_GetAttrString(array.get(), "flags"), doing);
      PyObject_SetAttrString(flags.get(), "writeable", Py_False);
    }
    return array;
  }

  // PyBUF_READ is what makes the array read-only, despite the cast
  PyRef memory = checked(
      PyMemoryView_FromMemory(
          const_cast<char *>(reinterpret_cast<const char *>(data)),
          static_cast<Py_ssize_t>(tensor.size() * sizeof(double)),
          writable ? PyBUF_WRITE : PyBUF_READ),
      doing);
  PyRef flat = checked(PyObject_CallMethod(numpy.get(), "frombuffer", "Os",
                                           memory.get(), "float64"),
                       doing);
  return checked(PyObject_CallMethod(flat.get(), "reshape", "(O)", shape.get()),
                 doing);
}

} // namespace

struct PythonBackend::Impl {
  PyRef numpy;
  PyRef model; ///< What the factory returned
};

PythonBackend::PythonBackend(const InferenceConfig &config)
    : InferenceBackend(config), m_impl(new Impl()) {
  const std::string module_name = m_config.get("python_module");
  const std::string factory_name =
      m_config.get("python_factory", "create_emulator");
  EMULATOR_INFER_REQUIRE(!module_name.empty(),
                         "The Python backend needs the python_module option.");

  if (Py_IsInitialized() == 0) {
    // 0: leave signal handling to the host model
    Py_InitializeEx(0);
  }
  GilGuard gil;
  FpeGuard no_fpe;

  prepend_sys_path(m_config.get("python_path"));
  m_impl->numpy = checked(PyImport_ImportModule("numpy"), "importing numpy");

  PyRef module = checked(PyImport_ImportModule(module_name.c_str()),
                         "importing '" + module_name +
                             "' (is its directory in python_path?)");
  PyRef factory =
      checked(PyObject_GetAttrString(module.get(), factory_name.c_str()),
              "looking up " + module_name + "." + factory_name);
  PyRef settings = config_dict(m_config);
  m_impl->model = checked(
      PyObject_CallFunctionObjArgs(factory.get(), settings.get(), nullptr),
      "calling " + module_name + "." + factory_name + "(config)");
}

PythonBackend::~PythonBackend() = default;

bool PythonBackend::infer(const TensorMap &inputs, TensorMap &outputs) {
  EMULATOR_INFER_REQUIRE(m_impl, "The Python backend was finalized.");

  GilGuard gil;
  FpeGuard no_fpe;

  PyRef in = checked(PyDict_New(), "building the inputs");
  for (const auto &tensor : inputs) {
    set_item(in, tensor.name(),
             as_numpy(m_impl->numpy, tensor, tensor.cdata(), false));
  }
  PyRef out = checked(PyDict_New(), "building the outputs");
  for (auto &tensor : outputs) {
    set_item(out, tensor.name(),
             as_numpy(m_impl->numpy, tensor, tensor.data(), true));
  }

  checked(PyObject_CallMethod(m_impl->model.get(), "infer", "OO", in.get(),
                              out.get()),
          "calling infer(inputs, outputs)");
  return true;
}

void PythonBackend::finalize() {
  if (!m_impl) {
    return;
  }
  // Dropping the model is what frees it, whether or not finalize() succeeds
  const std::unique_ptr<Impl> impl = std::move(m_impl);

  GilGuard gil;
  if (PyObject_HasAttrString(impl->model.get(), "finalize") != 0) {
    checked(PyObject_CallMethod(impl->model.get(), "finalize", nullptr),
            "calling finalize()");
  }
}

} // namespace inference
} // namespace emulator
