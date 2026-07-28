/**
 * @file python_interpreter.cpp
 * @brief Embedded CPython session management.
 */

#include "python_interpreter.hpp"

#include <cfenv>

// feenableexcept and friends are glibc extensions.  Everywhere else there is
// nothing to guard against, so FpeGuard becomes a no-op rather than an
// #error: trapping is what E3SM debug builds do on Linux.
#if defined(__GLIBC__)
#define EMULATOR_HAVE_FEENABLEEXCEPT 1
#endif

namespace emulator {
namespace inference {

PyInterpreter &PyInterpreter::instance() {
  static PyInterpreter singleton;
  return singleton;
}

void PyInterpreter::initialize() {
  if (m_customers == 0 && Py_IsInitialized() == 0) {
    // Signal handlers stay with the host process: Python must not install a
    // SIGINT handler in the middle of a batch job.
    Py_InitializeEx(0);
    EMULATOR_INFER_REQUIRE(Py_IsInitialized() != 0,
                           "Could not start the embedded Python interpreter.");
  }
  ++m_customers;
}

void PyInterpreter::finalize() {
  if (m_customers > 0) {
    --m_customers;
  }
  // The interpreter stays up on purpose; see the header.
}

void PyInterpreter::add_sys_path(const std::string &path) {
  if (path.empty()) {
    return;
  }
  PyObject *sys_path = PySys_GetObject("path"); // borrowed
  if (sys_path == nullptr || PyList_Check(sys_path) == 0) {
    return;
  }
  PyRef entry = py_string(path);
  const Py_ssize_t n = PyList_Size(sys_path);
  for (Py_ssize_t i = 0; i < n; ++i) {
    if (PyObject_RichCompareBool(PyList_GetItem(sys_path, i), entry.get(),
                                 Py_EQ) == 1) {
      return; // already there
    }
  }
  PyList_Insert(sys_path, 0, entry.get());
}

// ---------------------------------------------------------------------------

FpeGuard::FpeGuard() {
#ifdef EMULATOR_HAVE_FEENABLEEXCEPT
  m_saved_excepts = fegetexcept();
  if (m_saved_excepts > 0) {
    fedisableexcept(m_saved_excepts);
  }
#endif
}

FpeGuard::~FpeGuard() {
#ifdef EMULATOR_HAVE_FEENABLEEXCEPT
  if (m_saved_excepts > 0) {
    feclearexcept(FE_ALL_EXCEPT);
    feenableexcept(m_saved_excepts);
  }
#endif
}

// ---------------------------------------------------------------------------

PyRef py_string(const std::string &s) {
  return PyRef(
      PyUnicode_FromStringAndSize(s.data(), static_cast<Py_ssize_t>(s.size())));
}

std::string py_to_string(PyObject *obj) {
  if (obj == nullptr) {
    return "";
  }
  PyRef str(PyObject_Str(obj));
  if (!str) {
    PyErr_Clear();
    return "";
  }
  Py_ssize_t len = 0;
  const char *data = PyUnicode_AsUTF8AndSize(str.get(), &len);
  if (data == nullptr) {
    PyErr_Clear();
    return "";
  }
  return std::string(data, static_cast<std::size_t>(len));
}

std::string py_take_error() {
  if (PyErr_Occurred() == nullptr) {
    return "";
  }
  PyObject *type = nullptr;
  PyObject *value = nullptr;
  PyObject *traceback = nullptr;
  PyErr_Fetch(&type, &value, &traceback);
  PyErr_NormalizeException(&type, &value, &traceback);
  PyRef type_ref(type);
  PyRef value_ref(value);
  PyRef traceback_ref(traceback);

  std::string message = py_to_string(value_ref.get());
  if (message.empty()) {
    message = py_to_string(type_ref.get());
  }

  // Format the traceback the way Python would, when there is one: an error
  // inside a model is nearly useless without it.
  if (traceback_ref) {
    PyRef tb_module(PyImport_ImportModule("traceback"));
    PyRef lines =
        tb_module ? PyRef(PyObject_CallMethod(tb_module.get(),
                                              "format_exception", "OOO",
                                              type_ref.get(), value_ref.get(),
                                              traceback_ref.get()))
                  : PyRef();
    if (lines && PyList_Check(lines.get())) {
      std::string formatted;
      const Py_ssize_t n = PyList_Size(lines.get());
      for (Py_ssize_t i = 0; i < n; ++i) {
        formatted += py_to_string(PyList_GetItem(lines.get(), i));
      }
      if (!formatted.empty()) {
        message = formatted;
      }
    }
  }

  PyErr_Clear();
  return message;
}

void py_throw(const std::string &context) {
  const std::string details = py_take_error();
  throw InferenceError("Python error while " + context +
                       (details.empty() ? "." : ":\n" + details));
}

} // namespace inference
} // namespace emulator
