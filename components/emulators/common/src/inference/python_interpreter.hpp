/**
 * @file python_interpreter.hpp
 * @brief Embedded CPython session management and small RAII helpers.
 *
 * Only compiled when EMULATOR_ENABLE_PYTHON is on.  Everything here uses the
 * plain CPython C API so that embedding needs nothing beyond the Python
 * development headers -- no pybind11 build dependency, and numpy stays a
 * *runtime* dependency only.  EAMxx solves the same problem with pybind11 in
 * `share/core/eamxx_pysession.hpp`; the semantics here are deliberately the
 * same (reference-counted session, never finalize an interpreter we did not
 * start) so that the two can coexist in one process.
 */

#ifndef E3SM_EMULATOR_PYTHON_INTERPRETER_HPP
#define E3SM_EMULATOR_PYTHON_INTERPRETER_HPP

// Python.h must come before the standard headers it configures.
#include <Python.h>

#include <string>

#include "inference_error.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Reference-counted embedded Python interpreter.
 *
 * Process-wide with a customer count, because several emulators may want
 * Python in one process.  If the interpreter is already running when the
 * first customer arrives -- the process is a Python process, or EAMxx's
 * PySession got there first -- this attaches to it and never finalizes it.
 */
class PyInterpreter {
public:
  static PyInterpreter &instance();

  PyInterpreter(const PyInterpreter &) = delete;
  PyInterpreter &operator=(const PyInterpreter &) = delete;

  /// @brief Start (or attach to) the interpreter and add one customer.
  void initialize();

  /**
   * @brief Drop one customer.
   *
   * Deliberately does *not* stop the interpreter, even at zero customers:
   * CPython cannot be brought back up once a C extension has been imported
   * (numpy 2.x refuses with "cannot load module more than once per process",
   * and PyTorch behaves the same way), so a second emulator created later in
   * the same run would fail to import anything.  What actually holds memory
   * -- the model and its weights -- is released when the backend drops its
   * Python references.
   */
  void finalize();

  /// @brief How many backends currently hold the interpreter.
  ///
  /// Exposed so a test can pin that a failed initialization does not leak
  /// one, which is otherwise invisible until a much later run.
  int num_customers() const { return m_customers; }

  /**
   * @brief Prepend a directory to sys.path.  Requires the GIL.
   *
   * Prepending means a run directory can shadow an installed module, which
   * is what someone editing an emulator in place expects.  Duplicates are
   * not added.
   */
  void add_sys_path(const std::string &path);

private:
  PyInterpreter() = default;

  int m_customers = 0;
};

/**
 * @brief RAII GIL acquisition for the current thread.
 *
 * Every call into Python must be wrapped in one of these.  Safe when the
 * calling thread already holds the GIL.
 */
class PyGilGuard {
public:
  PyGilGuard() : m_state(PyGILState_Ensure()) {}
  ~PyGilGuard() { PyGILState_Release(m_state); }
  PyGilGuard(const PyGilGuard &) = delete;
  PyGilGuard &operator=(const PyGilGuard &) = delete;

private:
  PyGILState_STATE m_state;
};

/**
 * @brief RAII disable of floating-point exception traps.
 *
 * Importing numpy raises benign FPEs; with trapping enabled -- as an E3SM
 * debug build does -- the process dies inside the import.  EAMxx's PySession
 * wraps its imports the same way.
 */
class FpeGuard {
public:
  FpeGuard();
  ~FpeGuard();
  FpeGuard(const FpeGuard &) = delete;
  FpeGuard &operator=(const FpeGuard &) = delete;

private:
  int m_saved_excepts = 0;
};

/**
 * @brief Owning, move-only handle for a PyObject*.
 *
 * Construction *steals* a reference, which is what the C API returns almost
 * everywhere; PyRef::borrow() increments first, for borrowed references.
 */
class PyRef {
public:
  PyRef() = default;

  /// @brief Take ownership of a new reference.
  explicit PyRef(PyObject *obj) : m_obj(obj) {}

  /// @brief Add a reference to a borrowed object and own that one.
  static PyRef borrow(PyObject *obj) {
    Py_XINCREF(obj);
    return PyRef(obj);
  }

  ~PyRef() { Py_XDECREF(m_obj); }

  PyRef(const PyRef &) = delete;
  PyRef &operator=(const PyRef &) = delete;
  PyRef(PyRef &&other) noexcept : m_obj(other.m_obj) { other.m_obj = nullptr; }
  PyRef &operator=(PyRef &&other) noexcept {
    if (this != &other) {
      Py_XDECREF(m_obj);
      m_obj = other.m_obj;
      other.m_obj = nullptr;
    }
    return *this;
  }

  PyObject *get() const { return m_obj; }
  explicit operator bool() const { return m_obj != nullptr; }

  /// @brief Attribute lookup returning an owning handle (null if absent).
  PyRef attr(const char *name) const {
    if (m_obj == nullptr || PyObject_HasAttrString(m_obj, name) == 0) {
      return PyRef();
    }
    return PyRef(PyObject_GetAttrString(m_obj, name));
  }

private:
  PyObject *m_obj = nullptr;
};

/**
 * @brief Consume the pending Python exception and format it, with traceback.
 *
 * Returns an empty string if no exception is set.  Requires the GIL.
 */
std::string py_take_error();

/**
 * @brief Throw an InferenceError describing the pending Python exception.
 * @param context What we were attempting, e.g. "importing module 'foo'"
 */
[[noreturn]] void py_throw(const std::string &context);

/// @brief Convert a C++ string to a Python str (owning handle).
PyRef py_string(const std::string &s);

/// @brief Read any Python object as a C++ string (empty if that fails).
std::string py_to_string(PyObject *obj);

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_PYTHON_INTERPRETER_HPP
