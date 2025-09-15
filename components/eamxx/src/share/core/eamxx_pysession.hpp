#ifndef EAMXX_PY_SESSION_HPP
#define EAMXX_PY_SESSION_HPP

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

#include <ekat_assert.hpp>
#include <ekat_fpe.hpp>

#include <filesystem>
#include <string>
#include <any>

namespace scream {

/*
 * A singleton class holding a py interpreter context
 *
 * This class is needed b/c we must ensure that a py interpreter session
 * is active when we call any py code. To avoid having to do that every
 * time we call python (directly or indirectly), users can initialize
 * the PySession singleton before doing any py stuff, and finalizing it
 * when they are done. For instance, call PySession::get().initialize()
 * inside an atm proc constructor, and PySession::get().finalize() inside
 * the destructor is a good way to ensure the py interpreter session is
 * active during the lifetime of the atm process.
 * A few notes on the class:
 *   - it is impl-ed via singleton design
 *   - it keeps a ref count of how many times initialize() was called,
 *     so that it can release the py interpreter session only when the
 *     last customer has called finalize()
 *   - it offers convenience functions to add a path to sys.path, so that
 *     we can later import py modules from arbitrary locations
 *   - the safe_import method can be used to temporarily disable FPEs during
 *     the import operation (and re-enable them right after). This is needed
 *     for instance when importing numpy, as FPEs are generated inside the
 *     import operation, and cannot be avoided (but are benign).
 */

class PySession {
public:
  // Avoid accidentally initializing a COPY of the session
  PySession(const PySession&) = delete;

  void initialize ();
  void finalize ();

  bool is_initialized () const {return num_customers>0; }
  static PySession& get () {
    static PySession s;
    return s;
  }

  // These are needed in order to be able to import py modules from non-standard paths
  // Even the current path (which is usually automatically added by the py interpreted)
  // has to be manually added when calling py from C++.
  void add_path (const std::string& path);
  void add_curr_path ();

  // Imports a python module preventing spurious FPEs.All FPEs are DISABLED
  // during the import operation, and re-enabled afterwards.
  pybind11::module safe_import (const std::string& module_name) const;

private:

  PySession () = default;

  // Keep track of how many eamxx places are requesting py support
  int num_customers = 0;

  // Created inside initialize and destroyed inside finalize
  std::any py_guard;
};

// ==================== IMPLEMENTATION ===================== //

inline void PySession::initialize () {
  if (num_customers==0) {
    // Note: if Py interpreter is already inited, we ASSUME someone else
    // is handling the interpreter initialization/finalization
    if (not Py_IsInitialized()) {
      py_guard = std::make_shared<pybind11::scoped_interpreter>();
    }
  }
  ++num_customers;
}

inline void PySession::finalize () {
  EKAT_REQUIRE_MSG (num_customers>0,
      "Error! Invalid number of customers.\n"
      "  Did you call PySession::finalize() without calling PySession::initialize()?\n");

  --num_customers;
  if (num_customers==0) {
    py_guard.reset();
  }
}

inline void PySession::add_path (const std::string& path)
{
  EKAT_REQUIRE_MSG (is_initialized(),
      "Error! Cannot modify python's sys.path, since PySession was not initialized yet.\n");

  try {
    // Import the sys module
    pybind11::module sysModule = pybind11::module::import("sys");

    // Get the sys.path list
    pybind11::list sysPath = sysModule.attr("path");

    // Append the new path to sys.path
    sysPath.append(path);
  } catch (const pybind11::error_already_set& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    throw std::runtime_error("Could not modify sys.path. Aborting.");
  }
}

inline void PySession::add_curr_path ()
{
  auto curr_path = std::filesystem::current_path();
  add_path(curr_path.string());
}

inline pybind11::module PySession::safe_import (const std::string& module_name) const
{
  pybind11::module m;

  // Disable FPEs while loading the module, then immediately re-enable them
  auto fpes = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();
  try {
    m = pybind11::module::import(module_name.c_str());
  } catch (const pybind11::error_already_set& e) {
    std::cout << "[PySession::safe_import] Error! Python module import failed.\n"
                 " - module name: " + module_name + "\n"
                 " - pybind11 error: " + std::string(e.what()) + "\n"
                 "Did you forget to call PySession::add_path to add the module location to sys.path?\n";
    throw e;
  }
  ekat::enable_fpes(fpes);

  EKAT_REQUIRE_MSG (not m.is_none(),
      "Error! Could not import module '" + module_name + "'.\n");

  return m;
}

} // namespace scream

#endif // EAMXX_PY_SESSION_HPP
