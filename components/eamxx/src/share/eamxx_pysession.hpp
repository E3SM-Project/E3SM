#ifndef EAMXX_PY_SESSION_HPP
#define EAMXX_PY_SESSION_HPP

#include <pybind11/pybind11.h>

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

} // namespace scream

#endif // EAMXX_PY_SESSION_HPP
