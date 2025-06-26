#ifndef EAMXX_PY_SESSION_HPP
#define EAMXX_PY_SESSION_HPP

#include <string>

namespace scream {

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

private:

  PySession () = default;

  // Keep track of how many eamxx places are requesting py support
  int num_customers = 0;

  // If when we call initialize the 1st time (num_customers=0) we see
  // that py is already inited, we ASSUME its init/finalize is EXTERNALLY handled
  bool should_finalize = false;
};

} // namespace scream

#endif // EAMXX_PY_SESSION_HPP
