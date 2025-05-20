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
