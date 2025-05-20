#ifndef EAMXX_PY_SESSION_HPP
#define EAMXX_PY_SESSION_HPP

namespace scream {

struct PySession {

  void initialize ();
  void finalize ();

  static PySession& get () {
    static PySession s;
    return s;
  }

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
