#include <eamxx_pysession.hpp>

#include <ekat/ekat_assert.hpp>
#include <Python.h>

namespace scream {

void PySession::initialize () {
  if (num_customers==0) {
    if (not Py_IsInitialized()) {
      Py_Initialize();
      should_finalize = true;
    }
  }
  ++num_customers;
}
void PySession::finalize () {
  --num_customers;
  if (num_customers==0 and should_finalize) {
    Py_Finalize();
  }
  EKAT_REQUIRE_MSG (num_customers>=0,
      "Error! Negative number of customers.\n"
      "  Did you call PySession::finalize() without calling PySession::initialize()?\n");
}

} // namespace scream
