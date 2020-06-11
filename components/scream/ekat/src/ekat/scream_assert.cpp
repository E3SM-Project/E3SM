#include <iostream>

#include <mpi.h>

#include "scream_assert.hpp"
#include "scream_session.hpp"

namespace scream {
namespace error {

void runtime_check(bool cond, const std::string& message, int code) {
  if (!cond) {
    runtime_abort(message,code);
  }
}

void runtime_abort(const std::string& message, int code) {
  std::cerr << message << std::endl << "Exiting..." << std::endl;

  // Finalize scream (e.g., finalize kokkos);
  finalize_scream_session();

  // Check if mpi is active. If so, use MPI_Abort, otherwise, simply std::abort
  int flag;
  MPI_Initialized(&flag);
  if (flag!=0) {
    MPI_Abort(MPI_COMM_WORLD, code);
  } else {
    std::abort();
  }
}

} // namespace scream

void enable_fpes (const int mask) {
  // Make sure we don't throw because one of those exceptions
  // was already set, due to previous calculations
  feclearexcept(mask);

  feenableexcept(mask);
}

void disable_fpes (const int mask) {
  fedisableexcept(mask);
}

int get_enabled_fpes () {
  return fegetexcept();
}

void disable_all_fpes () {
  disable_fpes(FE_ALL_EXCEPT);
}

} // namespace error
