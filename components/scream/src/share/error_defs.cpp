#include <iostream>

#include <mpi.h>

#include "error_defs.hpp"
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
  finalize_scream_session();
  MPI_Abort(MPI_COMM_WORLD, code);
}

} // namespace scream 
} // namespace error
