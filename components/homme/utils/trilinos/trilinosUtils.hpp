#ifndef TRILINOS_UTILS_HPP
#define TRILINOS_UTILS_HPP

#include <iostream>
#include <string>

// check return values
void check_flag(void* flagvalue, const std::string message, int opt, int myid = -1);

extern "C" {
  // Fortran functions to call MPI finalize or abort in HOMME
  void haltrun();
  void abortrun();
}
#endif
