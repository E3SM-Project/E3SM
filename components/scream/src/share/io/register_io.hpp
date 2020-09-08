#ifndef SCREAM_REGISTER_IO_HPP
#define SCREAM_REGISTER_IO_HPP

#include "scream_scorpio_interface.hpp"

namespace scream {

  inline void register_io (const ekat::ParameterList &io_file_params) {

  printf("Registering SCORPIO output files...\n");
  io_file_params.print();
  printf(" -----\n");
  printf("Registering SCORPIO output files... DONE\n");

  }

}

#endif // SCREAM_REGISTER_IO_HPP
