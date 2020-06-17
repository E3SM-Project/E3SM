#ifndef SCREAM_SCORPIO_INTERFACE_HPP
#define SCREAM_SCORPIO_INTERFACE_HPP

#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"

using scream::Real;
using scream::Int;
namespace scream {
namespace scorpio {

  void eam_history_write();
  void eam_init_pio_1(const int mpicom, const int compid);
  void eam_init_pio_2();
//void grid_write_data_array(std::string filename, Real &hbuf, std::string varname);
//void grid_write_data_array(std::string filename, Int  &hbuf, std::string varname);

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
