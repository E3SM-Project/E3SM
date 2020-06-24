#ifndef SCREAM_SCORPIO_INTERFACE_HPP
#define SCREAM_SCORPIO_INTERFACE_HPP

#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"

using scream::Real;
using scream::Int;
namespace scream {
namespace scorpio {

  void eam_history_write();
  void eam_init_pio_subsystem(const int mpicom, const int compid);
  void register_outfile(const std::string& filename);
  void register_dimension(const std::string& filename,const std::string& shortname, const std::string& longname, const int length);
//  void eam_init_pio_1(const int mpicom, const int compid);
  void eam_init_pio_2();
//void grid_write_data_array(std::string filename, Real &hbuf, std::string varname);
//void grid_write_data_array(std::string filename, Int  &hbuf, std::string varname);

} // namespace scorpio
} // namespace scream

#endif // define SCREAM_SCORPIO_INTERFACE_HPP 
