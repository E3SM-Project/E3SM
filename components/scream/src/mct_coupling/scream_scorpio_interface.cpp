#include "scream_scorpio_interface.hpp"

#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"

#include <string>

using scream::Real;
using scream::Int;
extern "C" {

// Fortran routines to be called from C++
//void grid_write_data_array_c(char* filename, Real* hbuf, char* varname);
  void eam_history_write_c();
  void eam_init_pio_subsystem_c(const int mpicom, const int compid);
//  void register_outfile_c(const char **filename_f);
  void register_outfile_c(const std::string (&filename));
  void register_dimension_c(const char **filename_f, const int length); //const char **shortname_f, const char **longname_f, const int length);
//  void register_dimension_c(const std::string &filename,const int length); //, const std::string (&shortname), const std::string (&longname));
//  void eam_init_pio_1_c(const int mpicom, const int compid);
  void eam_init_pio_2_c();

} // extern C

namespace scream {
namespace scorpio {
/* ----------------------------------------------------------------- */
void convert_string_to_char(std::string str_in,char* str_out) {

  str_in.copy(str_out,str_in.size()+1);
  str_out[str_in.size()] = '\0';

}
/* ----------------------------------------------------------------- */
void eam_init_pio_subsystem(const int mpicom, const int compid) {
  eam_init_pio_subsystem_c(mpicom,compid);
}
/* ----------------------------------------------------------------- */
void register_outfile(const std::string& filename) {

  const char * filename_f = filename.c_str();
  //register_outfile_c(&filename_f);
  register_outfile_c(filename);
}
/* ----------------------------------------------------------------- */
void register_dimension(const std::string &filename, const std::string& shortname, const std::string& longname, const int length) {

  std::cout << "ASD c++ : " << filename << "\n" << std::flush;
  const char * filename_f = filename.c_str();
//  const char * shortname_f = shortname.c_str();
//  const char * longname_f = longname.c_str();
  register_dimension_c(&filename_f, length); //&shortname_f, &longname_f, length);
//  register_dimension_c(filename, length); //, shortname, longname);
}
/* ----------------------------------------------------------------- */
void eam_init_pio_2() {
  eam_init_pio_2_c(); 
}
/* ----------------------------------------------------------------- */
void eam_history_write() {
  eam_history_write_c(); 
}
/* ----------------------------------------------------------------- */
//void grid_write_data_array(std::string filename, Real &hbuf, std::string varname) {
//
//  char filename_f[256];
//  char varname_f[256];
//  convert_string_to_char(filename,filename_f);
//  convert_string_to_char(varname,varname_f);
//
//  // Pass the arguements from C++ to Fortran 
//  void grid_write_data_array_c(filename_f,hbuf,varname_f
//
//};
/* ----------------------------------------------------------------- */

} // namespace scorpio
} // namespace scream
