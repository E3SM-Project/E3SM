#include "scream_scorpio_interface.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"

using scream::Real;
using scream::Int;
extern "C" {

// Fortran routines to be called from C++
//void grid_write_data_array_c(char* filename, Real* hbuf, char* varname);
  void eam_history_write_c();
  void eam_init_pio_subsystem_c(const int mpicom, const int compid);
  void register_outfile_c(char* filename_f);
  void eam_init_pio_1_c(const int mpicom, const int compid);
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

  std::cout << "ASD Reg outfile in C\n" << std::flush;
  char filename_f[filename.size()+1];
  convert_string_to_char(filename,filename_f);
  std::cout << "ASD Reg outfile in C - Half\n" << std::flush;
  register_outfile_c(filename_f);
  std::cout << "ASD Reg outfile in C - Done\n" << std::flush;
}
/* ----------------------------------------------------------------- */
//void eam_init_pio_1(const int mpicom, const int compid) {
//  eam_init_pio_1_c(mpicom,compid); 
//}
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
