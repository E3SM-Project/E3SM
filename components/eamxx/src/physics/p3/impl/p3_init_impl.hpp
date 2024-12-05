#ifndef P3_INIT_IMPL_HPP
#define P3_INIT_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

extern "C" {
  void micro_p3_utils_init_c(scream::Real Cpair, scream::Real Rair, scream::Real RH2O, scream::Real RHO_H2O,
                 scream::Real MWH2O, scream::Real MWdry, scream::Real gravit, scream::Real LatVap, scream::Real LatIce,
                 scream::Real CpLiq, scream::Real Tmelt, scream::Real Pi, bool masterproc);
  void p3_init_c(const char** lookup_file_dir, int* info, const bool& write_tables);
}

namespace scream {
namespace p3 {

/*
 * Implementation of p3 init. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */
template <typename S, typename D>
void Functions<S,D>
::p3_init (const bool write_tables, const bool masterproc) {
  static bool is_init = false;
  if (!is_init) {
    using c = scream::physics::Constants<Real>;
    micro_p3_utils_init_c(c::Cpair, c::Rair, c::RH2O, c::RHO_H2O,
                          c::MWH2O, c::MWdry, c::gravit, c::LatVap, c::LatIce,
                          c::CpLiq, c::Tmelt, c::Pi, masterproc);
    static const char* dir = SCREAM_DATA_DIR "/tables";
    Int info;
    p3_init_c(&dir, &info, write_tables);
    EKAT_REQUIRE_MSG(info == 0, "p3_init_c returned info " << info);
    is_init = true;
  }
}

} // namespace p3
} // namespace scream

#endif
