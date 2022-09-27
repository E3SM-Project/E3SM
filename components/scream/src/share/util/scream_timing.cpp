#include "share/util/scream_timing.hpp"

#include <gptl.h>

namespace scream {

void init_gptl (bool& was_already_inited) {
#ifdef SCREAM_CIME_BUILD
  was_already_inited = true;
#else
  auto ierr = GPTLinitialize();
  was_already_inited = (ierr!=0);
#endif
}
void finalize_gptl () {
  GPTLfinalize();
}

void start_timer (const std::string& name) {
  GPTLstart(name.c_str());
}

void stop_timer (const std::string& name) {
  GPTLstop(name.c_str());
}

void write_timers_to_file (const ekat::Comm& comm, const std::string& fname) {
  GPTLpr_summary_file (comm.mpi_comm(),fname.c_str());
}

} // namespace scream
