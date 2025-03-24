#ifndef SCREAM_TIMING_HPP
#define SCREAM_TIMING_HPP

#include <ekat/mpi/ekat_comm.hpp>

#include <string>

namespace scream {

// The following simply wrap GPTL calls. We encourage using
// these (rather than raw GPTL calls), to make SCREAM insensitive
// to any future refactor that might change how we do timing.
void init_gptl (bool& was_already_inited);
void finalize_gptl ();
void start_timer (const std::string& name);
void stop_timer (const std::string& name);

void write_timers_to_file (const ekat::Comm& comm, const std::string& fname);

} // namespace scream

#endif // SCREAM_TIMING_HPP
