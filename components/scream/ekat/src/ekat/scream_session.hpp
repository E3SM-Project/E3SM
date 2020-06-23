#ifndef SCREAM_SESSION_HPP
#define SCREAM_SESSION_HPP

namespace scream {

void initialize_scream_session();
void initialize_scream_session(int argc, char **argv);
extern "C" {
void finalize_scream_session();
} // extern "C"
} // namespace scream

#endif // SCREAM_SESSION_HPP
