#ifndef INCLUDE_COMPOSE_CEDR_HPP
#define INCLUDE_COMPOSE_CEDR_HPP

#include "compose.hpp"
#include "compose_kokkos.hpp"

namespace homme {

#ifdef COMPOSE_TIMERS
struct Timer {
  Timer (const std::string& name_) : name("CEDR_" + name_) { GPTLstart(name.c_str()); }
  ~Timer () { Kokkos::fence(); GPTLstop(name.c_str()); }
private:
  const std::string name;
};
#else
struct Timer {
  Timer (const std::string&) {}
};
#endif

} // namespace homme

#endif
