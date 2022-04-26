#include "share/util/scream_utils.hpp"

#include <sys/resource.h>

namespace scream {

long long get_mem_usage (const MemoryUnits u) {
  struct rusage r_usage;
  getrusage(RUSAGE_SELF,&r_usage);

  long long mem = r_usage.ru_maxrss;

  switch (u) {
    case B  : mem *= 1000;              break;
    case KB :                           break;
    case MB : mem /= 1000;              break;
    case GB : mem /= 1000*1000;         break;
    case GiB: mem /= 1024;              // Fallthrough
    case MiB: mem /= 1024;              // Fallthrough
    case KiB: mem *= 1000; mem /= 1024; break;
    default:
      EKAT_ERROR_MSG ("Invalid choice for memory units: " + std::to_string(u) + "\n");
  }

  return mem;
}

} // namespace scream
