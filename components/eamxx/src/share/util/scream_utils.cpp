#include "share/util/scream_utils.hpp"

#if defined(SCREAM_ENABLE_STATM)
#include <stdio.h>
#elif defined(SCREAM_ENABLE_GETRUSAGE)
#include <sys/resource.h>
#endif

namespace scream {

long long get_mem_usage (const MemoryUnits u) {

  long long mem = -1;

#if defined(SCREAM_ENABLE_STATM)
  FILE* fh = fopen("/proc/self/statm","r");
  int size;
  EKAT_REQUIRE_MSG(fscanf(fh, "%d %lld", &size, &mem)==2,
      "Error! Something went wrong while probing file '/proc/self/statm'.\n");
  fclose(fh);
#elif defined(SCREAM_ENABLE_GETRUSAGE)
  struct rusage r_usage;
  getrusage(RUSAGE_SELF,&r_usage);
  mem = r_usage.ru_maxrss;
#else
  // Silence compiler warning
  (void) u;
  return mem;
#endif

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
