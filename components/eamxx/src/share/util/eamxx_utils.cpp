#include "share/util/eamxx_utils.hpp"
#include <glob.h>

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

std::vector<std::string> filename_glob(const std::vector<std::string>& patterns) {
  std::vector<std::string> all_files;
  for (const auto& pattern : patterns) {
      auto files = globloc(pattern);
      all_files.insert(all_files.end(), files.begin(), files.end());
  }
  return all_files;
}

std::vector<std::string> globloc(const std::string& pattern) {
  // glob struct resides on the stack
  glob_t glob_result;
  memset(&glob_result, 0, sizeof(glob_result));

  int return_value = ::glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
  if (return_value != 0) {
    globfree(&glob_result);
    EKAT_REQUIRE_MSG(return_value == 0, "glob() failed with return value " + std::to_string(return_value));
  }

  std::vector<std::string> filenames;
  for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
      filenames.push_back(std::string(glob_result.gl_pathv[i]));
  }

  globfree(&glob_result);
  return filenames;
}

} // namespace scream
