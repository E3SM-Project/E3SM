

# File diagnostic\_factory.cpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**diagnostics**](dir_8c59288332f7499e80b506ac6e5a3903.md) **>** [**diagnostic\_factory.cpp**](diagnostic__factory_8cpp.md)

[Go to the documentation of this file](diagnostic__factory_8cpp.md)


```C++


#include "diagnostic_factory.hpp"
#include "horiz_avg_diagnostic.hpp"
#include "vert_slice_diagnostic.hpp"
#include <regex>

namespace emulator {

namespace {
// Regex patterns for diagnostic parsing
const std::regex HORIZ_AVG_PATTERN(R"((.+)_(horiz_avg|global_mean)$)");
const std::regex VERT_SLICE_PATTERN(R"((.+)_at_lev(\d+)$)");
} // namespace

bool is_derived_diagnostic(const std::string &name) {
  return std::regex_match(name, HORIZ_AVG_PATTERN) ||
         std::regex_match(name, VERT_SLICE_PATTERN);
}

std::string get_base_field_name(const std::string &diag_name) {
  std::smatch match;

  if (std::regex_match(diag_name, match, HORIZ_AVG_PATTERN)) {
    return match[1].str();
  }

  if (std::regex_match(diag_name, match, VERT_SLICE_PATTERN)) {
    return match[1].str();
  }

  return diag_name; // Not a diagnostic pattern
}

std::unique_ptr<DerivedDiagnostic>
create_diagnostic(const std::string &diag_name,
                  const DiagnosticMetadata &metadata) {
  std::smatch match;

  // Check horizontal average pattern
  if (std::regex_match(diag_name, match, HORIZ_AVG_PATTERN)) {
    std::string field_name = match[1].str();
    return std::make_unique<HorizAvgDiagnostic>(
        field_name, metadata.area_weights, metadata.comm);
  }

  // Check vertical slice pattern
  if (std::regex_match(diag_name, match, VERT_SLICE_PATTERN)) {
    std::string field_name = match[1].str();
    int level_idx = std::stoi(match[2].str());
    return std::make_unique<VertSliceDiagnostic>(field_name, level_idx,
                                                 metadata.nlevs);
  }

  // Not a diagnostic pattern
  return nullptr;
}

} // namespace emulator
```


