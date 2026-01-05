

# File atm\_field\_data\_provider.cpp

[**File List**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**impl**](dir_6975f7b28201ba7a9e865ff30c48a340.md) **>** [**atm\_field\_data\_provider.cpp**](atm__field__data__provider_8cpp.md)

[Go to the documentation of this file](atm__field__data__provider_8cpp.md)


```C++


#include "atm_field_data_provider.hpp"
#include <algorithm>
#include <regex>

namespace emulator {
namespace impl {

AtmFieldDataProvider::AtmFieldDataProvider(AtmFieldManager &fields,
                                           int ncols_local)
    : m_fields(fields), m_ncols(ncols_local) {}

const std::vector<double> *
AtmFieldDataProvider::get_field(const std::string &name) const {
  // First try direct field access
  std::vector<double> *ptr = m_fields.get_field_ptr(name);
  if (ptr != nullptr) {
    return ptr;
  }

  // Check if it's a stacked field
  if (is_stacked_field(name)) {
    build_stacked_cache(name);
    auto it = m_stacked_cache.find(name);
    if (it != m_stacked_cache.end()) {
      return &it->second;
    }
  }

  return nullptr;
}

std::vector<std::string> AtmFieldDataProvider::get_field_names() const {
  if (m_field_names_cached) {
    return std::vector<std::string>(m_all_field_names.begin(),
                                    m_all_field_names.end());
  }

  // Collect all known field names from field manager
  // First add hardcoded fields by trying common names
  static const char *common_fields[] = {
      // Import fields
      "shf", "cflx", "lhf", "wsx", "wsy", "lwup", "asdir", "aldir", "asdif",
      "aldif", "ts", "sst", "snowhland", "snowhice", "tref", "qref", "u10",
      "u10withgusts", "icefrac", "ocnfrac", "lndfrac",
      // Export fields
      "zbot", "ubot", "vbot", "tbot", "ptem", "shum", "dens", "pbot", "pslv",
      "lwdn", "rainc", "rainl", "snowc", "snowl", "swndr", "swvdr", "swndf",
      "swvdf", "swnet"};

  for (const char *name : common_fields) {
    if (m_fields.get_field_ptr(name) != nullptr) {
      m_all_field_names.insert(name);
    }
  }

  // Add all dynamic fields
  for (const auto &pair : m_fields.dynamic_fields) {
    m_all_field_names.insert(pair.first);
  }

  // Add stacked field basenames
  for (const auto &pair : m_stacked_field_levels) {
    m_all_field_names.insert(pair.first + "_3d");
  }

  m_field_names_cached = true;
  return std::vector<std::string>(m_all_field_names.begin(),
                                  m_all_field_names.end());
}

int AtmFieldDataProvider::get_field_nlevs(const std::string &name) const {
  // Check stacked fields first
  auto it = m_stacked_field_levels.find(name);
  if (it != m_stacked_field_levels.end()) {
    return static_cast<int>(it->second.size());
  }

  // Check if it's a "_3d" stacked field
  if (name.size() > 3 && name.substr(name.size() - 3) == "_3d") {
    std::string basename = name.substr(0, name.size() - 3);
    auto it2 = m_stacked_field_levels.find(basename);
    if (it2 != m_stacked_field_levels.end()) {
      return static_cast<int>(it2->second.size());
    }
  }

  // Single slice from larger field
  const std::vector<double> *ptr = get_field(name);
  if (ptr != nullptr) {
    int total_size = static_cast<int>(ptr->size());
    if (total_size > m_ncols) {
      return total_size / m_ncols;
    }
  }

  return 1; // Default: 2D field
}

void AtmFieldDataProvider::detect_stacked_fields() {
  m_stacked_field_levels.clear();

  // Scan dynamic fields for slice patterns
  for (const auto &pair : m_fields.dynamic_fields) {
    std::string basename;
    int level_idx;

    if (parse_slice_pattern(pair.first, basename, level_idx)) {
      m_stacked_field_levels[basename].push_back(level_idx);
    }
  }

  // Sort levels for each basename
  for (auto &pair : m_stacked_field_levels) {
    std::sort(pair.second.begin(), pair.second.end());
  }

  // Invalidate field names cache
  m_field_names_cached = false;
}

bool AtmFieldDataProvider::is_stacked_field(const std::string &name) const {
  // Direct check
  if (m_stacked_field_levels.find(name) != m_stacked_field_levels.end()) {
    return true;
  }

  // Check "_3d" suffix
  if (name.size() > 3 && name.substr(name.size() - 3) == "_3d") {
    std::string basename = name.substr(0, name.size() - 3);
    return m_stacked_field_levels.find(basename) !=
           m_stacked_field_levels.end();
  }

  return false;
}

bool AtmFieldDataProvider::parse_slice_pattern(const std::string &name,
                                               std::string &basename,
                                               int &level_idx) const {
  // Pattern: basename_N where N is an integer
  static const std::regex pattern(R"((.+)_(\d+)$)");
  std::smatch match;

  if (std::regex_match(name, match, pattern)) {
    basename = match[1].str();
    level_idx = std::stoi(match[2].str());
    return true;
  }

  return false;
}

const std::vector<double> &
AtmFieldDataProvider::get_stacked_field(const std::string &basename) const {
  build_stacked_cache(basename);

  auto it = m_stacked_cache.find(basename);
  if (it != m_stacked_cache.end()) {
    return it->second;
  }

  // Return empty vector if not found
  static const std::vector<double> empty;
  return empty;
}

void AtmFieldDataProvider::build_stacked_cache(
    const std::string &basename) const {
  // Check if already cached
  if (m_stacked_cache.find(basename) != m_stacked_cache.end()) {
    return;
  }

  // Get actual basename (strip _3d if present)
  std::string actual_basename = basename;
  if (basename.size() > 3 && basename.substr(basename.size() - 3) == "_3d") {
    actual_basename = basename.substr(0, basename.size() - 3);
  }

  auto it = m_stacked_field_levels.find(actual_basename);
  if (it == m_stacked_field_levels.end()) {
    return; // Not a stacked field
  }

  const std::vector<int> &levels = it->second;
  int nlevs = static_cast<int>(levels.size());

  // Allocate stacked buffer: [nlevs, ncols] in row-major order
  std::vector<double> stacked(nlevs * m_ncols, 0.0);

  // Copy each level
  for (int lev = 0; lev < nlevs; ++lev) {
    std::string slice_name =
        actual_basename + "_" + std::to_string(levels[lev]);
    std::vector<double> *slice = m_fields.get_field_ptr(slice_name);

    if (slice != nullptr && slice->size() >= static_cast<size_t>(m_ncols)) {
      // Copy to stacked buffer
      std::copy(slice->begin(), slice->begin() + m_ncols,
                stacked.begin() + lev * m_ncols);
    }
  }

  m_stacked_cache[basename] = std::move(stacked);
}

} // namespace impl
} // namespace emulator
```


