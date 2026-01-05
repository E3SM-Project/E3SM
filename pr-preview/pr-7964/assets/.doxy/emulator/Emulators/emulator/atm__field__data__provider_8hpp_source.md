

# File atm\_field\_data\_provider.hpp

[**File List**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**impl**](dir_6975f7b28201ba7a9e865ff30c48a340.md) **>** [**atm\_field\_data\_provider.hpp**](atm__field__data__provider_8hpp.md)

[Go to the documentation of this file](atm__field__data__provider_8hpp.md)


```C++


#ifndef ATM_FIELD_DATA_PROVIDER_HPP
#define ATM_FIELD_DATA_PROVIDER_HPP

#include "../../../common/src/emulator_output_stream.hpp"
#include "atm_field_manager.hpp"
#include <regex>
#include <set>
#include <string>
#include <vector>

namespace emulator {
namespace impl {

class AtmFieldDataProvider : public FieldDataProvider {
public:
  AtmFieldDataProvider(AtmFieldManager &fields, int ncols_local);

  ~AtmFieldDataProvider() override = default;

  const std::vector<double> *get_field(const std::string &name) const override;

  std::vector<std::string> get_field_names() const override;

  int get_ncols() const override { return m_ncols; }

  int get_field_nlevs(const std::string &name) const override;

  void detect_stacked_fields();

  bool is_stacked_field(const std::string &name) const;

  const std::vector<double> &
  get_stacked_field(const std::string &basename) const;

private:
  AtmFieldManager &m_fields;
  int m_ncols;

  // Stacked field detection and caching
  // Maps basename → list of level indices found
  mutable std::map<std::string, std::vector<int>> m_stacked_field_levels;

  // Cache for stacked field data
  mutable std::map<std::string, std::vector<double>> m_stacked_cache;

  // Set of all known field names (for get_field_names)
  mutable std::set<std::string> m_all_field_names;
  mutable bool m_field_names_cached = false;

  bool parse_slice_pattern(const std::string &name, std::string &basename,
                           int &level_idx) const;

  void build_stacked_cache(const std::string &basename) const;
};

} // namespace impl
} // namespace emulator

#endif // ATM_FIELD_DATA_PROVIDER_HPP
```


