

# File coupling\_fields.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**coupling\_fields.hpp**](coupling__fields_8hpp.md)

[Go to the documentation of this file](coupling__fields_8hpp.md)


```C++


#ifndef COUPLING_FIELDS_HPP
#define COUPLING_FIELDS_HPP

#include <map>
#include <sstream>
#include <string>

namespace emulator {

class CouplingFieldsBase {
public:
  virtual ~CouplingFieldsBase() = default;

  virtual void initialize(const std::string &export_fields,
                          const std::string &import_fields) {
    parse_field_list(export_fields, export_map, num_exports);
    parse_field_list(import_fields, import_map, num_imports);
  }

  int get_export_index(const std::string &name) const {
    auto it = export_map.find(name);
    return (it != export_map.end()) ? it->second : -1;
  }

  int get_import_index(const std::string &name) const {
    auto it = import_map.find(name);
    return (it != import_map.end()) ? it->second : -1;
  }

  int num_exports = 0; 
  int num_imports = 0; 

protected:
  std::map<std::string, int> export_map; 
  std::map<std::string, int> import_map; 

  void parse_field_list(const std::string &fields,
                        std::map<std::string, int> &field_map, int &count) {
    std::istringstream ss(fields);
    std::string field;
    int idx = 0;
    while (std::getline(ss, field, ':')) {
      if (!field.empty()) {
        field_map[field] = idx++;
      }
    }
    count = idx;
  }
};

} // namespace emulator

#endif // COUPLING_FIELDS_HPP
```


