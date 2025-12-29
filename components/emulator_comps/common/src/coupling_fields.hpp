/**
 * @file coupling_fields.hpp
 * @brief Base class for managing coupling field indices between components.
 *
 * Provides utilities for parsing and looking up field indices used in
 * data exchange between coupled model components via MCT (Model Coupling
 * Toolkit).
 */

#ifndef COUPLING_FIELDS_HPP
#define COUPLING_FIELDS_HPP

#include <map>
#include <sstream>
#include <string>

namespace emulator {

/**
 * @brief Base utility class for coupling field index management.
 *
 * Parses colon-separated field lists (as provided by MCT) and maintains
 * a mapping from field names to their indices. Derived emulator classes
 * (EmulatorAtm, EmulatorOcn, etc.) use this to look up component-specific
 * field indices for import/export operations.
 *
 * ## Usage Example
 * ```cpp
 * CouplingFieldsBase fields;
 * fields.initialize("Sa_z:Sa_u:Sa_v", "Sx_t:Sx_avsdr");
 * int z_idx = fields.get_export_index("Sa_z");  // Returns 0
 * int t_idx = fields.get_import_index("Sx_t");  // Returns 0
 * ```
 *
 * @see EmulatorComp for how this is used in component implementations
 */
class CouplingFieldsBase {
public:
  virtual ~CouplingFieldsBase() = default;

  /**
   * @brief Initialize field indices from colon-separated field lists.
   *
   * Parses the export and import field strings and builds internal
   * lookup maps. Field strings use the MCT format: "field1:field2:field3".
   *
   * @param export_fields Colon-separated list of export (a2x) field names
   * @param import_fields Colon-separated list of import (x2a) field names
   */
  virtual void initialize(const std::string &export_fields,
                          const std::string &import_fields) {
    parse_field_list(export_fields, export_map, num_exports);
    parse_field_list(import_fields, import_map, num_imports);
  }

  /**
   * @brief Look up the index of an export field by name.
   * @param name Field name to look up
   * @return Field index (0-based), or -1 if not found
   */
  int get_export_index(const std::string &name) const {
    auto it = export_map.find(name);
    return (it != export_map.end()) ? it->second : -1;
  }

  /**
   * @brief Look up the index of an import field by name.
   * @param name Field name to look up
   * @return Field index (0-based), or -1 if not found
   */
  int get_import_index(const std::string &name) const {
    auto it = import_map.find(name);
    return (it != import_map.end()) ? it->second : -1;
  }

  int num_exports = 0; ///< Total number of export fields
  int num_imports = 0; ///< Total number of import fields

protected:
  std::map<std::string, int> export_map; ///< Export field name to index
  std::map<std::string, int> import_map; ///< Import field name to index

  /**
   * @brief Parse a colon-separated field list into a name->index map.
   * @param fields Input string of colon-separated field names
   * @param field_map Output map from field name to index
   * @param count Output total count of fields parsed
   */
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
