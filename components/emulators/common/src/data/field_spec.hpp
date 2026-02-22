/**
 * @file field_spec.hpp
 * @brief Field metadata specification for DataView fields.
 */

#ifndef FIELD_SPEC_HPP
#define FIELD_SPEC_HPP

#include <string>

namespace emulator {

/**
 * @brief Metadata describing a single field in a DataView.
 *
 * Carries naming, units, and global size information needed
 * for IO variable definitions and coupler field identification.
 *
 * Field names follow the MCT convention (e.g. "Sa_z", "Faxa_lwdn").
 */
struct FieldSpec {
  std::string name;        ///< Short field name (MCT-style, e.g. "Sa_z")
  std::string units;       ///< Physical units (e.g. "m", "K", "W/m2")
  std::string long_name;   ///< Human-readable description
  int global_size = 0;     ///< Total global grid points for this field

  FieldSpec() = default;

  FieldSpec(std::string name_, std::string units_ = "",
            std::string long_name_ = "", int global_size_ = 0)
      : name(std::move(name_)), units(std::move(units_)),
        long_name(std::move(long_name_)), global_size(global_size_) {}
};

} // namespace emulator

#endif // FIELD_SPEC_HPP
