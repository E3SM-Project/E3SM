#ifndef SCREAM_SCALAR_WRAPPER_HPP
#define SCREAM_SCALAR_WRAPPER_HPP

#include "share/util/eamxx_data_type.hpp"

namespace scream {

// Small struct that allows to pass around a scalar without having to template on the scalar type,
// and expose all implementations in hpp files.
// NOTE: when storing the constant, the type of the constant passed will be stored.
//       When retrieving the constant, we check that the requested template type is
//       such as to NOT cause loss of precision. E.g., if storing 1 (an int), we can
//       later retrieve either as int or Real, but if storing 1.0 (a double), we will
//       get an error if retrieving as an int, even if in this specific case the cast
//       would not cause loss of precision.
struct ScalarWrapper {
  DataType type;     // For type safety checks
  double value;   // Ensure any precision is handled

  ScalarWrapper () { type = DataType::Invalid; }

  template<typename ST>
  ScalarWrapper (const ST x) {
    set(x);
  }

  template<typename ST>
  void set (const ST x) {
    value = x;
    type = get_data_type<ST>();
  }

  template<typename ST>
  ST as () const {
    EKAT_REQUIRE_MSG (not is_narrowing_conversion(type,get_data_type<ST>()),
        "Error! Requested scalar type would cause a loss of precision.\n");

    return static_cast<ST>(value);
  }

  static ScalarWrapper one  () { return ScalarWrapper(1); }
  static ScalarWrapper zero () { return ScalarWrapper(0); }
};

} // namespace scream

#endif // SCREAM_SCALAR_WRAPPER_HPP
