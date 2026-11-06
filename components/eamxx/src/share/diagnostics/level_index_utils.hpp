#ifndef EAMXX_LEVEL_INDEX_UTILS_HPP
#define EAMXX_LEVEL_INDEX_UTILS_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_assert.hpp>

#include <string>

// Helpers shared by PressureLevelIndex/HeightLevelIndex and their consumers
// (FieldAtPressureLevel/FieldAtHeight).

namespace scream
{
namespace level_index_utils
{

// Convert a pressure value, expressed in the given units ('mb', 'hPa', 'Pa'), to Pa
inline Real pressure_in_Pa (const std::string& value, const std::string& units)
{
  EKAT_REQUIRE_MSG (units=="mb" or units=="hPa" or units=="Pa",
      "Error! Invalid pressure units.\n"
      " - input units: " + units + "\n"
      " - valid units: 'mb', 'hPa', 'Pa'\n");
  const Real v = std::stod(value);
  return (units=="Pa") ? v : v*100;
}

} // namespace level_index_utils
} // namespace scream

#endif // EAMXX_LEVEL_INDEX_UTILS_HPP
