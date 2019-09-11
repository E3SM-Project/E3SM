#ifndef SCREAM_GRID_UTILS_HPP
#define SCREAM_GRID_UTILS_HPP

#include <string>

namespace scream
{

enum class GridType {
  Undefined,
  SE_NodeBased,
  SE_CellBased
};

inline std::string e2str (const GridType type) {
  std::string str;
  switch (type) {
    case GridType::Undefined:
      str = "Undefined";
      break;
    case GridType::SE_NodeBased:
      str = "SE Physics";
      break;
    case GridType::SE_CellBased:
      str = "SE Dynamics";
      break;
    default:
      str = "INVALID";
      break;
  }

  return str;
}

} // namespace scream

#endif // SCREAM_GRID_UTILS_HPP
