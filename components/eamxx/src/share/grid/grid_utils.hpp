#ifndef SCREAM_GRID_UTILS_HPP
#define SCREAM_GRID_UTILS_HPP

#include <string>

namespace scream
{

/*
 * An enum for the type of grid.
 *
 * Possible choices, and their meaning:
 *
 *  - Undefined: a placeholder to spot uninitialized stuff
 *  - SE: the dofs are the gauss points (GP) of a Spectral Element mesh.
 *  - Point: the dofs are simply a range of gids, and there's no assumed connection between dofs.
 */

enum class GridType {
  Undefined,
  SE,         // Spectral Element
  Point       // Mesh-free set of points
};

inline std::string e2str (const GridType type) {
  std::string str;
  switch (type) {
    case GridType::Undefined:
      str = "Undefined";
      break;
    case GridType::SE:
      str = "SE";
      break;
    case GridType::Point:
      str = "Point";
      break;
    default:
      str = "INVALID";
      break;
  }

  return str;
}

} // namespace scream

#endif // SCREAM_GRID_UTILS_HPP
