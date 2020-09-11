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
 *  - SE_CellBased and SE_NodeBased: the dofs are the gauss points (GP) of a Spectral Element mesh.
 *    The two differ based on how dofs are indexed:
 *      - CellBased: a typical dof is indexed as (ie, igp, jgp). Each rank is guaranteed to own a
 *                   full element. Dofs on element edges/corners are repeated across elements.
 *      - NodeBased: a typical dof is indexed simply with its local id. Dofs are unique across
 *                   elements, which implies ranks may own only some of the GPs inside some elements.
 *  - MeshFree: the dof a simply a range of gids, and there's no assumed connection between any
 *              two dofs. The typical indexing is the same as a NodeBased (local id).
 */

enum class GridType {
  Undefined,
  SE_NodeBased,
  SE_CellBased,
  MeshFree
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
    case GridType::MeshFree:
      str = "Mesh Free";
      break;
    default:
      str = "INVALID";
      break;
  }

  return str;
}

} // namespace scream

#endif // SCREAM_GRID_UTILS_HPP
