#ifndef EAMXX_STRUCTURED_GRID_UTILS_HPP
#define EAMXX_STRUCTURED_GRID_UTILS_HPP

#include "share/grid/abstract_grid.hpp"

#include <string>
#include <utility>

namespace scream {

/*
 * Utilities to deal with "structured" point grids.
 *
 * A structured grid is a PointGrid whose dofs can be seen as an (ny,nx) array (x
 * being the fastest varying index), so that dof with gid g (0-based, after subtracting
 * the min gid) is at position (iy,ix) = (g/nx, g%nx). Structured grids come in two flavors:
 *  - lat-lon: the lat (lon) of each dof only depends on iy (ix). Here, lat/lon geo data
 *    are 1d coordinate arrays of length ny/nx, and the dims are named "lat"/"lon".
 *  - curvilinear: lat and lon of a dof depend on both iy and ix (e.g., a HRRR grid).
 *    Here, lat/lon geo data are regular fields over the dofs, and the dims are named "ny"/"nx".
 *
 * A grid advertises being structured by storing the geo data "structured_y" and "structured_x".
 * These are 1d fields (tag CMP) of length ny and nx. The names in their layouts are the names
 * of the dims, in the output file, that replace the COL dim. Their values are not meant to be
 * saved to file (for lat-lon grids they are the lat/lon values, which are saved as "lat"/"lon",
 * while for curvilinear grids they are just 1,...,ny and 1,...,nx).
 */

constexpr const char* structured_y_name = "structured_y";
constexpr const char* structured_x_name = "structured_x";

inline bool is_structured_grid (const AbstractGrid& grid)
{
  return grid.has_geometry_data(structured_y_name) and
         grid.has_geometry_data(structured_x_name);
}

// The dims that replace COL in the output file. Dim y is the slowest varying one, while x
// is the fastest. The names are the names of the dims in the file (e.g., "lat"/"lon" or "ny"/"nx").
struct StructuredGridDims {
  std::string y_name;
  std::string x_name;
  int ny;
  int nx;
};

inline StructuredGridDims
get_structured_grid_dims (const AbstractGrid& grid)
{
  auto get_dim = [&](const char* geo_name) {
    const auto& fl = grid.get_geometry_data(geo_name).get_header().get_identifier().get_layout();
    return std::make_pair(fl.names()[0],fl.dim(0));
  };
  const auto y = get_dim(structured_y_name);
  const auto x = get_dim(structured_x_name);
  return {y.first, x.first, y.second, x.second};
}

} // namespace scream

#endif // EAMXX_STRUCTURED_GRID_UTILS_HPP
