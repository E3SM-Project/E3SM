#ifndef SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP
#define SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP

#include "share/grid/abstract_grid.hpp"

namespace scream
{

/*
 * A base class for (horizontal) interpolation remappers
 *
 * This base class simply implements one method, common to all interpolation
 * remappers, which reads a map file, and grabs the sparse matrix triplets
 * that are needed.
 */

class HorizInterpRemapperBase
{
public:
  virtual ~HorizInterpRemapperBase () = default;

protected:

  enum class OwnedBy {
    Col,
    Row
  };

  struct Triplet {
    using gid_type = AbstractGrid::gid_type;
    // Note: unfortunately, C++17 does not support emplace-ing POD
    //       types as aggregates unless a ctor is declared. C++20 does though.
    Triplet () = default;
    Triplet(const gid_type rr, const gid_type cc, const Real ww)
      : row(rr), col(cc), w(ww) {}
    gid_type row;
    gid_type col;
    Real  w;
  };

  std::vector<Triplet>
  get_my_triplets (const std::string& map_file,
                   const ekat::Comm&  comm,
                   const std::shared_ptr<const AbstractGrid>& grid,
                   const OwnedBy owned_by) const;
};

} // namespace scream

#endif // SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP
