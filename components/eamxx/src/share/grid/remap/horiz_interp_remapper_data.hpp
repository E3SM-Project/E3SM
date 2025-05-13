#ifndef EAMXX_HORIZ_INTERP_REMAP_DATA_HPP
#define EAMXX_HORIZ_INTERP_REMAP_DATA_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/grid/remap/abstract_remapper.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <memory>
#include <map>
#include <string>

namespace scream {

enum class InterpType {
  Refine,
  Coarsen
};

// A small struct to hold horiz remap data, which can be shared across multiple horiz remappers
// NOTE: the client will call the build method, which will read the map file, and create the
//       CRS matrix data for online interpolation.
struct HorizRemapperData {
  using KT = KokkosTypes<DefaultDevice>;
  template<typename T>
  using view_1d = typename KT::template view_1d<T>;

  // The last argument specifies the base index for gids in the map file
  // For ncremap-type files, all indices are 1-based
  void build (const std::string& map_file,
              const std::shared_ptr<const AbstractGrid>& fine_grid,
              const ekat::Comm& comm,
              const InterpType type);

  // The coarse grid data
  std::shared_ptr<AbstractGrid> coarse_grid;
  std::shared_ptr<AbstractGrid> ov_coarse_grid;
  
  // The CRS matrix data for online interpolation
  view_1d<int>    row_offsets;
  view_1d<int>    col_lids;
  view_1d<Real>   weights;

  int num_customers = 0;
private:
  using gid_type = AbstractGrid::gid_type;

  InterpType                          type;
  std::shared_ptr<const AbstractGrid> fine_grid;
  ekat::Comm                          comm;

  struct Triplet {
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
  get_my_triplets (const std::string& map_file) const;

  void create_coarse_grids (const std::vector<Triplet>& triplets);

  // Not a const ref, since we'll sort the triplets according to
  // how row gids appear in the coarse grid
  void create_crs_matrix_structures (std::vector<Triplet>& triplets);
};

} // namespace scream

#endif // EAMXX_HORIZ_INTERP_REMAP_DATA_HPP
