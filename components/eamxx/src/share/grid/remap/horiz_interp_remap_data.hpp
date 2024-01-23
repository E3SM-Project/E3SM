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

// A small struct to hold horiz remap data, which can
// be shared across multiple horiz remappers
struct HorizRemapData {
  using KT = KokkosTypes<DefaultDevice>;
  template<typename T>
  using view_1d = typename KT::template view_1d<T>;

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

// A singleton struct, which will hold the horiz remappers data
struct HorizRemapDataRepo {
  static HorizRemapDataRepo& instance () {
    static HorizRemapDataRepo hrdr;
    return hrdr;
  }

  // If data for this map file is NOT present, proceed to build it
  const HorizRemapData&
  get_data (const std::string& map_file,
            const std::shared_ptr<const AbstractGrid>& fine_grid,
            const ekat::Comm& comm,
            const InterpType type);

  void release_data (const std::string& map_file);
private:
  HorizRemapDataRepo () = default;

  std::map<AbstractRemapper*,std::string> customer_to_data;
  std::map<std::string,HorizRemapData>    data_repo;
};

} // namespace scream

#endif // EAMXX_HORIZ_INTERP_REMAP_DATA_HPP
