#ifndef EAMXX_HORIZ_INTERP_REMAP_DATA_HPP
#define EAMXX_HORIZ_INTERP_REMAP_DATA_HPP

#include "share/grid/abstract_grid.hpp"

#include <memory>
#include <map>
#include <string>

namespace scream {

enum class InterpType {
  Refine,
  Coarsen
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

// A small struct to hold horiz remap data, which can be shared across multiple horiz remappers
struct HorizRemapperData {
public:
  using KT = KokkosTypes<DefaultDevice>;
  template<typename T>
  using view_1d = typename KT::template view_1d<T>;

  void build (const std::shared_ptr<const AbstractGrid>& grid,
              const std::string& map_file);

  // The CRS matrix data for online interpolation
  view_1d<int>    m_row_offsets;
  view_1d<int>    m_col_lids;
  view_1d<Real>   m_weights;

  // Allow to retrieve how the remap was built
  bool m_coarsening;      // The src grid is finer than the tgt grid
  bool m_built_from_src;  // The input grid was the src grid

  std::shared_ptr<const AbstractGrid> m_input_grid;     // The grid that was passed to build
  std::shared_ptr<AbstractGrid>       m_generated_grid; // The grid we generated during build

  // This will be an overlap version of either the src or tgt grid, depending on whether
  // the remap is fine->coarse or coarse->fine
  std::shared_ptr<AbstractGrid> m_overlap_grid;
private:

  // Read sparse matrix in triplets form. The triplets are split uniformly across ranks
  std::vector<Triplet>  read_mat_triplets (const std::string& map_file);

  // Gather the sparse matrix triplets that this rank needs to perform
  // its local part of the mat-vec product
  std::vector<Triplet>  get_my_triplets (const std::vector<Triplet>& triplets);

  void create_ov_grid (const std::vector<Triplet>& triplets);

  // Not a const ref, since we'll sort the triplets according to
  // how row gids appear in the coarse grid
  void create_crs_matrix_structures (std::vector<Triplet>& triplets);
};

// A small struct to hold horiz remap data, which can be shared across multiple horiz remappers
// NOTE: the client will call the build method, which will read the map file, and create the
//       CRS matrix data for online interpolation.
class HorizRemapperDataRepo {
public:

  static HorizRemapperDataRepo& instance () {
    static HorizRemapperDataRepo repo;
    return repo;
  };

  std::shared_ptr<const HorizRemapperData>
  get_data (const std::shared_ptr<const AbstractGrid>& grid,
            const std::string& map_file);

private:
  HorizRemapperDataRepo () = default;

  std::map<std::string,std::weak_ptr<HorizRemapperData>> m_repo;
};

} // namespace scream

#endif // EAMXX_HORIZ_INTERP_REMAP_DATA_HPP
