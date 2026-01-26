#ifndef SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP
#define SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP

#include "share/remap/abstract_remapper.hpp"
#include "share/remap/horiz_interp_remapper_data.hpp"

namespace scream
{

/*
 * A base class for (horizontal) interpolation remappers
 *
 * This base class simply implements one method, common to all interpolation
 * remappers, which reads a map file, and grabs the sparse matrix triplets
 * that are needed.
 */

class HorizInterpRemapperBase : public AbstractRemapper
{
protected:
  HorizInterpRemapperBase (const grid_ptr_type& grid,
                           const std::string& map_file);

  void registration_ends_impl () override;

  using KT = KokkosTypes<DefaultDevice>;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;

  void create_ov_fields ();

  void clean_up ();

  // Derived classes may elect to do different things (e.g., depending on interpolation direction
  // or MPI communication strategy)
  virtual void setup_mpi_data_structures () = 0;

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  template<int N>
  void local_mat_vec (const Field& f_src, const Field& f_tgt) const;

  // The fine and coarse grids. Depending on interpolation direction, they could be
  // respectively m_src_grid and m_tgt_grid or viceversa
  grid_ptr_type   m_fine_grid;
  grid_ptr_type   m_coarse_grid;

  // The intermediate grid, where some dofs are replicated across ranks,
  // to allow local calculation of mat-vec product.
  grid_ptr_type   m_overlap_grid;

  // Version of the fields on the intermediate overlapped grid
  std::vector<Field>    m_ov_fields;

  // Store whether each field needs remap. Only fields with COL dim do.
  // NOTE: use int and NOT bool, as vector<bool> is evil
  std::vector<int>    m_needs_remap;

  // We need to keep this (and not just its content) so that the weak_ptr in HorizRemapperDataRepo
  // does not expire. This allows other remappers that need the same data to reuse it rather than
  // have the repo re-create it anew.
  std::shared_ptr<const HorizRemapperData> m_remap_data;
};

} // namespace scream

#endif // SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP
