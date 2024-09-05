#ifndef SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP
#define SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP

#include "share/grid/remap/abstract_remapper.hpp"
#include "share/grid/remap/horiz_interp_remapper_data.hpp"

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
public:
  HorizInterpRemapperBase (const grid_ptr_type& fine_grid,
                           const std::string& map_file,
                           const InterpType type);

  ~HorizInterpRemapperBase ();

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override;
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override;

  bool compatible_layouts (const layout_type& src,
                           const layout_type& tgt) const override {
    // Same type of layout, and same sizes except for possibly the first one
    // Note: we can't do [src|tgt].size()/[src|tgt].dim(0)), since there may
    // be 0 src/tgt gids on some ranks, which means src/tgt.dim(0)=0.
    using namespace ShortFieldTagsNames;

    // Use congruence, since we don't really care about dimension names, only tags/extents
    return src.clone().strip_dim(COL).congruent(tgt.clone().strip_dim(COL));
  }

protected:

  FieldLayout create_layout (const FieldLayout& fl_in,
                             const grid_ptr_type& grid) const;

  const identifier_type& do_get_src_field_id (const int ifield) const override {
    return m_src_fields[ifield].get_header().get_identifier();
  }
  const identifier_type& do_get_tgt_field_id (const int ifield) const override {
    return m_tgt_fields[ifield].get_header().get_identifier();
  }
  const field_type& do_get_src_field (const int ifield) const override {
    return m_src_fields[ifield];
  }
  const field_type& do_get_tgt_field (const int ifield) const override {
    return m_tgt_fields[ifield];
  }

  void do_registration_begins () override { /* Nothing to do here */ }
  void do_register_field (const identifier_type& src, const identifier_type& tgt) override;
  void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) override;
  void do_registration_ends () override;

  void do_remap_bwd () override {
    EKAT_ERROR_MSG ("HorizInterpRemapperBase only supports fwd remapping.\n");
  }

  using KT = KokkosTypes<DefaultDevice>;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;

  void create_ov_fields ();

  void clean_up ();

  // Derived classes will do different things, depending on m_type and the
  // MPI strategy they use (P2P or RMA)
  virtual void setup_mpi_data_structures () = 0;

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  template<int N>
  void local_mat_vec (const Field& f_src, const Field& f_tgt) const;

  // The fine and coarse grids. Depending on m_type, they could be
  // respectively m_src_grid and m_tgt_grid or viceversa
  // Note: coarse grid is non-const, so that we can add geo data later.
  //       This helps with m_type=Coarsen, which is typically during
  //       model output, so that we can coarsen also geo data.
  grid_ptr_type   m_fine_grid;
  std::shared_ptr<AbstractGrid> m_coarse_grid;

  // An version of the coarse grid where this rank owns all the ids
  // needed for the local mat-vec product. Depending on m_type, this
  // can be on the src or tgt side.
  grid_ptr_type   m_ov_coarse_grid;

  // Source, target, and overlapped intermediate fields
  std::vector<Field>    m_src_fields;
  std::vector<Field>    m_ov_fields;
  std::vector<Field>    m_tgt_fields;

  // ----- Sparse matrix CRS representation ---- //
  view_1d<int>    m_row_offsets;
  view_1d<int>    m_col_lids;
  view_1d<Real>   m_weights;

  // Keep track of this, since we need to tell the remap data repo
  // we are releasing the data for our map file.
  std::string     m_map_file;

  InterpType      m_type;

  ekat::Comm      m_comm;

  static std::map<std::string,HorizRemapperData> s_remapper_data;
};

} // namespace scream

#endif // SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP
