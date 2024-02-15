#ifndef SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP
#define SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP

#include "share/grid/remap/abstract_remapper.hpp"

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
  enum class InterpType {
    Refine,
    Coarsen
  };

public:
  HorizInterpRemapperBase (const grid_ptr_type& fine_grid,
                           const std::string& map_file,
                           const InterpType type);

  virtual ~HorizInterpRemapperBase () = default;

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override;
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override;

  bool compatible_layouts (const layout_type& src,
                           const layout_type& tgt) const override {
    // Same type of layout, and same sizes except for possibly the first one
    // Note: we can't do [src|tgt].size()/[src|tgt].dim(0)), since there may
    // be 0 src/tgt gids on some ranks, which means src/tgt.dim(0)=0.
    using namespace ShortFieldTagsNames;

    return src.strip_dim(COL)==tgt.strip_dim(COL);
  }

protected:

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

  using gid_type = AbstractGrid::gid_type;
  using KT = KokkosTypes<DefaultDevice>;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;

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

  InterpType      m_type;

  ekat::Comm      m_comm;
};

} // namespace scream

#endif // SCREAM_HORIZ_INTERP_REMAPPER_BASE_HPP
