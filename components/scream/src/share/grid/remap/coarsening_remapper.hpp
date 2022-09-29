#ifndef SCREAM_COARSENING_REMAPPER_HPP
#define SCREAM_COARSENING_REMAPPER_HPP

#include "share/grid/remap/abstract_remapper.hpp"

#include "scream_config.h"

#include <mpi.h>

namespace scream
{

// An abstract interface for a remapper

// A remapper is basically a functor, that, given two fields,
// copies the first into the second, or viceversa. The copy must
// account for different layouts and/or different mpi distributions.
// This concept can be extended to remaps that involve interpolation,
// but as of now (07/2019) it is not the intent and projected use
// of this class in the scream framework
class CoarseningRemapper : public AbstractRemapper
{
public:

  CoarseningRemapper (const grid_ptr_type& src_grid,
                      // const grid_ptr_type& tgt_grid,
                      const std::string& map_file);

  ~CoarseningRemapper ();

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override;
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override;

protected:
  static grid_ptr_type create_tgt_grid (const std::string& map_file,
                                        const grid_ptr_type& src_grid);

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

  void do_remap_fwd () const override;

  void do_remap_bwd () const  override {
    EKAT_ERROR_MSG ("CoarseningRemapper only supports fwd remapping.\n");
  }

protected:

  using KT = KokkosTypes<DefaultDevice>;
  using gid_t = AbstractGrid::gid_type;

  template<int N>
  using RPack = ekat::Pack<Real,N>;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;
  template<typename T>
  using view_2d = typename KT::template view_2d<T>;

  void create_ov_tgt_fields ();
  void setup_mpi_data_structures ();

  view_1d<gid_t>::HostMirror
  get_my_triplets_gids (const std::string& map_file) const;

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  template<int N>
  void local_mat_vec (const Field& f_src, const Field& f_tgt) const;

protected:
  ekat::Comm            m_comm;

  static constexpr bool MpiOnDev = SCREAM_MPI_ON_DEVICE;

  // An "overlapped" tgt grid, that is a version of the tgt grid where
  // ranks own all rows that are affected by local dofs in their src grid
  grid_ptr_type         m_ov_tgt_grid;

  // std::vector<MPI_Datatype> m_col_dtypes;
  std::vector<std::map<int,MPI_Datatype>> m_local_dtypes;
  std::vector<std::map<int,MPI_Datatype>> m_remote_dtypes;
  std::vector<int> m_remote_pids;

  std::vector<MPI_Win>  m_tgt_fields_win;

  MPI_Group             m_rma_post_group;
  MPI_Group             m_rma_start_group;

  std::vector<Field>    m_src_fields;
  std::vector<Field>    m_ov_tgt_fields;
  std::vector<Field>    m_tgt_fields;

  // Sparse matrix representation in triplet form
  view_2d<int>    m_row_col_lids;
  view_1d<Real>   m_weights;
};

} // namespace scream

#endif // SCREAM_COARSENING_REMAPPER_HPP
