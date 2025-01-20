#ifndef EAMXX_IOP_REMAPPER_HPP
#define EAMXX_IOP_REMAPPER_HPP

#include "share/grid/remap/abstract_remapper.hpp"

namespace scream
{

/*
 *  A remapper for IOP-like needs
 *
 *  This remapper does the following operations:
 *    - extract col closest to given lat/lon coordinates
 *    - broadcast the closest col to ALL the model columns
 */

class IOPRemapper : public AbstractRemapper
{
public:
  IOPRemapper (const grid_ptr_type src_grid,
               const grid_ptr_type tgt_grid,
               const Real lat, const Real lon);

  ~IOPRemapper () = default;

protected:

  struct ClosestColInfo {
    // MPI rank which owns the closest column to the requested lat/lon
    int mpi_rank;
    // Local column index of on rank=mpi_rank (-1 on all other ranks)
    int col_lid = -1;
  };  

// CUDA requires the parent fcn of a KOKKOS_LAMBDA to have public access
#ifdef EAMXX_ENABLE_GPU
public:
#endif
  void setup_closest_col_info (const Real lat, const Real lon);
protected:

  void registration_ends_impl () override;

  void remap_fwd_impl () override;

  std::vector<Field>    m_single_col_fields;

  ClosestColInfo        m_closest_col_info;

  ekat::Comm            m_comm;
};

} // namespace scream

#endif // EAMXX_IOP_REMAPPER_HPP
