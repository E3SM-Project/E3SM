#ifndef EAMXX_VERTICAL_REMAPPER_MAM4_HPP
#define EAMXX_VERTICAL_REMAPPER_MAM4_HPP

// #include "share/grid/remap/abstract_remapper.hpp"
#include "share/grid/remap/vertical_remapper.hpp"
#include <mam4xx/mam4.hpp>

namespace scream
{

/*
 * A remapper to interpolate fields on a separate vertical grid
 */

class VerticalRemapperMAM4 : public VerticalRemapper
{
public:

  enum VertRemapType {
    None,
    MAM4_PSRef, // Reconstructs a reference 3d pressure from time-dep PS in input data
    MAM4_ZONAL,
    MAM4_ELEVATED_EMISSIONS,
  };

  VerticalRemapperMAM4 (const grid_ptr_type& src_grid,
                        const grid_ptr_type& tgt_grid,
                        const VertRemapType& vremp_type);

  ~VerticalRemapperMAM4 () = default;

  void set_source_pressure (const Field& p);
  void set_target_pressure (const Field& p);
  void set_source_pressure (const std::string& file_name);

protected:


  void remap_fwd_impl () override;

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void apply_vertical_interpolation (const Field& f_src, const Field& f_tgt,
                                     const Field& p_src, const Field& p_tgt) const;

protected:

  using KT = KokkosTypes<DefaultDevice>;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;
  template<typename T>
  using view_2d = typename KT::template view_2d<T>;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;
  VertRemapType         m_vremap_type;
};

} // namespace scream

#endif // EAMXX_VERTICAL_REMAPPER_MAM4_HPP
