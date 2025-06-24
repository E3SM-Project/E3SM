#ifndef EAMXX_VERTICAL_REMAPPER_MAM4_HPP
#define EAMXX_VERTICAL_REMAPPER_MAM4_HPP

#include "share/grid/remap/abstract_remapper.hpp"
#include <mam4xx/mam4.hpp>

namespace scream
{

/*
 * A remapper to interpolate fields on a separate vertical grid
 */

class VerticalRemapperMAM4 : public AbstractRemapper
{
public:

// Use fixed value as mask value
  // VerticalRemapperMAM4 (const grid_ptr_type& src_grid,
  //                   const std::string& map_file);

  VerticalRemapperMAM4 (const grid_ptr_type& src_grid,
                    const grid_ptr_type& tgt_grid);

  ~VerticalRemapperMAM4 () = default;

  void set_source_pressure (const Field& p);
  void set_target_pressure (const Field& p);

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

  ekat::Comm            m_comm;

  // Vertical profile fields, both for source and target
  Field                 m_src_pmid;
  // Field                 m_src_pint;
  Field                 m_tgt_pmid;
  // Field                 m_tgt_pint;

  // bool m_src_int_same_as_mid = false;
  // bool m_tgt_int_same_as_mid = false;
};

} // namespace scream

#endif // EAMXX_VERTICAL_REMAPPER_MAM4_HPP
