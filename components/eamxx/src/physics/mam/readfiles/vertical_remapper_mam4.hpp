#ifndef EAMXX_VERTICAL_REMAPPER_MAM4_HPP
#define EAMXX_VERTICAL_REMAPPER_MAM4_HPP

#include "share/grid/remap/vertical_remapper.hpp"

namespace scream
{

/*
* A remapper is used to interpolate fields on a separate vertical grid.
* This remapper has three types of vertical interpolation that
* are required by MAM4XX. We have ported the vertical
* routines from EAM, and in this remapper,
* we invoke the ported routines.
 */

class VerticalRemapperMAM4 : public VerticalRemapper
{
public:

  enum VertRemapType {
    None,
    MAM4_PSRef, // Reconstructs a reference 3d pressure from time-dep PS in input data
    MAM4_ZONAL, // For zonal-type files that are employed in LINOZ.
    MAM4_ELEVATED_EMISSIONS,// for vertical interpolation using altitude instead of pressure.
  };

  VerticalRemapperMAM4 (const grid_ptr_type& src_grid,
                        const grid_ptr_type& tgt_grid,
                        const VertRemapType& vremp_type);

  ~VerticalRemapperMAM4 () = default;

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
  VertRemapType         m_vremap_type;
};

} // namespace scream

#endif // EAMXX_VERTICAL_REMAPPER_MAM4_HPP
