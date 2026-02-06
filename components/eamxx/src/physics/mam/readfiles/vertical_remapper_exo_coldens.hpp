#ifndef EAMXX_VERTICAL_REMAPPER_EXO_COLDENS_MAM4_HPP
#define EAMXX_VERTICAL_REMAPPER_EXO_COLDENS_MAM4_HPP

#include "share/remap/vertical_remapper.hpp"

namespace scream
{

/*
* A remapper is used to interpolate fields on a separate vertical grid.
* This remapper has three types of vertical interpolation that
* are required by MAM4XX. We have ported the vertical
* routines from EAM, and in this remapper,
* we invoke the ported routines.
 */

class VerticalRemapperExoColdensMAM4 : public VerticalRemapper
{
public:


  VerticalRemapperExoColdensMAM4 (const grid_ptr_type& src_grid,
                        const grid_ptr_type& tgt_grid);

  ~VerticalRemapperExoColdensMAM4 () = default;

  void set_delta_pressure(const std::string& file_name, const Field& pint);

protected:

  int m_exo_ki;
  Real m_exo_delp;


  void remap_fwd_impl () override;

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void apply_vertical_interpolation (const Field& f_src, const Field& f_tgt) const;
};



} // namespace scream

#endif // EAMXX_VERTICAL_REMAPPER_MAM4_HPP
