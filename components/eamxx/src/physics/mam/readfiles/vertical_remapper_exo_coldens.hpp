#ifndef EAMXX_VERTICAL_REMAPPER_EXO_COLDENS_MAM4_HPP
#define EAMXX_VERTICAL_REMAPPER_EXO_COLDENS_MAM4_HPP

#include "share/remap/abstract_remapper.hpp"
namespace scream
/**
 * @class VerticalRemapperExoColdensMAM4
 * @brief Vertical remapper specialized for MAM4 Exo O3 column density field.
 */
{
class VerticalRemapperExoColdensMAM4 : public AbstractRemapper
{
public:


  VerticalRemapperExoColdensMAM4 (const grid_ptr_type& src_grid,
                        const grid_ptr_type& tgt_grid);

  ~VerticalRemapperExoColdensMAM4 () = default;

  void set_delta_pressure(const std::string& file_name, const Field& pint);

protected:

  int m_exo_ki;
  Real m_exo_delp;

  using KT = KokkosTypes<DefaultDevice>;


  void remap_fwd_impl () override;

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void apply_vertical_interpolation (const Field& f_src, const Field& f_tgt) const;
};



} // namespace scream

#endif // EAMXX_VERTICAL_REMAPPER_EXO_COLDENS_MAM4_HPP
