#ifndef EAMXX_VERTICAL_REMAPPER_ELEVATED_EMISSIONS_MAM4_HPP
#define EAMXX_VERTICAL_REMAPPER_ELEVATED_EMISSIONS_MAM4_HPP

#include "share/remap/abstract_remapper.hpp"

namespace scream
{

/**
 * @class VerticalRemapperElevatedEmissionsMAM4
 * @brief Vertical remapper specialized for MAM4 elevated emissions.
 *
 * Remaps emission species from altitude-based source levels (km) to model
 * interface height levels (z_mam4_int) using the MAM4xx rebin routine.
 */
class VerticalRemapperElevatedEmissionsMAM4 : public AbstractRemapper
{
public:

  VerticalRemapperElevatedEmissionsMAM4 (const grid_ptr_type& src_grid,
                                         const grid_ptr_type& tgt_grid);

  ~VerticalRemapperElevatedEmissionsMAM4 () = default;

  void set_source_interface_height (const std::string& file_name);
  void set_target_interface_height (const Field& z_iface);

protected:

  Field m_alt_int_src;
  Field m_z_iface_tgt;

  using KT = KokkosTypes<DefaultDevice>;

  void remap_fwd_impl () override;

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void apply_vertical_interpolation (const Field& f_src, const Field& f_tgt) const;
};

} // namespace scream

#endif // EAMXX_VERTICAL_REMAPPER_ELEVATED_EMISSIONS_MAM4_HPP
