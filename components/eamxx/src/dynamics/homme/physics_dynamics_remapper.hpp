#ifndef SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
#define SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP

#include "share/remap/abstract_remapper.hpp"
#include "share/remap/field_views_repo.hpp"

#include "share/core/eamxx_config.hpp"

#include <ekat_pack.hpp>

namespace Homme {
class BoundaryExchange;
}

namespace scream
{

// Performs remap from physics to dynamics grids
// Note: the opposite direction (dyn->phys) is not supported here. Use UniqueGridRemapper for that.
class PhysicsDynamicsRemapper : public AbstractRemapper
{
public:
  using grid_ptr_type   = typename AbstractRemapper::grid_ptr_type;

  using device_type     = typename Field::device_t;
  using KT              = KokkosTypes<device_type>;

  template<typename T,int N>
  using view_Nd = typename KT::template view_ND<T,N>;
  template<typename T>
  using view_1d = view_Nd<T,1>;

  using pack_type = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  PhysicsDynamicsRemapper (const grid_ptr_type& phys_grid,
                           const grid_ptr_type& dyn_grid);

  ~PhysicsDynamicsRemapper () = default;

  bool compatible_layouts (const FieldLayout& src,
                           const FieldLayout& tgt) const override {
    return src.type()==tgt.type();
  }

  bool is_valid_tgt_layout (const FieldLayout& layout) const override {
    // We don't want fields with TimeLevel in it. Just subview the fields instead
    return AbstractRemapper::is_valid_tgt_layout(layout) and
           not ekat::contains(layout.tags(),FieldTag::TimeLevel);
  }

protected:

  // Registration methods
  void registration_ends_impl () override;

  void setup_boundary_exchange ();

  grid_ptr_type     m_dyn_grid;
  grid_ptr_type     m_phys_grid;

  int m_num_phys_cols;
  typename Field::view_dev_t<const int**>  m_lid2elgp;

  std::shared_ptr<Homme::BoundaryExchange>  m_be;

  view_1d<int>  m_p2d;

#ifdef KOKKOS_ENABLE_CUDA
public:
  // This enum should be morally private, but CUDA complains that
  // it cannot be private/protected
#endif
  enum AllocPropType : int {
    PackAlloc      = 0,
    RealAlloc      = 2
  };

protected:

  FieldViewsRepo  m_phys_repo;
  FieldViewsRepo  m_dyn_repo;

  view_1d<Int> m_layout;
  view_1d<Int> m_pack_alloc_property;

  // Only meaningful for 3d fields. Set to -1 for 2d fields, in case wrongfully used
  view_1d<Int> m_num_levels;

  void initialize_device_variables();

  // Remap methods
  void remap_fwd_impl () override;

  // phys->dyn requires a halo-exchange. Since not all entries in dyn
  // are overwritten before the exchange, to avoid leftover garbage,
  // we need to set all entries of dyn to zero.
  template <typename ScalarT, typename MT>
  KOKKOS_FUNCTION
  void set_dyn_to_zero(const MT& team) const;

  template <typename MT>
  KOKKOS_FUNCTION
  void local_remap_fwd_2d (const MT& team) const;

  template <typename ScalarT, typename MT>
  KOKKOS_FUNCTION
  void local_remap_fwd_3d (const MT& team) const;

public:
  struct RemapFwdTag {};

  template<typename MT>
  KOKKOS_INLINE_FUNCTION
  void operator()(const RemapFwdTag&, const MT &team) const;
};

} // namespace scream

#endif // SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
