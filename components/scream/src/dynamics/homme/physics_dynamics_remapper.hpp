#ifndef SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
#define SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP

#include "share/scream_config.hpp"

#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"

#include "share/grid/remap/abstract_remapper.hpp"
#include "share/grid/se_grid.hpp"
#include "share/util/scream_utils.hpp"

#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_assert.hpp"

// Homme includes
#include "Context.hpp"
#include "HommexxEnums.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"
#include "Types.hpp"
#include "mpi/Connectivity.hpp"
#include "mpi/BoundaryExchange.hpp"
#include "mpi/MpiBuffersManager.hpp"

namespace scream
{

// Performs remap from physics to dynamics grids, and viceversa
class PhysicsDynamicsRemapper : public AbstractRemapper
{
public:
  using base_type       = AbstractRemapper;
  using field_type      = typename base_type::field_type;
  using identifier_type = typename base_type::identifier_type;
  using layout_type     = typename base_type::layout_type;
  using grid_type       = typename base_type::grid_type;
  using grid_ptr_type   = typename base_type::grid_ptr_type;

  using device_type     = typename field_type::device_t;

  using pack_type = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using small_pack_type = ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>;

  PhysicsDynamicsRemapper (const grid_ptr_type& phys_grid,
                           const grid_ptr_type& dyn_grid);

  ~PhysicsDynamicsRemapper () = default;

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override;
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override;

  bool compatible_layouts (const layout_type& src,
                           const layout_type& tgt) const override {
    auto tgt_tags = tgt.tags();
    ekat::erase(tgt_tags,FieldTag::TimeLevel);

    return get_layout_type(src.tags())==get_layout_type(tgt_tags);
  }

protected:

  // Getters
  const identifier_type& do_get_src_field_id (const int ifield) const override {
    return m_phys[ifield].get_header().get_identifier();
  }
  const identifier_type& do_get_tgt_field_id (const int ifield) const override {
    return m_dyn[ifield].get_header().get_identifier();
  }
  const field_type& do_get_src_field (const int ifield) const override {
    return m_phys[ifield];
  }
  const field_type& do_get_tgt_field (const int ifield) const override {
    return m_dyn[ifield];
  }

  // Registration methods
  void do_registration_begins () override {}
  void do_register_field (const identifier_type& src, const identifier_type& tgt) override;
  void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) override;
  void do_registration_ends () override;

  void setup_boundary_exchange ();

  std::vector<field_type>   m_phys;
  std::vector<field_type>   m_dyn;

  std::vector<bool>         m_is_state_field;

  grid_ptr_type     m_dyn_grid;
  grid_ptr_type     m_phys_grid;

  int m_num_phys_cols;
  typename grid_type::lid_to_idx_map_type    m_lid2elgp;

  std::shared_ptr<Homme::BoundaryExchange>  m_be[HOMMEXX_NUM_TIME_LEVELS];

  KokkosTypes<DefaultDevice>::view_1d<int>  m_p2d;

  template<typename DataType>
  ::Homme::ExecViewUnmanaged<DataType>
  getHommeView(const Field& f) {
    auto scream_view = f.template get_view<DataType>();
    return ::Homme::ExecViewUnmanaged<DataType>(scream_view.data(),scream_view.layout());
  }

#ifdef HOMMEXX_ENABLE_GPU
public:
  // These structs and function should be morally private, but CUDA complains that
  // they cannot be private/protected
#endif
  void create_p2d_map ();
  struct Pointer {
          Real* ptr = nullptr;
    const Real* cptr = nullptr;
  };

  enum AllocPropType : int {
    PackAlloc      = 0,
    SmallPackAlloc = 1,
    RealAlloc      = 2
  };

  struct Dims {
    int size;
    Kokkos::Array<int,6> dims;
  };

  struct RemapFwdTag {};
  struct RemapBwdTag {};

protected:

  KokkosTypes<DefaultDevice>::view_1d<Pointer>  phys_ptrs;
  KokkosTypes<DefaultDevice>::view_1d<Dims>     phys_dims;
  KokkosTypes<DefaultDevice>::view_1d<Int>      phys_layout;

  KokkosTypes<DefaultDevice>::view_1d<Pointer>  dyn_ptrs;
  KokkosTypes<DefaultDevice>::view_1d<Dims>     dyn_dims;
  KokkosTypes<DefaultDevice>::view_1d<Int>      dyn_layout;

  KokkosTypes<DefaultDevice>::view_1d<bool>     has_parent;
  KokkosTypes<DefaultDevice>::view_1d<Int>      pack_alloc_property;
  KokkosTypes<DefaultDevice>::view_1d<bool>     is_state_field_dev;
  KokkosTypes<DefaultDevice>::view<int>         states_tl_idx;
  // Only meaningful for 3d fields. Set to -1 for 2d fields, in case wrongfully used
  KokkosTypes<DefaultDevice>::view_1d<Int>      num_levels;

  void initialize_device_variables();

  template<typename ScalarT, typename AllocType>
  void compute_view_dims(const AllocType &alloc_prop, const std::vector<int> &field_dims, Dims &view_dims);

  // Remap methods
  void do_remap_fwd () const override;
  void do_remap_bwd () const override;

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

  template <typename MT>
  KOKKOS_FUNCTION
  void local_remap_bwd_2d (const MT& team) const;

  template <typename ScalarT, typename MT>
  KOKKOS_FUNCTION
  void local_remap_bwd_3d (const MT& team) const;

  template<typename ScalarT, int N, typename PtrType>
  KOKKOS_FUNCTION
  Unmanaged<typename KokkosTypes<device_type>::template view_ND<ScalarT,N>>
  reshape (PtrType ptr, const Dims& dims) const {
    using uview_nd = Unmanaged<typename KokkosTypes<device_type>::template view_ND<ScalarT,N>>;

    // Note: if ScalarT is const T, then we would be ok even if pointee_t is non const.
    //       However, likely, this is an error that will arise at run time, due to
    //       getting ptr instead of cptr out of the Pointer struct, and if the Pointer
    //       struct was set up from a read-only field, cptr will be nullptr.
    //       A static assert exposes this issue at compile time. To fix compilation errors,
    //       extract ptr from Pointer with non-const ScalarT, and cptr for const ScalarT.
    using pointee_t = typename std::remove_pointer<PtrType>::type;
    static_assert (std::is_const<ScalarT>::value == std::is_const<pointee_t>::value,
        "Error! Mismatching const qualifiers for input pointer type and ScalarT.\n");

    uview_nd ret_view;

    switch (dims.size) {
      case 1:
        ret_view = uview_nd(reinterpret_cast<ScalarT*>(ptr),
                            dims.dims[0]);
        break;
      case 2:
        ret_view = uview_nd(reinterpret_cast<ScalarT*>(ptr),
                            dims.dims[0],
                            dims.dims[1]);
        break;
      case 3:
        ret_view = uview_nd(reinterpret_cast<ScalarT*>(ptr),
                            dims.dims[0],
                            dims.dims[1],
                            dims.dims[2]);
        break;
      case 4:
        ret_view = uview_nd(reinterpret_cast<ScalarT*>(ptr),
                            dims.dims[0],
                            dims.dims[1],
                            dims.dims[2],
                            dims.dims[3]);
        break;
      case 5:
        ret_view = uview_nd(reinterpret_cast<ScalarT*>(ptr),
                            dims.dims[0],
                            dims.dims[1],
                            dims.dims[2],
                            dims.dims[3],
                            dims.dims[4]);
        break;
      case 6:
        ret_view = uview_nd(reinterpret_cast<ScalarT*>(ptr),
                            dims.dims[0],
                            dims.dims[1],
                            dims.dims[2],
                            dims.dims[3],
                            dims.dims[4],
                            dims.dims[5]);
        break;
      default:
        EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");

    }

    return ret_view;
  }

public:
  template<typename MT>
  KOKKOS_INLINE_FUNCTION
  void operator()(const RemapFwdTag&, const MT &team) const;
  template<typename MT>
  KOKKOS_INLINE_FUNCTION
  void operator()(const RemapBwdTag&, const MT &team) const;

};

} // namespace scream

#endif // SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
