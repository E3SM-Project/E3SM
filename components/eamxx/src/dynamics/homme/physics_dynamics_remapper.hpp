#ifndef SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
#define SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP

#include "share/scream_config.hpp"

#include "share/grid/remap/abstract_remapper.hpp"

#include "ekat/ekat_pack.hpp"

namespace Homme {
class BoundaryExchange;
}

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
  using KT              = KokkosTypes<device_type>;

  template<typename T,int N>
  using view_Nd = typename KT::template view_ND<T,N>;
  template<typename T>
  using view_1d = view_Nd<T,1>;

  using pack_type = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using small_pack_type = ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>;

  PhysicsDynamicsRemapper (const grid_ptr_type& phys_grid,
                           const grid_ptr_type& dyn_grid);

  ~PhysicsDynamicsRemapper () = default;

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override;
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override;

  bool compatible_layouts (const layout_type& src,
                           const layout_type& tgt) const override {
    return get_layout_type(src.tags())==get_layout_type(tgt.tags());
  }

protected:

  // Getters
  const identifier_type& do_get_src_field_id (const int ifield) const override {
    return m_phys_fields[ifield].get_header().get_identifier();
  }
  const identifier_type& do_get_tgt_field_id (const int ifield) const override {
    return m_dyn_fields[ifield].get_header().get_identifier();
  }
  const field_type& do_get_src_field (const int ifield) const override {
    return m_phys_fields[ifield];
  }
  const field_type& do_get_tgt_field (const int ifield) const override {
    return m_dyn_fields[ifield];
  }

  // Registration methods
  void do_registration_begins () override {}
  void do_register_field (const identifier_type& src, const identifier_type& tgt) override;
  void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) override;
  void do_registration_ends () override;

  void setup_boundary_exchange ();

  std::vector<field_type>   m_phys_fields;
  std::vector<field_type>   m_dyn_fields;

  grid_ptr_type     m_dyn_grid;
  grid_ptr_type     m_phys_grid;

  int m_num_phys_cols;
  typename Field::view_dev_t<const int**>  m_lid2elgp;

  std::shared_ptr<Homme::BoundaryExchange>  m_be;

  view_1d<int>  m_p2d;

#ifdef KOKKOS_ENABLE_CUDA
public:
  // These structs and function should be morally private, but CUDA complains that
  // they cannot be private/protected
#endif
  void create_p2d_map ();

  enum AllocPropType : int {
    PackAlloc      = 0,
    SmallPackAlloc = 1,
    RealAlloc      = 2
  };

  // A container structure to hold the physically-shaped views for the fields.
  // Notice that only one of the vNd will be set, while the others will be empty.
  template<typename T>
  struct ViewsContainer {
    view_Nd<T,1>   v1d;
    view_Nd<T,2>   v2d;
    view_Nd<T,3>   v3d;
    view_Nd<T,4>   v4d;
    view_Nd<T,5>   v5d;
  };

  // A views repo contains a view of ViewsContainer (one per field)
  // We have both const and non const, as well as device and host copies.
  struct ViewsRepo {
    using views_t  = view_1d<ViewsContainer<Real>>;
    using cviews_t = view_1d<ViewsContainer<const Real>>;
    using hviews_t = typename views_t::HostMirror;
    using hcviews_t = typename cviews_t::HostMirror;

    views_t   views;
    cviews_t  cviews;

    hviews_t  h_views;
    hcviews_t h_cviews;
  };

protected:

  ViewsRepo   m_phys_repo;
  ViewsRepo   m_dyn_repo;

  // NOTE: one could deduce from OldViewT whether the return type
  //       should have a const value type. But the code is a bit tedious,
  //       so we'll just force to call this as pack_view<const T>(v).
  template<typename NewValueT, typename OldViewT>
  KOKKOS_INLINE_FUNCTION
  view_Nd<NewValueT,OldViewT::Rank> pack_view (const OldViewT& v) const {
    constexpr int N = OldViewT::Rank;
    Kokkos::LayoutRight kl;
    for (int i=0; i<N; ++i) {
      kl.dimension[i] = v.extent(i);
    }
    using OldValueT = typename OldViewT::traits::value_type;
    auto pack_size = sizeof(NewValueT) / sizeof(OldValueT);
    kl.dimension[N-1] /= pack_size;
    auto ptr = reinterpret_cast<NewValueT*>(v.data());

    // Naively, you would return tmp directly. However, if v is a
    // subview of another view along the second extent, you need
    // to fix the stride, since the ctor from ptr+layout would have
    // no knowledge about that.
    view_Nd<NewValueT,N> tmp (ptr,kl);
    auto vm = tmp.impl_map();
    vm.m_impl_offset.m_stride = v.impl_map().stride_0() / pack_size;
    return view_Nd<NewValueT,N> (tmp.impl_track(),vm);
  }

  view_1d<Int> m_layout;
  view_1d<Int> m_pack_alloc_property;

  // Only meaningful for 3d fields. Set to -1 for 2d fields, in case wrongfully used
  view_1d<Int> m_num_levels;

  // List of phys/dyn fields that are subfields of other fields.
  // For each field, we store the last value of subview info, to check if
  // anything changed (only matters for 'dynamic' subfields).
  // Note: here "dynamic" is in the sense explained in Field::subfield
  std::map<int,SubviewInfo> m_subfield_info_dyn;
  std::map<int,SubviewInfo> m_subfield_info_phys;

  void initialize_device_variables();

  bool subfields_info_has_changed (const std::map<int,SubviewInfo>& subfield_info,
                                   const std::vector<field_type>& fields) const;
  void update_subfields_views (const std::map<int,SubviewInfo>& subfield_info,
                               const ViewsRepo& repo,
                               const std::vector<field_type>& fields) const;

  // Remap methods
  void do_remap_fwd () override;
  void do_remap_bwd () override;

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

public:
  struct RemapFwdTag {};
  struct RemapBwdTag {};

  template<typename MT>
  KOKKOS_INLINE_FUNCTION
  void operator()(const RemapFwdTag&, const MT &team) const;
  template<typename MT>
  KOKKOS_INLINE_FUNCTION
  void operator()(const RemapBwdTag&, const MT &team) const;
};

} // namespace scream

#endif // SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
