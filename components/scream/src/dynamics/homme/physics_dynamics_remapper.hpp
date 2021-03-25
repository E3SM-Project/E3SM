#ifndef SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
#define SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP

#include "share/scream_config.hpp"

#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"

#include "share/grid/remap/abstract_remapper.hpp"

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

using Homme::etoi;

// Performs remap from physics to dynamics grids, and viceversa
template<typename RealType>
class PhysicsDynamicsRemapper : public AbstractRemapper<RealType>
{
public:
  using real_type       = RealType;
  using base_type       = AbstractRemapper<real_type>;
  using field_type      = typename base_type::field_type;
  using identifier_type = typename base_type::identifier_type;
  using layout_type     = typename base_type::layout_type;
  using grid_ptr_type   = typename base_type::grid_ptr_type;

  using pack_type = ekat::Pack<RealType,SCREAM_PACK_SIZE>;
  using small_pack_type = ekat::Pack<RealType,SCREAM_SMALL_PACK_SIZE>;

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

  std::shared_ptr<Homme::BoundaryExchange>  m_be[HOMMEXX_NUM_TIME_LEVELS];

  KokkosTypes<DefaultDevice>::view_1d<int>  m_p2d;

  template<typename DataType>
  ::Homme::ExecViewUnmanaged<DataType>
  getHommeView(const Field<RealType>& f) {
    auto scream_view = f.template get_reshaped_view<DataType>();
    return ::Homme::ExecViewUnmanaged<DataType>(scream_view.data(),scream_view.layout());
  }

public:
  // These functions should be morally privade, but CUDA does not allow extended host-device lambda
  // to have private/protected access within the class

  struct Pointer {
    KOKKOS_FORCEINLINE_FUNCTION
    Real* get() { return ptr; }

    Real* ptr;
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

  KokkosTypes<DefaultDevice>::view_1d<Pointer>  phys_ptrs;
  KokkosTypes<DefaultDevice>::view_1d<Dims>     phys_dims;
  KokkosTypes<DefaultDevice>::view_1d<Int>      phys_layout;

  KokkosTypes<DefaultDevice>::view_1d<Pointer> dyn_ptrs;
  KokkosTypes<DefaultDevice>::view_1d<Dims>    dyn_dims;
  KokkosTypes<DefaultDevice>::view_1d<Int>     dyn_layout;

  KokkosTypes<DefaultDevice>::view_1d<bool>                   has_parent;
  KokkosTypes<DefaultDevice>::view_1d<Int>                    pack_alloc_property;
  KokkosTypes<DefaultDevice>::view_1d<bool>                   is_state_field_dev;
  KokkosTypes<DefaultDevice>::view_1d<Kokkos::pair<int,int>>  time_levels;

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

  template <typename ScalarT, typename MT, typename VT1, typename VT2>
  KOKKOS_FUNCTION
  void local_remap_fwd_2d (const MT& team, const int num_cols, const VT1 lid2elgp, const VT2 p2d) const;

  template <typename ScalarT, typename MT, typename VT1, typename VT2>
  KOKKOS_FUNCTION
  void local_remap_fwd_3d (const MT& team, const int num_cols, const VT1 lid2elgp, const VT2 p2d) const;

  template <typename ScalarT, typename MT, typename VT1, typename VT2>
  KOKKOS_FUNCTION
  void local_remap_bwd_2d (const MT& team, const int num_cols, const VT1 lid2elgp, const VT2 p2d) const;

  template <typename ScalarT, typename MT, typename VT1, typename VT2>
  KOKKOS_FUNCTION
  void local_remap_bwd_3d (const MT& team, const int num_cols, const VT1 lid2elgp, const VT2 p2d) const;

  void create_p2d_map ();
};

// ================= IMPLEMENTATION ================= //

template<typename RealType>
PhysicsDynamicsRemapper<RealType>::
PhysicsDynamicsRemapper (const grid_ptr_type& phys_grid,
                         const grid_ptr_type& dyn_grid)
 : base_type(phys_grid,dyn_grid)
{
  EKAT_REQUIRE_MSG(dyn_grid->type()==GridType::SE,     "Error! Input dynamics grid is not a SE grid.\n");
  EKAT_REQUIRE_MSG(phys_grid->type()==GridType::Point, "Error! Input physics grid is not a Point grid.\n");

  m_dyn_grid  = dyn_grid;
  m_phys_grid = phys_grid;

  // For each phys dofs, we find a corresponding dof in the dyn grid.
  // Notice that such dyn dof may not be unique (if phys dof is on an edge
  // of a SE element), but we don't care. We just need to find a match.
  // The BoundaryExchange already takes care of syncing all shared dyn dofs.
  create_p2d_map ();
}

template<typename RealType>
FieldLayout PhysicsDynamicsRemapper<RealType>::
create_src_layout (const FieldLayout& tgt_layout) const {
  using namespace ShortFieldTagsNames;

  auto tags = tgt_layout.tags();
  auto dims = tgt_layout.dims();

  // Note down the position of the first 'GaussPoint' tag.
  auto it_pos = ekat::find(tags,GP);
  EKAT_REQUIRE_MSG (it_pos!=tags.end(),
      "Error! Did not find the tag 'GaussPoint' in the dynamics layout.\n");
  int pos = std::distance(tags.begin(),it_pos);

  // We replace 'Element' with 'Column'. The number of columns is taken from the src grid.
  tags[0] = COL;
  dims[0] = this->m_src_grid->get_num_local_dofs();

  // Delete GP tags/dims
  ekat::erase(tags,GP);
  ekat::erase(tags,GP);
  dims.erase(dims.begin()+pos);
  dims.erase(dims.begin()+pos);

  // If the tgt layout contains the TimeLevel tag, we slice it off.
  auto it_tl = ekat::find(tags,TL);
  if (it_tl!=tags.end()) {
    pos = std::distance(tags.begin(),it_tl);
    tags.erase(tags.begin()+pos);
    dims.erase(dims.begin()+pos);
  }

  return FieldLayout(tags,dims);
}

template<typename RealType>
FieldLayout PhysicsDynamicsRemapper<RealType>::
create_tgt_layout (const FieldLayout& src_layout) const {
  using namespace ShortFieldTagsNames;

  auto tags = src_layout.tags();
  auto dims = src_layout.dims();

  // Replace COL with EL, and num_cols with num_elems
  tags[0] = EL;
  dims[0] = this->m_tgt_grid->get_num_local_dofs() / (HOMMEXX_NP*HOMMEXX_NP);

  // For position of GP and NP, it's easier to switch between 2d and 3d
  auto lt = get_layout_type(src_layout.tags());
  switch (lt) {
    case LayoutType::Scalar2D:
    case LayoutType::Vector2D:
      // Simple: GP/NP are at the end.
      // Push back GP/NP twice
      tags.push_back(GP);
      tags.push_back(GP);
      dims.push_back(HOMMEXX_NP);
      dims.push_back(HOMMEXX_NP);
      break;
    case LayoutType::Scalar3D:
    case LayoutType::Vector3D:
      {
        // Replace last tag/tim with GP/NP, then push back GP/NP and LEV/nvl
        tags.back() = GP;
        dims.back() = HOMMEXX_NP;

        tags.push_back(GP);
        dims.push_back(HOMMEXX_NP);

        tags.push_back(src_layout.tags().back()); // LEV or ILEV
        dims.push_back(src_layout.dims().back()); // nlev or nlev+1
        break;
      }
    default:
      EKAT_ERROR_MSG("Error! Unrecognized layout type.\n");
  }

  return FieldLayout(tags,dims);
}

template<typename RealType>
void PhysicsDynamicsRemapper<RealType>::
do_register_field (const identifier_type& src, const identifier_type& tgt)
{
  m_phys.push_back(field_type(src));
  m_dyn.push_back(field_type(tgt));
}

template<typename RealType>
void PhysicsDynamicsRemapper<RealType>::
do_bind_field (const int ifield, const field_type& src, const field_type& tgt)
{
  const auto& tgt_layout = tgt.get_header().get_identifier().get_layout();
  const auto& tgt_tags = tgt_layout.tags();
  const auto& tgt_dims = tgt_layout.dims();

  const bool has_time_level  = ekat::contains(tgt_tags,FieldTag::TimeLevel);
  if (has_time_level) {
    const bool valid_tl_dim = (tgt_dims[1]==HOMMEXX_NUM_TIME_LEVELS);
    EKAT_REQUIRE_MSG (valid_tl_dim,
        "Error! Field has the TimeLevel tag, but it does not appear to be 'state'.");
  }

  m_is_state_field.push_back(has_time_level);
  m_phys[ifield] = src;
  m_dyn[ifield] = tgt;

  // If this was the last field to be bound, we can setup the BE and
  // precompute fields needed on device during remapper
  if (this->m_state==RepoState::Closed &&
      (this->m_num_bound_fields+1)==this->m_num_registered_fields) {
    setup_boundary_exchange ();
    initialize_device_variables();
  }
}

template<typename RealType>
void PhysicsDynamicsRemapper<RealType>::
do_registration_ends ()
{
  // If we have all fields allocated, we can setup the BE and
  // precompute fields needed on device during remapper
  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    setup_boundary_exchange ();
    initialize_device_variables();
  }
}

// Helper function to compute the dimensions of
// field views.
template<typename RealType>
template<typename ScalarT, typename AllocType>
void PhysicsDynamicsRemapper<RealType>::
compute_view_dims (const AllocType& alloc_prop, const std::vector<int>& field_dims, Dims& view_dims)
{
  int num_values = alloc_prop.get_alloc_size()/sizeof(ScalarT);
  int N = field_dims.size();

  view_dims.size = N;
  for (int i=0; i<N-1; ++i) {
    view_dims.dims[i] = field_dims[i];
    num_values /= field_dims[i];
  }
  view_dims.dims[N-1] = num_values;
}

template<typename RealType>
void PhysicsDynamicsRemapper<RealType>::
initialize_device_variables()
{
  const int num_fields = m_phys.size();

  phys_ptrs   = decltype(phys_ptrs)      ("phys_ptrs",   num_fields);
  phys_dims   = decltype(phys_dims)      ("phys_dims",   num_fields);
  phys_layout = decltype(phys_layout)    ("phys_layout", num_fields);

  dyn_ptrs   = decltype(dyn_ptrs)   ("dyn_ptrs",   num_fields);
  dyn_dims   = decltype(dyn_dims)   ("dyn_dims",   num_fields);
  dyn_layout = decltype(dyn_layout) ("dyn_layout", num_fields);

  has_parent          = decltype(has_parent)          ("has_parent",               num_fields);
  pack_alloc_property = decltype(pack_alloc_property) ("phys_pack_alloc_property", num_fields);
  is_state_field_dev  = decltype(is_state_field_dev)  ("is_state_field_dev",       num_fields);
  time_levels         = decltype(time_levels)         ("time_levels",              1);

  auto h_phys_ptrs   = Kokkos::create_mirror_view(phys_ptrs);
  auto h_phys_layout = Kokkos::create_mirror_view(phys_layout);
  auto h_phys_dims   = Kokkos::create_mirror_view(phys_dims);

  auto h_dyn_ptrs   = Kokkos::create_mirror_view(dyn_ptrs);
  auto h_dyn_layout = Kokkos::create_mirror_view(dyn_layout);
  auto h_dyn_dims   = Kokkos::create_mirror_view(dyn_dims);

  auto h_has_parent          = Kokkos::create_mirror_view(has_parent);
  auto h_pack_alloc_property = Kokkos::create_mirror_view(pack_alloc_property);
  auto h_is_state_field_dev  = Kokkos::create_mirror_view(is_state_field_dev);
  auto h_time_levels         = Kokkos::create_mirror_view(time_levels);

  const auto& tl = Homme::Context::singleton().get<Homme::TimeLevel>();

  for (int i=0; i<num_fields; ++i) {
    const auto& phys = m_phys[i];
    const auto& dyn  = m_dyn[i];

    const auto& ph = phys.get_header();
    const auto& dh = dyn.get_header();

    const auto& phys_dim = ph.get_identifier().get_layout().dims();
    const auto& dyn_dim  = dh.get_identifier().get_layout().dims();

    // If field has a parent, then its view has been subviewed. We
    // do not want to remmap subviews, only the view of the parent,
    // so we store this and no other info for these fields, then
    // skip these fields during the remap.
    if (ph.get_parent().lock() != nullptr) {
      EKAT_REQUIRE_MSG(dh.get_parent().lock() != nullptr,
                       "Error! If physics field has parent,"
                       "dynamics field must also have parent.");
      h_has_parent(i) = true;
      continue;
    }
    else {
      h_has_parent(i) = false;
    }

    // Store view pointers
    h_phys_ptrs(i).ptr = phys.get_view().data();
    h_dyn_ptrs(i).ptr  = dyn.get_view().data();

    // Store phys layout
    const auto phys_lt = get_layout_type(ph.get_identifier().get_layout().tags());
    h_phys_layout(i) = etoi(phys_lt);

    // Store allocation properties
    const auto& phys_alloc_prop = ph.get_alloc_properties();
    const auto& dyn_alloc_prop  = dh.get_alloc_properties();
    if (phys_alloc_prop.template is_compatible<pack_type>() &&
        dyn_alloc_prop.template  is_compatible<pack_type>()) {
      h_pack_alloc_property(i) = AllocPropType::PackAlloc;

      // Store dimensions of phys/dyn view
      compute_view_dims<pack_type>(phys_alloc_prop, phys_dim, h_phys_dims(i));
      compute_view_dims<pack_type>(dyn_alloc_prop,  dyn_dim,  h_dyn_dims(i));
    } else if (phys_alloc_prop.template is_compatible<small_pack_type>() &&
               dyn_alloc_prop.template  is_compatible<small_pack_type>()) {
      h_pack_alloc_property(i) = AllocPropType::SmallPackAlloc;

      // Store dimensions of phys/dyn view
      compute_view_dims<small_pack_type>(phys_alloc_prop, phys_dim, h_phys_dims(i));
      compute_view_dims<small_pack_type>(dyn_alloc_prop,  dyn_dim,  h_dyn_dims(i));
    } else {
      h_pack_alloc_property(i) = AllocPropType::RealAlloc;

      // Store dimensions of phys/dyn view
      compute_view_dims<Real>(phys_alloc_prop, phys_dim, h_phys_dims(i));
      compute_view_dims<Real>(dyn_alloc_prop,  dyn_dim,  h_dyn_dims(i));
    }

    // Store time levels
    h_is_state_field_dev(i) = m_is_state_field[i];
    h_time_levels(0).first  = tl.n0;
    h_time_levels(0).second = tl.np1;
  }

  Kokkos::deep_copy(phys_ptrs,   h_phys_ptrs);
  Kokkos::deep_copy(phys_layout, h_phys_layout);
  Kokkos::deep_copy(phys_dims,   h_phys_dims);

  Kokkos::deep_copy(dyn_ptrs,   h_dyn_ptrs);
  Kokkos::deep_copy(dyn_layout, h_dyn_layout);
  Kokkos::deep_copy(dyn_dims,   h_dyn_dims);

  Kokkos::deep_copy(has_parent,          h_has_parent);
  Kokkos::deep_copy(time_levels,         h_time_levels);
  Kokkos::deep_copy(pack_alloc_property, h_pack_alloc_property);
  Kokkos::deep_copy(is_state_field_dev,  h_is_state_field_dev);
}

template<typename RealType>
template <typename ScalarT, typename MT>
KOKKOS_FUNCTION
void PhysicsDynamicsRemapper<RealType>::
set_dyn_to_zero(const MT& team) const
{
  const int i = team.league_rank();
  const int itl = time_levels(0).first;
  const auto& dim_d = dyn_dims(i).dims;

  switch (phys_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    {
      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                             dim_d[0],
                                                                             dim_d[1],
                                                                             dim_d[2],
                                                                             dim_d[3]);
        auto v = ekat::subview_1(dyn,itl);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, v.size()), [&](const int& k) {
          int k0   = k%dim_d[0];
          int ktmp = k/dim_d[0];
          int k1   = ktmp%dim_d[2];
          int k2   = ktmp/dim_d[2];
          v(k0, k1, k2) = 0;
        });
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*[NP][NP]>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                                  dim_d[0]);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, dyn.size()), [&](const int& k) {
          int k0   = k%dim_d[0];
          int ktmp = k/dim_d[0];
          int k1   = ktmp%dim_d[1];
          int k2   = ktmp/dim_d[1];
          dyn(k0, k1, k2) = 0;
        });
      }
      break;
    }
    case etoi(LayoutType::Vector2D):
    {
      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                              dim_d[0],
                                                                              dim_d[1],
                                                                              dim_d[2],
                                                                              dim_d[3],
                                                                              dim_d[4]);
        auto v = ekat::subview_1(dyn,itl);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, v.size()), [&](const int& k) {
          int k0   = k%dim_d[0];
          int ktmp = k/dim_d[0];
          int k1   = ktmp%dim_d[2];
          ktmp     = ktmp/dim_d[2];
          int k2   = ktmp%dim_d[3];
          int k3   = ktmp/dim_d[3];
          v(k0, k1, k2, k3) = 0;
        });
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT**[NP][NP]>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                                   dim_d[0],
                                                                                   dim_d[1]);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, dyn.size()), [&](const int& k) {
          int k0   = k%dim_d[0];
          int ktmp = k/dim_d[0];
          int k1   = ktmp%dim_d[1];
          ktmp     = ktmp/dim_d[1];
          int k2   = ktmp%dim_d[2];
          int k3   = ktmp/dim_d[2];
          dyn(k0, k1, k2, k3) = 0;
        });
      }
      break;
    }
    case etoi(LayoutType::Scalar3D):
    {
      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                              dim_d[0],
                                                                              dim_d[1],
                                                                              dim_d[2],
                                                                              dim_d[3],
                                                                              dim_d[4]);
        auto v = ekat::subview_1(dyn,itl);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, v.size()), [&](const int& k) {
          int k0   = k%dim_d[0];
          int ktmp = k/dim_d[0];
          int k1   = ktmp%dim_d[2];
          ktmp     = ktmp/dim_d[2];
          int k2   = ktmp%dim_d[3];
          int k3   = ktmp/dim_d[3];
          v(k0, k1, k2, k3) = 0;
        });

      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                             dim_d[0],
                                                                             dim_d[1],
                                                                             dim_d[2],
                                                                             dim_d[3]);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, dyn.size()), [&](const int& k) {
          int k0   = k%dim_d[0];
          int ktmp = k/dim_d[0];
          int k1   = ktmp%dim_d[1];
          ktmp     = ktmp/dim_d[1];
          int k2   = ktmp%dim_d[2];
          int k3   = ktmp/dim_d[2];
          dyn(k0, k1, k2, k3) = 0;
        });
      }
      break;
    }
    case etoi(LayoutType::Vector3D):
    {
      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT******>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                               dim_d[0],
                                                                               dim_d[1],
                                                                               dim_d[2],
                                                                               dim_d[3],
                                                                               dim_d[4],
                                                                               dim_d[5]);
        auto v = ekat::subview_1(dyn,itl);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, v.size()), [&](const int& k) {
          int k0   = k%dim_d[0];
          int ktmp = k/dim_d[0];
          int k1   = ktmp%dim_d[2];
          ktmp     = ktmp/dim_d[2];
          int k2   = ktmp%dim_d[3];
          ktmp     = ktmp/dim_d[3];
          int k3   = ktmp%dim_d[4];
          int k4   = ktmp/dim_d[4];
          v(k0, k1, k2, k3, k4) = 0;
        });
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                              dim_d[0],
                                                                              dim_d[1],
                                                                              dim_d[2],
                                                                              dim_d[3],
                                                                              dim_d[4]);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, dyn.size()), [&](const int& k) {
          int k0   = k%dim_d[0];
          int ktmp = k/dim_d[0];
          int k1   = ktmp%dim_d[1];
          ktmp     = ktmp/dim_d[1];
          int k2   = ktmp%dim_d[2];
          ktmp     = ktmp/dim_d[2];
          int k3   = ktmp%dim_d[3];
          int k4   = ktmp/dim_d[3];
          dyn(k0, k1, k2, k3, k4) = 0;
        });
      }
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}


template<typename RealType>
void PhysicsDynamicsRemapper<RealType>::
do_remap_fwd() const
{
  using KT = KokkosTypes<DefaultDevice>;

  const auto& tl = Homme::Context::singleton().get<Homme::TimeLevel>();

  const int num_fields = m_phys.size();
  const int num_cols   = m_phys_grid->get_num_local_dofs();
  const auto lid2elgp  = m_dyn_grid->get_lid_to_idx_map();
  const auto p2d       = m_p2d;

  const auto field_loop = KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();

    if (has_parent(i)) return;

    switch (phys_layout(i)) {
      case etoi(LayoutType::Scalar2D):
      case etoi(LayoutType::Vector2D):
      {
        if (pack_alloc_property(i) == AllocPropType::PackAlloc) {
          set_dyn_to_zero<pack_type>(team);
          team.team_barrier();

          local_remap_fwd_2d<pack_type>(team, num_cols, lid2elgp, p2d);
        } else if (pack_alloc_property(i) == AllocPropType::SmallPackAlloc) {
          set_dyn_to_zero<small_pack_type>(team);
          team.team_barrier();

          local_remap_fwd_2d<small_pack_type>(team, num_cols, lid2elgp, p2d);
        } else {
          set_dyn_to_zero<Real>(team);
          team.team_barrier();

          local_remap_fwd_2d<Real>(team, num_cols, lid2elgp, p2d);
        }
        break;
      }
      case etoi(LayoutType::Scalar3D):
      case etoi(LayoutType::Vector3D):
      {
        if (pack_alloc_property(i) == AllocPropType::PackAlloc) {
          set_dyn_to_zero<pack_type>(team);
          team.team_barrier();

          local_remap_fwd_3d<pack_type>(team, num_cols, lid2elgp, p2d);
        } else if (pack_alloc_property(i) == AllocPropType::SmallPackAlloc) {
          set_dyn_to_zero<small_pack_type>(team);
          team.team_barrier();

          local_remap_fwd_3d<small_pack_type>(team, num_cols, lid2elgp, p2d);
        } else {
          set_dyn_to_zero<Real>(team);
          team.team_barrier();

          local_remap_fwd_3d<Real>(team, num_cols, lid2elgp, p2d);
        }
        break;
      }
      default:
        EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
    }
  };

  const auto concurrency = KT::ExeSpace::concurrency();
#ifdef KOKKOS_ENABLE_CUDA
  const int team_size = std::min(1024, std::min(128*num_cols,32*(concurrency/num_fields+31)/32));
#else
  const int team_size = (concurrency<num_fields ? 1 : concurrency/num_fields);
#endif
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_team_policy_force_team_size(num_fields, team_size);
  Kokkos::parallel_for(policy, field_loop);
  Kokkos::fence();

  // Exchange only the current time levels
  m_be[tl.n0]->exchange();
}

template<typename RealType>
void PhysicsDynamicsRemapper<RealType>::
do_remap_bwd() const
{
  using KT = KokkosTypes<DefaultDevice>;

  const int num_fields = m_phys.size();
  const int num_cols   = m_phys_grid->get_num_local_dofs();
  const auto lid2elgp  = m_dyn_grid->get_lid_to_idx_map();
  const auto p2d       = m_p2d;

  const auto field_loop = KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();

    if (has_parent(i)) return;

    switch (phys_layout(i)) {
      case etoi(LayoutType::Scalar2D):
      case etoi(LayoutType::Vector2D):
      {
        if (pack_alloc_property(i) == AllocPropType::PackAlloc) {
          local_remap_bwd_2d<pack_type>(team, num_cols, lid2elgp, p2d);
        } else if (pack_alloc_property(i) == AllocPropType::SmallPackAlloc) {
          local_remap_bwd_2d<small_pack_type>(team, num_cols, lid2elgp, p2d);
        } else {
          local_remap_bwd_2d<Real>(team, num_cols, lid2elgp, p2d);
        }
        break;
      }
      case etoi(LayoutType::Scalar3D):
      case etoi(LayoutType::Vector3D):
      {
        if (pack_alloc_property(i) == AllocPropType::PackAlloc) {
          local_remap_bwd_3d<pack_type>(team, num_cols, lid2elgp, p2d);
        } else if (pack_alloc_property(i) == AllocPropType::SmallPackAlloc) {
          local_remap_bwd_3d<small_pack_type>(team, num_cols, lid2elgp, p2d);
        } else {
          local_remap_bwd_3d<Real>(team, num_cols, lid2elgp, p2d);
        }
        break;
      }
      default:
        EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
    }
  };

  const auto concurrency = KT::ExeSpace::concurrency();
#ifdef KOKKOS_ENABLE_CUDA
  const int team_size = std::min(1024, std::min(128*num_cols,32*(concurrency/nfields+31)/32));
#else
  const int team_size = (concurrency<num_fields ? 1 : concurrency/num_fields);
#endif
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_team_policy_force_team_size(num_fields, team_size);
  Kokkos::parallel_for(policy, field_loop);
  Kokkos::fence();
}

template<typename RealType>
void PhysicsDynamicsRemapper<RealType>::
setup_boundary_exchange () {
  // TODO: should we check that the BE was not already setup? I don't see this happening, but even if it does,
  //       there should be no side effect if we re-build it. We waste some time, sure, but should yield a correct result.

  auto& c = Homme::Context::singleton();

  using Scalar = Homme::Scalar;

  int num_2d = 0;
  int num_3d_mid = 0;
  int num_3d_int = 0;
  const int num_fields = m_dyn.size();
  for (int i=0; i<num_fields; ++i) {
    const auto& layout = m_dyn[i].get_header().get_identifier().get_layout();
    const auto lt = get_layout_type(layout.tags());
    switch (lt) {
      case LayoutType::Scalar2D:
        ++num_2d;
        break;
      case LayoutType::Vector2D:
        if (m_is_state_field[i]) {
          // A scalar state: we only exchange the timelevel we remapped
          num_2d += 1;
        } else {
          // A vector field (not a state): remap all components
          num_2d += layout.dim(1);
        }
        break;
      case LayoutType::Tensor2D:
        if (m_is_state_field[i]) {
          // A vector state: we only exchange the timelevel we remapped
          num_2d += layout.dim(2);
        } else {
          // A tensor field (not a state): remap all components
          num_2d += layout.dim(1)*layout.dim(2);
        }
        break;
      case LayoutType::Scalar3D:
        if (layout.dims().back()==HOMMEXX_NUM_PHYSICAL_LEV) {
          ++num_3d_mid;
        } else if (layout.dims().back()==HOMMEXX_NUM_INTERFACE_LEV) {
          ++num_3d_int;
        } else {
          EKAT_ERROR_MSG ("Error! Unexpected vertical level extent.\n");
        }
        break;
      case LayoutType::Vector3D:
        if (m_is_state_field[i]) {
          // A scalar state: we only exchange the timelevel we remapped
          if (layout.dims().back()==HOMMEXX_NUM_PHYSICAL_LEV) {
            ++num_3d_mid;
          } else if (layout.dims().back()==HOMMEXX_NUM_INTERFACE_LEV) {
            ++num_3d_int;
          } else {
            EKAT_ERROR_MSG ("Error! Unexpected vertical level extent.\n");
          }
        } else {
          // A vector field (not a state): remap all components
          if (layout.dims().back()==HOMMEXX_NUM_PHYSICAL_LEV) {
            num_3d_mid += layout.dim(1);
          } else if (layout.dims().back()==HOMMEXX_NUM_INTERFACE_LEV) {
            num_3d_int += layout.dim(1);
          } else {
            EKAT_ERROR_MSG ("Error! Unexpected vertical level extent.\n");
          }
        }
        break;
      case LayoutType::Tensor3D:
        if (m_is_state_field[i]) {
          auto num_slices = layout.dim(2);
          // A vector state: we only exchange the timelevel we remapped
          if (layout.dims().back()==HOMMEXX_NUM_PHYSICAL_LEV) {
            num_3d_mid += num_slices;
          } else if (layout.dims().back()==HOMMEXX_NUM_INTERFACE_LEV) {
            num_3d_int += num_slices;
          } else {
            EKAT_ERROR_MSG ("Error! Unexpected vertical level extent.\n");
          }
        } else {
          // A tensor field (not a state): remap all components
          if (layout.dims().back()==HOMMEXX_NUM_PHYSICAL_LEV) {
            num_3d_mid += layout.dim(1)*layout.dim(2);
          } else if (layout.dims().back()==HOMMEXX_NUM_INTERFACE_LEV) {
            num_3d_int += layout.dim(1)*layout.dim(2);
          } else {
            EKAT_ERROR_MSG ("Error! Unexpected vertical level extent.\n");
          }
        }
        break;
    default:
      ekat::error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
    }
  }

  // Make sure stuff is created in the context first
  c.create_if_not_there<Homme::MpiBuffersManagerMap>();

  auto bm   = c.get<Homme::MpiBuffersManagerMap>()[Homme::MPI_EXCHANGE];
  auto conn = c.get_ptr<Homme::Connectivity>();

  EKAT_REQUIRE_MSG (bm, "Error! Homme's MpiBuffersManager shared pointer is null.\n");
  EKAT_REQUIRE_MSG (conn, "Error! Homme's Connectivity shared pointer is null.\n");

  constexpr int NLEV = HOMMEXX_NUM_LEV;
  constexpr int NINT = HOMMEXX_NUM_LEV_P;
  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  for (int it=0; it<HOMMEXX_NUM_TIME_LEVELS; ++it) {
    m_be[it]= std::make_shared<Homme::BoundaryExchange>(conn,bm);

    auto be = m_be[it];
    be->set_num_fields(0,num_2d,num_3d_mid,num_3d_int);

    // If some fields are already bound, set them in the bd exchange
    for (int i=0; i<num_fields; ++i) {
      const auto& layout = m_dyn[i].get_header().get_identifier().get_layout();
      const auto& dims = layout.dims();
      const auto lt = get_layout_type(layout.tags());
      switch (lt) {
        case LayoutType::Scalar2D:
          be->register_field(getHommeView<Real*[NP][NP]>(m_dyn[i]));
          break;
        case LayoutType::Vector2D:
          if (m_is_state_field[i]) {
            // A scalar state: exchange one time level only
            be->register_field(getHommeView<Real**[NP][NP]>(m_dyn[i]),1,it);
          } else {
            // A 2d vector (not a state): exchange all slices
            be->register_field(getHommeView<Real**[NP][NP]>(m_dyn[i]),dims[1],0);
          }
          break;
        case LayoutType::Tensor2D:
          if (m_is_state_field[i]) {
            // A vector state: exchange one time level only
            be->register_field(getHommeView<Real***[NP][NP]>(m_dyn[i]),it,dims[2],0);
          } else {
            // A 2d tensor (not a state): exchange all slices
            for (int idim=0; idim<dims[1]; ++idim) {
              // Homme::BoundaryExchange only exchange one slice of the outer dim at a time,
              // so loop on the outer dim and register each slice individually.
              be->register_field(getHommeView<Real***[NP][NP]>(m_dyn[i]),idim,dims[2],0);
            }
          }
          break;
        case LayoutType::Scalar3D:
          if (dims.back()==HOMMEXX_NUM_PHYSICAL_LEV) {
            be->register_field(getHommeView<Scalar*[NP][NP][NLEV]>(m_dyn[i]));
          } else {
            be->register_field(getHommeView<Scalar*[NP][NP][NINT]>(m_dyn[i]));
          }
          break;
        case LayoutType::Vector3D:
          if (m_is_state_field[i]) {
            // A state: exchange one time level only
            if (dims.back()==HOMMEXX_NUM_PHYSICAL_LEV) {
              be->register_field(getHommeView<Scalar*[NTL][NP][NP][NLEV]>(m_dyn[i]),1,it);
            } else {
              be->register_field(getHommeView<Scalar*[NTL][NP][NP][NINT]>(m_dyn[i]),1,it);
            }
          } else {
            // Not a state: exchange all slices
            if (dims.back()==HOMMEXX_NUM_PHYSICAL_LEV) {
              be->register_field(getHommeView<Scalar**[NP][NP][NLEV]>(m_dyn[i]),dims[1],0);
            } else {
              be->register_field(getHommeView<Scalar**[NP][NP][NINT]>(m_dyn[i]),dims[1],0);
            }
          }
          break;
        case LayoutType::Tensor3D:
          if (m_is_state_field[i]) {
            // This must either be v.
            auto num_slices = layout.dim(2);
            auto slice = it;
            be->register_field(getHommeView<Scalar***[NP][NP][NLEV]>(m_dyn[i]),slice,num_slices,0);
          } else {
            EKAT_ERROR_MSG ("Error! We were not expected a rank-6 fields in homme that is not a state.\n");
          }
          break;
      default:
        ekat::error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
      }
    }
    be->registration_completed();
  }
}

template<typename RealType>
template <typename ScalarT, typename MT, typename VT1, typename VT2>
KOKKOS_FUNCTION
void PhysicsDynamicsRemapper<RealType>::
local_remap_fwd_2d (const MT& team, const int num_cols, const VT1 lid2elgp, const VT2 p2d) const
{
  const int i = team.league_rank();

  const auto& dim_p = phys_dims(i).dims;
  const auto& dim_d = dyn_dims(i).dims;

  switch (phys_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    {
      auto phys = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*>> (reinterpret_cast<ScalarT*>(phys_ptrs(i).get()),
                                                                         dim_p[0]);

      const auto tr = Kokkos::TeamThreadRange(team, num_cols);
      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT**[NP][NP]>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                                   dim_d[0],
                                                                                   dim_d[1]);

        const auto f = [&] (const int icol) {
          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          dyn(elgp[0],time_levels(0).first,elgp[1],elgp[2]) = phys(icol);
        };
        Kokkos::parallel_for(tr, f);
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*[NP][NP]>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                                  dim_d[0]);

        const auto tr = Kokkos::TeamThreadRange(team, num_cols);
        const auto f = [&] (const int icol) {
          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          dyn(elgp[0],elgp[1],elgp[2]) = phys(icol);
        };
        Kokkos::parallel_for(tr, f);
      }
      break;
    }
    case etoi(LayoutType::Vector2D):
    {
      auto phys = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT**>> (reinterpret_cast<ScalarT*>(phys_ptrs(i).get()),
                                                                          dim_p[0],
                                                                          dim_p[1]);

      const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_p[1]);
      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT***[NP][NP]>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                                    dim_d[0],
                                                                                    dim_d[1],
                                                                                    dim_d[2]);

        const auto f = [&] (const int idx) {
          const int icol = idx/dim_p[1];
          const int idim = idx%dim_p[1];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          dyn(elgp[0],time_levels(0).first,idim,elgp[1],elgp[2]) = phys(icol,idim);
        };
        Kokkos::parallel_for(tr, f);
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT**[NP][NP]>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                                   dim_d[0],
                                                                                   dim_d[1]);

        const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_p[1]);
        const auto f = [&] (const int idx) {
          const int icol = idx/dim_p[1];
          const int idim = idx%dim_p[1];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          dyn(elgp[0],idim,elgp[1],elgp[2]) = phys(icol,idim);
        };
        Kokkos::parallel_for(tr, f);
      }
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template<typename RealType>
template <typename ScalarT, typename MT, typename VT1, typename VT2>
KOKKOS_FUNCTION
void PhysicsDynamicsRemapper<RealType>::
local_remap_fwd_3d (const MT& team, const int num_cols, const VT1 lid2elgp, const VT2 p2d) const
{
  const int i = team.league_rank();
  const auto& dim_p = phys_dims(i).dims;
  const auto& dim_d = dyn_dims(i).dims;

  switch (phys_layout(i)) {
    case etoi(LayoutType::Scalar3D):
    {
      auto phys = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT**>> (reinterpret_cast<ScalarT*>(phys_ptrs(i).get()),
                                                                          dim_p[0],
                                                                          dim_p[1]);

      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                              dim_d[0],
                                                                              dim_d[1],
                                                                              dim_d[2],
                                                                              dim_d[3],
                                                                              dim_d[4]);

        const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_d[4]);
        const auto f = [&] (const int idx) {
          const int icol = idx/dim_d[4];
          const int ilev = idx%dim_d[4];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          dyn(elgp[0],time_levels(0).first,elgp[1],elgp[2],ilev) = phys(icol,ilev);
        };
        Kokkos::parallel_for(tr, f);
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                             dim_d[0],
                                                                             dim_d[1],
                                                                             dim_d[2],
                                                                             dim_d[3]);

        const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_d[3]);
        const auto f = [&] (const int idx) {
          const int icol = idx/dim_d[3];
          const int ilev = idx%dim_d[3];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          dyn(elgp[0],elgp[1],elgp[2],ilev) = phys(icol,ilev);
        };
        Kokkos::parallel_for(tr, f);
      }
      break;
    }
    case etoi(LayoutType::Vector3D):
    {
      auto phys = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT***>> (reinterpret_cast<ScalarT*>(phys_ptrs(i).get()),
                                                                           dim_p[0],
                                                                           dim_p[1],
                                                                           dim_p[2]);

      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT******>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                               dim_d[0],
                                                                               dim_d[1],
                                                                               dim_d[2],
                                                                               dim_d[3],
                                                                               dim_d[4],
                                                                               dim_d[5]);

        const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_p[1]*dim_d[5]);
        const auto f = [&] (const int idx) {
          const int icol =  idx/(dim_p[1]*dim_d[5]);
          const int idim = (idx/dim_d[5])%dim_p[1];
          const int ilev =  idx%dim_d[5];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          dyn(elgp[0],time_levels(0).first,idim,elgp[1],elgp[2],ilev) = phys(icol,idim,ilev);
        };
        Kokkos::parallel_for(tr, f);
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                              dim_d[0],
                                                                              dim_d[1],
                                                                              dim_d[2],
                                                                              dim_d[3],
                                                                              dim_d[4]);

        const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_p[1]*dim_d[4]);
        const auto f = [&] (const int idx) {
          const int icol =  idx/(dim_p[1]*dim_d[4]);
          const int idim = (idx/dim_d[4])%dim_p[1];
          const int ilev =  idx%dim_d[4];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          dyn(elgp[0],idim,elgp[1],elgp[2],ilev) = phys(icol,idim,ilev);
        };
        Kokkos::parallel_for(tr, f);
      }
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template<typename RealType>
template <typename ScalarT, typename MT, typename VT1, typename VT2>
KOKKOS_FUNCTION
void PhysicsDynamicsRemapper<RealType>::
local_remap_bwd_2d (const MT& team, const int num_cols, const VT1 lid2elgp, const VT2 p2d) const
{
  const int i = team.league_rank();
  const auto& dim_p = phys_dims(i).dims;
  const auto& dim_d = dyn_dims(i).dims;

  switch (phys_dims(i).size) {
    case 1:
    {
      auto phys = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*>> (reinterpret_cast<ScalarT*>(phys_ptrs(i).get()),
                                                                         dim_p[0]);

      const auto tr = Kokkos::TeamThreadRange(team, num_cols);
      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT**[NP][NP]>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                                   dim_d[0],
                                                                                   dim_d[1]);

        const auto f = [&] (const int icol) {
          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          phys(icol) = dyn(elgp[0],time_levels(0).second,elgp[1],elgp[2]);
        };
        Kokkos::parallel_for(tr, f);
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*[NP][NP]>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                                  dim_d[0]);

        const auto f = [&] (const int icol) {
          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          phys(icol) = dyn(elgp[0],elgp[1],elgp[2]);
        };
        Kokkos::parallel_for(tr, f);
      }
      break;
    }
    case 2:
    {
      auto phys = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT**>> (reinterpret_cast<ScalarT*>(phys_ptrs(i).get()),
                                                                          dim_p[0],
                                                                          dim_p[1]);

      const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_p[1]);
      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT***[NP][NP]>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                                    dim_d[0],
                                                                                    dim_d[1],
                                                                                    dim_d[2]);

        const auto f = [&] (const int idx) {
          const int icol = idx/dim_p[1];
          const int idim = idx%dim_p[1];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          phys(icol,idim) = dyn(elgp[0],time_levels(0).second,idim,elgp[1],elgp[2]);
        };
        Kokkos::parallel_for(tr, f);
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT**[NP][NP]>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                                   dim_d[0],
                                                                                   dim_d[1]);

        const auto f = [&] (const int idx) {
          const int icol = idx/dim_p[1];
          const int idim = idx%dim_p[1];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          phys(icol,idim) = dyn(elgp[0],idim,elgp[1],elgp[2]);
        };
        Kokkos::parallel_for(tr, f);
      }
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template<typename RealType>
template <typename ScalarT, typename MT, typename VT1, typename VT2>
KOKKOS_FUNCTION
void PhysicsDynamicsRemapper<RealType>::
local_remap_bwd_3d (const MT& team, const int num_cols, const VT1 lid2elgp, const VT2 p2d) const
{
  const int i = team.league_rank();
  const auto& dim_p = phys_dims(i).dims;
  const auto& dim_d = dyn_dims(i).dims;

  switch (phys_dims(i).size) {
    case 2:
    {
      auto phys = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT**>> (reinterpret_cast<ScalarT*>(phys_ptrs(i).get()),
                                                                          dim_p[0],
                                                                          dim_p[1]);

      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                              dim_d[0],
                                                                              dim_d[1],
                                                                              dim_d[2],
                                                                              dim_d[3],
                                                                              dim_d[4]);

        const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_p[1]);
        const auto f = [&] (const int idx) {
          const int icol = idx/dim_p[1];
          const int ilev = idx%dim_p[1];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          phys(icol,ilev) = dyn(elgp[0],time_levels(0).second,elgp[1],elgp[2],ilev);
        };
        Kokkos::parallel_for(tr, f);
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                             dim_d[0],
                                                                             dim_d[1],
                                                                             dim_d[2],
                                                                             dim_d[3]);

        const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_p[1]);
        const auto f = [&] (const int idx) {
          const int icol = idx/dim_p[1];
          const int ilev = idx%dim_p[1];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          phys(icol,ilev) = dyn(elgp[0],elgp[1],elgp[2],ilev);
        };
        Kokkos::parallel_for(tr, f);
      }
      break;
    }
    case 3:
    {
      auto phys = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT***>> (reinterpret_cast<ScalarT*>(phys_ptrs(i).get()),
                                                                           dim_p[0],
                                                                           dim_p[1],
                                                                           dim_p[2]);

      if (is_state_field_dev(i)) {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT******>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                               dim_d[0],
                                                                               dim_d[1],
                                                                               dim_d[2],
                                                                               dim_d[3],
                                                                               dim_d[4],
                                                                               dim_d[5]);

        const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_p[1]*dim_p[2]);
        const auto f = [&] (const int idx) {
          const int icol =  idx/(dim_p[1]*dim_p[2]);
          const int idim = (idx/dim_p[2])%dim_p[1];
          const int ilev =  idx%dim_p[2];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          phys(icol,idim,ilev) = dyn(elgp[0],time_levels(0).second,idim,elgp[1],elgp[2],ilev);
        };
        Kokkos::parallel_for(tr, f);
      } else {
        auto dyn = Unmanaged<KokkosTypes<DefaultDevice>::view<ScalarT*****>> (reinterpret_cast<ScalarT*>(dyn_ptrs(i).get()),
                                                                              dim_d[0],
                                                                              dim_d[1],
                                                                              dim_d[2],
                                                                              dim_d[3],
                                                                              dim_d[4]);

        const auto tr = Kokkos::TeamThreadRange(team, num_cols*dim_p[1]*dim_p[2]);
        const auto f = [&] (const int idx) {
          const int icol =  idx/(dim_p[1]*dim_p[2]);
          const int idim = (idx/dim_p[2])%dim_p[1];
          const int ilev =  idx%dim_p[2];

          const auto& elgp = Kokkos::subview(lid2elgp,p2d(icol),Kokkos::ALL());
          phys(icol,idim,ilev) = dyn(elgp[0],idim,elgp[1],elgp[2],ilev);
        };
        Kokkos::parallel_for(tr, f);
      }
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template<typename RealType>
void PhysicsDynamicsRemapper<RealType>::
create_p2d_map () {
  auto num_phys_dofs = m_phys_grid->get_num_local_dofs();
  auto num_dyn_dofs  = m_dyn_grid->get_num_local_dofs();

  auto dyn_gids  = m_dyn_grid->get_dofs_gids();
  auto phys_gids = m_phys_grid->get_dofs_gids();

  auto policy = KokkosTypes<DefaultDevice>::RangePolicy(0,num_phys_dofs);
  m_p2d = decltype(m_p2d) ("",num_phys_dofs);
  auto p2d = m_p2d;

  Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idof){
    auto gid = phys_gids(idof);
    bool found = false;
    for (int i=0; i<num_dyn_dofs; ++i) {
      if (dyn_gids(i)==gid) {
        p2d(idof) = i;
        found = true;
        break;
      }
    }
    EKAT_KERNEL_ASSERT_MSG (found, "Error! Physics grid gid not found in the dynamics grid.\n");
  });
}

} // namespace scream

#endif // SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
