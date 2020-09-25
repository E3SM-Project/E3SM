#ifndef SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
#define SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP

#include "share/scream_config.hpp"

#include "dynamics/homme/hommexx_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"

#include "share/grid/remap/abstract_remapper.hpp"
#include "share/grid/se_grid.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_assert.hpp"

// Homme includes
#include "Context.hpp"
#include "TimeLevel.hpp"
#include "Types.hpp"
#include "mpi/MpiContext.hpp"
#include "mpi/Connectivity.hpp"
#include "mpi/BoundaryExchange.hpp"
#include "mpi/BuffersManager.hpp"

namespace scream
{

template<typename DataType,typename ScalarType,typename DeviceType>
::Homme::ExecViewUnmanaged<DataType>
getHommeView(const Field<ScalarType,DeviceType>& f) {
  auto scream_view = f.template get_reshaped_view<DataType>();
  return ::Homme::ExecViewUnmanaged<DataType>(scream_view.data(),scream_view.layout());
}

// Performs remap from physics to dynamics grids, and viceversa
template<typename ScalarType, typename DeviceType>
class PhysicsDynamicsRemapper : public AbstractRemapper<ScalarType,DeviceType>
{
public:
  using scalar_type     = ScalarType;
  using device_type     = DeviceType;
  using base_type       = AbstractRemapper<scalar_type,device_type>;
  using field_type      = typename base_type::field_type;
  using identifier_type = typename base_type::identifier_type;
  using layout_type     = typename base_type::layout_type;
  using grid_ptr_type   = typename base_type::grid_ptr_type;
  using kt              = KokkosTypes<DeviceType>;

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
  void do_unregister_field (const int ifield) override;
  void do_registration_ends () override;

  // Remap methods
  void do_remap_fwd () const override;
  void do_remap_bwd () const override;

  void setup_boundary_exchange ();

  std::vector<field_type>   m_phys;
  std::vector<field_type>   m_dyn;

  std::vector<bool>         m_is_state_field;
  std::vector<bool>         m_is_tracer_field;

  std::shared_ptr<const SEGrid>   m_phys_grid;

  std::shared_ptr<Homme::BoundaryExchange>  m_be[HOMMEXX_NUM_TIME_LEVELS][HOMMEXX_Q_NUM_TIME_LEVELS];

public:
  // These functions should be morally privade, but CUDA does not allow extended host-device lambda
  // to have private/protected access within the class
  void local_remap_fwd_2d (const field_type& phys, const field_type& dyn, const int itl) const;
  template<typename ScalarT>
  void local_remap_fwd_3d_impl (const field_type& phys, const field_type& dyn, const int itl) const;

  void remap_bwd_2d (const field_type& phys, const field_type& dyn, const int itl) const;
  template<typename ScalarT>
  void remap_bwd_3d_impl (const field_type& phys, const field_type& dyn, const int itl) const;
};

// ================= IMPLEMENTATION ================= //

template<typename ScalarType, typename DeviceType>
PhysicsDynamicsRemapper<ScalarType,DeviceType>::
PhysicsDynamicsRemapper (const grid_ptr_type& phys_grid,
                         const grid_ptr_type& dyn_grid)
 : base_type(phys_grid,dyn_grid)
{
  EKAT_REQUIRE_MSG(dyn_grid->type()==GridType::SE_CellBased,  "Error! Input dynamics grid is not a Dynamics grid.\n");
  EKAT_REQUIRE_MSG(phys_grid->type()==GridType::SE_NodeBased, "Error! Input physics grid is not a Physics grid.\n");

  m_phys_grid = std::dynamic_pointer_cast<const SEGrid>(phys_grid);
}

template<typename ScalarType, typename DeviceType>
FieldLayout PhysicsDynamicsRemapper<ScalarType,DeviceType>::
create_src_layout (const FieldLayout& tgt_layout) const {
  namespace SFTN = ShortFieldTagsNames;

  auto tags = tgt_layout.tags();
  auto dims = tgt_layout.dims();

  // Note down the position of the first 'GaussPoint' tag.
  int pos = std::distance(tags.cbegin(),ekat::find(tags,SFTN::GP));

  // We replace 'Element' with 'Column'. The number of columns is taken from the src grid.
  tags[0] = SFTN::COL;
  dims[0] = this->m_src_grid->get_num_local_dofs();

  // Delete GP tags/dims
  ekat::erase(tags,SFTN::GP);
  ekat::erase(tags,SFTN::GP);
  dims.erase(dims.begin()+pos);
  dims.erase(dims.begin()+pos);

  // If the tgt layout contains the TimeLevel tag, we slice it off.
  auto it_tl = ekat::find(tags,SFTN::TL);
  if (it_tl!=tags.end()) {
    pos = std::distance(tags.cbegin(),it_tl);
    tags.erase(tags.begin()+pos);
    dims.erase(dims.begin()+pos);
  }

  return FieldLayout(tags,dims);
}

template<typename ScalarType, typename DeviceType>
FieldLayout PhysicsDynamicsRemapper<ScalarType,DeviceType>::
create_tgt_layout (const FieldLayout& src_layout) const {
  namespace SFTN = ShortFieldTagsNames;

  auto tags = src_layout.tags();
  auto dims = src_layout.dims();

  // Replace COL with EL, and num_cols with num_elems
  tags[0] = SFTN::EL;
  dims[0] = this->m_tgt_grid->get_num_local_dofs() / (HOMMEXX_NP*HOMMEXX_NP);

  // For position of GP and NP, it's easier to switch between 2d and 3d
  auto lt = get_layout_type(tags);
  switch (lt) {
    case LayoutType::Scalar2D:
    case LayoutType::Vector2D:
    case LayoutType::Tensor2D:
      // Simple: GP/NP are at the end.
      // Push back GP/NP twice
      tags.push_back(SFTN::GP);
      tags.push_back(SFTN::GP);
      dims.push_back(HOMMEXX_NP);
      dims.push_back(HOMMEXX_NP);
      break;
    case LayoutType::Scalar3D:
    case LayoutType::Vector3D:
    case LayoutType::Tensor3D:
      {
        // Replace last tag/tim with GP/NP, then push back GP/NP and VL/nvl

        // Note down num levels
        const int nvl = dims.back();
        tags.back() = SFTN::GP;
        dims.back() = HOMMEXX_NP;

        tags.push_back(SFTN::GP);
        dims.push_back(HOMMEXX_NP);

        tags.push_back(SFTN::VL);
        dims.push_back(nvl);
        break;
      }
    default:
      EKAT_ERROR_MSG("Error! Unrecognized layout type.\n");
  }

  return FieldLayout(tags,dims);
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_register_field (const identifier_type& src, const identifier_type& tgt)
{
  m_phys.push_back(field_type(src));
  m_dyn.push_back(field_type(tgt));
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_bind_field (const int ifield, const field_type& src, const field_type& tgt)
{
  const auto& tgt_layout = tgt.get_header().get_identifier().get_layout();
  const auto& tgt_tags = tgt_layout.tags();
  const auto& tgt_dims = tgt_layout.dims();

  const bool has_time_level  = ekat::contains(tgt_tags,FieldTag::TimeLevel);
  if (has_time_level) {
    const auto& data = tgt.get_header().get_extra_data();

    const bool is_tracer = data.find("Is Tracer State")!=data.end() &&
                           ekat::any_cast<bool>(data.at("Is Tracer State"));
    const bool valid_tl_dim = (tgt_dims[1]==HOMMEXX_NUM_TIME_LEVELS) || (tgt_dims[1]==HOMMEXX_Q_NUM_TIME_LEVELS);
    EKAT_REQUIRE_MSG (valid_tl_dim, "Error! Field has the TimeLevel tag, but it does not appear to be either a 'state' or 'tracer state'.\n");
    m_is_tracer_field.push_back(is_tracer);
  } else {
    m_is_tracer_field.push_back(false);
  }
  m_is_state_field.push_back(has_time_level);
  m_phys[ifield] = src;
  m_dyn[ifield] = tgt;

  // If this was the last field to be bound, we can setup the BE
  if (this->m_state==RepoState::Closed &&
      (this->m_num_bound_fields+1)==this->m_num_registered_fields) {
    setup_boundary_exchange ();
  }
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_unregister_field (const int ifield)
{
  m_phys.erase(m_phys.begin()+ifield);
  m_dyn.erase(m_dyn.begin()+ifield);
  m_is_state_field.erase(m_is_state_field.begin()+ifield);
  m_is_tracer_field.erase(m_is_tracer_field.begin()+ifield);

  // If unregistering this field makes all fields bound, we can setup the BE
  if (this->m_state==RepoState::Closed &&
      (this->m_num_bound_fields==(this->m_num_registered_fields+1))) {
    setup_boundary_exchange ();
  }
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_registration_ends ()
{
  // If we have all fields allocated, we can setup the BE
  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    setup_boundary_exchange ();
  }
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_remap_fwd() const
{
  using pack_type = ekat::Pack<ScalarType,SCREAM_PACK_SIZE>;
  using small_pack_type = ekat::Pack<ScalarType,SCREAM_SMALL_PACK_SIZE>;

  const auto& tl = Homme::Context::singleton().get_time_level();

  const int num_fields = m_phys.size();
  for (int i=0; i<num_fields; ++i) {
    const auto& phys = m_phys[i];
    const auto& dyn  = m_dyn[i];

    const int itl = m_is_state_field[i] ? (m_is_tracer_field[i] ? tl.np1_qdp : tl.np1) : -1;

    // phys->dyn requires a halo-exchange. Since not all entries in dyn
    // are overwritten before the exchange, to avoid leftover garbage,
    // we need to set all entries of dyn to zero.
    Kokkos::deep_copy(dyn.get_view(),0.0);

    const auto& ph = phys.get_header();
    const auto& dh = dyn.get_header();

    const auto& phys_alloc_prop = ph.get_alloc_properties();
    const auto& dyn_alloc_prop  = dh.get_alloc_properties();

    const bool phys_pack_alloc = phys_alloc_prop.template is_allocation_compatible_with_value_type<pack_type>();
    const bool dyn_pack_alloc  = dyn_alloc_prop.template  is_allocation_compatible_with_value_type<pack_type>();
    const bool phys_small_pack_alloc = phys_alloc_prop.template is_allocation_compatible_with_value_type<small_pack_type>();
    const bool dyn_small_pack_alloc  = dyn_alloc_prop.template  is_allocation_compatible_with_value_type<small_pack_type>();

    const auto phys_lt = get_layout_type(ph.get_identifier().get_layout().tags());
    switch (phys_lt) {
      case LayoutType::Scalar2D:
      case LayoutType::Vector2D:
      case LayoutType::Tensor2D:
        local_remap_fwd_2d(phys,dyn,itl);
        break;
      case LayoutType::Scalar3D:
      case LayoutType::Vector3D:
      case LayoutType::Tensor3D:
        if (phys_pack_alloc && dyn_pack_alloc) {
          local_remap_fwd_3d_impl<pack_type>(phys,dyn,itl);
        } else if (phys_small_pack_alloc && dyn_small_pack_alloc) {
          local_remap_fwd_3d_impl<small_pack_type>(phys,dyn,itl);
        } else {
          local_remap_fwd_3d_impl<Real>(phys,dyn,itl);
        }
        break;
      default:
        ekat::error::runtime_abort("Error! Unhandled case in switch statement.\n");
    }
  }

  // Exchange only the current time levels
  m_be[tl.np1][tl.np1_qdp]->exchange();
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_remap_bwd() const {
  using pack_type = ekat::Pack<ScalarType,SCREAM_PACK_SIZE>;
  using small_pack_type = ekat::Pack<ScalarType,SCREAM_SMALL_PACK_SIZE>;

  const auto& tl = Homme::Context::singleton().get_time_level();

  const int num_fields = m_dyn.size();
  for (int i=0; i<num_fields; ++i) {
    const auto& phys = m_phys[i];
    const auto& dyn  = m_dyn[i];

    const auto& ph = phys.get_header();
    const auto& dh = dyn.get_header();

    const auto& phys_alloc_prop = ph.get_alloc_properties();
    const auto& dyn_alloc_prop  = dh.get_alloc_properties();

    const bool phys_pack_alloc = phys_alloc_prop.template is_allocation_compatible_with_value_type<pack_type>();
    const bool dyn_pack_alloc  = dyn_alloc_prop.template  is_allocation_compatible_with_value_type<pack_type>();
    const bool phys_small_pack_alloc = phys_alloc_prop.template is_allocation_compatible_with_value_type<small_pack_type>();
    const bool dyn_small_pack_alloc  = dyn_alloc_prop.template  is_allocation_compatible_with_value_type<small_pack_type>();

    const int itl = m_is_state_field[i] ? (m_is_tracer_field[i] ? tl.np1_qdp : tl.np1) : -1;

    const LayoutType lt = get_layout_type(ph.get_identifier().get_layout().tags());
    switch (lt) {
      case LayoutType::Scalar2D:
      case LayoutType::Vector2D:
      case LayoutType::Tensor2D:
        remap_bwd_2d(phys,dyn,itl);
        break;
      case LayoutType::Scalar3D:
      case LayoutType::Vector3D:
      case LayoutType::Tensor3D:
        if (phys_pack_alloc && dyn_pack_alloc) {
          remap_bwd_3d_impl<pack_type>(phys,dyn,itl);
        } else if (phys_small_pack_alloc && dyn_small_pack_alloc) {
          remap_bwd_3d_impl<small_pack_type>(phys,dyn,itl);
        } else {
          remap_bwd_3d_impl<Real>(phys,dyn,itl);
        }
        break;
      default:
        ekat::error::runtime_abort("Error! Unhandled case in switch statement.\n");
    }
  }
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
setup_boundary_exchange () {
  // TODO: should we check that the BE was not already setup? I don't see this happening, but even if it does,
  //       there should be no side effect if we re-build it. We waste some time, sure, but should yield a correct result.

  using Scalar = Homme::Scalar;

  int num_2d = 0;
  int num_3d = 0;
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
          // A state: we only exchange the timelevel we remapped
          num_2d += 1;
        } else {
          // Not a state: remap all slices
          num_2d += layout.dim(1);
        }
        break;
      case LayoutType::Tensor2D:
        if (m_is_state_field[i]) {
          // A state: we only exchange the timelevel we remapped
          num_2d += layout.dim(2);
        } else {
          // Not a state: remap all slices
          num_2d += layout.dim(1)*layout.dim(2);
        }
        break;
      case LayoutType::Scalar3D:
        ++num_3d;
        break;
      case LayoutType::Vector3D:
        if (m_is_state_field[i]) {
          // A state: we only exchange the timelevel we remapped
          num_3d += 1;
        } else {
          // Not a state: remap all slices
          num_3d += layout.dim(1);
        }
        break;
      case LayoutType::Tensor3D:
        if (m_is_state_field[i]) {
          // A state: we only exchange the timelevel we remapped
          num_3d += layout.dim(2);
        } else {
          // Not a state: remap all slices
          num_3d += layout.dim(1)*layout.dim(2);
        }
        break;
    default:
      ekat::error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
    }
  }

  auto bm   = Homme::MpiContext::singleton().get_buffers_manager(Homme::MPI_EXCHANGE);
  auto conn = Homme::MpiContext::singleton().get_connectivity();
  constexpr int NL= HOMMEXX_NUM_LEV;
  for (int it=0; it<HOMMEXX_NUM_TIME_LEVELS; ++it) {
    for (int itq=0; itq<HOMMEXX_Q_NUM_TIME_LEVELS; ++itq) {
      m_be[it][itq] = std::make_shared<Homme::BoundaryExchange>(conn,bm);

      auto be = m_be[it][itq];
      be->set_num_fields(0,num_2d,num_3d);

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
              if (m_is_tracer_field[i]) {
                // A tracer state: exchange one time level only
                be->register_field(getHommeView<Real**[NP][NP]>(m_dyn[i]),1,itq);
              } else {
                // A state: exchange one time level only
                be->register_field(getHommeView<Real**[NP][NP]>(m_dyn[i]),1,it);
              }
            } else {
              // Not a state: exchange all slices
              be->register_field(getHommeView<Real**[NP][NP]>(m_dyn[i]),dims[1],0);
            }
            break;
          case LayoutType::Tensor2D:
            if (m_is_state_field[i]) {
              if (m_is_tracer_field[i]) {
                // A tracer state: exchange one time level only
                be->register_field(getHommeView<Real***[NP][NP]>(m_dyn[i]),itq,dims[2],0);
              } else {
                // A state: exchange one time level only
                be->register_field(getHommeView<Real***[NP][NP]>(m_dyn[i]),it,dims[2],0);
              }
            } else {
              // Not a state: exchange all slices
              for (int idim=0; idim<dims[1]; ++idim) {
                // Homme::BoundaryExchange only exchange one slice of the outer dim at a time,
                // so loop on the outer dim and register each slice individually.
                be->register_field(getHommeView<Real***[NP][NP]>(m_dyn[i]),idim,dims[2],0);
              }
            }
            break;
          case LayoutType::Scalar3D:
            be->register_field(getHommeView<Scalar*[NP][NP][HOMMEXX_NUM_LEV]>(m_dyn[i]));
            break;
          case LayoutType::Vector3D:
            if (m_is_state_field[i]) {
              if (m_is_tracer_field[i]) {
                // A tracer state: exchange one time level only
                be->register_field(getHommeView<Scalar**[NP][NP][NL]>(m_dyn[i]),1,0);
              } else {
                // A state: exchange one time level only
                be->register_field(getHommeView<Scalar**[NP][NP][NL]>(m_dyn[i]),1,0);
              }
            } else {
              // Not a state: exchange all slices
              be->register_field(getHommeView<Scalar**[NP][NP][NL]>(m_dyn[i]),dims[1],0);
            }
            break;
          case LayoutType::Tensor3D:
            if (m_is_state_field[i]) {
              if (m_is_tracer_field[i]) {
                // A tracer state: exchange one time level only
                be->register_field(getHommeView<Scalar***[NP][NP][NL]>(m_dyn[i]),itq,dims[2],0);
              } else {
                // A state: exchange one time level only
                be->register_field(getHommeView<Scalar***[NP][NP][NL]>(m_dyn[i]),it,dims[2],0);
              }
            } else {
              // Not a state: exchange all slices
              for (int idim=0; idim<dims[1]; ++idim) {
                // Homme::BoundaryExchange only exchange one slice of the outer dim at a time,
                // so loop on the outer dim and register each slice individually.
                be->register_field(getHommeView<Scalar***[NP][NP][NL]>(m_dyn[i]),idim,dims[2],0);
              }
            }
            break;
        default:
          ekat::error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
        }
      }
      be->registration_completed();
    }
  }
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
local_remap_fwd_2d(const field_type& phys_field, const field_type& dyn_field, const int itl) const
{
  using RangePolicy = Kokkos::RangePolicy<typename kt::ExeSpace>;

  auto p2d = m_phys_grid->get_dofs_map();
  const int num_cols = p2d.extent_int(0);

  const auto& phys_layout   = phys_field.get_header().get_identifier().get_layout();
  const auto& phys_dims = phys_layout.dims();
  switch (get_layout_type(phys_layout.tags())) {
    case LayoutType::Scalar2D:
    {
      auto phys = phys_field.template get_reshaped_view<Real*>();
      if (itl>=0) {
        auto dyn  = dyn_field.template get_reshaped_view<Real**[NP][NP]>();
        Kokkos::parallel_for(RangePolicy(0,num_cols),
                             KOKKOS_LAMBDA(const int icol) {
          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          dyn(elgp[0],itl,elgp[1],elgp[2]) = phys(icol);
        });
      } else {
        auto dyn = dyn_field.template get_reshaped_view<Real*[NP][NP]>();
        Kokkos::parallel_for(RangePolicy(0,num_cols),
                             KOKKOS_LAMBDA(const int icol) {
          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          dyn(elgp[0],elgp[1],elgp[2]) = phys(icol);
        });
      }
      break;
    }
    case LayoutType::Vector2D:
    {
      auto phys = phys_field.template get_reshaped_view<Real**>();
      if (itl>=0) {
        auto dyn = dyn_field.template get_reshaped_view<Real***[NP][NP]>();
        const int dim = phys_dims[1];
        Kokkos::parallel_for(RangePolicy(0,num_cols*dim),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol = idx / dim;
          const int idim = idx % dim;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          dyn(elgp[0],itl,idim,elgp[1],elgp[2]) = phys(icol,idim);
        });
      } else {
        auto dyn = dyn_field.template get_reshaped_view<Real**[NP][NP]>();
        const int dim = phys_dims[1];
        Kokkos::parallel_for(RangePolicy(0,num_cols*dim),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol = idx / dim;
          const int idim = idx % dim;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          dyn(elgp[0],idim,elgp[1],elgp[2]) = phys(icol,idim);
        });
      }
      break;
    }
    case LayoutType::Tensor2D:
    {
      auto phys = phys_field.template get_reshaped_view<Real***>();
      auto dyn  = dyn_field.template get_reshaped_view<Real***[NP][NP]>();
      const int dim1 = phys_dims[1];
      const int dim2 = phys_dims[2];
      Kokkos::pair<int,int> ordering;
      if (phys_layout.tag(1)==dyn_field.get_header().get_identifier().get_layout().tag(1)) {
        ordering.first=0;
        ordering.second=1;
      } else {
        ordering.first=1;
        ordering.second=0;
      }
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim1*dim2),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / (dim1*dim2);
        const int dims [2] = { (idx/dim2)%dim1 , idx%dim2 };

        const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
        dyn(elgp[0],dims[ordering.first],dims[ordering.second],elgp[1],elgp[2]) = phys(icol,dims[0],dims[1]);
      });
      break;
    }
    default:
      ekat::error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
  }
  Kokkos::fence();
}

template<typename ScalarType, typename DeviceType>
template<typename ScalarT>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
local_remap_fwd_3d_impl(const field_type& phys_field, const field_type& dyn_field, const int itl) const {
  using RangePolicy = Kokkos::RangePolicy<typename kt::ExeSpace>;

  auto p2d = m_phys_grid->get_dofs_map();
  const int num_cols = p2d.extent_int(0);

  const auto& phys_layout   = phys_field.get_header().get_identifier().get_layout();
  const auto& phys_dims = phys_layout.dims();

  constexpr int pack_size = sizeof(ScalarT) / sizeof(Real);
  const int NumVerticalLevels = (phys_dims.back() + pack_size - 1) / pack_size;
  switch (get_layout_type(phys_layout.tags())) {
    case LayoutType::Scalar3D:
    {
      auto phys = phys_field.template get_reshaped_view<ScalarT**>();
      if (itl>=0) {
        auto dyn  = dyn_field.template get_reshaped_view<ScalarT*****>();
        Kokkos::parallel_for(RangePolicy(0,num_cols*NumVerticalLevels),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol = idx / NumVerticalLevels;
          const int ilev = idx % NumVerticalLevels;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          dyn(elgp[0],itl,elgp[1],elgp[2],ilev) = phys(icol,ilev);
        });
      } else {
        auto dyn  = dyn_field.template get_reshaped_view<ScalarT****>();
        Kokkos::parallel_for(RangePolicy(0,num_cols*NumVerticalLevels),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol = idx / NumVerticalLevels;
          const int ilev = idx % NumVerticalLevels;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          dyn(elgp[0],elgp[1],elgp[2],ilev) = phys(icol,ilev);
        });
      }
      break;
    }
    case LayoutType::Vector3D:
    {
      auto phys = phys_field.template get_reshaped_view<ScalarT***>();
      const int dim = phys_dims[1];
      if (itl>=0) {
        auto dyn  = dyn_field.template get_reshaped_view<ScalarT******>();
        Kokkos::parallel_for(RangePolicy(0,num_cols*dim*NumVerticalLevels),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol =  idx / (dim*NumVerticalLevels);
          const int idim = (idx / NumVerticalLevels) % dim;
          const int ilev =  idx % NumVerticalLevels;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          dyn(elgp[0],itl,idim,elgp[1],elgp[2],ilev) = phys(icol,idim,ilev);
        });
      } else {
        auto dyn  = dyn_field.template get_reshaped_view<ScalarT*****>();
        Kokkos::parallel_for(RangePolicy(0,num_cols*dim*NumVerticalLevels),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol =  idx / (dim*NumVerticalLevels);
          const int idim = (idx / NumVerticalLevels) % dim;
          const int ilev =  idx % NumVerticalLevels;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          dyn(elgp[0],idim,elgp[1],elgp[2],ilev) = phys(icol,idim,ilev);
        });
      }
      break;
    }
    case LayoutType::Tensor3D:
    {
      auto phys = phys_field.template get_reshaped_view<ScalarT****>();
      auto dyn  = dyn_field.template get_reshaped_view<ScalarT******>();
      const int dim1 = phys_dims[1];
      const int dim2 = phys_dims[2];
      Kokkos::pair<int,int> ordering;
      if (phys_layout.tag(1)==dyn_field.get_header().get_identifier().get_layout().tag(1)) {
        ordering.first=0;
        ordering.second=1;
      } else {
        ordering.first=1;
        ordering.second=0;
      }
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim1*dim2*NumVerticalLevels),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol =  idx / (dim1*dim2*NumVerticalLevels);
        const int dims [2] = { (idx/NumVerticalLevels)%dim1 , (idx/NumVerticalLevels)%dim2 };
        const int ilev =  idx % NumVerticalLevels;

        const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
        dyn(elgp[0],dims[ordering.first],dims[ordering.second],elgp[1],elgp[1],ilev) = phys(icol,dims[0],dims[1],ilev);
      });
      break;
    }
    default:
      ekat::error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
  }
  Kokkos::fence();
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
remap_bwd_2d(const field_type& phys_field, const field_type& dyn_field, const int itl) const {
  using RangePolicy = Kokkos::RangePolicy<typename kt::ExeSpace>;

  auto p2d = m_phys_grid->get_dofs_map();
  const int num_cols = p2d.extent_int(0);

  const auto& phys_layout = phys_field.get_header().get_identifier().get_layout();
  const auto& dyn_layout  = dyn_field.get_header().get_identifier().get_layout();

  const auto& phys_dims = phys_layout.dims();
  switch (phys_dims.size()) {
    case 1:
    {
      auto phys = phys_field.template get_reshaped_view<Real*>();
      if (itl>=0) {
        auto dyn  = dyn_field.template get_reshaped_view<Real**[NP][NP]>();
        Kokkos::parallel_for(RangePolicy(0,num_cols),
                             KOKKOS_LAMBDA(const int icol) {
          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          phys(icol) = dyn(elgp[0],itl,elgp[1],elgp[2]);
        });
      } else {
        auto dyn  = dyn_field.template get_reshaped_view<Real*[NP][NP]>();
        Kokkos::parallel_for(RangePolicy(0,num_cols),
                             KOKKOS_LAMBDA(const int icol) {
          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          phys(icol) = dyn(elgp[0],elgp[1],elgp[2]);
        });
      }
      break;
    }
    case 2:
    {
      auto phys = phys_field.template get_reshaped_view<Real**>();
      const int dim = phys_dims[1];
      if (itl>=0) {
        auto dyn  = dyn_field.template get_reshaped_view<Real***[NP][NP]>();
        Kokkos::parallel_for(RangePolicy(0,num_cols*dim),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol = idx / dim;
          const int idim = idx % dim;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          phys(icol,idim) = dyn(elgp[0],itl,idim,elgp[1],elgp[2]);
        });
      } else {
        auto dyn  = dyn_field.template get_reshaped_view<Real**[NP][NP]>();
        Kokkos::parallel_for(RangePolicy(0,num_cols*dim),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol = idx / dim;
          const int idim = idx % dim;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          phys(icol,idim) = dyn(elgp[0],idim,elgp[1],elgp[2]);
        });
      }
      break;
    }
    case 3:
    {
      auto phys = phys_field.template get_reshaped_view<Real***>();
      auto dyn  = dyn_field.template get_reshaped_view<Real***[NP][NP]>();
      const int dim1 = phys_dims[1];
      const int dim2 = phys_dims[2];
      Kokkos::pair<int,int> ordering;
      if (phys_layout.tag(1)==dyn_layout.tag(1)) {
        ordering.first=0;
        ordering.second=1;
      } else {
        ordering.first=1;
        ordering.second=0;
      }
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim1*dim2),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / (dim1*dim2);
        const int dims [2] = { (idx/dim2)%dim1 , idx%dim2 };

        const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
        phys(icol,dims[ordering.first],dims[ordering.second]) = dyn(elgp[0],dims[0],dims[1],elgp[1],elgp[2]);
      });
      break;
    }
    default:
      ekat::error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
  }
  Kokkos::fence();
}

template<typename ScalarType, typename DeviceType>
template<typename ScalarT>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
remap_bwd_3d_impl(const field_type& phys_field, const field_type& dyn_field, const int itl) const {
  using RangePolicy = Kokkos::RangePolicy<typename kt::ExeSpace>;

  auto p2d = m_phys_grid->get_dofs_map();
  const int num_cols = p2d.extent_int(0);

  const auto& phys_layout = phys_field.get_header().get_identifier().get_layout();
  const auto& dyn_layout  = dyn_field.get_header().get_identifier().get_layout();

  const auto& phys_dims = phys_layout.dims();

  constexpr int pack_size = sizeof(ScalarT) / sizeof(Real);
  const int NumVerticalLevels = (phys_dims.back() + pack_size - 1) / pack_size;
  switch (phys_dims.size()) {
    case 2:
    {
      auto phys = phys_field.template get_reshaped_view<ScalarT**>();
      if (itl>=0) {
        auto dyn = dyn_field.template get_reshaped_view<ScalarT*****>();
        Kokkos::parallel_for(RangePolicy(0,num_cols*NumVerticalLevels),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol = idx / NumVerticalLevels;
          const int ilev = idx % NumVerticalLevels;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          phys(icol,ilev) = dyn(elgp[0],itl,elgp[1],elgp[2],ilev);
        });
      } else {
        auto dyn = dyn_field.template get_reshaped_view<ScalarT****>();
        Kokkos::parallel_for(RangePolicy(0,num_cols*NumVerticalLevels),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol = idx / NumVerticalLevels;
          const int ilev = idx % NumVerticalLevels;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          phys(icol,ilev) = dyn(elgp[0],elgp[1],elgp[2],ilev);
        });
      }
      break;
    }
    case 3:
    {
      auto phys = phys_field.template get_reshaped_view<ScalarT***>();
      const int dim = phys_dims[1];
      if (itl>=0) {
        auto dyn = dyn_field.template get_reshaped_view<ScalarT******>();
        Kokkos::parallel_for(RangePolicy(0,num_cols*dim*NumVerticalLevels),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol =  idx / (dim*NumVerticalLevels);
          const int idim = (idx / NumVerticalLevels) % dim;
          const int ilev =  idx % NumVerticalLevels;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          phys(icol,idim,ilev) = dyn(elgp[0],itl,idim,elgp[1],elgp[2],ilev);
        });
      } else {
        auto dyn = dyn_field.template get_reshaped_view<ScalarT*****>();
        Kokkos::parallel_for(RangePolicy(0,num_cols*dim*NumVerticalLevels),
                             KOKKOS_LAMBDA(const int idx) {
          const int icol =  idx / (dim*NumVerticalLevels);
          const int idim = (idx / NumVerticalLevels) % dim;
          const int ilev =  idx % NumVerticalLevels;

          const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
          phys(icol,idim,ilev) = dyn(elgp[0],idim,elgp[1],elgp[2],ilev);
        });
      }
      break;
    }
    case 4:
    {
      auto phys = phys_field.template get_reshaped_view<ScalarT****>();
      auto dyn  = dyn_field.template get_reshaped_view<ScalarT******>();
      const int dim1 = phys_dims[1];
      const int dim2 = phys_dims[2];
      Kokkos::pair<int,int> ordering;
      if (phys_layout.tag(1)==dyn_layout.tag(1)) {
        ordering.first=0;
        ordering.second=1;
      } else {
        ordering.first=1;
        ordering.second=0;
      }
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim1*dim2*NumVerticalLevels),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / (dim1*dim2*NumVerticalLevels);
        const int dims [2] = { (idx / (dim2*NumVerticalLevels)) % dim1 , (idx / NumVerticalLevels) % dim2 };
        const int ilev =  idx % NumVerticalLevels;

        const auto& elgp = Kokkos::subview(p2d,icol,Kokkos::ALL());
        phys(icol,dims[ordering.first],dims[ordering.second],ilev) = dyn(elgp[0],dims[0],dims[1],elgp[1],elgp[2],ilev);
      });
      break;
    }
    default:
      ekat::error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
  }
  Kokkos::fence();
}

} // namespace scream

#endif // SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
