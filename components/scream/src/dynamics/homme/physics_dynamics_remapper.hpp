#ifndef SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
#define SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP

#include "share/field/remap/abstract_remapper.hpp"
#include "share/field/field_tag.hpp"
#include "share/scream_pack.hpp"
#include "share/field/remap/remap_utils.hpp"

// Homme includes
#include "Types.hpp"
#include "mpi/MpiContext.hpp"
#include "mpi/Connectivity.hpp"
#include "mpi/BoundaryExchange.hpp"
#include "mpi/BuffersManager.hpp"

namespace scream
{

template<>
struct util::ScalarProperties<Homme::Scalar> {
  using scalar_type = Homme::Real;
  static constexpr bool is_pack = true;
};

template<>
struct util::TypeName<Homme::Scalar> {
  static std::string name () {
    return "Homme::Scalar";
  }
};

template<int N>
struct IntN {
  Int idx[N];
};

template<typename DataType,typename ScalarType,typename DeviceType>
Homme::ExecViewUnmanaged<DataType>
getHommeView(const Field<ScalarType,DeviceType>& f) {
  auto scream_view = f.template get_reshaped_view<DataType>();
  return Homme::ExecViewUnmanaged<DataType>(scream_view.data(),scream_view.layout());
}

// Performs remap from physics to dynamics grids, and viceversa
// Note: this class *ONLY* performs local remap. This means that,
//       when going from physics to dynamics layout, this remap
//       must be followed by a halo exchange
template<typename ScalarType, typename DeviceType>
class PhysicsDynamicsRemapper : public AbstractRemapper<ScalarType,DeviceType>
{
public:
  using base_type       = AbstractRemapper<ScalarType,DeviceType>;
  using field_type      = typename base_type::field_type;
  using identifier_type = typename base_type::identifier_type;
  using layout_type     = typename base_type::layout_type;
  using kt              = KokkosTypes<DeviceType>;
  using p2d_map_type    = typename kt::template view<IntN<3>*>;

  PhysicsDynamicsRemapper (const p2d_map_type& phys_to_dyn_map)
   : m_phys_to_dyn_map (phys_to_dyn_map)
  {
    // Nothing to do here
  }

  ~PhysicsDynamicsRemapper () = default;

protected:

  const layout_type& do_get_src_layout (const int ifield) const {
    return m_phys[ifield].get_header().get_identifier().get_layout();
  }
  const layout_type& do_get_tgt_layout (const int ifield) const {
    return m_dyn[ifield].get_header().get_identifier().get_layout();
  }

  void do_registration_start () override;
  void do_register_field (const field_type& src, const field_type& tgt) override;
  void do_registration_complete () override;

  void do_remap_fwd () const override;
  void do_remap_bwd () const override;

  void local_remap_fwd_2d (const field_type& src, const field_type& tgt, const LayoutType lt) const;
  template<typename ScalarT>
  void local_remap_fwd_3d_impl (const field_type& src, const field_type& tgt, const LayoutType lt) const;

  void remap_bwd_2d (const field_type& src, const field_type& tgt, const LayoutType lt) const;
  template<typename ScalarT>
  void remap_bwd_3d_impl (const field_type& src, const field_type& tgt, const LayoutType lt) const;

  void setup_boundary_exchange ();

  p2d_map_type              m_phys_to_dyn_map;

  std::vector<field_type>   m_phys;
  std::vector<field_type>   m_dyn;

  std::shared_ptr<Homme::BoundaryExchange>  m_be;
};

// ================= IMPLEMENTATION ================= //

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_registration_start ()
{
  m_phys.reserve(this->get_num_fields());
  m_dyn.reserve(this->get_num_fields());
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_register_field (const field_type& src, const field_type& tgt)
{
  error::runtime_check(static_cast<int>(m_phys.size())<this->get_num_fields(),
                       "Error! You already registered " + std::to_string(m_phys.size()) + " fields. \n"
                       "       Did you call 'set_num_fields' with the wrong input (" + std::to_string(this->get_num_fields()) + ")?\n");
  error::runtime_check(src.is_allocated(), "Error! Physics field is not yet allocated.\n");
  error::runtime_check(tgt.is_allocated(), "Error! Dynamics field is not yet allocated.\n");


  const auto phys = get_layout_specs(src.get_header().get_identifier().get_layout()).first;
  const auto dyn  = get_layout_specs(tgt.get_header().get_identifier().get_layout()).first;
  const auto plt = get_layout_specs(src.get_header().get_identifier().get_layout()).second;
  const auto dlt = get_layout_specs(tgt.get_header().get_identifier().get_layout()).second;

  error::runtime_check(phys==PhysDyn::Phys, "Error! Source field does not have a Physics layout.\n");
  error::runtime_check(dyn==PhysDyn::Dyn, "Error! Target field does not have a Dynamics layout.\n");
  error::runtime_check(plt!=LayoutType::Tensor2D, "Error! PhysicsDynamicsRemapper does not support 2d tensors yet.\n");
  error::runtime_check(dlt==plt, "Error! Source and target layouts do not match.\n");

  m_phys.push_back(src);
  m_dyn.push_back(tgt);
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_registration_complete ()
{
  using Scalar = Homme::Scalar;
  constexpr int NUM_LEV = Homme::NUM_LEV;

  auto bm   = Homme::MpiContext::singleton().get_buffers_manager(Homme::MPI_EXCHANGE);
  auto conn = Homme::MpiContext::singleton().get_connectivity();
  m_be = std::make_shared<Homme::BoundaryExchange>(conn,bm);

  int num_2d = 0;
  int num_3d = 0;
  for (int i=0; i<this->get_num_fields(); ++i) {
    auto layout = m_dyn[i].get_header().get_identifier().get_layout();
    auto lt = get_layout_specs(layout).second;
    switch (lt) {
      case LayoutType::Scalar2D:
        ++num_2d;
        break;
      case LayoutType::Vector2D:
        num_2d += layout.dim(1);
        break;
      case LayoutType::Tensor2D:
        num_2d += layout.dim(1)*layout.dim(2);
        break;
      case LayoutType::Scalar3D:
        ++num_3d;
        break;
      case LayoutType::Vector3D:
        num_3d += layout.dim(1);
        break;
      case LayoutType::Tensor3D:
        num_3d += layout.dim(1)*layout.dim(2);
        break;
    default:
      error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
    }
  }

  m_be->set_num_fields(0,num_2d,num_3d);
  for (int i=0; i<this->get_num_fields(); ++i) {
    auto dims = m_dyn[i].get_header().get_identifier().get_layout().dims();
    auto lt = get_layout_specs(m_dyn[i].get_header().get_identifier().get_layout()).second;
    switch (lt) {
      case LayoutType::Scalar2D:
        m_be->register_field(getHommeView<Real*[NP][NP]>(m_dyn[i]));
        break;
      case LayoutType::Vector2D:
        m_be->register_field(getHommeView<Real**[NP][NP]>(m_dyn[i]),dims[1],0);
        break;
      case LayoutType::Tensor2D:
        for (int idim=0; idim<dims[1]; ++idim) {
          // Homme::BoundaryExchange only exchange one slice of the outer dim at a time,
          // so loop on the outer dim and register each slice individually.
          m_be->register_field(getHommeView<Real***[NP][NP]>(m_dyn[i]),idim,dims[2],0);
        }
        break;
      case LayoutType::Scalar3D:
        m_be->register_field(getHommeView<Scalar*[NP][NP][NUM_LEV]>(m_dyn[i]));
        break;
      case LayoutType::Vector3D:
        m_be->register_field(getHommeView<Scalar**[NP][NP][NUM_LEV]>(m_dyn[i]),dims[1],0);
        break;
      case LayoutType::Tensor3D:
        for (int idim=0; idim<dims[1]; ++idim) {
          // Homme::BoundaryExchange only exchange one slice of the outer dim at a time,
          // so loop on the outer dim and register each slice individually.
          m_be->register_field(getHommeView<Scalar***[NP][NP][NUM_LEV]>(m_dyn[i]),idim,dims[2],0);
        }
        break;
    default:
      error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
    }
  }
  m_be->registration_completed();
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_remap_fwd() const {
  for (int i=0; i<this->get_num_fields(); ++i) {
    const auto& phys = m_phys[i];
    const auto& dyn  = m_dyn[i];

    const auto lt = get_layout_specs(phys.get_header().get_identifier().get_layout()).second;
    const auto& tgt_alloc_prop = dyn.get_header().get_alloc_properties();
    const auto& src_alloc_prop = phys.get_header().get_alloc_properties();
    using pack_type = pack::Pack<ScalarType,SCREAM_PACK_SIZE>;
    using small_pack_type = pack::Pack<ScalarType,SCREAM_SMALL_PACK_SIZE>;
    switch (lt) {
      case LayoutType::Scalar2D:
      case LayoutType::Vector2D:
      case LayoutType::Tensor2D:
        local_remap_fwd_2d(phys,dyn,lt);
        break;
      case LayoutType::Scalar3D:
      case LayoutType::Vector3D:
      case LayoutType::Tensor3D:
        if (src_alloc_prop.template is_allocation_compatible_with_value_type<pack_type>() &&
            tgt_alloc_prop.template is_allocation_compatible_with_value_type<pack_type>())
        {
          local_remap_fwd_3d_impl<pack_type>(phys,dyn,lt);
        } else if (src_alloc_prop.template is_allocation_compatible_with_value_type<small_pack_type>() &&
                   tgt_alloc_prop.template is_allocation_compatible_with_value_type<small_pack_type>())
        {
          local_remap_fwd_3d_impl<small_pack_type>(phys,dyn,lt);
        } else {
          local_remap_fwd_3d_impl<Real>(phys,dyn,lt);
        }
        break;
      default:
        error::runtime_abort("Error! Unhandled case in switch statement.\n");
    }
  }
  m_be->exchange();
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
do_remap_bwd() const {
  for (int i=0; i<this->get_num_fields(); ++i) {
    const auto& phys = m_phys[i];
    const auto& dyn  = m_dyn[i];

    const auto lt = get_layout_specs(phys.get_header().get_identifier().get_layout()).second;
    const auto& tgt_alloc_prop = phys.get_header().get_alloc_properties();
    const auto& src_alloc_prop = dyn.get_header().get_alloc_properties();
    using pack_type = pack::Pack<ScalarType,SCREAM_PACK_SIZE>;
    using small_pack_type = pack::Pack<ScalarType,SCREAM_SMALL_PACK_SIZE>;
    switch (lt) {
      case LayoutType::Scalar2D:
      case LayoutType::Vector2D:
      case LayoutType::Tensor2D:
        remap_bwd_2d(dyn,phys,lt);
        break;
      case LayoutType::Scalar3D:
      case LayoutType::Vector3D:
      case LayoutType::Tensor3D:
        if (src_alloc_prop.template is_allocation_compatible_with_value_type<pack_type>() &&
            tgt_alloc_prop.template is_allocation_compatible_with_value_type<pack_type>())
        {
          remap_bwd_3d_impl<pack_type>(dyn,phys,lt);
        } else if (src_alloc_prop.template is_allocation_compatible_with_value_type<small_pack_type>() &&
                   tgt_alloc_prop.template is_allocation_compatible_with_value_type<small_pack_type>())
        {
          remap_bwd_3d_impl<small_pack_type>(dyn,phys,lt);
        } else {
          remap_bwd_3d_impl<Real>(dyn,phys,lt);
        }
        break;
      default:
        error::runtime_abort("Error! Unhandled case in switch statement.\n");
    }
  }
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
local_remap_fwd_2d(const field_type& src_field, const field_type& tgt_field, const LayoutType lt) const
{
  using RangePolicy = Kokkos::RangePolicy<typename kt::ExeSpace>;

  auto p2d = m_phys_to_dyn_map;
  const int num_cols = p2d.extent_int(0);

  switch (lt) {
    case LayoutType::Scalar2D:
    {
      auto phys = src_field.template get_reshaped_view<Real*>();
      auto dyn  = tgt_field.template get_reshaped_view<Real*[NP][NP]>();
      Kokkos::parallel_for(RangePolicy(0,num_cols),
                           KOKKOS_LAMBDA(const int icol) {
        const auto& elgp = p2d(icol).idx;
        dyn(elgp[0],elgp[1],elgp[2]) = phys(icol);
      });
      break;
    }
    case LayoutType::Vector2D:
    {
      auto phys = src_field.template get_reshaped_view<Real**>();
      auto dyn  = tgt_field.template get_reshaped_view<Real**[NP][NP]>();
      const int dim = phys.extent_int(1);
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / dim;
        const int idim = idx % dim;

        const auto& elgp = p2d(icol).idx;
        dyn(elgp[0],idim,elgp[1],elgp[2]) = phys(icol,idim);
      });
      break;
    }
    case LayoutType::Tensor2D:
    {
      auto phys = src_field.template get_reshaped_view<Real***>();
      auto dyn  = tgt_field.template get_reshaped_view<Real***[NP][NP]>();
      const int dim1 = phys.extent_int(1);
      const int dim2 = phys.extent_int(2);
      int ordering[2];
      if (src_field.get_header().get_identifier().get_layout().tag(1)==tgt_field.get_header().get_identifier().get_layout().tag(1)) {
        ordering[0]=0;
        ordering[0]=1;
      } else {
        ordering[0]=1;
        ordering[0]=0;
      }
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim1*dim2),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / (dim1*dim2);
        const int dims [2] = { (idx/dim2)%dim1 , idx%dim2 };

        const auto& elgp = p2d(icol).idx;
        dyn(elgp[0],dims[ordering[0]],dims[ordering[1]],elgp[1],elgp[2]) = phys(icol,dims[0],dims[1]);
      });
      break;
    }
    default:
      error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
  }
}

template<typename ScalarType, typename DeviceType>
template<typename ScalarT>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
local_remap_fwd_3d_impl(const field_type& src_field, const field_type& tgt_field, const LayoutType lt) const {
  using RangePolicy = Kokkos::RangePolicy<typename kt::ExeSpace>;

  auto p2d = m_phys_to_dyn_map;
  const int num_cols = p2d.extent_int(0);

  switch (lt) {
    case LayoutType::Scalar3D:
    {
      auto phys = src_field.template get_reshaped_view<ScalarT**>();
      auto dyn  = tgt_field.template get_reshaped_view<ScalarT****>();
      const int NumVerticalLevels = dyn.extent_int(3);
      Kokkos::parallel_for(RangePolicy(0,num_cols*NumVerticalLevels),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / NumVerticalLevels;
        const int ilev = idx % NumVerticalLevels;

        const auto& elgp = p2d(icol).idx;
        dyn(elgp[0],elgp[1],elgp[2],ilev) = phys(icol,ilev);
      });
      break;
    }
    case LayoutType::Vector3D:
    {
      auto phys = src_field.template get_reshaped_view<ScalarT***>();
      auto dyn  = tgt_field.template get_reshaped_view<ScalarT*****>();
      const int dim = phys.extent_int(1);
      const int NumVerticalLevels = dyn.extent_int(4);
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim*NumVerticalLevels),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol =  idx / (dim*NumVerticalLevels);
        const int idim = (idx / NumVerticalLevels) % dim;
        const int ilev =  idx % NumVerticalLevels;

        const auto& elgp = p2d(icol).idx;
        dyn(elgp[0],idim,elgp[1],elgp[2],ilev) = phys(icol,idim,ilev);
      });
      break;
    }
    case LayoutType::Tensor3D:
    {
      auto phys = src_field.template get_reshaped_view<ScalarT****>();
      auto dyn  = tgt_field.template get_reshaped_view<ScalarT******>();
      const int dim1 = phys.extent_int(1);
      const int dim2 = phys.extent_int(2);
      const int NumVerticalLevels = dyn.extent_int(5);
      int ordering[2];
      if (src_field.get_header().get_identifier().get_layout().tag(1)==tgt_field.get_header().get_identifier().get_layout().tag(1)) {
        ordering[0]=0;
        ordering[0]=1;
      } else {
        ordering[0]=1;
        ordering[0]=0;
      }
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim1*dim2*NumVerticalLevels),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol =  idx / (dim1*dim2*NumVerticalLevels);
        const int dims [2] = { (idx/NumVerticalLevels)%dim1 , (idx/NumVerticalLevels)%dim2 };
        const int ilev =  idx % NumVerticalLevels;

        const auto& elgp = p2d(icol).idx;
        dyn(elgp[0],dims[ordering[0]],dims[ordering[1]],elgp[1],elgp[1],ilev) = phys(icol,dims[0],dims[1],ilev);
      });
      break;
    }
    default:
      error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
  }
}

template<typename ScalarType, typename DeviceType>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
remap_bwd_2d(const field_type& src_field, const field_type& tgt_field, const LayoutType lt) const {
  using RangePolicy = Kokkos::RangePolicy<typename kt::ExeSpace>;

  auto p2d = m_phys_to_dyn_map;
  const int num_cols = p2d.extent_int(0);

  switch (lt) {
    case LayoutType::Scalar2D:
    {
      auto dyn  = src_field.template get_reshaped_view<Real*[NP][NP]>();
      auto phys = tgt_field.template get_reshaped_view<Real*>();
      Kokkos::parallel_for(RangePolicy(0,num_cols),
                           KOKKOS_LAMBDA(const int icol) {
        const auto& elgp = p2d(icol).idx;
        phys(icol) = dyn(elgp[0],elgp[1],elgp[2]);
      });
      break;
    }
    case LayoutType::Vector2D:
    {
      auto dyn  = src_field.template get_reshaped_view<Real**[NP][NP]>();
      auto phys = tgt_field.template get_reshaped_view<Real**>();
      const int dim = dyn.extent_int(1);
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / dim;
        const int idim = idx % dim;

        const auto& elgp = p2d(icol).idx;
        phys(icol,idim) = dyn(elgp[0],idim,elgp[1],elgp[2]);
      });
      break;
    }
    case LayoutType::Tensor2D:
    {
      auto dyn  = src_field.template get_reshaped_view<Real***[NP][NP]>();
      auto phys = tgt_field.template get_reshaped_view<Real***>();
      const int dim1 = dyn.extent_int(1);
      const int dim2 = dyn.extent_int(2);
      int ordering[2];
      if (src_field.get_header().get_identifier().get_layout().tag(1)==tgt_field.get_header().get_identifier().get_layout().tag(1)) {
        ordering[0]=0;
        ordering[0]=1;
      } else {
        ordering[0]=1;
        ordering[0]=0;
      }
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim1*dim2),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / (dim1*dim2);
        const int dims [2] = { (idx/dim2)%dim1 , idx%dim2 };

        const auto& elgp = p2d(icol).idx;
        phys(icol,dims[ordering[0]],dims[ordering[1]]) = dyn(elgp[0],dims[0],dims[1],elgp[1],elgp[2]);
      });
      break;
    }
    default:
      error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
  }
}

template<typename ScalarType, typename DeviceType>
template<typename ScalarT>
void PhysicsDynamicsRemapper<ScalarType,DeviceType>::
remap_bwd_3d_impl(const field_type& src_field, const field_type& tgt_field, const LayoutType lt) const {
  using RangePolicy = Kokkos::RangePolicy<typename kt::ExeSpace>;

  auto p2d = m_phys_to_dyn_map;
  const int num_cols = p2d.extent_int(0);

  switch (lt) {
    case LayoutType::Scalar3D:
    {
      auto dyn  = src_field.template get_reshaped_view<ScalarT****>();
      auto phys = tgt_field.template get_reshaped_view<ScalarT**>();
      const int NumVerticalLevels = dyn.extent_int(3);
      Kokkos::parallel_for(RangePolicy(0,num_cols*NumVerticalLevels),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / NumVerticalLevels;
        const int ilev = idx % NumVerticalLevels;

        const auto& elgp = p2d(icol).idx;
        phys(icol,ilev) = dyn(elgp[0],elgp[1],elgp[2],ilev);
      });
      break;
    }
    case LayoutType::Vector3D:
    {
      auto dyn  = src_field.template get_reshaped_view<ScalarT*****>();
      auto phys = tgt_field.template get_reshaped_view<ScalarT***>();
      const int dim = dyn.extent_int(1);
      const int NumVerticalLevels = dyn.extent_int(4);
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim*NumVerticalLevels),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol =  idx / (dim*NumVerticalLevels);
        const int idim = (idx / NumVerticalLevels) % dim;
        const int ilev =  idx % NumVerticalLevels;

        const auto& elgp = p2d(icol).idx;
        phys(icol,idim,ilev) = dyn(elgp[0],idim,elgp[1],elgp[2],ilev);
      });
      break;
    }
    case LayoutType::Tensor3D:
    {
      auto dyn  = src_field.template get_reshaped_view<ScalarT******>();
      auto phys = tgt_field.template get_reshaped_view<ScalarT****>();
      const int dim1 = dyn.extent_int(1);
      const int dim2 = dyn.extent_int(2);
      const int NumVerticalLevels = dyn.extent_int(5);
      int ordering [2];
      if (src_field.get_header().get_identifier().get_layout().tag(1)==tgt_field.get_header().get_identifier().get_layout().tag(1)) {
        ordering[0]=0;
        ordering[0]=1;
      } else {
        ordering[0]=1;
        ordering[0]=0;
      }
      Kokkos::parallel_for(RangePolicy(0,num_cols*dim1*dim2*NumVerticalLevels),
                           KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / (dim1*dim2*NumVerticalLevels);
        const int dims [2] = { (idx / (dim2*NumVerticalLevels)) % dim1 , (idx / NumVerticalLevels) % dim2 };
        const int ilev =  idx % NumVerticalLevels;

        const auto& elgp = p2d(icol).idx;
        phys(icol,dims[ordering[0]],dims[ordering[1]],ilev) = dyn(elgp[0],dims[0],dims[1],elgp[1],elgp[2],ilev);
      });
      break;
    }
    default:
      error::runtime_abort("Error! Invalid layout. This is an internal error. Please, contact developers\n");
  }
}

} // namespace scream

#endif // SCREAM_PHYSICS_DYNAMICS_REMAPPER_HPP
