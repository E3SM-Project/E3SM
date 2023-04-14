#include "dynamics/homme/physics_dynamics_remapper.hpp"

// Homme includes
#include "Types.hpp"
#include "Context.hpp"
#include "HommexxEnums.hpp"
#include "mpi/Connectivity.hpp"
#include "mpi/MpiBuffersManager.hpp"
#include "mpi/BoundaryExchange.hpp"

#include "dynamics/homme/homme_dynamics_helpers.hpp"
#include "dynamics/homme/homme_dimensions.hpp"
#include "share/grid/se_grid.hpp"
#include "share/util/scream_utils.hpp"

#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_assert.hpp"

namespace {

// Utility function to convert a view from a Field into a Homme-compatible view
template<typename DataType>
::Homme::ExecViewUnmanaged<DataType>
getHommeView(const scream::Field& f) {
  auto p = f.get_header().get_parent().lock();
  auto scream_view = f.template get_view<DataType>();
  using homme_view_t = ::Homme::ExecViewUnmanaged<DataType>;
  if (p!=nullptr) {
    // Need to fix the mapping stride, so that it can correctly map the subfield.
    homme_view_t tmp(scream_view.data(),scream_view.layout());
    auto vm = tmp.impl_map();
    vm.m_impl_offset.m_stride = scream_view.impl_map().stride_0();
    return homme_view_t(scream_view.impl_track(),vm);
  } else {
    return homme_view_t(scream_view.data(),scream_view.layout());
  }
}

} // anonymous namespace

namespace scream
{

PhysicsDynamicsRemapper::
PhysicsDynamicsRemapper (const grid_ptr_type& phys_grid,
                         const grid_ptr_type& dyn_grid)
 : base_type(phys_grid,dyn_grid)
{
  EKAT_REQUIRE_MSG(dyn_grid->type()==GridType::SE,     "Error! Input dynamics grid is not a SE grid.\n");
  EKAT_REQUIRE_MSG(phys_grid->type()==GridType::Point, "Error! Input physics grid is not a Point grid.\n");

  m_dyn_grid  = dyn_grid;
  m_phys_grid = phys_grid;

  m_num_phys_cols = phys_grid->get_num_local_dofs();
  m_lid2elgp      = m_dyn_grid->get_lid_to_idx_map().get_view<const int**>();

  // For each phys dofs, we find a corresponding dof in the dyn grid.
  // Notice that such dyn dof may not be unique (if phys dof is on an edge
  // of a SE element), but we don't care. We just need to find a match.
  // The BoundaryExchange already takes care of syncing all shared dyn dofs.
  create_p2d_map ();
}

FieldLayout PhysicsDynamicsRemapper::
create_src_layout (const FieldLayout& tgt_layout) const {
  using namespace ShortFieldTagsNames;

  auto tags = tgt_layout.tags();
  auto dims = tgt_layout.dims();

  EKAT_REQUIRE_MSG (!ekat::contains(tags,FieldTag::TimeLevel),
      "Error! Cannot remap to a Field with the TimeLevel tag.\n"
      "       Please, provide the proper subfield instead.\n");

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

  return FieldLayout(tags,dims);
}

FieldLayout PhysicsDynamicsRemapper::
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

void PhysicsDynamicsRemapper::
do_register_field (const identifier_type& src, const identifier_type& tgt)
{
  m_phys_fields.push_back(field_type(src));
  m_dyn_fields.push_back(field_type(tgt));

  EKAT_REQUIRE_MSG (src.data_type()==field_valid_data_types().at<Real>(),
      "Error! PD remapper only works for Real data type.\n");
  EKAT_REQUIRE_MSG (tgt.data_type()==field_valid_data_types().at<Real>(),
      "Error! PD remapper only works for Real data type.\n");
}

void PhysicsDynamicsRemapper::
do_bind_field (const int ifield, const field_type& src, const field_type& tgt)
{
  const auto& tgt_layout = tgt.get_header().get_identifier().get_layout();
  const auto& tgt_tags = tgt_layout.tags();

  EKAT_REQUIRE_MSG (!ekat::contains(tgt_tags,FieldTag::TimeLevel),
      "Error! Cannot remap to a Field with the TimeLevel tag.\n"
      "       Please, provide the proper subfield instead.\n");

  m_phys_fields[ifield] = src;
  m_dyn_fields[ifield] = tgt;

  // If this was the last field to be bound, we can setup the BE and
  // precompute fields needed on device during remapper
  if (this->m_state==RepoState::Closed &&
      (this->m_num_bound_fields+1)==this->m_num_registered_fields) {
    setup_boundary_exchange ();
    initialize_device_variables();
  }
}

void PhysicsDynamicsRemapper::
do_registration_ends ()
{
  // If we have all fields allocated, we can setup the BE and
  // precompute fields needed on device during remapper
  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    setup_boundary_exchange ();
    initialize_device_variables();
  }
}

void PhysicsDynamicsRemapper::
initialize_device_variables()
{
  m_layout              = decltype(m_layout)              ("layout", this->m_num_fields);
  m_pack_alloc_property = decltype(m_pack_alloc_property) ("pack_alloc_property", this->m_num_fields);
  m_num_levels          = decltype(m_num_levels)          ("num_physical_levels", this->m_num_fields);

  for (auto which : {'P','D'}) {
    auto& repo   = which=='P' ? m_phys_repo : m_dyn_repo;
    auto& fields = which=='P' ? m_phys_fields: m_dyn_fields;

    repo.views  = ViewsRepo::views_t ("views", this->m_num_fields);
    repo.cviews = ViewsRepo::cviews_t("cviews",this->m_num_fields);

    repo.h_views  = Kokkos::create_mirror_view(repo.views);
    repo.h_cviews = Kokkos::create_mirror_view(repo.cviews);

    auto get_view = [&] (const int i, const Field& f) {
      const auto rank = f.get_header().get_identifier().get_layout().rank();
      switch (rank) {
        case 1: repo.h_cviews[i].v1d = f.get_view<const Real*>();      break;
        case 2: repo.h_cviews[i].v2d = f.get_view<const Real**>();     break;
        case 3: repo.h_cviews[i].v3d = f.get_view<const Real***>();    break;
        case 4: repo.h_cviews[i].v4d = f.get_view<const Real****>();   break;
        case 5: repo.h_cviews[i].v5d = f.get_view<const Real*****>();  break;
      }
      if (not f.is_read_only()) {
        switch (rank) {
          case 1: repo.h_views[i].v1d = f.get_view<Real*>();      break;
          case 2: repo.h_views[i].v2d = f.get_view<Real**>();     break;
          case 3: repo.h_views[i].v3d = f.get_view<Real***>();    break;
          case 4: repo.h_views[i].v4d = f.get_view<Real****>();   break;
          case 5: repo.h_views[i].v5d = f.get_view<Real*****>();  break;
        }
      }
    };
    for (int i=0; i<this->m_num_fields; ++i) {
      get_view(i,fields[i]);
    }
    Kokkos::deep_copy(repo.views,  repo.h_views);
    Kokkos::deep_copy(repo.cviews, repo.h_cviews);
  }

  auto h_layout              = Kokkos::create_mirror_view(m_layout);
  auto h_pack_alloc_property = Kokkos::create_mirror_view(m_pack_alloc_property);
  auto h_num_levels          = Kokkos::create_mirror_view(m_num_levels);

  // Some info that is the same for both dyn and phys
  for (int i=0; i<this->m_num_fields; ++i) {
    const auto& phys = m_phys_fields[i];
    const auto& dyn  = m_dyn_fields[i];

    const auto& ph = phys.get_header();
    const auto& dh = dyn.get_header();

    // A dynamic subfield will need some special treatment at runtime
    // Namely, we'll need to re-extract the view every time,
    // since the subview info may have changed
    if (ph.get_parent().lock()) {
      EKAT_REQUIRE_MSG (ph.get_parent().lock()->get_parent().lock()==nullptr,
          "Error! We do not support remapping of subfields of other subfields.\n");
      m_subfield_info_phys[i] = ph.get_alloc_properties().get_subview_info();
    }
    if (dh.get_parent().lock()) {
      EKAT_REQUIRE_MSG (dh.get_parent().lock()->get_parent().lock()==nullptr,
          "Error! We do not support remapping of subfields of other subfields.\n");
      m_subfield_info_dyn[i] = dh.get_alloc_properties().get_subview_info();
    }

    const auto& pl = ph.get_identifier().get_layout();

    auto lt = get_layout_type(pl.tags());
    h_layout(i) = etoi(lt);

    const bool is_field_3d = lt==LayoutType::Scalar3D || lt==LayoutType::Vector3D;
    h_num_levels(i) = is_field_3d ? pl.dims().back() : -1;

    const auto& pap = ph.get_alloc_properties();
    const auto& dap = dh.get_alloc_properties();
    if (is_field_3d &&
        pap.template is_compatible<pack_type>() &&
        dap.template is_compatible<pack_type>()) {
      h_pack_alloc_property(i) = AllocPropType::PackAlloc;
    } else if (is_field_3d &&
               pap.template is_compatible<small_pack_type>() &&
               dap.template is_compatible<small_pack_type>()) {
      h_pack_alloc_property(i) = AllocPropType::SmallPackAlloc;
    } else {
      h_pack_alloc_property(i) = AllocPropType::RealAlloc;
    }
  }
  Kokkos::deep_copy(m_layout,              h_layout             );
  Kokkos::deep_copy(m_pack_alloc_property, h_pack_alloc_property);
  Kokkos::deep_copy(m_num_levels,          h_num_levels         );
}

bool PhysicsDynamicsRemapper::
subfields_info_has_changed (const std::map<int,SubviewInfo>& subfield_info,
                            const std::vector<field_type>& fields) const
{
  for (const auto& it : subfield_info) {
    const auto& f = fields[it.first];
    const auto& info_old = it.second;
    const auto& info_new = f.get_header().get_alloc_properties().get_subview_info();
    if ( not(info_old==info_new) ) {
      return true;
    }
  }
  return false;
}

void PhysicsDynamicsRemapper::
update_subfields_views (const std::map<int,SubviewInfo>& subfield_info,
                        const ViewsRepo& repo,
                        const std::vector<field_type>& fields) const
{
  auto get_view = [&] (const int i, const Field& f) {
    const auto rank = f.get_header().get_identifier().get_layout().rank();
    switch (rank) {
      case 1: repo.h_cviews[i].v1d = f.get_view<const Real*>();      break;
      case 2: repo.h_cviews[i].v2d = f.get_view<const Real**>();     break;
      case 3: repo.h_cviews[i].v3d = f.get_view<const Real***>();    break;
      case 4: repo.h_cviews[i].v4d = f.get_view<const Real****>();   break;
      case 5: repo.h_cviews[i].v5d = f.get_view<const Real*****>();  break;
    }
    if (not f.is_read_only()) {
      switch (rank) {
        case 1: repo.h_views[i].v1d = f.get_view<Real*>();      break;
        case 2: repo.h_views[i].v2d = f.get_view<Real**>();     break;
        case 3: repo.h_views[i].v3d = f.get_view<Real***>();    break;
        case 4: repo.h_views[i].v4d = f.get_view<Real****>();   break;
        case 5: repo.h_views[i].v5d = f.get_view<Real*****>();  break;
      }
    }
  };

  for (const auto& it : subfield_info) {
    const auto& f = fields[it.first];
    if ( not(it.second==f.get_header().get_alloc_properties().get_subview_info()) ){
      get_view(it.first,fields[it.first]);
    }
  }
  Kokkos::deep_copy(repo.views,  repo.h_views);
  Kokkos::deep_copy(repo.cviews, repo.h_cviews);
}

template <typename ScalarT, typename MT>
KOKKOS_FUNCTION
void PhysicsDynamicsRemapper::
set_dyn_to_zero(const MT& team) const
{
  const int i = team.league_rank();

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    {
      auto v = m_dyn_repo.views[i].v3d;
      int dim1 = v.extent(1);
      int dim2 = v.extent(2);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, v.size()), [&](const int& k) {
        int k0   = (k / dim2) / dim1;
        int k1   = (k / dim2) % dim1;
        int k2   =  k % dim2;
        v(k0, k1, k2) = 0;
      });
      break;
    }
    case etoi(LayoutType::Vector2D): // Fallthrough
    case etoi(LayoutType::Scalar3D):
    {
      auto v = m_dyn_repo.views[i].v4d;
      int dim1 = v.extent(1);
      int dim2 = v.extent(2);
      int dim3 = v.extent(3);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, v.size()), [&](const int& k) {
        int k0   = ((k / dim3) / dim2) / dim1;
        int k1   = ((k / dim3) / dim2) % dim1;
        int k2   =  (k / dim3) % dim2;
        int k3   =   k % dim3;
        v(k0, k1, k2, k3) = 0;
      });
      break;
    }
    case etoi(LayoutType::Vector3D):
    {
      auto v = m_dyn_repo.views[i].v5d;
      int dim1 = v.extent(1);
      int dim2 = v.extent(2);
      int dim3 = v.extent(3);
      int dim4 = v.extent(4);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, v.size()), [&](const int& k) {
        int k0   = (((k / dim4) / dim3) / dim2) / dim1;
        int k1   = (((k / dim4) / dim3) / dim2) % dim1;
        int k2   =  ((k / dim4) / dim3) % dim2;
        int k3   =   (k / dim4) % dim3;
        int k4   =    k % dim4;
        v(k0, k1, k2, k3, k4) = 0;
      });
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

void PhysicsDynamicsRemapper::
do_remap_fwd()
{
  // When remapping from phys to dyn, we need to perform a BEX
  // on the dyn fields. The BEX is 'static', meaning that the stored
  // views cannot be changed after it's been setup. Therefore,
  // if subfield info has changed for some dyn fields, we're toast.
  EKAT_REQUIRE_MSG (
      not subfields_info_has_changed(m_subfield_info_dyn,m_dyn_fields),
      "Error! P->D remapping is not supported if the some of the dyn fields\n"
      "       are subfields whose subview info has changed since setup.\n"
      "       Note: see field.hpp and field_alloc_prop.hpp for an explanation\n"
      "       of what a subfield and subview info).\n");

  // Check if we need to update the views for subfields on phys grid
  update_subfields_views(m_subfield_info_phys,m_phys_repo,m_phys_fields);

  using TeamPolicy = typename KT::TeamTagPolicy<RemapFwdTag>;

  const auto concurrency = KT::ExeSpace::concurrency();
#ifdef KOKKOS_ENABLE_CUDA
#ifdef KOKKOS_ENABLE_DEBUG
  const int team_size = std::min(256, std::min(128*m_num_phys_cols,32*(concurrency/this->m_num_fields+31)/32));
#else
  const int team_size = std::min(1024, std::min(128*m_num_phys_cols,32*(concurrency/this->m_num_fields+31)/32));
#endif
#endif

#ifdef KOKKOS_ENABLE_HIP
  const int team_size = std::min(256, std::min(128*m_num_phys_cols,32*(concurrency/this->m_num_fields+31)/32));
#endif

//should exclude above cases of CUDA and HIP
#ifndef EAMXX_ENABLE_GPU
  const int team_size = (concurrency<this->m_num_fields ? 1 : concurrency/this->m_num_fields);
#endif


  // TeamPolicy over this->m_num_fields
  const TeamPolicy policy(this->m_num_fields,team_size);
  Kokkos::parallel_for(policy, *this);
  Kokkos::fence();

  // Exchange element halo
  m_be->exchange();
}

void PhysicsDynamicsRemapper::
do_remap_bwd()
{
  // Check if we need to update the views for subfields
  update_subfields_views(m_subfield_info_dyn,m_dyn_repo,m_dyn_fields);
  update_subfields_views(m_subfield_info_phys,m_phys_repo,m_phys_fields);

  using TeamPolicy = typename KT::TeamTagPolicy<RemapBwdTag>;

  const auto concurrency = KT::ExeSpace::concurrency();
#ifdef KOKKOS_ENABLE_CUDA
  const int num_levs  = m_phys_grid->get_num_vertical_levels();
  const int team_size = std::min(128,32*(int)ceil(((Real)num_levs)/32));
#else
  const int team_size = (concurrency<this->m_num_fields*m_num_phys_cols ? 1 : concurrency/(this->m_num_fields*m_num_phys_cols));
#endif

  // TeamPolicy over m_num_phys_cols*this->m_num_fields. Unlike do_remap_fwd,
  // here we do not require setting dyn=0, allowing us to extend
  // the TeamPolicy
  const TeamPolicy policy(this->m_num_fields*m_num_phys_cols,team_size);
  Kokkos::parallel_for(policy, *this);
  Kokkos::fence();
}

void PhysicsDynamicsRemapper::
setup_boundary_exchange () {

  auto& c = Homme::Context::singleton();

  using Scalar = Homme::Scalar;

  int num_2d = 0;
  int num_3d_mid = 0;
  int num_3d_int = 0;
  for (int i=0; i<this->m_num_fields; ++i) {
    const auto& layout = m_dyn_fields[i].get_header().get_identifier().get_layout();
    const auto lt = get_layout_type(layout.tags());
    switch (lt) {
      case LayoutType::Scalar2D:
        ++num_2d;
        break;
      case LayoutType::Vector2D:
        num_2d += layout.dim(1);
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
        // A vector field (not a state): remap all components
        if (layout.dims().back()==HOMMEXX_NUM_PHYSICAL_LEV) {
          num_3d_mid += layout.dim(1);
        } else if (layout.dims().back()==HOMMEXX_NUM_INTERFACE_LEV) {
          num_3d_int += layout.dim(1);
        } else {
          EKAT_ERROR_MSG ("Error! Unexpected vertical level extent.\n");
        }
        break;
    default:
      EKAT_ERROR_MSG("Error! Invalid layout. This is an internal error. Please, contact developers\n");
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
  m_be = std::make_shared<Homme::BoundaryExchange>(conn,bm);
  m_be->set_num_fields(0,num_2d,num_3d_mid,num_3d_int);

  // If some fields are already bound, set them in the bd exchange
  for (int i=0; i<this->m_num_fields; ++i) {
    const auto& layout = m_dyn_fields[i].get_header().get_identifier().get_layout();
    const auto& dims = layout.dims();
    const auto lt = get_layout_type(layout.tags());
    switch (lt) {
      case LayoutType::Scalar2D:
        m_be->register_field(getHommeView<Real*[NP][NP]>(m_dyn_fields[i]));
        break;
      case LayoutType::Vector2D:
        m_be->register_field(getHommeView<Real**[NP][NP]>(m_dyn_fields[i]),dims[1],0);
        break;
      case LayoutType::Scalar3D:
        if (dims.back()==HOMMEXX_NUM_PHYSICAL_LEV) {
          m_be->register_field(getHommeView<Scalar*[NP][NP][NLEV]>(m_dyn_fields[i]));
        } else {
          m_be->register_field(getHommeView<Scalar*[NP][NP][NINT]>(m_dyn_fields[i]));
        }
        break;
      case LayoutType::Vector3D:
        if (dims.back()==HOMMEXX_NUM_PHYSICAL_LEV) {
          m_be->register_field(getHommeView<Scalar**[NP][NP][NLEV]>(m_dyn_fields[i]),dims[1],0);
        } else {
          m_be->register_field(getHommeView<Scalar**[NP][NP][NINT]>(m_dyn_fields[i]),dims[1],0);
        }
        break;
    default:
      EKAT_ERROR_MSG("Error! Invalid layout. This is an internal error. Please, contact developers\n");
    }
  }
  m_be->registration_completed();
}

template <typename MT>
KOKKOS_FUNCTION
void PhysicsDynamicsRemapper::
local_remap_fwd_2d (const MT& team) const
{
  const int i = team.league_rank();

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    {
      auto phys = m_phys_repo.cviews[i].v1d;
      auto dyn  = m_dyn_repo.views[i].v3d;

      const auto tr = Kokkos::TeamVectorRange(team, m_num_phys_cols);
      const auto f = [&] (const int icol) {
        const auto& elgp = Kokkos::subview(m_lid2elgp,m_p2d(icol),Kokkos::ALL());
        dyn(elgp[0],elgp[1],elgp[2]) = phys(icol);
      };
      Kokkos::parallel_for(tr, f);
      break;
    }
    case etoi(LayoutType::Vector2D):
    {
      auto phys = m_phys_repo.cviews[i].v2d;
      auto dyn  = m_dyn_repo.views[i].v4d;

      const int vec_dim = phys.extent(1);
      const auto tr = Kokkos::TeamVectorRange(team, m_num_phys_cols*vec_dim);
      const auto f = [&] (const int idx) {
        const int icol = idx / vec_dim;
        const int idim = idx % vec_dim;

        const auto& elgp = Kokkos::subview(m_lid2elgp,m_p2d(icol),Kokkos::ALL());
        dyn(elgp[0],idim,elgp[1],elgp[2]) = phys(icol,idim);
      };
      Kokkos::parallel_for(tr, f);
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template <typename ScalarT, typename MT>
KOKKOS_FUNCTION
void PhysicsDynamicsRemapper::
local_remap_fwd_3d (const MT& team) const
{
  const int i = team.league_rank();

  constexpr int PackSize = sizeof(ScalarT) / sizeof(Real);
  using PI = ekat::PackInfo<PackSize>;
  const int num_packs = PI::num_packs(m_num_levels(i));

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar3D):
    {
      auto phys = pack_view<const ScalarT>(m_phys_repo.cviews[i].v2d);
      auto dyn  = pack_view<      ScalarT>(m_dyn_repo.views[i].v4d);

      const auto tr = Kokkos::TeamVectorRange(team, m_num_phys_cols*num_packs);
      const auto f = [&] (const int idx) {
        const int icol = idx / num_packs;
        const int ilev = idx % num_packs;

        const auto& elgp = Kokkos::subview(m_lid2elgp,m_p2d(icol),Kokkos::ALL());
        dyn(elgp[0],elgp[1],elgp[2],ilev) = phys(icol,ilev);
      };
      Kokkos::parallel_for(tr, f);
      break;
    }
    case etoi(LayoutType::Vector3D):
    {
      auto phys = pack_view<const ScalarT>(m_phys_repo.cviews[i].v3d);
      auto dyn  = pack_view<      ScalarT>(m_dyn_repo.views[i].v5d);
      const int vec_dim = phys.extent(1);

      const auto tr = Kokkos::TeamVectorRange(team, m_num_phys_cols*vec_dim*num_packs);
      const auto f = [&] (const int idx) {
        const int icol = (idx / num_packs) / vec_dim;
        const int idim = (idx / num_packs) % vec_dim;
        const int ilev =  idx % num_packs;

        const auto& elgp = Kokkos::subview(m_lid2elgp,m_p2d(icol),Kokkos::ALL());
        dyn(elgp[0],idim,elgp[1],elgp[2],ilev) = phys(icol,idim,ilev);
      };
      Kokkos::parallel_for(tr, f);
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template <typename MT>
KOKKOS_FUNCTION
void PhysicsDynamicsRemapper::
local_remap_bwd_2d (const MT& team) const
{
  const int rank = team.league_rank();
  const int i    = rank % this->m_num_fields;
  const int icol = rank / this->m_num_fields;
  const auto& elgp = Kokkos::subview(m_lid2elgp,m_p2d(icol),Kokkos::ALL());

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    {
      auto phys = m_phys_repo.views[i].v1d;
      auto dyn  = m_dyn_repo.cviews[i].v3d;

      phys(icol) = dyn(elgp[0],elgp[1],elgp[2]);
      break;
    }
    case etoi(LayoutType::Vector2D):
    {
      auto phys = m_phys_repo.views[i].v2d;
      auto dyn  = m_dyn_repo.cviews[i].v4d;
      const int vec_dim = phys.extent(1);

      const auto tr = Kokkos::TeamVectorRange(team, vec_dim);
      const auto f = [&] (const int idim) {
        phys(icol,idim) = dyn(elgp[0],idim,elgp[1],elgp[2]);
      };
      Kokkos::parallel_for(tr, f);
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template <typename ScalarT, typename MT>
KOKKOS_FUNCTION
void PhysicsDynamicsRemapper::
local_remap_bwd_3d (const MT& team) const
{
  const int rank = team.league_rank();
  const int i    = rank % this->m_num_fields;
  const int icol = rank / this->m_num_fields;
  const auto& elgp = Kokkos::subview(m_lid2elgp,m_p2d(icol),Kokkos::ALL());

  constexpr int PackSize = sizeof(ScalarT) / sizeof(Real);
  using PI = ekat::PackInfo<PackSize>;
  const int num_packs = PI::num_packs(m_num_levels(i));

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar3D):
    {
      auto phys = pack_view<      ScalarT>(m_phys_repo.views[i].v2d);
      auto dyn  = pack_view<const ScalarT>(m_dyn_repo.cviews[i].v4d);

      const auto tr = Kokkos::TeamVectorRange(team, num_packs);
      const auto f = [&] (const int ilev) {
        phys(icol,ilev) = dyn(elgp[0],elgp[1],elgp[2],ilev);
      };
      Kokkos::parallel_for(tr, f);
      break;
    }
    case etoi(LayoutType::Vector3D):
    {
      auto phys = pack_view<      ScalarT>(m_phys_repo.views[i].v3d);
      auto dyn  = pack_view<const ScalarT>(m_dyn_repo.cviews[i].v5d);
      const int vec_dim = phys.extent(1);

      const auto tr = Kokkos::TeamVectorRange(team, vec_dim*num_packs);
      const auto f = [&] (const int idx) {
        const int idim = idx / num_packs;
        const int ilev = idx % num_packs;
        phys(icol,idim,ilev) = dyn(elgp[0],idim,elgp[1],elgp[2],ilev);
      };
      Kokkos::parallel_for(tr, f);
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

void PhysicsDynamicsRemapper::
create_p2d_map () {
  auto num_phys_dofs = m_phys_grid->get_num_local_dofs();
  auto num_dyn_dofs  = m_dyn_grid->get_num_local_dofs();

  auto se_dyn = std::dynamic_pointer_cast<const SEGrid>(m_dyn_grid);
  EKAT_REQUIRE_MSG(se_dyn, "Error! Something went wrong casting dyn grid to a SEGrid.\n");

  using gid_t = AbstractGrid::gid_type;
  auto dyn_gids  = se_dyn->get_cg_dofs_gids().get_view<const gid_t*>();
  auto phys_gids = m_phys_grid->get_dofs_gids().get_view<const gid_t*>();

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
    (void)found;
  });
}

template<typename MT>
KOKKOS_INLINE_FUNCTION
void PhysicsDynamicsRemapper::
operator()(const RemapFwdTag&, const MT& team) const
{
  const int i = team.league_rank();

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    case etoi(LayoutType::Vector2D):
      set_dyn_to_zero<Real>(team);
      team.team_barrier();

      local_remap_fwd_2d(team);
      break;
    case etoi(LayoutType::Scalar3D):
    case etoi(LayoutType::Vector3D):
      if (m_pack_alloc_property(i) == AllocPropType::PackAlloc) {
        set_dyn_to_zero<pack_type>(team);
        team.team_barrier();

        local_remap_fwd_3d<pack_type>(team);
      } else if (m_pack_alloc_property(i) == AllocPropType::SmallPackAlloc) {
        set_dyn_to_zero<small_pack_type>(team);
        team.team_barrier();

        local_remap_fwd_3d<small_pack_type>(team);
      } else {
        set_dyn_to_zero<Real>(team);
        team.team_barrier();

        local_remap_fwd_3d<Real>(team);
      }
      break;
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template<typename MT>
KOKKOS_INLINE_FUNCTION
void PhysicsDynamicsRemapper::
operator()(const RemapBwdTag&, const MT& team) const
{
  const int rank = team.league_rank();
  const int i = rank % this->m_num_fields;

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    case etoi(LayoutType::Vector2D):
      local_remap_bwd_2d(team);
      break;
    case etoi(LayoutType::Scalar3D):
    case etoi(LayoutType::Vector3D):
      if (m_pack_alloc_property(i) == AllocPropType::PackAlloc) {
        local_remap_bwd_3d<pack_type>(team);
      } else if (m_pack_alloc_property(i) == AllocPropType::SmallPackAlloc) {
        local_remap_bwd_3d<small_pack_type>(team);
      } else {
        local_remap_bwd_3d<Real>(team);
      }
      break;
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

} // namespace scream
