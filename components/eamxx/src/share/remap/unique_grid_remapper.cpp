#include "share/remap/unique_grid_remapper.hpp"

#include "share/util/eamxx_utils.hpp"
#include "share/util/eamxx_bound_functor.hpp"
#include "share/util/eamxx_parallel_for_nd.hpp"
#include "share/field/field_visit.hpp"

#include <ekat_pack_utils.hpp>
#include <ekat_assert.hpp>

namespace scream
{

KokkosTypes<DefaultDevice>::view_1d<int>
compute_unique_grid_dof_map (const AbstractGrid& grid,
                             const AbstractGrid& unique_grid)
{
  using gid_type = AbstractGrid::gid_type;

  const int num_unique_dofs = unique_grid.get_num_local_dofs();
  const int num_dofs        = grid.get_num_local_dofs();

  auto gids        = grid.get_unique_grid_gids().get_view<const gid_type*>();
  auto unique_gids = unique_grid.get_dofs_gids().get_view<const gid_type*>();

  auto policy = KokkosTypes<DefaultDevice>::RangePolicy(0,num_unique_dofs);
  KokkosTypes<DefaultDevice>::view_1d<int> u2g ("",num_unique_dofs);

  Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idof){
    auto gid = unique_gids(idof);
    bool found = false;
    for (int i=0; i<num_dofs; ++i) {
      if (gids(i)==gid) {
        u2g(idof) = i;
        found = true;
        break;
      }
    }
    EKAT_KERNEL_ASSERT_MSG (found, "Error! Unique grid gid not found among the source grid dofs.\n");
    (void)found;
  });

  return u2g;
}

UniqueGridRemapper::
UniqueGridRemapper (const grid_ptr_type& grid)
 : AbstractRemapper(grid,grid->get_unique_grid())
{
  const auto unique_grid = m_tgt_grid;

  EKAT_REQUIRE_MSG (grid!=unique_grid,
      "Error! The source grid has no replicated dofs (it is its own unique grid). Use an IdentityRemapper.\n"
      " - grid name: " + grid->name() + "\n");
  EKAT_REQUIRE_MSG (grid->get_num_vertical_levels()==unique_grid->get_num_vertical_levels(),
      "Error! The source grid and its unique grid have a different number of vertical levels.\n"
      " - grid name       : " + grid->name() + "\n"
      " - unique grid name: " + unique_grid->name() + "\n");

  // Fields are addressed with K indices per dof. We only support K=1 or K=3 on the source
  // grid (e.g., (col) or (elem,gp,gp)), and K=1 on the target grid.
  auto get_num_idx = [](const AbstractGrid& g) {
    return g.get_lid_to_idx_map().get_header().get_identifier().get_layout().dim(1);
  };
  m_src_K = get_num_idx(*grid);
  EKAT_REQUIRE_MSG (m_src_K==1 or m_src_K==3,
      "Error! Only grids whose dofs are addressed with 1 or 3 indices are supported as source.\n"
      " - grid name: " + grid->name() + "\n"
      " - number of indices per dof: " + std::to_string(m_src_K) + "\n");
  EKAT_REQUIRE_MSG (get_num_idx(*unique_grid)==1,
      "Error! The target grid must have dofs addressed with a single index.\n"
      " - grid name: " + unique_grid->name() + "\n"
      " - number of indices per dof: " + std::to_string(get_num_idx(*unique_grid)) + "\n");

  m_bwd_allowed = false;

  m_num_unique_cols = unique_grid->get_num_local_dofs();
  m_src_lid2idx     = grid->get_lid_to_idx_map().get_view<const int**>();
  m_unique2src      = compute_unique_grid_dof_map(*grid,*unique_grid);
}

bool UniqueGridRemapper::
is_valid_src_layout (const FieldLayout& layout) const
{
  // We only support the layouts handled by the remap kernels.
  // Also, we don't want fields with TimeLevel in it. Just subview the fields instead
  const auto lt = layout.type();
  const bool supported = lt==LayoutType::Scalar2D or lt==LayoutType::Vector2D or
                         lt==LayoutType::Scalar3D or lt==LayoutType::Vector3D;
  return supported and AbstractRemapper::is_valid_src_layout(layout) and
         not ekat::contains(layout.tags(),FieldTag::TimeLevel);
}

void UniqueGridRemapper::
registration_ends_impl ()
{
  // Vector fields are remapped one component at a time. So, for each of them, we store
  // a (dynamic) subfield for its components, and re-target it at every remap. This way, the
  // kernels only deal with scalar fields, whose layout is (col[,lev]) on the unique grid
  // (so the source grid is addressed with the dof indices, followed by [lev]).
  // NOTE: if a vector field is itself a (dynamic) subfield, nothing special is needed: the views
  //       of a subfield are extracted by recursively subviewing the parent, at every call.
  m_is_vec.clear();
  m_ncomp.clear();
  m_src_scalar.clear();
  m_tgt_scalar.clear();
  for (int i=0; i<m_num_fields; ++i) {
    const auto& tl = m_tgt_fields[i].get_header().get_identifier().get_layout();
    const bool is_vec = tl.is_vector_layout();
    m_is_vec.push_back(is_vec);
    m_ncomp.push_back(is_vec ? tl.get_vector_dim() : 1);
    m_src_scalar.push_back(is_vec ? m_src_fields[i].get_component(0,true) : m_src_fields[i]);
    m_tgt_scalar.push_back(is_vec ? m_tgt_fields[i].get_component(0,true) : m_tgt_fields[i]);
  }
}

void UniqueGridRemapper::
set_component (const int i, const int c)
{
  if (m_is_vec[i]) {
    m_src_scalar[i].get_header().get_alloc_properties().reset_subview_idx(c);
    m_tgt_scalar[i].get_header().get_alloc_properties().reset_subview_idx(c);
  }
}

namespace {

// Copy kernel: tgt(icol[,ilev]) = src(<indices of the source dof>[,ilev])
//
// The source field is addressed with the dof itself (if K==1), or with the three indices
// (e,g1,g2) in the lid->idx map (if K==3). The latter are in the same order as the horizontal
// tags in the grid's 2d scalar layout.
//
// NOTE: kernels going through parallel_for_nd must be functors with a templated operator()
template<int K>
struct UniqueGridCopyKernel {
  template<typename TgtView, typename SrcView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const TgtView& tgt, const SrcView& src, const int icol, const Idx... rest) const {
    const int dof = unique2src(icol);
    if constexpr (K==1) {
      tgt(icol,rest...) = src(dof,rest...);
    } else {
      tgt(icol,rest...) = src(lid2idx(dof,0),lid2idx(dof,1),lid2idx(dof,2),rest...);
    }
  }

  UniqueGridRemapper::KT::view_ND<int,1>  unique2src;
  Field::view_dev_t<const int**>          lid2idx;
};

// Launch the copy kernel for the scalar fields. N is the rank of the target field
template<int K, int N>
void remap_scalar_fwd (const Field& src, const Field& tgt,
                       const UniqueGridCopyKernel<K>& kernel)
{
  using details::get_view_variant;
  constexpr int NS = K==1 ? N : N+2;

  const auto& dims = tgt.get_header().get_identifier().get_layout().dims();
  std::visit([&](const auto& t, const auto& s) {
    using exec_space = typename std::decay_t<decltype(t)>::traits::execution_space;
    parallel_for_nd<N,exec_space>(dims,bind_functor(kernel,t,s));
  }, get_view_variant<N,Real>(tgt), get_view_variant<NS,const Real>(src));
}

template<int K>
void remap_scalar_fwd (const Field& src, const Field& tgt,
                       const UniqueGridCopyKernel<K>& kernel)
{
  switch (tgt.rank()) {
    case 1: remap_scalar_fwd<K,1>(src,tgt,kernel); break;
    case 2: remap_scalar_fwd<K,2>(src,tgt,kernel); break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in UniqueGridRemapper.\n"
          " - field: " + tgt.name() + "\n"
          " - rank : " + std::to_string(tgt.rank()) + "\n");
  }
}

} // anonymous namespace

void UniqueGridRemapper::
remap_fwd_impl ()
{
  for (int i=0; i<m_num_fields; ++i) {
    for (int c=0; c<m_ncomp[i]; ++c) {
      set_component(i,c);
      if (m_src_K==1) {
        remap_scalar_fwd<1>(m_src_scalar[i],m_tgt_scalar[i],UniqueGridCopyKernel<1>{m_unique2src,m_src_lid2idx});
      } else {
        remap_scalar_fwd<3>(m_src_scalar[i],m_tgt_scalar[i],UniqueGridCopyKernel<3>{m_unique2src,m_src_lid2idx});
      }
    }
  }

  // The kernels are async, but remap is blocking
  Kokkos::fence();
}

} // namespace scream
