#include "share/remap/unique_grid_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/util/eamxx_timing.hpp"
#include "share/util/eamxx_bound_functor.hpp"
#include "share/util/eamxx_parallel_for_nd.hpp"
#include "share/field/field_visit.hpp"

#include <ekat_pack_utils.hpp>
#include <ekat_assert.hpp>

#include <algorithm>
#include <map>

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

  m_num_unique_cols = unique_grid->get_num_local_dofs();
  m_src_lid2idx     = grid->get_lid_to_idx_map().get_view<const int**>();
  m_unique2src      = compute_unique_grid_dof_map(*grid,*unique_grid);
}

UniqueGridRemapper::MpiRequests::
~MpiRequests ()
{
  // Free MPI requests (unless MPI is already finalized, in which case there's nothing to do)
  int finalized;
  MPI_Finalized(&finalized);
  if (not finalized) {
    for (auto& req : send) {
      MPI_Request_free(&req);
    }
    for (auto& req : recv) {
      MPI_Request_free(&req);
    }
  }
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


// Pack kernel: copy entry k of column icol of a (scalar) unique grid field into the send buffer.
// All fields on the same column are spliced together, and all columns for the same pid are contiguous.
// Within a pid segment, the data is grouped by field, and within a field, by component.
struct UniqueGridPackKernel {
  template<typename TgtView>
  KOKKOS_INLINE_FUNCTION
  void operator() (const TgtView& tgt, const int idx, const int k) const {
    const int pid        = exp_pids(idx);
    const int pid_offset = exp_offsets(pid);
    const int ncols_pid  = exp_offsets(pid+1) - pid_offset;
    const int icol       = exp_lids(idx);
    const int offset = pid_offset*total_col_size
                     + ncols_pid*field_col_offset
                     + (idx-pid_offset)*col_size
                     + comp*comp_size;

    if constexpr (TgtView::rank()==1)
      buf(offset+k) = tgt(icol);    // k==0
    else
      buf(offset+k) = tgt(icol,k);
  }

  UniqueGridRemapper::KT::view_ND<int,1>   exp_pids, exp_lids, exp_offsets;
  UniqueGridRemapper::KT::view_ND<Real,1>  buf;
  int  total_col_size, field_col_offset, col_size, comp, comp_size;
};

// Unpack kernel: copy entry k of column u of the recv buffer to all the source dofs with such gid
template<int K>
struct UniqueGridUnpackKernel {
  template<typename SrcView>
  KOKKOS_INLINE_FUNCTION
  void operator() (const SrcView& src, const int idx, const int k) const {
    const int pid        = imp_pids(idx);
    const int pid_offset = imp_offsets(pid);
    const int ncols_pid  = imp_offsets(pid+1) - pid_offset;
    const int u          = imp_lids(idx);
    const int offset = pid_offset*total_col_size
                     + ncols_pid*field_col_offset
                     + (idx-pid_offset)*col_size
                     + comp*comp_size;

    // The same unique column may correspond to multiple source dofs
    const Real val = buf(offset+k);
    constexpr int NIDX = K==1 ? 1 : 3;      // indices addressing the dof
    constexpr bool has_lev = SrcView::rank()>NIDX;
    for (int d=u2src_offsets(u); d<u2src_offsets(u+1); ++d) {
      const int dof = u2src_dofs(d);
      if constexpr (K==1) {
        if constexpr (has_lev) src(dof,k) = val; else src(dof) = val;
      } else {
        const int e = lid2idx(dof,0), g1 = lid2idx(dof,1), g2 = lid2idx(dof,2);
        if constexpr (has_lev) src(e,g1,g2,k) = val; else src(e,g1,g2) = val;
      }
    }
  }

  UniqueGridRemapper::KT::view_ND<int,1>   imp_pids, imp_lids, imp_offsets;
  UniqueGridRemapper::KT::view_ND<int,1>   u2src_offsets, u2src_dofs;
  Field::view_dev_t<const int**>           lid2idx;
  UniqueGridRemapper::KT::view_ND<Real,1>  buf;
  int  total_col_size, field_col_offset, col_size, comp, comp_size;
};

template<int K, int N>
void unpack_scalar_impl (const Field& src, const UniqueGridUnpackKernel<K>& kernel, const int num_imports)
{
  constexpr int NS = K==1 ? N : N+2;
  std::visit([&](const auto& s) {
    using exec_space = typename std::decay_t<decltype(s)>::traits::execution_space;
    parallel_for_nd<2,exec_space>({num_imports,kernel.comp_size},bind_functor(kernel,s));
  }, details::get_view_variant<NS,Real>(src));
}

template<int K>
void unpack_scalar_impl (const Field& src, const Field& tgt, const UniqueGridUnpackKernel<K>& kernel,
                         const int num_imports)
{
  switch (tgt.rank()) {
    case 1: unpack_scalar_impl<K,1>(src,kernel,num_imports); break;
    case 2: unpack_scalar_impl<K,2>(src,kernel,num_imports); break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in UniqueGridRemapper.
"
          " - field: " + tgt.name() + "
"
          " - rank : " + std::to_string(tgt.rank()) + "
");
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

void UniqueGridRemapper::
remap_bwd_impl ()
{
  if (m_num_fields==0) {
    return;
  }

  if (not m_bwd_setup_done) {
    setup_bwd();
  }

  // Post the recvs first, so that they are ready when sends are fired
  if (not m_requests->recv.empty()) {
    check_mpi_call(MPI_Startall(m_requests->recv.size(),m_requests->recv.data()),
                   "[UniqueGridRemapper] starting persistent recv requests.\n");
  }

  pack_and_send();
  recv_and_unpack();

  // Wait for all sends to be completed
  if (not m_requests->send.empty()) {
    check_mpi_call(MPI_Waitall(m_requests->send.size(),m_requests->send.data(), MPI_STATUSES_IGNORE),
                   "[UniqueGridRemapper] waiting on persistent send requests.\n");
  }
}

void UniqueGridRemapper::
setup_bwd ()
{
  using gid_type = AbstractGrid::gid_type;

  const auto& comm = m_src_grid->get_comm();
  const int nranks = comm.size();

  // Group the local source dofs by unique gid. Each unique gid is a "column" of the overlap grid
  auto unique_gids_f = m_src_grid->get_unique_grid_gids();
  unique_gids_f.sync_to_host();
  auto src_gids_h = unique_gids_f.get_view<const gid_type*,Host>();

  const int num_src_dofs = m_src_grid->get_num_local_dofs();
  std::map<gid_type,std::vector<int>> gid2dofs;
  for (int i=0; i<num_src_dofs; ++i) {
    gid2dofs[src_gids_h(i)].push_back(i);
  }
  const int num_ov_cols = gid2dofs.size();

  auto ov_grid = std::make_shared<PointGrid>("unique_grid_overlap",num_ov_cols,m_tgt_grid->get_num_vertical_levels(),comm);
  auto ov_gids = ov_grid->get_dofs_gids();
  auto ov_gids_h = ov_gids.get_view<gid_type*,Host>();

  m_u2src_offsets = view_1d<int>("u2src_offsets",num_ov_cols+1);
  m_u2src_dofs    = view_1d<int>("u2src_dofs",num_src_dofs);
  auto offsets_h = Kokkos::create_mirror_view(m_u2src_offsets);
  auto dofs_h    = Kokkos::create_mirror_view(m_u2src_dofs);
  int u = 0;
  int pos = 0;
  offsets_h(0) = 0;
  for (const auto& [gid,dofs] : gid2dofs) {
    ov_gids_h(u) = gid;
    for (auto dof : dofs) {
      dofs_h(pos++) = dof;
    }
    offsets_h(++u) = pos;
  }
  ov_gids.sync_to_dev();
  Kokkos::deep_copy(m_u2src_offsets,offsets_h);
  Kokkos::deep_copy(m_u2src_dofs,dofs_h);

  // Figure out who owns the unique grid columns we need. Every source gid must be found among the unique grid gids.
  GridImportExport imp_exp (m_tgt_grid,ov_grid);
  m_num_imports = imp_exp.num_imports();
  m_num_exports = imp_exp.num_exports();
  m_imp_pids    = imp_exp.import_pids();
  m_imp_lids    = imp_exp.import_lids();
  m_imp_offsets = imp_exp.import_pid_offsets();
  m_exp_pids    = imp_exp.export_pids();
  m_exp_lids    = imp_exp.export_lids();
  m_exp_offsets = imp_exp.export_pid_offsets();

  // Compute the size of a column of each field, as well as of all fields spliced together
  using namespace ShortFieldTagsNames;
  m_f_col_size.resize(m_num_fields);
  m_f_col_offset.resize(m_num_fields);
  m_total_col_size = 0;
  for (int i=0; i<m_num_fields; ++i) {
    const auto& fl = m_tgt_fields[i].get_header().get_identifier().get_layout();
    m_f_col_size[i]   = fl.clone().strip_dim(COL).size();
    m_f_col_offset[i] = m_total_col_size;
    m_total_col_size += m_f_col_size[i];
  }

  // Create the send/recv buffers
  m_send_buffer = view_1d<Real>("UniqueGridRemapper::send_buf",m_num_exports*m_total_col_size);
  m_recv_buffer = view_1d<Real>("UniqueGridRemapper::recv_buf",m_num_imports*m_total_col_size);
  m_mpi_send_buffer = Kokkos::create_mirror_view(typename decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(typename decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  // Create the persistent requests
  m_requests = std::make_shared<MpiRequests>();
  auto imp_offsets_h = imp_exp.import_pid_offsets_h();
  auto exp_offsets_h = imp_exp.export_pid_offsets_h();
  const auto mpi_comm = comm.mpi_comm();
  const auto mpi_real = ekat::get_mpi_type<Real>();
  for (int pid=0; pid<nranks; ++pid) {
    const int ncols_send = exp_offsets_h(pid+1)-exp_offsets_h(pid);
    if (ncols_send>0) {
      auto send_ptr = m_mpi_send_buffer.data() + exp_offsets_h(pid)*m_total_col_size;
      auto& req = m_requests->send.emplace_back();
      check_mpi_call(MPI_Send_init(send_ptr,ncols_send*m_total_col_size,mpi_real,pid,0,mpi_comm,&req),
                     "[UniqueGridRemapper] creating persistent send request.\n");
    }
    const int ncols_recv = imp_offsets_h(pid+1)-imp_offsets_h(pid);
    if (ncols_recv>0) {
      auto recv_ptr = m_mpi_recv_buffer.data() + imp_offsets_h(pid)*m_total_col_size;
      auto& req = m_requests->recv.emplace_back();
      check_mpi_call(MPI_Recv_init(recv_ptr,ncols_recv*m_total_col_size,mpi_real,pid,0,mpi_comm,&req),
                     "[UniqueGridRemapper] creating persistent recv request.\n");
    }
  }

  m_bwd_setup_done = true;
}

void UniqueGridRemapper::
setup_bwd ()
{
  using gid_type = AbstractGrid::gid_type;

  const auto& comm = m_src_grid->get_comm();
  const int nranks = comm.size();

  // Group the local source dofs by unique gid. Each unique gid is a "column" of the overlap grid
  auto unique_gids_f = m_src_grid->get_unique_grid_gids();
  unique_gids_f.sync_to_host();
  auto src_gids_h = unique_gids_f.get_view<const gid_type*,Host>();

  const int num_src_dofs = m_src_grid->get_num_local_dofs();
  std::map<gid_type,std::vector<int>> gid2dofs;
  for (int i=0; i<num_src_dofs; ++i) {
    gid2dofs[src_gids_h(i)].push_back(i);
  }
  const int num_ov_cols = gid2dofs.size();

  auto ov_grid = std::make_shared<PointGrid>("unique_grid_overlap",num_ov_cols,m_tgt_grid->get_num_vertical_levels(),comm);
  auto ov_gids = ov_grid->get_dofs_gids();
  auto ov_gids_h = ov_gids.get_view<gid_type*,Host>();

  m_u2src_offsets = view_1d<int>("u2src_offsets",num_ov_cols+1);
  m_u2src_dofs    = view_1d<int>("u2src_dofs",num_src_dofs);
  auto offsets_h = Kokkos::create_mirror_view(m_u2src_offsets);
  auto dofs_h    = Kokkos::create_mirror_view(m_u2src_dofs);
  int u = 0;
  int pos = 0;
  offsets_h(0) = 0;
  for (const auto& [gid,dofs] : gid2dofs) {
    ov_gids_h(u) = gid;
    for (auto dof : dofs) {
      dofs_h(pos++) = dof;
    }
    offsets_h(++u) = pos;
  }
  ov_gids.sync_to_dev();
  Kokkos::deep_copy(m_u2src_offsets,offsets_h);
  Kokkos::deep_copy(m_u2src_dofs,dofs_h);

  // Figure out who owns the unique grid columns we need. Every source gid must be found among the unique grid gids.
  GridImportExport imp_exp (m_tgt_grid,ov_grid);
  m_num_imports = imp_exp.num_imports();
  m_num_exports = imp_exp.num_exports();
  m_imp_pids    = imp_exp.import_pids();
  m_imp_lids    = imp_exp.import_lids();
  m_imp_offsets = imp_exp.import_pid_offsets();
  m_exp_pids    = imp_exp.export_pids();
  m_exp_lids    = imp_exp.export_lids();
  m_exp_offsets = imp_exp.export_pid_offsets();

  // Compute the size of a column of each field, as well as of all fields spliced together
  using namespace ShortFieldTagsNames;
  m_f_col_size.resize(m_num_fields);
  m_f_col_offset.resize(m_num_fields);
  m_total_col_size = 0;
  for (int i=0; i<m_num_fields; ++i) {
    const auto& fl = m_tgt_fields[i].get_header().get_identifier().get_layout();
    m_f_col_size[i]   = fl.clone().strip_dim(COL).size();
    m_f_col_offset[i] = m_total_col_size;
    m_total_col_size += m_f_col_size[i];
  }

  // Create the send/recv buffers
  m_send_buffer = view_1d<Real>("UniqueGridRemapper::send_buf",m_num_exports*m_total_col_size);
  m_recv_buffer = view_1d<Real>("UniqueGridRemapper::recv_buf",m_num_imports*m_total_col_size);
  m_mpi_send_buffer = Kokkos::create_mirror_view(typename decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(typename decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  // Create the persistent requests
  m_requests = std::make_shared<MpiRequests>();
  auto imp_offsets_h = imp_exp.import_pid_offsets_h();
  auto exp_offsets_h = imp_exp.export_pid_offsets_h();
  const auto mpi_comm = comm.mpi_comm();
  const auto mpi_real = ekat::get_mpi_type<Real>();
  for (int pid=0; pid<nranks; ++pid) {
    const int ncols_send = exp_offsets_h(pid+1)-exp_offsets_h(pid);
    if (ncols_send>0) {
      auto send_ptr = m_mpi_send_buffer.data() + exp_offsets_h(pid)*m_total_col_size;
      auto& req = m_requests->send.emplace_back();
      check_mpi_call(MPI_Send_init(send_ptr,ncols_send*m_total_col_size,mpi_real,pid,0,mpi_comm,&req),
                     "[UniqueGridRemapper] creating persistent send request.\n");
    }
    const int ncols_recv = imp_offsets_h(pid+1)-imp_offsets_h(pid);
    if (ncols_recv>0) {
      auto recv_ptr = m_mpi_recv_buffer.data() + imp_offsets_h(pid)*m_total_col_size;
      auto& req = m_requests->recv.emplace_back();
      check_mpi_call(MPI_Recv_init(recv_ptr,ncols_recv*m_total_col_size,mpi_real,pid,0,mpi_comm,&req),
                     "[UniqueGridRemapper] creating persistent recv request.\n");
    }
  }

  m_bwd_setup_done = true;
}

void UniqueGridRemapper::
pack_and_send ()
{
  if (m_timers_enabled)
    start_timer(name()+" pack");

  for (int i=0; i<m_num_fields; ++i) {
    for (int c=0; c<m_ncomp[i]; ++c) {
      set_component(i,c);
      pack_scalar(i,c);
    }
  }
  Kokkos::fence();

  // If MPI does not use dev pointers, we need to deep copy from dev to host
  if (not MpiOnDev) {
    Kokkos::deep_copy(m_mpi_send_buffer,m_send_buffer);
  }

  if (m_timers_enabled)
    stop_timer(name()+" pack");

  if (not m_requests->send.empty()) {
    check_mpi_call(MPI_Startall(m_requests->send.size(),m_requests->send.data()),
                   "[UniqueGridRemapper] starting persistent send requests.\n");
  }
}

void UniqueGridRemapper::
recv_and_unpack ()
{
  if (not m_requests->recv.empty()) {
    check_mpi_call(MPI_Waitall(m_requests->recv.size(),m_requests->recv.data(), MPI_STATUSES_IGNORE),
                   "[UniqueGridRemapper] waiting on persistent recv requests.\n");
  }

  if (m_timers_enabled)
    start_timer(name()+" unpack");

  // If MPI does not use dev pointers, we need to deep copy from host to dev
  if (not MpiOnDev) {
    Kokkos::deep_copy(m_recv_buffer,m_mpi_recv_buffer);
  }

  for (int i=0; i<m_num_fields; ++i) {
    for (int c=0; c<m_ncomp[i]; ++c) {
      set_component(i,c);
      unpack_scalar(i,c);
    }
  }
  Kokkos::fence();

  if (m_timers_enabled)
    stop_timer(name()+" unpack");
}

void UniqueGridRemapper::
pack_scalar (const int i, const int c)
{
  if (m_num_exports==0)
    return;

  const auto& tgt = m_tgt_scalar[i];
  const int comp_size = m_f_col_size[i] / m_ncomp[i];
  UniqueGridPackKernel kernel {
    m_exp_pids, m_exp_lids, m_exp_offsets, m_send_buffer,
    m_total_col_size, m_f_col_offset[i], m_f_col_size[i], c, comp_size
  };

  auto launch = [&](auto N_c) {
    constexpr int N = decltype(N_c)::value;
    std::visit([&](const auto& t) {
      using exec_space = typename std::decay_t<decltype(t)>::traits::execution_space;
      parallel_for_nd<2,exec_space>({m_num_exports,comp_size},bind_functor(kernel,t));
    }, details::get_view_variant<N,const Real>(tgt));
  };
  switch (tgt.rank()) {
    case 1: launch(std::integral_constant<int,1>()); break;
    case 2: launch(std::integral_constant<int,2>()); break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in UniqueGridRemapper.
"
          " - field: " + tgt.name() + "
"
          " - rank : " + std::to_string(tgt.rank()) + "
");
  }
}

void UniqueGridRemapper::
unpack_scalar (const int i, const int c)
{
  if (m_num_imports==0)
    return;

  const int comp_size = m_f_col_size[i] / m_ncomp[i];
  if (m_src_K==1) {
    UniqueGridUnpackKernel<1> kernel {
      m_imp_pids, m_imp_lids, m_imp_offsets, m_u2src_offsets, m_u2src_dofs, m_src_lid2idx, m_recv_buffer,
      m_total_col_size, m_f_col_offset[i], m_f_col_size[i], c, comp_size
    };
    unpack_scalar_impl<1>(m_src_scalar[i],m_tgt_scalar[i],kernel,m_num_imports);
  } else {
    UniqueGridUnpackKernel<3> kernel {
      m_imp_pids, m_imp_lids, m_imp_offsets, m_u2src_offsets, m_u2src_dofs, m_src_lid2idx, m_recv_buffer,
      m_total_col_size, m_f_col_offset[i], m_f_col_size[i], c, comp_size
    };
    unpack_scalar_impl<3>(m_src_scalar[i],m_tgt_scalar[i],kernel,m_num_imports);
  }
}

} // namespace scream
