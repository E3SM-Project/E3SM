#include "share/remap/unique_grid_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/util/eamxx_timing.hpp"
#include "share/field/field_visit.hpp"

#include <ekat_pack_utils.hpp>
#include <ekat_assert.hpp>

#include <algorithm>
#include <iostream>
#include <limits>
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

  // One entry for each (field,component). All the kernels process all entries at once
  m_entry_field.clear();
  m_entry_comp.clear();
  for (int i=0; i<m_num_fields; ++i) {
    for (int c=0; c<m_ncomp[i]; ++c) {
      m_entry_field.push_back(i);
      m_entry_comp.push_back(c);
    }
  }
  m_num_entries = m_entry_field.size();

  m_entries   = KT::view_ND<UniqueGridEntryDesc,1>("entries",m_num_entries);
  m_entries_h = Kokkos::create_mirror_view(m_entries);
  for (int e=0; e<m_num_entries; ++e) {
    const auto& sl = m_tgt_scalar[m_entry_field[e]].get_header().get_identifier().get_layout();
    m_entries_h(e) = UniqueGridEntryDesc();
    m_entries_h(e).nlev = sl.rank()==2 ? sl.dim(1) : 1;
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

using ExeSpace   = Field::device_t::execution_space;
using TeamPolicy = Kokkos::TeamPolicy<ExeSpace>;
using Member     = typename TeamPolicy::member_type;

template<typename T>
using kview_1d = UniqueGridRemapper::KT::view_ND<T,1>;

// Get the data pointer and the strides of the (rank<=4) view of a field. The views are extracted
// at every call, so that any (dynamic) subview info is up to date.
// NOTE: we always ask for a view with const value type, which any field can provide.
//       Writing through the pointer is the responsibility of the caller.
template<int N>
void get_ptr_and_strides (const Field& f, Real*& ptr, long long* strides)
{
  std::visit([&](const auto& v) {
    ptr = const_cast<Real*>(v.data());
    for (int i=0; i<N; ++i) {
      strides[i] = v.stride(i);
    }
  }, details::get_view_variant<N,const Real>(f));
}

void get_ptr_and_strides (const Field& f, Real*& ptr, long long* strides)
{
  switch (f.rank()) {
    case 1: get_ptr_and_strides<1>(f,ptr,strides); break;
    case 2: get_ptr_and_strides<2>(f,ptr,strides); break;
    case 3: get_ptr_and_strides<3>(f,ptr,strides); break;
    case 4: get_ptr_and_strides<4>(f,ptr,strides); break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in UniqueGridRemapper.\n"
          " - field: " + f.name() + "\n"
          " - rank : " + std::to_string(f.rank()) + "\n");
  }
}

// The offset of the source grid dof 'dof', for the entry described by d
template<int K>
KOKKOS_INLINE_FUNCTION
long long src_dof_offset (const UniqueGridEntryDesc& d,
                          const Field::view_dev_t<const int**>& lid2idx,
                          const int dof)
{
  if constexpr (K==1) {
    return dof*d.src_s[0];
  } else {
    return static_cast<long long>(lid2idx(dof,0))*d.src_s[0]
         + static_cast<long long>(lid2idx(dof,1))*d.src_s[1]
         + static_cast<long long>(lid2idx(dof,2))*d.src_s[2];
  }
}

// Forward kernel: for each entry, copy the unique grid column icol from one of its replicas on the source grid
//   tgt(icol[,ilev]) = src(<indices of the source dof>[,ilev])
// The source field is addressed with the dof itself (if K==1), or with the three indices
// (e,g1,g2) in the lid->idx map (if K==3). The latter are in the same order as the horizontal
// tags in the grid's 2d scalar layout.
template<int K>
struct UniqueGridFwdKernel {
  KOKKOS_INLINE_FUNCTION
  void operator() (const Member& team) const {
    const int rank = team.league_rank();
    const int e    = rank % num_entries;
    const int icol = rank / num_entries;

    const auto& d = entries(e);
    const int dof = unique2src(icol);

    Real* t = d.tgt_ptr + icol*d.tgt_s[0];
    const Real* s = d.src_ptr + src_dof_offset<K>(d,lid2idx,dof);
    const long long ts = d.tgt_s[1];
    const long long ss = d.src_s[3];
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,d.nlev),
                         [&](const int k) {
      t[k*ts] = s[k*ss];
    });
  }

  int                                      num_entries;
  kview_1d<UniqueGridEntryDesc>             entries;
  kview_1d<int>                            unique2src;
  Field::view_dev_t<const int**>           lid2idx;
};

// Pack kernel: copy a column of the unique grid fields (one team for each (entry,export)) into the send buffer.
// All fields on the same column are spliced together, and all columns for the same pid are contiguous.
// Within a pid segment, the data is grouped by field, and within a field, by column and component.
struct UniqueGridPackKernel {
  KOKKOS_INLINE_FUNCTION
  void operator() (const Member& team) const {
    const int rank = team.league_rank();
    const int e    = rank % num_entries;
    const int idx  = rank / num_entries;

    const auto& d = entries(e);
    const int pid        = exp_pids(idx);
    const int pid_offset = exp_offsets(pid);
    const int ncols_pid  = exp_offsets(pid+1) - pid_offset;
    const int icol       = exp_lids(idx);
    const long long offset = static_cast<long long>(pid_offset)*total_col_size
                           + static_cast<long long>(ncols_pid)*d.field_off
                           + static_cast<long long>(idx-pid_offset)*d.col_size
                           + d.comp_off;

    const Real* t = d.tgt_ptr + icol*d.tgt_s[0];
    const long long ts = d.tgt_s[1];
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,d.nlev),
                         [&](const int k) {
      buf(offset+k) = t[k*ts];
    });
  }

  int                            num_entries;
  int                            total_col_size;
  kview_1d<UniqueGridEntryDesc>   entries;
  kview_1d<int>                  exp_pids, exp_lids, exp_offsets;
  kview_1d<Real>                 buf;
};

// Unpack kernel: copy a column from the recv buffer (one team for each (entry,import)) to all the
// source dofs with such gid. The imports processed are
//    idx = first + j,   with j in [0,count),  and then idx += skip_cnt if idx>=skip_beg
// (this allows to process the imports from this rank, and the other ones, separately)
template<int K>
struct UniqueGridUnpackKernel {
  KOKKOS_INLINE_FUNCTION
  void operator() (const Member& team) const {
    const int rank = team.league_rank();
    const int e    = rank % num_entries;
    int idx = first + rank / num_entries;
    if (idx>=skip_beg) {
      idx += skip_cnt;
    }

    const auto& d = entries(e);
    const int pid        = imp_pids(idx);
    const int pid_offset = imp_offsets(pid);
    const int ncols_pid  = imp_offsets(pid+1) - pid_offset;
    const int u          = imp_lids(idx);
    const long long offset = static_cast<long long>(pid_offset)*total_col_size
                           + static_cast<long long>(ncols_pid)*d.field_off
                           + static_cast<long long>(idx-pid_offset)*d.col_size
                           + d.comp_off;

    // The same unique column may correspond to multiple source dofs
    const int beg = u2src_offsets(u);
    const int end = u2src_offsets(u+1);
    const long long ss = d.src_s[3];
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,d.nlev),
                         [&](const int k) {
      const Real val = buf(offset+k);
      for (int dd=beg; dd<end; ++dd) {
        Real* s = d.src_ptr + src_dof_offset<K>(d,lid2idx,u2src_dofs(dd));
        s[k*ss] = val;
      }
    });
  }

  int                            num_entries;
  int                            total_col_size;
  int                            first, skip_beg, skip_cnt;
  kview_1d<UniqueGridEntryDesc>   entries;
  kview_1d<int>                  imp_pids, imp_lids, imp_offsets;
  kview_1d<int>                  u2src_offsets, u2src_dofs;
  Field::view_dev_t<const int**> lid2idx;
  kview_1d<Real>                 buf;
};

} // anonymous namespace

void UniqueGridRemapper::
update_entries ()
{
  if (m_timers_enabled)
    start_timer(name()+" update_entries");

  for (int e=0; e<m_num_entries; ++e) {
    const int i = m_entry_field[e];
    const int c = m_entry_comp[e];
    set_component(i,c);

    auto& d = m_entries_h(e);
    long long s[4] = {0,0,0,0};

    // Unique grid field: (col[,lev])
    get_ptr_and_strides(m_tgt_scalar[i],d.tgt_ptr,s);
    d.tgt_s[0] = s[0];
    d.tgt_s[1] = m_tgt_scalar[i].rank()==2 ? s[1] : 0;

    // Source field: (dof[,lev]) or (e,g1,g2[,lev])
    const int rank = m_src_scalar[i].rank();
    get_ptr_and_strides(m_src_scalar[i],d.src_ptr,s);
    if (m_src_K==1) {
      d.src_s[0] = s[0];
      d.src_s[3] = rank==2 ? s[1] : 0;
    } else {
      d.src_s[0] = s[0];
      d.src_s[1] = s[1];
      d.src_s[2] = s[2];
      d.src_s[3] = rank==4 ? s[3] : 0;
    }
  }
  Kokkos::deep_copy(m_entries,m_entries_h);

  if (m_timers_enabled)
    stop_timer(name()+" update_entries");
}

void UniqueGridRemapper::
remap_fwd_impl ()
{
  if (m_num_entries==0) {
    return;
  }

  update_entries();

  if (m_timers_enabled)
    start_timer(name()+" fwd_kernel");

  const TeamPolicy policy(m_num_entries*m_num_unique_cols,Kokkos::AUTO);
  if (m_src_K==1) {
    UniqueGridFwdKernel<1> kernel {m_num_entries,m_entries,m_unique2src,m_src_lid2idx};
    Kokkos::parallel_for("UniqueGridRemapper::fwd",policy,kernel);
  } else {
    UniqueGridFwdKernel<3> kernel {m_num_entries,m_entries,m_unique2src,m_src_lid2idx};
    Kokkos::parallel_for("UniqueGridRemapper::fwd",policy,kernel);
  }

  // The kernel is async, but remap is blocking
  ExeSpace().fence();

  if (m_timers_enabled)
    stop_timer(name()+" fwd_kernel");
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

  update_entries();

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
  const int myrank = comm.rank();

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

  auto imp_offsets_h = imp_exp.import_pid_offsets_h();
  auto exp_offsets_h = imp_exp.export_pid_offsets_h();

  // The columns owned by this rank (that we also need) are not sent through MPI
  m_self_imp_beg = imp_offsets_h(myrank);
  m_self_exp_beg = exp_offsets_h(myrank);
  m_self_cnt     = imp_offsets_h(myrank+1) - imp_offsets_h(myrank);
  EKAT_REQUIRE_MSG (m_self_cnt==exp_offsets_h(myrank+1)-exp_offsets_h(myrank),
      "Error! The number of columns this rank imports from itself does not match the number it exports to itself.\n");

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

  // Where each entry sits in a column's chunk of the buffers
  for (int e=0; e<m_num_entries; ++e) {
    const int i = m_entry_field[e];
    const int c = m_entry_comp[e];
    auto& d = m_entries_h(e);
    d.field_off = m_f_col_offset[i];
    d.col_size  = m_f_col_size[i];
    d.comp_off  = c*d.nlev;
  }

  // Create the send/recv buffers
  m_send_buffer = view_1d<Real>("UniqueGridRemapper::send_buf",m_num_exports*m_total_col_size);
  m_recv_buffer = view_1d<Real>("UniqueGridRemapper::recv_buf",m_num_imports*m_total_col_size);
  m_mpi_send_buffer = Kokkos::create_mirror_view(typename decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(typename decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  // Create the persistent requests (none for this rank)
  m_requests = std::make_shared<MpiRequests>();
  const auto mpi_comm = comm.mpi_comm();
  const auto mpi_real = ekat::get_mpi_type<Real>();
  int num_send_pids = 0;
  int num_recv_pids = 0;
  for (int pid=0; pid<nranks; ++pid) {
    if (pid==myrank) {
      continue;
    }
    const int ncols_send = exp_offsets_h(pid+1)-exp_offsets_h(pid);
    if (ncols_send>0) {
      auto send_ptr = m_mpi_send_buffer.data() + exp_offsets_h(pid)*m_total_col_size;
      auto& req = m_requests->send.emplace_back();
      check_mpi_call(MPI_Send_init(send_ptr,ncols_send*m_total_col_size,mpi_real,pid,0,mpi_comm,&req),
                     "[UniqueGridRemapper] creating persistent send request.\n");
      ++num_send_pids;
    }
    const int ncols_recv = imp_offsets_h(pid+1)-imp_offsets_h(pid);
    if (ncols_recv>0) {
      auto recv_ptr = m_mpi_recv_buffer.data() + imp_offsets_h(pid)*m_total_col_size;
      auto& req = m_requests->recv.emplace_back();
      check_mpi_call(MPI_Recv_init(recv_ptr,ncols_recv*m_total_col_size,mpi_real,pid,0,mpi_comm,&req),
                     "[UniqueGridRemapper] creating persistent recv request.\n");
      ++num_recv_pids;
    }
  }

  // Diagnostics: how balanced is the communication pattern?
  if (m_timers_enabled) {
    auto report = [&](const std::string& what, const int v) {
      int mn = v, mx = v, sm = v;
      comm.all_reduce(&mn,1,MPI_MIN);
      comm.all_reduce(&mx,1,MPI_MAX);
      comm.all_reduce(&sm,1,MPI_SUM);
      if (comm.am_i_root()) {
        std::cout << "[UniqueGridRemapper] " << name() << ": " << what
                  << " (min/max/avg over ranks): " << mn << " / " << mx << " / " << double(sm)/nranks << "\n";
      }
    };
    report("remote columns imported", m_num_imports-m_self_cnt);
    report("remote columns exported", m_num_exports-m_self_cnt);
    report("local columns (no comm) ", m_self_cnt);
    report("neighbors we send to    ", num_send_pids);
    report("neighbors we recv from  ", num_recv_pids);
  }

  Kokkos::deep_copy(m_entries,m_entries_h);
  m_bwd_setup_done = true;
}

void UniqueGridRemapper::
launch_pack ()
{
  if (m_num_exports==0 or m_num_entries==0) {
    return;
  }

  UniqueGridPackKernel kernel {m_num_entries,m_total_col_size,m_entries,m_exp_pids,m_exp_lids,m_exp_offsets,m_send_buffer};
  const TeamPolicy policy(m_num_entries*m_num_exports,Kokkos::AUTO);
  Kokkos::parallel_for("UniqueGridRemapper::pack",policy,kernel);
}

void UniqueGridRemapper::
launch_unpack (const int first, const int count, const int skip_beg, const int skip_cnt)
{
  if (count<=0 or m_num_entries==0) {
    return;
  }

  const TeamPolicy policy(m_num_entries*count,Kokkos::AUTO);
  if (m_src_K==1) {
    UniqueGridUnpackKernel<1> kernel {m_num_entries,m_total_col_size,first,skip_beg,skip_cnt,m_entries,
                                      m_imp_pids,m_imp_lids,m_imp_offsets,m_u2src_offsets,m_u2src_dofs,
                                      m_src_lid2idx,m_recv_buffer};
    Kokkos::parallel_for("UniqueGridRemapper::unpack",policy,kernel);
  } else {
    UniqueGridUnpackKernel<3> kernel {m_num_entries,m_total_col_size,first,skip_beg,skip_cnt,m_entries,
                                      m_imp_pids,m_imp_lids,m_imp_offsets,m_u2src_offsets,m_u2src_dofs,
                                      m_src_lid2idx,m_recv_buffer};
    Kokkos::parallel_for("UniqueGridRemapper::unpack",policy,kernel);
  }
}

void UniqueGridRemapper::
pack_and_send ()
{
  if (m_timers_enabled)
    start_timer(name()+" pack");

  launch_pack();
  ExeSpace().fence();

  // The columns that this rank owns and also needs don't go through MPI: copy them to the recv buffer.
  if (m_self_cnt>0) {
    const size_t n   = static_cast<size_t>(m_self_cnt)*m_total_col_size;
    const size_t sb  = static_cast<size_t>(m_self_exp_beg)*m_total_col_size;
    const size_t rb  = static_cast<size_t>(m_self_imp_beg)*m_total_col_size;
    Kokkos::deep_copy(Kokkos::subview(m_recv_buffer,Kokkos::make_pair(rb,rb+n)),
                      Kokkos::subview(m_send_buffer,Kokkos::make_pair(sb,sb+n)));
  }

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
  // While the messages are in flight, unpack the columns that did not need to be sent
  if (m_timers_enabled)
    start_timer(name()+" unpack_local");
  launch_unpack(m_self_imp_beg,m_self_cnt,std::numeric_limits<int>::max(),0);
  if (m_timers_enabled)
    stop_timer(name()+" unpack_local");

  if (m_timers_enabled)
    start_timer(name()+" wait");
  if (not m_requests->recv.empty()) {
    check_mpi_call(MPI_Waitall(m_requests->recv.size(),m_requests->recv.data(), MPI_STATUSES_IGNORE),
                   "[UniqueGridRemapper] waiting on persistent recv requests.\n");
  }
  if (m_timers_enabled)
    stop_timer(name()+" wait");

  if (m_timers_enabled)
    start_timer(name()+" unpack");

  // If MPI does not use dev pointers, we need to deep copy from host to dev. We must not overwrite
  // the columns of this rank (they did not go through MPI), so we copy the two ranges around them.
  if (not MpiOnDev) {
    const size_t ntot = static_cast<size_t>(m_num_imports)*m_total_col_size;
    const size_t b    = static_cast<size_t>(m_self_imp_beg)*m_total_col_size;
    const size_t e    = b + static_cast<size_t>(m_self_cnt)*m_total_col_size;
    if (b>0) {
      Kokkos::deep_copy(Kokkos::subview(m_recv_buffer,Kokkos::make_pair(size_t(0),b)),
                        Kokkos::subview(m_mpi_recv_buffer,Kokkos::make_pair(size_t(0),b)));
    }
    if (e<ntot) {
      Kokkos::deep_copy(Kokkos::subview(m_recv_buffer,Kokkos::make_pair(e,ntot)),
                        Kokkos::subview(m_mpi_recv_buffer,Kokkos::make_pair(e,ntot)));
    }
  }

  launch_unpack(0,m_num_imports-m_self_cnt,m_self_imp_beg,m_self_cnt);
  ExeSpace().fence();

  if (m_timers_enabled)
    stop_timer(name()+" unpack");
}

} // namespace scream
