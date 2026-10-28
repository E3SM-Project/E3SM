#include "share/remap/unique_grid_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/util/eamxx_timing.hpp"

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
  // We don't want fields with TimeLevel in it. Just subview the fields instead
  return AbstractRemapper::is_valid_src_layout(layout) and
         not ekat::contains(layout.tags(),FieldTag::TimeLevel);
}

void UniqueGridRemapper::
registration_ends_impl ()
{
  initialize_device_variables();
}

void UniqueGridRemapper::
initialize_device_variables ()
{
  m_src_repo.init(m_src_fields);
  m_tgt_repo.init(m_tgt_fields);

  m_layout              = decltype(m_layout)              ("layout", m_num_fields);
  m_pack_alloc_property = decltype(m_pack_alloc_property) ("pack_alloc_property", m_num_fields);
  m_num_levels          = decltype(m_num_levels)          ("num_physical_levels", m_num_fields);
  m_vec_dim             = decltype(m_vec_dim)             ("vec_dim", m_num_fields);

  auto h_layout              = Kokkos::create_mirror_view(m_layout);
  auto h_pack_alloc_property = Kokkos::create_mirror_view(m_pack_alloc_property);
  auto h_num_levels          = Kokkos::create_mirror_view(m_num_levels);
  auto h_vec_dim             = Kokkos::create_mirror_view(m_vec_dim);

  for (int i=0; i<m_num_fields; ++i) {
    const auto& sh = m_src_fields[i].get_header();
    const auto& th = m_tgt_fields[i].get_header();

    const auto& sl = sh.get_identifier().get_layout();

    auto lt = sl.type();
    h_layout(i) = etoi(lt);

    const bool is_field_3d = lt==LayoutType::Scalar3D || lt==LayoutType::Vector3D;
    h_num_levels(i) = is_field_3d ? sl.dims().back() : -1;
    h_vec_dim(i)    = sl.is_vector_layout() ? sl.get_vector_dim() : -1;

    const auto& sap = sh.get_alloc_properties();
    const auto& tap = th.get_alloc_properties();
    if (is_field_3d &&
        sap.template is_compatible<pack_type>() &&
        tap.template is_compatible<pack_type>()) {
      h_pack_alloc_property(i) = AllocPropType::PackAlloc;
    } else {
      h_pack_alloc_property(i) = AllocPropType::RealAlloc;
    }
  }
  Kokkos::deep_copy(m_layout,              h_layout             );
  Kokkos::deep_copy(m_pack_alloc_property, h_pack_alloc_property);
  Kokkos::deep_copy(m_num_levels,          h_num_levels         );
  Kokkos::deep_copy(m_vec_dim,             h_vec_dim            );
}

void UniqueGridRemapper::
remap_fwd_impl ()
{
  // Check if we need to update the views for subfields
  m_src_repo.update_subfields(m_src_fields);
  m_tgt_repo.update_subfields(m_tgt_fields);

  // TeamPolicy over m_num_unique_cols*m_num_fields
  using TeamPolicy = typename KT::TeamTagPolicy<RemapTag>;
  const TeamPolicy policy(m_num_fields*m_num_unique_cols,Kokkos::AUTO);
  Kokkos::parallel_for(policy, *this);
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

  // Check if we need to update the views for subfields
  m_src_repo.update_subfields(m_src_fields);
  m_tgt_repo.update_subfields(m_tgt_fields);

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
  m_f_col_size   = view_1d<int>("f_col_size",m_num_fields);
  m_f_col_offset = view_1d<int>("f_col_offset",m_num_fields);
  auto f_col_size_h   = Kokkos::create_mirror_view(m_f_col_size);
  auto f_col_offset_h = Kokkos::create_mirror_view(m_f_col_offset);
  m_total_col_size = 0;
  for (int i=0; i<m_num_fields; ++i) {
    const auto& fl = m_tgt_fields[i].get_header().get_identifier().get_layout();
    f_col_size_h(i)   = fl.clone().strip_dim(COL).size();
    f_col_offset_h(i) = m_total_col_size;
    m_total_col_size += f_col_size_h(i);
  }
  Kokkos::deep_copy(m_f_col_size,f_col_size_h);
  Kokkos::deep_copy(m_f_col_offset,f_col_offset_h);

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

  using TeamPolicy = typename KT::TeamTagPolicy<PackTag>;
  const TeamPolicy policy(m_num_fields*m_num_exports,Kokkos::AUTO);
  Kokkos::parallel_for(policy, *this);
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

  using TeamPolicy = typename KT::TeamTagPolicy<UnpackTag>;
  const TeamPolicy policy(m_num_fields*m_num_imports,Kokkos::AUTO);
  Kokkos::parallel_for(policy, *this);
  Kokkos::fence();

  if (m_timers_enabled)
    stop_timer(name()+" unpack");
}

// NOTE: in the kernels below, the source field is addressed with the dof 'lid' (if m_src_K==1),
//       or with the three indices (e,g1,g2) in m_src_lid2idx (if m_src_K==3). In the latter case,
//       the indices are in the same order as the horizontal tags in the grid's 2d scalar layout.
//       Any other dim in the layout is located by the layout type. E.g., for SE grids,
//       a vector dim sits between the first and second index.

template <typename MT>
KOKKOS_FUNCTION
void UniqueGridRemapper::
local_remap_2d (const MT& team) const
{
  const int rank = team.league_rank();
  const int i    = rank % m_num_fields;
  const int icol = rank / m_num_fields;
  const int lid  = m_unique2src(icol);
  const bool K3  = m_src_K==3;

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    {
      auto u = m_tgt_repo.views[i].v1d;
      if (K3) {
        auto s = m_src_repo.cviews[i].v3d;
        u(icol) = s(m_src_lid2idx(lid,0),m_src_lid2idx(lid,1),m_src_lid2idx(lid,2));
      } else {
        auto s = m_src_repo.cviews[i].v1d;
        u(icol) = s(lid);
      }
      break;
    }
    case etoi(LayoutType::Vector2D):
    {
      auto u = m_tgt_repo.views[i].v2d;
      const int vec_dim = m_vec_dim(i);

      const auto tr = Kokkos::TeamVectorRange(team, vec_dim);
      if (K3) {
        auto s = m_src_repo.cviews[i].v4d;
        const int e = m_src_lid2idx(lid,0), g1 = m_src_lid2idx(lid,1), g2 = m_src_lid2idx(lid,2);
        Kokkos::parallel_for(tr, [&] (const int idim) {
          u(icol,idim) = s(e,idim,g1,g2);
        });
      } else {
        auto s = m_src_repo.cviews[i].v2d;
        Kokkos::parallel_for(tr, [&] (const int idim) {
          u(icol,idim) = s(lid,idim);
        });
      }
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template <typename ScalarT, typename MT>
KOKKOS_FUNCTION
void UniqueGridRemapper::
local_remap_3d (const MT& team) const
{
  const int rank = team.league_rank();
  const int i    = rank % m_num_fields;
  const int icol = rank / m_num_fields;
  const int lid  = m_unique2src(icol);
  const bool K3  = m_src_K==3;

  constexpr int PackSize = sizeof(ScalarT) / sizeof(Real);
  using PI = ekat::PackInfo<PackSize>;
  const int num_packs = PI::num_packs(m_num_levels(i));

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar3D):
    {
      auto u = pack_view<      ScalarT>(m_tgt_repo.views[i].v2d);

      const auto tr = Kokkos::TeamVectorRange(team, num_packs);
      if (K3) {
        auto s = pack_view<const ScalarT>(m_src_repo.cviews[i].v4d);
        const int e = m_src_lid2idx(lid,0), g1 = m_src_lid2idx(lid,1), g2 = m_src_lid2idx(lid,2);
        Kokkos::parallel_for(tr, [&] (const int ilev) {
          u(icol,ilev) = s(e,g1,g2,ilev);
        });
      } else {
        auto s = pack_view<const ScalarT>(m_src_repo.cviews[i].v2d);
        Kokkos::parallel_for(tr, [&] (const int ilev) {
          u(icol,ilev) = s(lid,ilev);
        });
      }
      break;
    }
    case etoi(LayoutType::Vector3D):
    {
      auto u = pack_view<      ScalarT>(m_tgt_repo.views[i].v3d);
      const int vec_dim = m_vec_dim(i);

      const auto tr = Kokkos::TeamVectorRange(team, vec_dim*num_packs);
      if (K3) {
        auto s = pack_view<const ScalarT>(m_src_repo.cviews[i].v5d);
        const int e = m_src_lid2idx(lid,0), g1 = m_src_lid2idx(lid,1), g2 = m_src_lid2idx(lid,2);
        Kokkos::parallel_for(tr, [&] (const int idx) {
          const int idim = idx / num_packs;
          const int ilev = idx % num_packs;
          u(icol,idim,ilev) = s(e,idim,g1,g2,ilev);
        });
      } else {
        auto s = pack_view<const ScalarT>(m_src_repo.cviews[i].v3d);
        Kokkos::parallel_for(tr, [&] (const int idx) {
          const int idim = idx / num_packs;
          const int ilev = idx % num_packs;
          u(icol,idim,ilev) = s(lid,idim,ilev);
        });
      }
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

// Note: for each column, the data of a field is stored as
//   - Scalar2D: [1 entry]
//   - Vector2D: [cmp]
//   - Scalar3D: [lev]
//   - Vector3D: [cmp][lev]
// All fields on the same column are spliced together, and all columns for the same pid are contiguous.
// Within a pid segment, the data is grouped by field.
template<typename MT>
KOKKOS_INLINE_FUNCTION
void UniqueGridRemapper::
operator()(const PackTag&, const MT& team) const
{
  const int rank = team.league_rank();
  const int i    = rank % m_num_fields;
  const int idx  = rank / m_num_fields;

  const int pid        = m_exp_pids(idx);
  const int pid_offset = m_exp_offsets(pid);
  const int ncols_pid  = m_exp_offsets(pid+1) - pid_offset;
  const int icol       = m_exp_lids(idx);

  const int col_size = m_f_col_size(i);
  const int offset = pid_offset*m_total_col_size
                   + ncols_pid*m_f_col_offset(i)
                   + (idx-pid_offset)*col_size;

  auto buf = m_send_buffer;
  const auto tr = Kokkos::TeamVectorRange(team,col_size);
  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    {
      auto cg = m_tgt_repo.cviews[i].v1d;
      Kokkos::parallel_for(tr,[&](const int k) { buf(offset+k) = cg(icol); });
      break;
    }
    case etoi(LayoutType::Vector2D):
    case etoi(LayoutType::Scalar3D):
    {
      auto cg = m_tgt_repo.cviews[i].v2d;
      Kokkos::parallel_for(tr,[&](const int k) { buf(offset+k) = cg(icol,k); });
      break;
    }
    case etoi(LayoutType::Vector3D):
    {
      auto cg = m_tgt_repo.cviews[i].v3d;
      const int nlev = m_num_levels(i);
      Kokkos::parallel_for(tr,[&](const int k) { buf(offset+k) = cg(icol,k/nlev,k%nlev); });
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template<typename MT>
KOKKOS_INLINE_FUNCTION
void UniqueGridRemapper::
operator()(const UnpackTag&, const MT& team) const
{
  const int rank = team.league_rank();
  const int i    = rank % m_num_fields;
  const int idx  = rank / m_num_fields;

  const int pid        = m_imp_pids(idx);
  const int pid_offset = m_imp_offsets(pid);
  const int ncols_pid  = m_imp_offsets(pid+1) - pid_offset;
  const int u          = m_imp_lids(idx);

  const int col_size = m_f_col_size(i);
  const int offset = pid_offset*m_total_col_size
                   + ncols_pid*m_f_col_offset(i)
                   + (idx-pid_offset)*col_size;

  // The same unique column may correspond to multiple source dofs
  const int beg = m_u2src_offsets(u);
  const int end = m_u2src_offsets(u+1);
  const bool K3 = m_src_K==3;

  auto buf = m_recv_buffer;
  const auto tr = Kokkos::TeamVectorRange(team,col_size);
  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    {
      if (K3) {
        auto s = m_src_repo.views[i].v3d;
        Kokkos::parallel_for(tr,[&](const int k) {
          for (int d=beg; d<end; ++d) {
            const int dof = m_u2src_dofs(d);
            s(m_src_lid2idx(dof,0),m_src_lid2idx(dof,1),m_src_lid2idx(dof,2)) = buf(offset+k);
          }
        });
      } else {
        auto s = m_src_repo.views[i].v1d;
        Kokkos::parallel_for(tr,[&](const int k) {
          for (int d=beg; d<end; ++d) {
            s(m_u2src_dofs(d)) = buf(offset+k);
          }
        });
      }
      break;
    }
    case etoi(LayoutType::Vector2D):
    {
      if (K3) {
        auto s = m_src_repo.views[i].v4d;
        Kokkos::parallel_for(tr,[&](const int k) {
          for (int d=beg; d<end; ++d) {
            const int dof = m_u2src_dofs(d);
            s(m_src_lid2idx(dof,0),k,m_src_lid2idx(dof,1),m_src_lid2idx(dof,2)) = buf(offset+k);
          }
        });
      } else {
        auto s = m_src_repo.views[i].v2d;
        Kokkos::parallel_for(tr,[&](const int k) {
          for (int d=beg; d<end; ++d) {
            s(m_u2src_dofs(d),k) = buf(offset+k);
          }
        });
      }
      break;
    }
    case etoi(LayoutType::Scalar3D):
    {
      if (K3) {
        auto s = m_src_repo.views[i].v4d;
        Kokkos::parallel_for(tr,[&](const int k) {
          for (int d=beg; d<end; ++d) {
            const int dof = m_u2src_dofs(d);
            s(m_src_lid2idx(dof,0),m_src_lid2idx(dof,1),m_src_lid2idx(dof,2),k) = buf(offset+k);
          }
        });
      } else {
        auto s = m_src_repo.views[i].v2d;
        Kokkos::parallel_for(tr,[&](const int k) {
          for (int d=beg; d<end; ++d) {
            s(m_u2src_dofs(d),k) = buf(offset+k);
          }
        });
      }
      break;
    }
    case etoi(LayoutType::Vector3D):
    {
      const int nlev = m_num_levels(i);
      if (K3) {
        auto s = m_src_repo.views[i].v5d;
        Kokkos::parallel_for(tr,[&](const int k) {
          for (int d=beg; d<end; ++d) {
            const int dof = m_u2src_dofs(d);
            s(m_src_lid2idx(dof,0),k/nlev,m_src_lid2idx(dof,1),m_src_lid2idx(dof,2),k%nlev) = buf(offset+k);
          }
        });
      } else {
        auto s = m_src_repo.views[i].v3d;
        Kokkos::parallel_for(tr,[&](const int k) {
          for (int d=beg; d<end; ++d) {
            s(m_u2src_dofs(d),k/nlev,k%nlev) = buf(offset+k);
          }
        });
      }
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template<typename MT>
KOKKOS_INLINE_FUNCTION
void UniqueGridRemapper::
operator()(const RemapTag&, const MT& team) const
{
  const int rank = team.league_rank();
  const int i = rank % m_num_fields;

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    case etoi(LayoutType::Vector2D):
      local_remap_2d(team);
      break;
    case etoi(LayoutType::Scalar3D):
    case etoi(LayoutType::Vector3D):
      if (m_pack_alloc_property(i) == AllocPropType::PackAlloc) {
        local_remap_3d<pack_type>(team);
      } else {
        local_remap_3d<Real>(team);
      }
      break;
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

} // namespace scream
