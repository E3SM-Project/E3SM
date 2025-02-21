#include "refining_remapper_p2p.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/util/eamxx_utils.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>

#include <numeric>

namespace scream
{

RefiningRemapperP2P::
RefiningRemapperP2P (const grid_ptr_type& tgt_grid,
                     const std::string& map_file)
 : HorizInterpRemapperBase(tgt_grid,map_file,InterpType::Refine)
{
  // Nothing to do here
}

RefiningRemapperP2P::
~RefiningRemapperP2P ()
{
  clean_up();
}

void RefiningRemapperP2P::remap_fwd_impl ()
{
  // Fire the recv requests right away, so that if some other ranks
  // is done packing before us, we can start receiving their data
  if (not m_recv_req.empty()) {
    check_mpi_call(MPI_Startall(m_recv_req.size(),m_recv_req.data()),
                   "[RefiningRemapperP2P] starting persistent recv requests.\n");
  }

  // Do P2P communications
  pack_and_send ();
  recv_and_unpack ();

  // Perform local-mat vec
  // Helpef function, to establish if a field can be handled with packs
  auto can_pack_field = [](const Field& f) {
    const auto& ap = f.get_header().get_alloc_properties();
    return (ap.get_last_extent() % SCREAM_PACK_SIZE) == 0;
  };

  // Loop over each field, perform mat-vec
  constexpr auto COL = ShortFieldTagsNames::COL;
  for (int i=0; i<m_num_fields; ++i) {
    auto& f_tgt = m_tgt_fields[i];

    // It's ok to register fields that do not have the COL tag
    // These fields are simply copied from src to tgt.
    if (not f_tgt.get_header().get_identifier().get_layout().has_tag(COL)) {
      f_tgt.deep_copy(m_src_fields[i]);
      continue;
    }

    // Perform the local mat-vec. Recall that in these y=Ax products,
    // x is the overlapped src field, and y is the tgt field.
    const auto& f_ov = m_ov_fields[i];

    // If possible, dispatch kernel with SCREAM_PACK_SIZE
    if (can_pack_field(f_ov) and can_pack_field(f_tgt)) {
      local_mat_vec<SCREAM_PACK_SIZE>(f_ov,f_tgt);
    } else {
      local_mat_vec<1>(f_ov,f_tgt);
    }
  }

  // Wait for all sends to be completed
  if (not m_send_req.empty()) {
    check_mpi_call(MPI_Waitall(m_send_req.size(),m_send_req.data(), MPI_STATUSES_IGNORE),
                   "[RefiningRemapperP2P] waiting on persistent send requests.\n");
  }
}

void RefiningRemapperP2P::setup_mpi_data_structures ()
{
  using namespace ShortFieldTagsNames;

  const int nranks = m_comm.size();

  // Get cumulative col size of each field (to be used to compute offsets)
  m_fields_col_sizes_scan_sum.resize(m_num_fields+1,0);
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f = m_src_fields[i];
    const auto& fl = f.get_header().get_identifier().get_layout();

    // Fields without COL tag are nor remapped, so consider their col size as 0
    auto col_size = fl.has_tag(COL) ? fl.clone().strip_dim(COL).size() : 0;
    m_fields_col_sizes_scan_sum[i+1] = m_fields_col_sizes_scan_sum[i] + col_size;
  }
  auto total_col_size = m_fields_col_sizes_scan_sum.back();

  // ----------- Compute RECV metadata -------------- //

  // Figure out where ov_src cols are received from
  const int ncols_recv = m_ov_coarse_grid->get_num_local_dofs();

  m_imp_exp = std::make_shared<GridImportExport>(m_src_grid,m_ov_coarse_grid);

  // We can now compute the offset of each pid in the recv buffer
  m_pids_recv_offsets = view_1d<int>("",nranks+1);
  auto ncols_recv_h = m_imp_exp->num_imports_per_pid_h();
  auto pids_recv_offsets_h = Kokkos::create_mirror_view(m_pids_recv_offsets);
  pids_recv_offsets_h[0] = 0;
  for (int pid=0; pid<nranks; ++pid) {
    pids_recv_offsets_h(pid+1) = pids_recv_offsets_h(pid)
                               + ncols_recv_h(pid);
  }
  Kokkos::deep_copy(m_pids_recv_offsets,pids_recv_offsets_h);

  // Create the recv buffer(s)
  auto recv_buf_size = ncols_recv*total_col_size;
  m_recv_buffer = decltype(m_recv_buffer)("RefiningRemapperP2P::recv_buf",recv_buf_size);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  // ----------- Compute SEND metadata -------------- //

  m_pids_send_offsets = view_1d<int>("",nranks+1);
  auto pids_send_offsets_h = Kokkos::create_mirror_view(m_pids_send_offsets);
  auto ncols_send_h = m_imp_exp->num_exports_per_pid_h();
  pids_send_offsets_h[0] = 0;
  for (int pid=0; pid<nranks; ++pid) {
    pids_send_offsets_h(pid+1) = pids_send_offsets_h(pid)
                               + ncols_send_h(pid);
  }
  Kokkos::deep_copy(m_pids_send_offsets,pids_send_offsets_h);

  // Create the send buffer(s)
  auto send_buf_size = pids_send_offsets_h(nranks)*total_col_size;
  m_send_buffer = decltype(m_send_buffer)("RefiningRemapperP2P::send_buf",send_buf_size);
  m_mpi_send_buffer = Kokkos::create_mirror_view(decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);

  // ----------- Create Requests ------------ //

  const auto mpi_comm = m_comm.mpi_comm();
  const auto mpi_real = ekat::get_mpi_type<Real>();
  for (int pid=0; pid<nranks; ++pid) {
    // Send request
    if (ncols_send_h(pid)>0) {
      auto send_ptr = m_mpi_send_buffer.data() + pids_send_offsets_h(pid)*total_col_size;
      auto send_count = ncols_send_h(pid)*total_col_size;
      auto& req = m_send_req.emplace_back();
      MPI_Send_init (send_ptr, send_count, mpi_real, pid,
                     0, mpi_comm, &req);
    }
    // Recv request
    if (ncols_recv_h(pid)>0) {
      auto recv_ptr = m_mpi_recv_buffer.data() + pids_recv_offsets_h(pid)*total_col_size;
      auto recv_count = ncols_recv_h(pid)*total_col_size;
      auto& req = m_recv_req.emplace_back();
      MPI_Recv_init (recv_ptr, recv_count, mpi_real, pid,
                     0, mpi_comm, &req);
    }
  }
}

void RefiningRemapperP2P::pack_and_send ()
{
  using RangePolicy = typename KT::RangePolicy;
  using TeamMember  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  constexpr auto COL = ShortFieldTagsNames::COL;

  auto export_pids = m_imp_exp->export_pids();
  auto export_lids = m_imp_exp->export_lids();
  auto ncols_send  = m_imp_exp->num_exports_per_pid();
  auto pids_send_offsets = m_pids_send_offsets;
  auto send_buf = m_send_buffer;
  const int num_exports = export_pids.size();
  const int total_col_size = m_fields_col_sizes_scan_sum.back();
  for (int ifield=0; ifield<m_num_fields; ++ifield) {
    const auto& f = m_src_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    if (not fl.has_tag(COL)) {
      // No need to process this field. We'll deep copy src->tgt later
      continue;
    }
    const auto f_col_sizes_scan_sum = m_fields_col_sizes_scan_sum[ifield];
    switch (fl.rank()) {
      case 1:
      {
        const auto v = f.get_strided_view<const Real*>();
        auto pack = KOKKOS_LAMBDA(const int iexp) {
          auto pid = export_pids(iexp);
          auto icol = export_lids(iexp);
          auto pid_offset = pids_send_offsets(pid);
          auto pos_within_pid = iexp - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_send(pid)*f_col_sizes_scan_sum
                      + pos_within_pid;
          send_buf(offset) = v(icol);
        };
        Kokkos::parallel_for(RangePolicy(0,num_exports),pack);
        break;
      }
      case 2:
      {
        const auto v = f.get_view<const Real**>();
        const int dim1 = fl.dim(1);
        auto policy = ESU::get_default_team_policy(num_exports,dim1);
        auto pack = KOKKOS_LAMBDA(const TeamMember& team) {
          const int iexp = team.league_rank();
          const int icol = export_lids(iexp);
          const int pid  = export_pids(iexp);
          auto pid_offset = pids_send_offsets(pid);
          auto pos_within_pid = iexp - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_send(pid)*f_col_sizes_scan_sum
                      + pos_within_pid*dim1;
          auto col_pack = [&](const int& k) {
            send_buf(offset+k) = v(icol,k);
          };
          auto tvr = Kokkos::TeamVectorRange(team,dim1);
          Kokkos::parallel_for(tvr,col_pack);
        };
        Kokkos::parallel_for(policy,pack);
        break;
      }
      case 3:
      {
        const auto v = f.get_view<const Real***>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        const int f_col_size = dim1*dim2;
        auto policy = ESU::get_default_team_policy(num_exports,dim1*dim2);
        auto pack = KOKKOS_LAMBDA(const TeamMember& team) {
          const int iexp = team.league_rank();
          const int icol = export_lids(iexp);
          const int pid  = export_pids(iexp);
          auto pid_offset = pids_send_offsets(pid);
          auto pos_within_pid = iexp - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_send(pid)*f_col_sizes_scan_sum
                      + pos_within_pid*f_col_size;
          auto col_pack = [&](const int& idx) {
            const int j = idx / dim2;
            const int k = idx % dim2;
            send_buf(offset+idx) = v(icol,j,k);
          };
          auto tvr = Kokkos::TeamVectorRange(team,f_col_size);
          Kokkos::parallel_for(tvr,col_pack);
        };
        Kokkos::parallel_for(policy,pack);
        break;
      }
      case 4:
      {
        const auto v = f.get_view<const Real****>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        const int dim3 = fl.dim(3);
        const int f_col_size = dim1*dim2*dim3;
        auto policy = ESU::get_default_team_policy(num_exports,dim1*dim2*dim3);
        auto pack = KOKKOS_LAMBDA(const TeamMember& team) {
          const int iexp = team.league_rank();
          const int icol = export_lids(iexp);
          const int pid  = export_pids(iexp);
          auto pid_offset = pids_send_offsets(pid);
          auto pos_within_pid = iexp - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_send(pid)*f_col_sizes_scan_sum
                      + pos_within_pid*f_col_size;
          auto col_pack = [&](const int& idx) {
            const int j = (idx / dim3) / dim2;
            const int k = (idx / dim3) % dim2;
            const int l =  idx % dim3;
            send_buf(offset+idx) = v(icol,j,k,l);
          };
          auto tvr = Kokkos::TeamVectorRange(team,f_col_size);
          Kokkos::parallel_for(tvr,col_pack);
        };
        Kokkos::parallel_for(policy,pack);
        break;
      }
      default:
        EKAT_ERROR_MSG ("Unexpected field rank in RefiningRemapperP2P::pack.\n"
            "  - MPI rank  : " + std::to_string(m_comm.rank()) + "\n"
            "  - field name: " + f.name() + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }

  // Wait for all threads to be done packing
  Kokkos::fence();

  // If MPI does not use dev pointers, we need to deep copy from dev to host
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_mpi_send_buffer,m_send_buffer);
  }

  if (not m_send_req.empty()) {
    check_mpi_call(MPI_Startall(m_send_req.size(),m_send_req.data()),
                   "[RefiningRemapperP2P] start persistent send requests.\n");
  }
}

void RefiningRemapperP2P::recv_and_unpack ()
{
  if (not m_recv_req.empty()) {
    check_mpi_call(MPI_Waitall(m_recv_req.size(),m_recv_req.data(), MPI_STATUSES_IGNORE),
                   "[RefiningRemapperP2P] waiting on persistent recv requests.\n");
  }
  // If MPI does not use dev pointers, we need to deep copy from host to dev
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_recv_buffer,m_mpi_recv_buffer);
  }

  using RangePolicy = typename KT::RangePolicy;
  using TeamMember  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  constexpr auto COL = ShortFieldTagsNames::COL;

  auto import_pids = m_imp_exp->import_pids();
  auto import_lids = m_imp_exp->import_lids();
  auto ncols_recv  = m_imp_exp->num_imports_per_pid();
  auto pids_recv_offsets = m_pids_recv_offsets;
  auto recv_buf = m_recv_buffer;
  const int num_imports = import_pids.size();
  const int total_col_size = m_fields_col_sizes_scan_sum.back();
  for (int ifield=0; ifield<m_num_fields; ++ifield) {
          auto& f  = m_ov_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    if (not fl.has_tag(COL)) {
      // No need to process this field. We'll deep copy src->tgt later
      continue;
    }
    const auto f_col_sizes_scan_sum = m_fields_col_sizes_scan_sum[ifield];
    switch (fl.rank()) {
      case 1:
      {
        auto v = f.get_view<Real*>();
        auto unpack = KOKKOS_LAMBDA (const int idx) {
          const int pid  = import_pids(idx);
          const int icol = import_lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_recv(pid)*f_col_sizes_scan_sum
                      + pos_within_pid;
          v(icol) = recv_buf(offset);
        };
        Kokkos::parallel_for(RangePolicy(0,num_imports),unpack);
        break;
      }
      case 2:
      {
        auto v = f.get_view<Real**>();
        const int dim1 = fl.dim(1);
        auto policy = ESU::get_default_team_policy(num_imports,dim1);
        auto unpack = KOKKOS_LAMBDA (const TeamMember& team) {
          const int idx  = team.league_rank();
          const int pid  = import_pids(idx);
          const int icol = import_lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_recv(pid)*f_col_sizes_scan_sum
                      + pos_within_pid*dim1;
          auto col_unpack = [&](const int& k) {
            v(icol,k) = recv_buf(offset+k);
          };
          auto tvr = Kokkos::TeamVectorRange(team,dim1);
          Kokkos::parallel_for(tvr,col_unpack);
        };
        Kokkos::parallel_for(policy,unpack);
        break;
      }
      case 3:
      {
        auto v = f.get_view<Real***>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        const int f_col_size = dim1*dim2;
        auto policy = ESU::get_default_team_policy(num_imports,dim1*dim2);
        auto unpack = KOKKOS_LAMBDA (const TeamMember& team) {
          const int idx  = team.league_rank();
          const int pid  = import_pids(idx);
          const int icol = import_lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_recv(pid)*f_col_sizes_scan_sum
                      + pos_within_pid*f_col_size;
          auto col_unpack = [&](const int& idx) {
            const int j = idx / dim2;
            const int k = idx % dim2;
            v(icol,j,k) = recv_buf(offset+idx);
          };
          auto tvr = Kokkos::TeamVectorRange(team,f_col_size);
          Kokkos::parallel_for(tvr,col_unpack);
        };
        Kokkos::parallel_for(policy,unpack);
        break;
      }
      case 4:
      {
        auto v = f.get_view<Real****>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        const int dim3 = fl.dim(3);
        const int f_col_size = dim1*dim2*dim3;
        auto policy = ESU::get_default_team_policy(num_imports,dim1*dim2*dim3);
        auto unpack = KOKKOS_LAMBDA (const TeamMember& team) {
          const int idx  = team.league_rank();
          const int pid  = import_pids(idx);
          const int icol = import_lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_recv(pid)*f_col_sizes_scan_sum
                      + pos_within_pid*f_col_size;
          auto col_unpack = [&](const int& idx) {
            const int j = (idx / dim3) / dim2;
            const int k = (idx / dim3) % dim2;
            const int l =  idx % dim3;
            v(icol,j,k,l) = recv_buf(offset+idx);
          };
          auto tvr = Kokkos::TeamVectorRange(team,f_col_size);
          Kokkos::parallel_for(tvr,col_unpack);
        };
        Kokkos::parallel_for(policy,unpack);
        break;
      }
      default:
        EKAT_ERROR_MSG ("Unexpected field rank in RefiningRemapperP2P::unpack.\n"
            "  - MPI rank  : " + std::to_string(m_comm.rank()) + "\n"
            "  - field name: " + f.name() + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }
}

void RefiningRemapperP2P::clean_up ()
{
  // Clear all MPI related structures
  m_send_buffer         = view_1d<Real>();
  m_recv_buffer         = view_1d<Real>();
  m_mpi_send_buffer     = mpi_view_1d<Real>();
  m_mpi_recv_buffer     = mpi_view_1d<Real>();
  m_send_req.clear();
  m_recv_req.clear();
  m_imp_exp = nullptr;

  HorizInterpRemapperBase::clean_up();
}

} // namespace scream
