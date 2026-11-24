#include "share/remap/mpi_remapper.hpp"

#include "share/util/eamxx_utils.hpp"


#include <limits>

namespace scream
{

MpiRemapper::MpiRequests::
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

void MpiRemapper::
setup_mpi (const std::shared_ptr<const GridImportExport>& imp_exp,
           const CommDirection dir)
{
  using namespace ShortFieldTagsNames;

  m_imp_exp = imp_exp;
  m_send_is_export = dir==CommDirection::UniqueToOverlap;

  m_field_col_offset.assign(m_num_fields+1,0);
  for (int i=0; i<m_num_fields; ++i) {
    const auto& fl = m_tgt_fields[i].get_header().get_identifier().get_layout();
    m_field_col_offset[i+1] = m_field_col_offset[i] + (fl.has_tag(COL) ? fl.clone().strip_dim(COL).size() : 0);
  }
  const size_t col_size = total_col_size();
  if (col_size==0) {
    // None of the registered fields has columns to exchange
    return;
  }

  const auto& comm = m_src_grid->get_comm();
  const auto send_offsets = send_offsets_h();
  const auto recv_offsets = recv_offsets_h();

  // NOTE: sizes and offsets are computed with size_t, since the product of two ints can overflow
  m_send_buffer = view_1d<Real>("MpiRemapper::send_buf",static_cast<size_t>(num_send_cols())*col_size);
  m_recv_buffer = view_1d<Real>("MpiRemapper::recv_buf",static_cast<size_t>(num_recv_cols())*col_size);
  m_mpi_send_buffer = Kokkos::create_mirror_view(typename decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(typename decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  m_requests = std::make_shared<MpiRequests>();
  const auto mpi_comm = comm.mpi_comm();
  const auto mpi_real = ekat::get_mpi_type<Real>();
  for (int pid=0; pid<comm.size(); ++pid) {
    for (const bool send : {true,false}) {
      const auto& offsets = send ? send_offsets : recv_offsets;
      const int ncols = offsets(pid+1)-offsets(pid);
      if (ncols==0) {
        continue;
      }
      const size_t count = static_cast<size_t>(ncols)*col_size;
      EKAT_REQUIRE_MSG (count<=static_cast<size_t>(std::numeric_limits<int>::max()),
          "Error! The message to " + std::string(send ? "send to" : "receive from") + " a single rank is too large for MPI's int count.\n"
          " - remapper: " + name() + "\n"
          " - number of entries: " + std::to_string(count) + "\n");
      const size_t pos = static_cast<size_t>(offsets(pid))*col_size;
      if (send) {
        auto& req = m_requests->send.emplace_back();
        check_mpi_call(MPI_Send_init(m_mpi_send_buffer.data()+pos,static_cast<int>(count),mpi_real,pid,0,mpi_comm,&req),
                       "[" + name() + "] creating persistent send request.\n");
      } else {
        auto& req = m_requests->recv.emplace_back();
        check_mpi_call(MPI_Recv_init(m_mpi_recv_buffer.data()+pos,static_cast<int>(count),mpi_real,pid,0,mpi_comm,&req),
                       "[" + name() + "] creating persistent recv request.\n");
      }
    }
  }
}

void MpiRemapper::
clean_mpi ()
{
  m_requests.reset();
  m_imp_exp.reset();
  m_field_col_offset.clear();
  m_send_buffer     = view_1d<Real>();
  m_recv_buffer     = view_1d<Real>();
  m_mpi_send_buffer = mpi_view_1d<Real>();
  m_mpi_recv_buffer = mpi_view_1d<Real>();
}

void MpiRemapper::
start (std::vector<MPI_Request>& reqs, const std::string& what) const
{
  if (not reqs.empty()) {
    check_mpi_call(MPI_Startall(reqs.size(),reqs.data()),
                   "[" + name() + "] starting persistent " + what + " requests.\n");
  }
}

void MpiRemapper::
wait (std::vector<MPI_Request>& reqs, const std::string& what) const
{
  if (not reqs.empty()) {
    check_mpi_call(MPI_Waitall(reqs.size(),reqs.data(),MPI_STATUSES_IGNORE),
                   "[" + name() + "] waiting on persistent " + what + " requests.\n");
  }
}

void MpiRemapper::start_sends () { if (m_requests) start(m_requests->send,"send"); }
void MpiRemapper::start_recvs () { if (m_requests) start(m_requests->recv,"recv"); }
void MpiRemapper::wait_sends  () { if (m_requests) wait(m_requests->send,"send"); }
void MpiRemapper::wait_recvs  () { if (m_requests) wait(m_requests->recv,"recv"); }

void MpiRemapper::
copy_send_to_mpi ()
{
  if (not MpiOnDev) {
    Kokkos::deep_copy(m_mpi_send_buffer,m_send_buffer);
  }
}

void MpiRemapper::
copy_recv_from_mpi ()
{
  if (not MpiOnDev) {
    Kokkos::deep_copy(m_recv_buffer,m_mpi_recv_buffer);
  }
}

} // namespace scream
