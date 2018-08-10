/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "MpiContext.hpp"

#include "BuffersManager.hpp"
#include "Comm.hpp"
#include "Connectivity.hpp"

namespace Homme {

MpiContext::MpiContext() {}

MpiContext::~MpiContext() {}

Comm& MpiContext::get_comm() {
  if ( ! comm_) {
    comm_.reset(new Comm());
  }
  return *comm_;
}

void MpiContext::create_comm(const int f_comm) {
  // You should NOT create a C MPI_Comm from F90 twice during the same execution
  assert (!comm_);

  MPI_Comm c_comm = MPI_Comm_f2c(f_comm);
  comm_.reset(new Comm(c_comm));
}

std::shared_ptr<BuffersManager> MpiContext::get_buffers_manager(short int exchange_type) {
  if ( ! buffers_managers_) {
    buffers_managers_.reset(new BMMap());
  }

  if (!(*buffers_managers_)[exchange_type]) {
    (*buffers_managers_)[exchange_type] = std::make_shared<BuffersManager>(get_connectivity());
  }
  return (*buffers_managers_)[exchange_type];
}

std::shared_ptr<Connectivity> MpiContext::get_connectivity() {
  if ( ! connectivity_) connectivity_.reset(new Connectivity());
  return connectivity_;
}

void MpiContext::clear() {
  comm_ = nullptr;
  connectivity_ = nullptr;
  buffers_managers_ = nullptr;
}

MpiContext& MpiContext::singleton() {
  static MpiContext c;
  return c;
}

void MpiContext::finalize_singleton() {
  singleton().clear();
}

} // namespace Homme
