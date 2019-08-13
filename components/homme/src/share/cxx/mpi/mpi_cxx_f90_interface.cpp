/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "MpiContext.hpp"
#include "Comm.hpp"
#include "Connectivity.hpp"
#include "BoundaryExchange.hpp"

#include <map>

namespace Homme
{

extern "C"
{

void reset_cxx_comm (const MPI_Fint& f_comm)
{
  // f_comm must be a valid Fortran handle to a communicator
  MPI_Comm c_comm = MPI_Comm_f2c(f_comm);
  MpiContext::singleton().get_comm().reset_mpi_comm(c_comm);
}

void init_connectivity (const int& num_local_elems)
{
  Connectivity& connectivity = *MpiContext::singleton().get_connectivity();
  connectivity.set_num_elements(num_local_elems);
  connectivity.set_comm(MpiContext::singleton().get_comm());
}

void add_connection (const int& first_elem_lid,  const int& first_elem_gid,  const int& first_elem_pos,  const int& first_elem_pid,
                     const int& second_elem_lid, const int& second_elem_gid, const int& second_elem_pos, const int& second_elem_pid)
{
  // Check that F90 is in base 1
  if (first_elem_lid<=0  || first_elem_gid<=0  || first_elem_pos<=0  || first_elem_pid<=0 ||
      second_elem_lid<=0 || second_elem_gid<=0 || second_elem_pos<=0 || second_elem_pid<=0)
  {
    std::cout << "ERROR! We were assuming F90 indices started at 1, but it appears there is an exception.\n";
    std::abort();
  }

  // S, N, W, E is unpack order, so we want elem_pos to reflect
  // that. Convert here.
  const auto convert = [=] (int fpos) -> int {
    return fpos <= 4 ? (((fpos-1) + 2) % 4) : (fpos-1);
  };
  const int fep = convert(first_elem_pos), sep = convert(second_elem_pos);

  Connectivity& connectivity = *MpiContext::singleton().get_connectivity();
  connectivity.add_connection(first_elem_lid-1, first_elem_gid-1, fep, first_elem_pid-1,
                              second_elem_lid-1,second_elem_gid-1,sep,second_elem_pid-1);
}

void finalize_connectivity ()
{
  Connectivity& connectivity = *MpiContext::singleton().get_connectivity();

  connectivity.finalize();
}

} // extern "C"

} // namespace Homme
