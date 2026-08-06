/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Context.hpp"
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
  if (!Context::singleton().has<Comm>()) {
    Context::singleton().create<Comm>();
  }
  Context::singleton().get<Comm>().reset_mpi_comm(c_comm);
}

void init_connectivity (const int& num_local_elems, const int& max_corner_elems)
{
  Connectivity& connectivity = Context::singleton().create<Connectivity>();
  connectivity.set_num_elements(num_local_elems);
  connectivity.set_max_corner_elements(max_corner_elems);
  connectivity.set_comm(Context::singleton().get<Comm>());
}

// Extract dir (in 0:7) and dir_idx (in 0:max_corner_elements-1) from
//     fpos = (dir+1)*max_corner_elems + dir_idx.
static void convert (const int max_corner_elements, const int fpos,
                     std::uint8_t& dir, std::uint8_t& dir_idx)
{
  assert(fpos > 0);
  assert(fpos < 9*max_corner_elements);
  dir = fpos / max_corner_elements - 1;
  if (dir < 4) {
    // edge_mod::edgeVunpack_nlyr establishes the S, N, W, E unpack order for
    // element edges. Return the position in this order for this edge.
    dir = (dir + 2) % 4;
  }
  // In the unpack phase, corners are unpacked after edges. For the cubed-sphere
  // grid, BFB accumulation of corners is trivial since there is at most one
  // corner-attached element per corner. For an RRM grid, multiple elements can
  // be associated with a corner, and the order must follow the convention
  // established in gridgraph_mod.F90 and schedule_mod.F90.
  dir_idx = fpos % max_corner_elements;
  assert(dir_idx == 0 || dir >= 4);
}

void add_connection (const int& e1_lid, const int& e1_gid, const int& e1_pos, const int& e1_pid,
                     const int& e2_lid, const int& e2_gid, const int& e2_pos, const int& e2_pid)
{
  // Check that F90 is in base 1
  Errors::runtime_check(e1_lid >= 1 && e1_gid >= 1 && e1_pos >= 1 && e1_pid >= 1 &&
                        e2_lid >= 1 && e2_gid >= 1 && e2_pos >= 1 && e2_pid >= 1,
                        "add_connection: F90 indices should start at 1");

  Connectivity& connectivity = Context::singleton().get<Connectivity>();
  const auto max_corner_elements = connectivity.get_max_corner_elements();
  assert(max_corner_elements >= 1);

  std::uint8_t e1_d, e1_didx, e2_d, e2_didx;
  convert(max_corner_elements, e1_pos, e1_d, e1_didx);
  convert(max_corner_elements, e2_pos, e2_d, e2_didx);

  connectivity.add_connection(e1_lid-1, e1_gid-1, e1_d, e1_didx, e1_pid-1,
                              e2_lid-1, e2_gid-1, e2_d, e2_didx, e2_pid-1);
}

void finalize_connectivity ()
{
  Connectivity& connectivity = Context::singleton().get<Connectivity>();

  connectivity.finalize();
}

} // extern "C"

} // namespace Homme
