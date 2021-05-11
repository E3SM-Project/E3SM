/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Comm.hpp"

#include <assert.h>
#include "Config.hpp"
#include "ErrorDefs.hpp"

namespace Homme
{

Comm::Comm()
{
  check_mpi_inited();
  reset_mpi_comm (MPI_COMM_SELF);
}

Comm::Comm(MPI_Comm mpi_comm)
{
  check_mpi_inited();
  reset_mpi_comm (mpi_comm);
}

void Comm::reset_mpi_comm (MPI_Comm new_mpi_comm)
{
  m_mpi_comm = new_mpi_comm;

  MPI_Comm_size(m_mpi_comm,&m_size);
  MPI_Comm_rank(m_mpi_comm,&m_rank);

#ifndef NDEBUG
  MPI_Comm_set_errhandler(m_mpi_comm,MPI_ERRORS_RETURN);
#endif
}

void Comm::check_mpi_inited () const
{
  int flag;
  MPI_Initialized (&flag);
  assert (flag!=0);
}

} // namespace Homme
