/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "FunctorsBuffersManager.hpp"
#include "ErrorDefs.hpp"
#include "utilities/SubviewUtils.hpp"

namespace Homme {

FunctorsBuffersManager::FunctorsBuffersManager ()
 : m_num_2d_scalar_bufs (0)
 , m_num_2d_vector_bufs (0)
 , m_num_3d_scalar_mid_bufs (0)
 , m_num_3d_vector_mid_bufs (0)
 , m_num_3d_scalar_int_bufs (0)
 , m_num_3d_vector_int_bufs (0)
 , m_num_concurrent_teams   (-1)
 , m_allocated (false)
{
  // Nothing to be done here
}

void FunctorsBuffersManager::request_concurrency (const int num_concurrent_teams) {
  Errors::runtime_check (!m_allocated, "Error! Cannot request buffers after allocation happened.\n");

  m_num_concurrent_teams = std::max(m_num_concurrent_teams,num_concurrent_teams);
}

void FunctorsBuffersManager::request_2d_buffers (const int num_scalar_bufs, const int num_vector_bufs) {
  Errors::runtime_check (!m_allocated, "Error! Cannot request buffers after allocation happened.\n");
  m_num_2d_scalar_bufs = std::max(m_num_2d_scalar_bufs,num_scalar_bufs);
  m_num_2d_vector_bufs = std::max(m_num_2d_scalar_bufs,num_vector_bufs);
}

void FunctorsBuffersManager::request_3d_midpoint_buffers (const int num_scalar_bufs, const int num_vector_bufs) {
  Errors::runtime_check (!m_allocated, "Error! Cannot request buffers after allocation happened.\n");

  m_num_3d_scalar_mid_bufs = std::max(m_num_3d_scalar_mid_bufs,num_scalar_bufs);
  m_num_3d_vector_mid_bufs = std::max(m_num_3d_vector_mid_bufs,num_vector_bufs);
}

void FunctorsBuffersManager::request_3d_interface_buffers (const int num_scalar_bufs, const int num_vector_bufs) {
  m_num_3d_scalar_int_bufs = std::max(m_num_3d_scalar_int_bufs,num_scalar_bufs);
  m_num_3d_vector_int_bufs = std::max(m_num_3d_vector_int_bufs,num_vector_bufs);
}

void FunctorsBuffersManager:: allocate () {
  Errors::runtime_check(!m_allocated, "Error! Cannot call 'allocate' more than once.\n");

  m_scalar_2d_buffers = ExecViewManaged<Real**   [NP][NP]>("Scalar 2d buffers", m_num_2d_scalar_bufs, m_num_concurrent_teams);
  m_vector_2d_buffers = ExecViewManaged<Real**[2][NP][NP]>("Vector 2d buffers", m_num_2d_vector_bufs, m_num_concurrent_teams);

  m_scalar_3d_midpoint_buffers = ExecViewManaged<Scalar**   [NP][NP][NUM_LEV]>("Scalar 2d buffers", m_num_3d_scalar_mid_bufs, m_num_concurrent_teams);
  m_vector_3d_midpoint_buffers = ExecViewManaged<Scalar**[2][NP][NP][NUM_LEV]>("Vector 2d buffers", m_num_3d_vector_mid_bufs, m_num_concurrent_teams);

  m_scalar_3d_interface_buffers = ExecViewManaged<Scalar**   [NP][NP][NUM_LEV_P]>("Scalar 2d buffers", m_num_3d_scalar_int_bufs, m_num_concurrent_teams);
  m_vector_3d_interface_buffers = ExecViewManaged<Scalar**[2][NP][NP][NUM_LEV_P]>("Scalar 2d buffers", m_num_3d_vector_int_bufs, m_num_concurrent_teams);

  m_allocated = true;
}

ExecViewUnmanaged<Real *[NP][NP]>
FunctorsBuffersManager::get_2d_scalar_buffer (const int ibuf) const
{
  Errors::runtime_check(m_allocated, "Error! Buffer index out of bounds.\n");
  Errors::runtime_check(ibuf>=0 && ibuf<m_num_3d_scalar_mid_bufs, "Error! Buffer index out of bounds.\n");

  return Homme::subview(m_scalar_2d_buffers,ibuf);
}

ExecViewUnmanaged<Real *[2][NP][NP]>
FunctorsBuffersManager::get_2d_vector_buffer (const int ibuf) const
{
  Errors::runtime_check(m_allocated, "Error! Buffer index out of bounds.\n");
  Errors::runtime_check(ibuf>=0 && ibuf<m_num_3d_vector_mid_bufs, "Error! Buffer index out of bounds.\n");

  return Homme::subview(m_vector_2d_buffers,ibuf);
}

ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV]>
FunctorsBuffersManager::get_3d_scalar_midpoint_buffer (const int ibuf) const
{
  Errors::runtime_check(m_allocated, "Error! Buffer index out of bounds.\n");
  Errors::runtime_check(ibuf>=0 && ibuf<m_num_3d_scalar_mid_bufs, "Error! Buffer index out of bounds.\n");

  return Homme::subview(m_scalar_3d_midpoint_buffers,ibuf);
}

ExecViewUnmanaged<Scalar *[2][NP][NP][NUM_LEV]>
FunctorsBuffersManager::get_3d_vector_midpoint_buffer (const int ibuf) const
{
  Errors::runtime_check(m_allocated, "Error! Buffer index out of bounds.\n");
  Errors::runtime_check(ibuf>=0 && ibuf<m_num_3d_vector_mid_bufs, "Error! Buffer index out of bounds.\n");

  return Homme::subview(m_vector_3d_midpoint_buffers,ibuf);
}

ExecViewUnmanaged<Scalar *   [NP][NP][NUM_LEV_P]>
FunctorsBuffersManager::get_3d_scalar_interface_buffer(const int ibuf) const
{
  Errors::runtime_check(m_allocated, "Error! Buffer index out of bounds.\n");
  Errors::runtime_check(ibuf>=0 && ibuf<m_num_3d_scalar_int_bufs, "Error! Buffer index out of bounds.\n");

  return Homme::subview(m_scalar_3d_interface_buffers,ibuf);
}

ExecViewUnmanaged<Scalar *[2][NP][NP][NUM_LEV_P]>
FunctorsBuffersManager::get_3d_vector_interface_buffer(const int ibuf) const
{
  Errors::runtime_check(m_allocated, "Error! Buffer index out of bounds.\n");
  Errors::runtime_check(ibuf>=0 && ibuf<m_num_3d_vector_int_bufs, "Error! Buffer index out of bounds.\n");

  return Homme::subview(m_vector_3d_interface_buffers,ibuf);
}

} // namespace Homme
