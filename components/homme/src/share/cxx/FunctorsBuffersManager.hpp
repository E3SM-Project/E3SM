/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_FUNCTORS_BUFFERS_MANAGER_HPP
#define HOMMEXX_FUNCTORS_BUFFERS_MANAGER_HPP

#include "Types.hpp"

namespace Homme {

// Elements-dependent buffers.
struct FunctorsBuffersManager {

  FunctorsBuffersManager();

  // For each buffer id, request a number of concurrent buffers. This means that there
  // will ba eup to num_concurrent_buffers teams requesting the same buffer at the same time.
  void request_concurrency (const int num_concurrent_teams);

  void request_2d_buffers           (const int num_scalar_bufs, const int num_vector_bufs);
  void request_3d_midpoint_buffers  (const int num_scalar_bufs, const int num_vector_bufs);
  void request_3d_interface_buffers (const int num_scalar_bufs, const int num_vector_bufs);

  void allocate ();
  bool are_buffers_allocated () const { return m_allocated; }

  ExecViewUnmanaged<Real *   [NP][NP]> get_2d_scalar_buffer (const int ibuf) const;
  ExecViewUnmanaged<Real *[2][NP][NP]> get_2d_vector_buffer (const int ibuf) const;

  ExecViewUnmanaged<Scalar *   [NP][NP][NUM_LEV]>   get_3d_scalar_midpoint_buffer (const int ibuf) const;
  ExecViewUnmanaged<Scalar *[2][NP][NP][NUM_LEV]>   get_3d_vector_midpoint_buffer (const int ibuf) const;

  ExecViewUnmanaged<Scalar *   [NP][NP][NUM_LEV_P]> get_3d_scalar_interface_buffer (const int ibuf) const;
  ExecViewUnmanaged<Scalar *[2][NP][NP][NUM_LEV_P]> get_3d_vector_interface_buffer (const int ibuf) const;

protected:

  // The runtime dimensions are the number of buffers and the number of concurrent teams
  ExecViewManaged<Real**    [NP][NP]> m_scalar_2d_buffers;
  ExecViewManaged<Real** [2][NP][NP]> m_vector_2d_buffers;

  ExecViewManaged<Scalar**    [NP][NP][NUM_LEV]>   m_scalar_3d_midpoint_buffers;
  ExecViewManaged<Scalar** [2][NP][NP][NUM_LEV]>   m_vector_3d_midpoint_buffers;
  ExecViewManaged<Scalar**    [NP][NP][NUM_LEV_P]> m_scalar_3d_interface_buffers;
  ExecViewManaged<Scalar** [2][NP][NP][NUM_LEV_P]> m_vector_3d_interface_buffers;

  int   m_num_2d_scalar_bufs;
  int   m_num_2d_vector_bufs;
  int   m_num_3d_scalar_mid_bufs;
  int   m_num_3d_vector_mid_bufs;
  int   m_num_3d_scalar_int_bufs;
  int   m_num_3d_vector_int_bufs;

  int   m_num_concurrent_teams;
  bool  m_allocated;
};

} // Homme

#endif // HOMMEXX_FUNCTORS_BUFFERS_MANAGER_HPP
