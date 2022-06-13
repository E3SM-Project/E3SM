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
  ~FunctorsBuffersManager() = default;

  void request_size (const int num_doubles);

  Real* get_memory () const { return m_buffer.data(); }

  int allocated_size () const { return m_size; }

  void allocate ();

  // Reset m_size and allocate using input data.
  // This is useful if the user wants to allocate
  // buffer data outside of Homme.
  void allocate (Real* data, const int new_size);

  bool allocated () const { return m_allocated; }

protected:

  ExecViewManaged<Real*> m_buffer;
  int   m_size;
  bool  m_allocated;

private:

  void generate_random_data();
};

} // Homme

#endif // HOMMEXX_FUNCTORS_BUFFERS_MANAGER_HPP
