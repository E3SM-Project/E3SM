#ifndef SCREAM_ATM_BUFFERS_MANAGER_HPP
#define SCREAM_ATM_BUFFERS_MANAGER_HPP

#include "share/eamxx_types.hpp"
#include "ekat/ekat_assert.hpp"

namespace scream {

// Struct which allows for the allocation of a single
// memory buffer for all ATM processes.
struct ATMBufferManager {

  template <typename S>
  using view_1d = typename KokkosTypes<DefaultDevice>::template view_1d<S>;

  ATMBufferManager()
  {
    m_size      = 0;
    m_allocated = false;
  }

  ~ATMBufferManager() = default;

  // Each ATM process should request the number of bytes
  // needed for local variables. Since no two process runs at
  // the same time, the total allocation will be the maximum
  // of each request.
  void request_bytes (const size_t num_bytes) {
    ekat::error::runtime_check(num_bytes%sizeof(Real)==0,
                               "Error! Must request number of bytes which is divisible by sizeof(Real).\n");

    const size_t num_reals = num_bytes/sizeof(Real);
    m_size = std::max(num_reals, m_size);
  }

  Real* get_memory () const { return m_buffer.data(); }

  size_t allocated_bytes () const { return m_size*sizeof(Real); }

  void allocate () {
    ekat::error::runtime_check(!m_allocated, "Error! Cannot call 'allocate' more than once.\n");

    m_buffer = view_1d<Real>("",m_size);
    m_allocated = true;
  }

  bool allocated () const { return m_allocated; }

protected:

  view_1d<Real> m_buffer;
  size_t        m_size;
  bool          m_allocated;
};

} // scream

#endif // SCREAM_ATM_BUFFERS_MANAGER_HPP
