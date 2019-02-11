#ifndef HOMMEXX_DIAGNOSTICS_BASE_HPP
#define HOMMEXX_DIAGNOSTICS_BASE_HPP

#include <Types.hpp>

namespace Homme
{

class DiagnosticsBase
{
public:
  void init (const int num_elems,
             F90Ptr& elem_state_q_ptr, F90Ptr& elem_accum_qvar_ptr,
             F90Ptr& elem_accum_qmass_ptr, F90Ptr& elem_accum_q1mass_ptr);

  void prim_diag_scalars (const bool before_advance, const int ivar);

  virtual void prim_energy_halftimes (const bool before_advance, const int ivar) = 0;

  HostViewUnmanaged<Real*[QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]> h_Q;

protected:

  HostViewUnmanaged<Real*[4][QSIZE_D][NP][NP]>        h_Qvar;
  HostViewUnmanaged<Real*[4][QSIZE_D][NP][NP]>        h_Qmass;
  HostViewUnmanaged<Real*   [QSIZE_D][NP][NP]>        h_Q1mass;

  int m_num_elems;
};

} // namespace Homme

#endif // HOMMEXX_DIAGNOSTICS_BASE_HPP
