#ifndef HOMMEXX_DIAGNOSTICS_HPP
#define HOMMEXX_DIAGNOSTICS_HPP

#include <Types.hpp>

namespace Homme
{

class Diagnostics
{
public:
  void init (const int num_elems, F90Ptr& elem_state_q_ptr, F90Ptr& elem_accum_qvar_ptr, F90Ptr& elem_accum_qmass_ptr,
             F90Ptr& elem_accum_q1mass_ptr, F90Ptr& elem_accum_iener_ptr, F90Ptr& elem_accum_iener_wet_ptr,
             F90Ptr& elem_accum_kener_ptr, F90Ptr& elem_accum_pener_ptr);

  void prim_diag_scalars (const bool before_advance, const int ivar);
  void prim_energy_halftimes (const bool before_advance, const int ivar);

  HostViewUnmanaged<Real*[QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]> h_Q;

private:

  HostViewUnmanaged<Real*[4][QSIZE_D][NP][NP]>        h_Qvar;
  HostViewUnmanaged<Real*[4][QSIZE_D][NP][NP]>        h_Qmass;
  HostViewUnmanaged<Real*   [QSIZE_D][NP][NP]>        h_Q1mass;
  HostViewUnmanaged<Real*[4][NP][NP]>                 h_IEner;
  HostViewUnmanaged<Real*[4][NP][NP]>                 h_IEner_wet;
  HostViewUnmanaged<Real*[4][NP][NP]>                 h_KEner;
  HostViewUnmanaged<Real*[4][NP][NP]>                 h_PEner;

  int m_num_elems;
};

} // namespace Homme

#endif // HOMMEXX_DIAGNOSTICS_HPP

