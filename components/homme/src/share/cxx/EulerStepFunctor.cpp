/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "EulerStepFunctorImpl.hpp"

namespace Homme {

EulerStepFunctor
::EulerStepFunctor ()
  : is_setup(true) {
  p_ = std::make_shared<EulerStepFunctorImpl>();
}

// This constructor is useful for using buffer functionality without
// having all other Functor information available.
// If this constructor is used, the setup() function must be called
// before using any other EulerStepFunctor functions.
EulerStepFunctor
::EulerStepFunctor (const int num_elems)
  : is_setup(false)
{
  p_ = std::make_shared<EulerStepFunctorImpl>(num_elems);
}

void EulerStepFunctor::setup()
{
  // Sanity check
  assert (!is_setup);

  p_->setup();
  is_setup = true;
}

void EulerStepFunctor::reset (const SimulationParams& params) {
  p_->reset(params);
}

int EulerStepFunctor::requested_buffer_size () const {
  return p_->requested_buffer_size();
}
void EulerStepFunctor::init_buffers (const FunctorsBuffersManager& fbm) {
  p_->init_buffers(fbm);
}

void EulerStepFunctor::init_boundary_exchanges () {
  // The Functor needs to be fully setup to use this function
  assert (is_setup);

  p_->init_boundary_exchanges();
}

void EulerStepFunctor::precompute_divdp () {
  // The Functor needs to be fully setup to use this function
  assert (is_setup);

  p_->precompute_divdp();
}

void EulerStepFunctor
::euler_step (const int np1_qdp, const int n0_qdp, const Real dt,
              const Real rhs_multiplier, const DSSOption DSSopt) {
  // The Functor needs to be fully setup to use this function
  assert (is_setup);

  p_->euler_step(np1_qdp, n0_qdp, dt, rhs_multiplier, DSSopt);
}

void EulerStepFunctor
::qdp_time_avg (const int n0_qdp, const int np1_qdp) {
  // The Functor needs to be fully setup to use this function
  assert (is_setup);

  p_->qdp_time_avg(n0_qdp, np1_qdp);
}

} // namespace Homme
