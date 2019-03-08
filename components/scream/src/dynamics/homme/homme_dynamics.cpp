#include "homme_dynamics.hpp"

namespace scream
{

HommeDynamics::HommeDynamics (const ParameterList& /* params */) {
}

void HommeDynamics::initialize (const Comm& comm)
{
  m_dynamics_comm = comm;
}

void HommeDynamics::run (/* what inputs? */)
{

}

void HommeDynamics::finalize (/* what inputs? */)
{

}

void HommeDynamics::register_fields (FieldRepository<Real, device_type>& /* field_repo */) const {
  // register in/out fields in the repo
}

} // namespace scream
