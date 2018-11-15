#ifndef SCREAM_ATMOSPHERE_DYNAMICS_HPP
#define SCREAM_ATMOSPHERE_DYNAMICS_HPP

#include <share/atmosphere_process.hpp>

#include <string>

namespace scream
{

/*
 *  The class responsible to handle the atmosphere dynamics
 *
 *  The AD should store exactly ONE instance of this class stored
 *  in its list of subcomponents (the AD should make sure of this).
 *
 *  Note: for now, Scream is only going to accommodate HOMME as a dynamics
 *  dycore.
 */

class AtmosphereDynamics : public AtmosphereProcess
{
public:

  // Constructor(s)
  AtmosphereDynamics (const ParameterList& params);

  // The type of the subcomponent (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Dynamics; }

  // The name of the subcomponent 
  std::string name () const { return "dynamics"; }

  // The communicator used by the dynamics
  const Comm& get_comm () const { return m_dynamics_comm; }

  // These are the three main interfaces:
  void initialize (const Comm& comm);
  void run        (/* what inputs? */);
  void finalize   (/* what inputs? */);

  // Get the set of required/computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_computed_fields; }

protected:

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real*, ExecMemSpace, MemoryManaged>& /*f*/) { /* impl */ }
  void set_computed_field_impl (const Field<      Real*, ExecMemSpace, MemoryManaged>& /*f*/) { /* impl */ }

  std::set<FieldIdentifier> m_required_fields;
  std::set<FieldIdentifier> m_computed_fields;

  Comm      m_dynamics_comm;

};

AtmosphereProcess* create_atmosphere_dynamics (const ParameterList& params);

} // namespace scream

#endif // SCREAM_ATMOSPHERE_DYNAMICS_HPP
