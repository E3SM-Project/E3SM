#ifndef SCREAM_ATMOSPHERE_DYNAMICS_HPP
#define SCREAM_ATMOSPHERE_DYNAMICS_HPP

#include <share/atmosphere_process.hpp>

#include <string>
#include <list>

namespace scream
{

/*
 *  The class responsible to handle the atmosphere dynamics
 *
 *  The AD should only exactly ONE instance of this class stored
 *  in its list of subcomponents (the AD should make sure of this).
 *
 *  For now, Scream is only going to accommodate HOMME as a dynamics
 *  dycore.
 */

class AtmosphereDynamics : public AtmosphereProcess
{
public:

  // The type of the subcomponent (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Dynamics; }

  // The name of the subcomponent 
  std::string name () const { return "dynamics"; }

  // The communicator used by the dynamics
  const Comm& get_comm () const { return m_dynamics_comm; }

  // These are the three main interfaces:
  void initialize (/* what inputs? */);
  void run        (/* what inputs? */);
  void finalize   (/* what inputs? */);

  // These two methods allow the driver to figure out what subcomponent need
  // a given field and what subcomponent updates a given field.
  const std::list<std::string>&  get_required_fields () const { return m_input_fields;  }
  const std::list<std::string>&  get_computed_fields () const { return m_output_fields; }

protected:

  std::list<std::string> m_input_fields;
  std::list<std::string> m_output_fields;

  Comm      m_dynamics_comm;

};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_DYNAMICS_HPP
