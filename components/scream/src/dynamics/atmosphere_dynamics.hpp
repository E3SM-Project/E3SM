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
  AtmosphereDynamics () = default;

  // The type of the subcomponent (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Dynamics; }

  // The name of the subcomponent 
  std::string name () const { return "Scream::Dynamics"; }

  // These are the three main interfaces:
  void initialize (/* what inputs? */);
  void run        (/* what inputs? */);
  void finalize   (/* what inputs? */);

  // These two methods allow the driver to figure out what subcomponent need
  // a given field and what subcomponent updates a given field.
  const std::list<std::string>&  get_required_fields () const;
  const std::list<std::string>&  get_computed_fields () const;

protected:

};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_DYNAMICS_HPP
