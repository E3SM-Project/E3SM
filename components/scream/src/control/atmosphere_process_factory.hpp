#ifndef SCREAM_ATMOSPHERE_PROCESS_FACTORY_HPP
#define SCREAM_ATMOSPHERE_PROCESS_FACTORY_HPP

#include "share/atmosphere_process.hpp"
#include "share/parameter_list.hpp"

namespace scream {

/*
 *  Tiny struct to create atmosphere processes.
 *  
 *  The struct only contains one static method to create an atmosphere
 *  process (AP) given its name and a parameter list, nothing else.
 *  We could have simply put the body of the create function in place
 *  to any call of 'create', given its limited size. The only reason
 *  why we introduced this struct is to isolate concepts in the code.
 *  In particular, the only place (so far) that creates APs is the
 *  AtmosphereProcessGroup class (APG). By putting the code for the
 *  APs creation in a separate translation unit, we relieve the APG
 *  class from the need to know what APs are available. It is not the
 *  APG concern to know details of downstream classes.
 *  
 *  If you are a fan of the SOLID acronym in programming, this goes
 *  in the direction suggested by the 'D' in that acronym: dependency
 *  inversion principle.
 */

struct AtmosphereProcessFactory {
  static AtmosphereProcess* create (const std::string& name, const ParameterList& params);
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_FACTORY
