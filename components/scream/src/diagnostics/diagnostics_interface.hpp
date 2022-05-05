#ifndef EAMXX_DIAGNOSTICS_INTERFACE_HPP
#define EAMXX_DIAGNOSTICS_INTERFACE_HPP

/* 
 * It is important to add each new diagnostic to
 * this list of includes.
 */
#include "diagnostics/potential_temperature.hpp"
#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {

inline void register_diagnostics () {
  // Note, we use the atmosphere process factory to register diagnostics
  // because each diagnostic is a sub-class of an atmopshere process
  auto& proc_factory = AtmosphereProcessFactory::instance();
 
  proc_factory.register_product("Potential Temperature",&create_atmosphere_process<PotentialTemperatureDiagnostic>); 
}

class DiagnosticInterface
{
public:
  // Constructor
  DiagnosticInterface ();

  virtual ~DiagnosticInterface () = default;

  void add_diagnostic(const std::string name);
  void init();
  Field get_diagnostic(const std::string name);
  void finalize();
  
protected:
//  std::map<std::string,AtmosphereDiagnostic> m_diagnostics;

}; // DiagnosticInterface

} //scream

#endif //EAMXX_DIAGNOSTICS_INTERFACE_HPP

