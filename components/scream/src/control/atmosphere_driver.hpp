#ifndef SCREAM_ATMOSPHERE_DRIVER_HPP
#define SCREAM_ATMOSPHERE_DRIVER_HPP

#include <share/field_repository.hpp>
#include <share/scream_types.hpp>

#include <list>
#include <memory>

namespace scream {

// Forward declarations
class AtmosphereProcess;

namespace control {

/*
 *  The driver for the atmosphere component.
 *  
 *  This class (AD) is responsible to keep track of the different atmosphere
 *  subcomponents (ASC) (parametrizations and dynamics). The AD is responsible for ensuring
 *  that all the ASCs required by the current test case are correctly initialized,
 *  that they are called in the correct order, and that they do not leak memory at the end.
 *  
 *  The AD is also responsible for handling the FieldManager(s) of the Atmosphere.
 *  It must keep track of the different fields, and of the ASCs that require
 *  each field, ensuring all dependencies are met (in a DAG fashion).
 */

class AtmosphereDriver
{
public:
  template<typename VT>
  using field_repo_type = FieldRepository<VT>;

  using atm_process_type = AtmosphereProcess;

  // The initialization method should:
  //   1) create all the subcomponents needed, given the current simulation parameters
  //   2) initialize all the subcomponents 
  //   3) initialize the field manager(s)
  // The subcomponents are stored in the order requested by the user, so that when
  // going through the list, and calling their run method, they will be called in the
  // correct order. There should ALWAYS be a component that handles the dynamics. We should
  // make sure of that.
  void initialize ( /* inputs? */ );

  // The run method is responsible for advancing the atmosphere component by one atm time step
  // Inside here you should find calls to the run method of each subcomponent, including parametrizations
  // and dynamics (HOMME).
  void run ( /* inputs? */ );

  // Clean up the driver (includes cleaning up the parametrizations and the fm's);
  void finalize ( /* inputs */ );

protected:

  field_repo_type<ExecViewManaged<Real*>>         m_device_field_repo;
  field_repo_type<HostViewUnmanaged<Real*>>       m_host_field_repo;

  std::list<std::shared_ptr<atm_process_type>>    m_atm_processes;
};

int driver_stub();

}  // namespace control 
}  // namespace scream

#endif // SCREAM_ATMOSPHERE_DRIVER_HPP
