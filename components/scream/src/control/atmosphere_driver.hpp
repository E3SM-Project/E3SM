#ifndef SCREAM_ATMOSPHERE_DRIVER_HPP
#define SCREAM_ATMOSPHERE_DRIVER_HPP

#include "share/field/field_repository.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/scream_types.hpp"

#include "share/util/scream_time_stamp.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include <memory>

namespace scream {

// Forward declarations
class AtmosphereProcess;
class AtmosphereProcessGroup;

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

  // The initialization method should:
  //   1) create all the subcomponents needed, given the current simulation parameters
  //   2) initialize all the subcomponents
  //   3) initialize the field manager(s)
  // The subcomponents are stored in the order requested by the user, so that when
  // going through the list, and calling their run method, they will be called in the
  // correct order. There should ALWAYS be a component that handles the dynamics. We should
  // make sure of that.
  void initialize (const ekat::Comm& atm_comm,
                   const ekat::ParameterList& params,
                   const util::TimeStamp& t0 /*, inputs? */ );

  // The run method is responsible for advancing the atmosphere component by one atm time step
  // Inside here you should find calls to the run method of each subcomponent, including parameterizations
  // and dynamics (HOMME).
  void run (const Real dt);

  // Clean up the driver (includes cleaning up the parameterizations and the fm's);
  void finalize ( /* inputs */ );

  const FieldRepository<Real>& get_field_repo () const { return *m_field_repo; }
#ifdef SCREAM_DEBUG
  const FieldRepository<Real>& get_bkp_field_repo () const { return m_bkp_field_repo; }
#endif

  // Get atmosphere time stamp
  const util::TimeStamp& get_atm_time_stamp () const { return m_current_ts; }
protected:

  void register_groups ();
  void init_atm_inputs ();
  void inspect_atm_dag ();
#ifdef SCREAM_DEBUG
  void create_bkp_field_repo ();
#endif

  std::shared_ptr<FieldRepository<Real> >  m_field_repo;
#ifdef SCREAM_DEBUG
  FieldRepository<Real>                    m_bkp_field_repo;
#endif
  ekat::WeakPtrSet<FieldInitializer>                  m_field_initializers;

  std::shared_ptr<AtmosphereProcessGroup>             m_atm_process_group;

  std::shared_ptr<GridsManager>                       m_grids_manager;

  ekat::ParameterList                                 m_atm_params;

  // This are the time stamps of the start and end of the time step.
  util::TimeStamp                       m_old_ts;
  util::TimeStamp                       m_current_ts;

  // This is the comm containing all (and only) the processes assigned to the atmosphere
  ekat::Comm   m_atm_comm;
};

}  // namespace control
}  // namespace scream

#endif // SCREAM_ATMOSPHERE_DRIVER_HPP
