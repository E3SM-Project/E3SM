#ifndef SCREAM_ATMOSPHERE_DRIVER_HPP
#define SCREAM_ATMOSPHERE_DRIVER_HPP

#include "control/surface_coupling.hpp"

#include "share/field/field_repository.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/scream_types.hpp"
#include "share/io/output_manager.hpp"
#include "share/io/scorpio_input.hpp"

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

  AtmosphereDriver () = default;
  AtmosphereDriver (const ekat::Comm& atm_comm,
                    const ekat::ParameterList& params);

  // The default dtor is fine.
  ~AtmosphereDriver () = default;

  // ---- Begin initialization methods ---- //

  // Set comm for the whole atmosphere
  void set_comm (const ekat::Comm& atm_comm);

  // Set AD params
  void set_params (const ekat::ParameterList& params);

  // Create atm processes, without initializing them
  void create_atm_processes ();

  // Create needed grids, based on processes needs.
  void create_grids ();

  // Create fields as requested by all processes
  void create_fields ();

  // Sets a pre-built SurfaceCoupling object in the driver (for CIME runs only)
  void set_surface_coupling (const std::shared_ptr<SurfaceCoupling>& sc);

  // Load initial conditions for atm inputs
  void initialize_fields (const util::TimeStamp& t0);

  // Initialie I/O structures for output
  void initialize_output_manager ();

  // Call 'initialize' on all atm procs
  void initialize_atm_procs ();

  // Complete any leftover initialization task (e.g., some debug stuff)
  void finish_setup ();

  // ---- End of initialization methods ---- //

  // A wrapper of all of the above (except setting SurfaceCoupling),
  // which is handy for scream standalone tests.
  void initialize (const ekat::Comm& atm_comm,
                   const ekat::ParameterList& params,
                   const util::TimeStamp& t0);

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

  const std::shared_ptr<SurfaceCoupling>& get_surface_coupling () const { return m_surface_coupling; }

  // Get atmosphere time stamp
  const util::TimeStamp& get_atm_time_stamp () const { return m_current_ts; }

  const std::shared_ptr<GridsManager>& get_grids_manager () const { return m_grids_manager; }

protected:

  void register_groups ();
#ifdef SCREAM_DEBUG
  void create_bkp_field_repo ();
#endif

  std::shared_ptr<FieldRepository<Real> >  m_field_repo;
#ifdef SCREAM_DEBUG
  FieldRepository<Real>                    m_bkp_field_repo;
#endif

  std::shared_ptr<AtmosphereProcessGroup>             m_atm_process_group;

  std::shared_ptr<GridsManager>                       m_grids_manager;

  ekat::ParameterList                                 m_atm_params;

  OutputManager                                       m_output_manager;

  // Surface coupling stuff
  std::shared_ptr<SurfaceCoupling>            m_surface_coupling;

  // This are the time stamps of the start and end of the time step.
  util::TimeStamp                       m_current_ts;

  // This is the comm containing all (and only) the processes assigned to the atmosphere
  ekat::Comm   m_atm_comm;

  // Some status flags, used to make sure we call the init functions in the right order
  static constexpr int s_comm_set       =   1;
  static constexpr int s_params_set     =   2;
  static constexpr int s_procs_created  =   4;
  static constexpr int s_grids_created  =   8;
  static constexpr int s_fields_created =  16;
  static constexpr int s_sc_set         =  32;
  static constexpr int s_output_inited  =  64;
  static constexpr int s_fields_inited  = 128;
  static constexpr int s_procs_inited   = 256;

  // Lazy version to ensure s_atm_inited & flag is true for every flag,
  // even if someone adds new flags later on
  static constexpr int s_atm_inited     =  ~0;

  // Utility function to check the ad status
  void check_ad_status (const int flag, const bool must_be_set = true);


  // Current ad initialization status
  int m_ad_status = 0;
};

}  // namespace control
}  // namespace scream

#endif // SCREAM_ATMOSPHERE_DRIVER_HPP
